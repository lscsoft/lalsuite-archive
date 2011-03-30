#!/usr/bin/env python

# ==============================================================================
# Preamble
# ==============================================================================

import os,sys,re,operator,math
from StringIO import StringIO
from glue.segments import segment, segmentlist
from glue.ligolw import ligolw,lsctables,table,utils
from glue.ligolw.utils import segments as ligolw_segments
from glue.segmentdb import query_engine,segmentdb_utils
from pylal import llwapp
from pylal.dq.dqFrameUtils import grab_data
from pylal.dq.dqTriggerUtils import def_get_time

from scipy.stats import poisson

LIGOTimeGPS = lsctables.LIGOTimeGPS

# Some boilerplate to make segmentlists picklable
import copy_reg
copy_reg.pickle(type(segment(0,1)), lambda x:(segment,(x[0],x[1])))
copy_reg.pickle(type(segmentlist([])), lambda x:(segmentlist,([y for y in x],)))

from glue import git_version

__author__  = "Andrew P Lundgren <aplundgr@syr.edu>, Duncan Macleod <duncan.macleod@astro.cf.ac.uk>"
__version__ = "git id %s" % git_version.id
__date__    = git_version.date

"""
This module provides useful segment and veto tools for data quality investigations.
"""

# ==============================================================================
# Function to load segments from an xml file
# ==============================================================================

def fromsegmentxml(file):

  """
    Read a glue.segments.segmentlist from the file object file containing an
    xml segment table.
  """

  xmldoc,digest = utils.load_fileobj(file,gz=file.name.endswith(".gz"))
  seg_table = table.get_table(xmldoc,lsctables.SegmentTable.tableName)

  segs = segmentlist()
  for seg in seg_table:
    segs.append(segment(seg.start_time,seg.end_time))
  segs = segs.coalesce()

  xmldoc.unlink()

  return segs

# ==============================================================================
# Write to segment xml file
# ==============================================================================

def tosegmentxml(file,segs):

  """
    Write the glue.segments.segmentlist object segs to file object file in xml
    format with appropriate tables.
  """

  #== generate empty document
  xmldoc = ligolw.Document()
  xmldoc.appendChild(ligolw.LIGO_LW())
  xmldoc.childNodes[-1].appendChild(lsctables.New(lsctables.ProcessTable))
  xmldoc.childNodes[-1].appendChild(lsctables.New(lsctables.ProcessParamsTable))

  #== append process to table
  process = llwapp.append_process(xmldoc,program='pylal.dq.dqSegmentUtils',\
                                         version=__version__,\
                                         cvs_repository = 'lscsoft',\
                                         cvs_entry_time = __date__)

  gpssegs = segmentlist()
  for seg in segs:
    gpssegs.append(segment(LIGOTimeGPS(seg[0]),LIGOTimeGPS(seg[1])))

  #== append segs and seg definer
  segments_tables = ligolw_segments.LigolwSegments(xmldoc)
  segments_tables.segment_lists.append(ligolw_segments.\
                                       LigolwSegmentList(active=gpssegs))
  #== finalise
  segments_tables.coalesce()
  segments_tables.optimize()
  segments_tables.finalize(process)
  llwapp.set_process_end_time(process)

  #== write file
  utils.write_fileobj(xmldoc,file,gz=False)

# ==============================================================================
# Function to load segments from a csv file
# ==============================================================================

def fromsegmentcsv(csvfile):

  """
    Read a glue.segments.segmentlist object from the file object file containin
    a comma separated list of segments.
  """

  def CSVLineToSeg(line):
    tstart, tend = map(int, line.split(','))
    return segment(tstart, tend)

  segs = segmentlist([CSVLineToSeg(line) for line in csvfile])

  return segs.coalesce()

# ==============================================================================
# Function to parse a segment list for CBC analysable segments
# ==============================================================================

def CBCAnalyzableSegs(seglist):

  """
    Remove any segments shorter than 2064 seconds from seglist because ihope
    won't analyze them.
  """

  return segmentlist([seg for seg in seglist if abs(seg) >= 2064])

# ==============================================================================
# Function to pad a list of segments given start and end paddings
# ==============================================================================

def pad_segmentlist(seglist, start_pad, end_pad):

  """
    Given a veto segmentlist, start pad, and end pad, pads and coalesces the
    segments. Signs of start and end pad are disregarded - the segment is always
    expanded outward.
  """

  padded = lambda seg: segment(seg[0] - abs(start_pad), seg[1] + abs(end_pad))

  seglist = segmentlist([padded(seg) for seg in seglist])

  return seglist.coalesce()

# ==============================================================================
# Function to crop a list of segments
# ==============================================================================

def crop_segmentlist(seglist, end_chop = 30):

  """
    Given a segmentlist and time to chop, removes time from the end of each
    segment (defaults to 30 seconds).
  """

  chopped = lambda seg: segment(seg[0], max(seg[0], seg[1] - end_chop))

  seglist = segmentlist([chopped(seg) for seg in seglist])

  return seglist.coalesce()

# =============================================================================
# Function to return segments in given gps time range
# =============================================================================

def grab_segments(start,end,flag,segment_url='https://segdb.ligo.caltech.edu'):

  """
    Returns a segmentlist containing the segments during which the given flag
    was active in the given period.
  """

  # set times
  start = int( math.floor(start) )
  end   = int( math.ceil(end) )

  # set query engine
  connection        = segmentdb_utils.setup_database(segment_url)
  engine            = query_engine.LdbdQueryEngine(connection)

  # format flag name
  spec = flag.split(':')
  if len(spec) < 2 or len(spec) > 3:
    raise AttributeError, "Included segements must be of the form "+\
                          "ifo:name:version or ifo:name:*"

  ifo     = spec[0]
  name    = spec[1]

  if len(spec) is 3 and spec[2] is not '*':
    version = int(spec[2])
    if version < 1:
      raise AttributeError, "Segment version numbers must be greater than zero"
  else:
    version = '*'

  # expand segment definer
  segdefs = segmentdb_utils.expand_version_number(engine,(ifo,name,version,\
                                                          start,end,0,0))

  # query segs
  segs = segmentdb_utils.query_segments(engine, 'segment', segdefs)
  segs = reduce(operator.or_, segs).coalesce()

  return segs

# =============================================================================
# Function to generate segments for given ifos in period
# =============================================================================

def coinc_segments(start,end,ifos):

  """
    Returns a set of tuple of dictobjects giving singl,double,triple etc time
    science segments for each combination of the contents of the list ifos.
  """

  #== first, construct doubles and triples lists
  doubles=[]
  triples=[]
  for ifo_1 in ifos:
    for ifo_2 in ifos:
      if ifos.index(ifo_2)>ifos.index(ifo_1):
        doubles.append(ifo_1+ifo_2)
        for ifo_3 in ifos:
          if ifos.index(ifo_3)>ifos.index(ifo_2):
            triples.append(ifo_1+ifo_2+ifo_3)

  segments={}
  double_segments={}
  triple_segments={}
  #== grab science data for each ifo
  science_flag = {'G1':'G1:GEO-SCIENCE',\
                  'H1':'H1:DMT-SCIENCE',\
                  'H2':'H2:DMT-SCIENCE',\
                  'L1':'L1:DMT-SCIENCE',\
                  'V1':'V1:ITF_SCIENCEMODE'}
  for ifo in ifos:
    segments[ifo] = grab_segments(start,\
                                  end,\
                                  science_flag[ifo])
  #== grab double-coincidence segments
  for double in doubles:
    ifo_1 = double[0:2]
    ifo_2 = double[2:4]
    double_segments[double] = segments[ifo_1] & segments[ifo_2]

  #== grab triple-coincidence segments
  for triple in triples:
    ifo_1=triple[0:2]
    ifo_2=triple[2:4]
    ifo_3=triple[4:6]
    #== grab concident segments
    triple_segments[triple] = segments[ifo_1] & segments[ifo_2] \
                                              & segments[ifo_3]

  for double in doubles:
    for triple in triples:
      double_segments[double] = double_segments[double]-triple_segments[triple]

  if triples:
    return segments,double_segments,triple_segments
  elif doubles:
    return segments,double_segments
  else:
    return segments

# =============================================================================
# Function to calculate duty cycle and analysable time given segmentlist
# =============================================================================

def duty_cycle(seglist,cbc=False):
  science_time=0
  if cbc:  analysable_time=0
  for segment in seglist:
    seg_length = segment[1]-segment[0]
    science_time+=seg_length
    if cbc and seg_length >= 2048:  analysable_time+=seg_length

  if cbc:  return science_time,analysable_time
  else:  return science_time

# ==============================================================================
# Dump flags from segment database
# ==============================================================================

def dump_flags(ifos=None,segment_url=None,match=None,unmatch=None,latest=False,\
               squery="select ifos,name,version from segment_definer"):

  """
    Returns the list of all flags defined in the database.

    Keyword rguments:
      ifo : [ str | list ]
        list of ifos to query, or str for single ifo
      segment_url : str 
        url of segment database, defaults to contents of S6_SEGMENT_SERVER
        environment variable
      match : [ str | compiled regular expression ]
        regular expression to search against returned flag names, e.g, 'UPV'
      unmatch : str
        regular expression to negatively search against returned flag names
      squery : str
        SQL format query to grab information from the segment database
  """

  if isinstance(ifos,str):
    ifos = [ifos]

  # get url
  if not segment_url:
    segment_url = os.getenv('S6_SEGMENT_SERVER')

  # open connection to LDBD(W)Server
  myClient = segmentdb_utils.setup_database(segment_url)

  reply = StringIO(myClient.query(squery))
  xmldoc,digest = utils.load_fileobj(reply)
  seg_def_table = table.get_table(xmldoc,lsctables.SegmentDefTable.tableName)

  # sort table by ifo,name and version
  seg_def_table.sort(key=lambda flag: (flag.ifos[0],flag.name,\
                                       flag.version),reverse=True)

  flags = lsctables.New(type(seg_def_table))

  for row in seg_def_table:

    # test re match
    if match and not re.search(match,row.name):  continue

    # test re unmatch
    if unmatch and re.search(unmatch,row.name):  continue

    # only append latest versions of multiple flags
    flatest=True
    if latest:
      # get all flags with same ifo and name
      vflags = [f for f in flags if row.name==f.name and\
                row.get_ifos()==f.get_ifos()]
      # if later version in list, move on
      for f in vflags:
        if f.version>=row.version:
          flatest=False
          break
    if not flatest:
      continue

    # append those flags matching ifos requirement
    for ifo in ifos:
      if ifo in row.get_ifos():
        flags.append(row)
        break

  return flags

# ==============================================================================
# Calculate Poisson safety
# ==============================================================================

def poisson_safety( segs, injTable, livetime, returnall=False ):

  """
    Calculates the safety probability of a given segments.segmentlist segs
    based on the number of injections vetoed relative to random chance
    according to Poisson statistics.

    Arguments:

      segs : glue.segments.segmentlist
        list of segments to be tested
      injTable : [ SimBurstTable | SimInspiralTable | SimRingdownTable ]
        table of injections
      livetime : [ float ]
        livetime of search

    Keyword arguments : 

      returnall : [ False | True ]
        return the tuple (numbervetoed, numberexpected, probability), default:
        False returns probability

    """

  deadtime = segs.__abs__()

  get_time = def_get_time( injTable.tableName )
  injvetoed = len([ inj for inj in injTable if get_time(inj) in segs ])

  injexp = len(injTable) * float(deadtime) / float(livetime)

  prob = 1 - poisson.cdf( injvetoed-1, injexp )

  if returnall:
    return injvetoed,injexp,prob
  else:
    return prob

