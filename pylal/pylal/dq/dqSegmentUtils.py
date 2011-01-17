#!/usr/bin/env python

# ==============================================================================
# Preamble
# ==============================================================================

import os,sys,shlex,subprocess,operator,tempfile
from glue.segments import segment, segmentlist
from glue.ligolw import lsctables,table,utils
from glue.segmentdb import query_engine,segmentdb_utils

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

# =============================================================================
# Function to execute shell command and get output
# =============================================================================

def make_external_call(cmd,shell=False):

  """
    Execute shell command and capture standard output and errors. Does not
    support complex commands with pipes, e.g. `echo ${VAR} | grep insp` will
    fail. Returns tuple "(stdout,stderr)".
  """

  args = shlex.split(str(cmd))
  p = subprocess.Popen(args,shell=shell,\
                       stdout=subprocess.PIPE,stderr=subprocess.PIPE)
  p_out, p_err = p.communicate()
  if p.returncode != 0:
    raise ValueError, "Command %s failed. Stderr Output: \n%s" %( cmd, p_err)

  return p_out, p_err

# ==============================================================================
# Function to load segments from an xml file
# ==============================================================================

def fromsegmentxml(file):

  """
    Read a glue.segments.segmentlist from the file object file containing an
    xml segment table.
  """

  xmldoc,digest = utils.load_fileobj(file)
  seg_table = table.get_table(xmldoc,lsctables.SegmentTable.tableName)

  segs = segmentlist()
  for seg in seg_table:
    segs.append(segment(seg.start_time,seg.end_time))
  segs = segs.coalesce()

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

def grab_segments(start,end,flag):

  """
    Returns a segmentlist containing the segments during which the given flag
    was active in the given period.
  """

  # set times
  start = int(start)
  end   = int(end)

  # set query engine
  database_location = os.environ['S6_SEGMENT_SERVER']
  connection        = segmentdb_utils.setup_database(database_location)
  engine            = query_engine.LdbdQueryEngine(connection)

  # format flag name
  spec = flag.split(':')
  if len(spec) < 2 or len(spec) > 3:
    print >>sys.stderr, "Included segements must be of the form ifo:name:version or ifo:name:*"
    sys.exit(1)

  ifo     = spec[0]
  name    = spec[1]
  if len(spec) is 3 and spec[2] is not '*':
    version = int(spec[2])
    if version < 1:
      print >>sys.stderr, "Segment version numbers must be greater than zero"
      sys.exit(1)
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

def dump_flags(start=None,end=None,ifo=None,segment_url=None,\
               squery="select ifos,name,version from segment_definer"):

  """
    Returns the list of all flags defined in the database.

    Keyword rguments:
      start : [ float | int | LIGOTimeGPS ]
        GPS start time of requested period
      end : [ float | int | LIGOTimeGPS ]
        GPS end time of requested period
      ifo : [ "G1" | "H1" | "H2" | "L1" | "V1" ]
      segment_url : string 
        url of segment database, defaults to contents of S6_SEGMENT_SERVER
        environment variable
      squery : string
        SQL format query to grab information from the segment database
  """

  # get url
  if not segment_url:
    segment_url = os.getenv('S6_SEGMENT_SERVER')

  # open connection to LDBD(W)Server
  myClient = segmentdb_utils.setup_database(segment_url)

  tmp = tempfile.TemporaryFile()
  tmp.write(myClient.query(squery))
  tmp.seek(0)
  xmldoc,digest = utils.load_fileobj(tmp)
  tmp.close()
  seg_def_table = table.get_table(xmldoc,lsctables.SegmentDefTable.tableName)

  flags = []
  for line in seg_def_table:
    if ifo and line.ifo!=ifo:  continue
    flag = ':'.join([line.ifos.rstrip(),line.name,str(line.version)])
    flags.append(flag)

  return flags

