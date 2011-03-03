#!/usr/bin/env python

# =============================================================================
# Preamble
# =============================================================================

from __future__ import division
import re,os,sys,numpy,subprocess,datetime,shlex,urlparse,glob
from socket import getfqdn

from glue import segments,git_version
from glue.lal import Cache as LALCache
from glue.lal import CacheEntry as LALCacheEntry

from pylal import Fr,date
from pylal.xlal.datatypes.ligotimegps import LIGOTimeGPS

__author__ = "Duncan Macleod <duncan.macleod@ligo.org>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date

"""
This module provides frame manipulation routines for use in data quality investigations, and cache manipulation routines.
"""

# =============================================================================
# Execute shell command and get output
# =============================================================================

def make_external_call(cmd,shell=False):

  """
    Execute shell command and capture standard output and errors. 
    Returns tuple "(stdout,stderr)".
  """

  if shell:
    args=str(cmd)
  else: 
    args = shlex.split(str(cmd))

  p = subprocess.Popen(args,shell=shell,\
                       stdout=subprocess.PIPE,stderr=subprocess.PIPE)
  p_out, p_err = p.communicate()

  return p_out, p_err

# ==============================================================================
# Grab data from frames
# ==============================================================================

def grab_data(start,end,ifo,channel,type,\
              nds=False,verbose=False,dmt=False):

  """
    This function will return the frame data for the given channel of the given
    type in the given [start,end] time range and will construct a gps time
    vector to go with it. The nds option is not yet supported,
    and the dmt option will return data for dmt channels in frames not found by 
    ligo_data_find.

    >>>grab_data(960000000,960000001,'H1','IFO-SV_STATE_VECTOR','H1_RDS_R_L3')
    ([960000000.0,960000001.0,960000002.0,960000003.0,960000004.0,960000005.0],
     [15.0, 14.125, 13.0, 13.0, 13.0, 13.0])
  """

  time = []
  data = []

  # find FrCheck
  frcheck,err = make_external_call('which FrCheck')
  if err or not frcheck:
    raise ValueError, "FrCheck not found."
  frcheck = frcheck.replace('\n','')

  # generate framecache
  if verbose:
    print >>sys.stdout
    print >>sys.stdout, "Generating framecache..."
    sys.stdout.flush()

  if not dmt:
    cache = generate_cache(start,end,ifo[0:1],type)
  else:
    cache = dmt_cache(start,end,ifo[0:1],type)

  # loop over frames in cache
  for frame in cache:

    # check frame file exists
    if not os.path.isfile(frame.path()):  continue

    # check for Segmentation fault
    segtest = subprocess.Popen([frcheck, "-i", frame.path()],\
                                 stdout=subprocess.PIPE )
    if os.waitpid(segtest.pid,0)[1]==11:  
      print >>sys.stderr, "Warning. Segmentation fault detected with command:"
      print >>sys.stderr, "FrCheck -i "+frame.path()
      continue
    segtest.stdout.close()

    # try to extract data from frame
    frame_data,data_start,_,dt,_,_ = Fr.frgetvect1d(frame.path(),\
                                                      ':'.join([ifo,channel]))
    if frame_data==[]:
      print >>sys.stderr, "No data for %s:%s in %s" % (ifo,channel,frame)
      continue

    # construct time array
    frame_length = float(dt)*len(frame_data)
    frame_time = data_start+dt*numpy.arange(len(frame_data))

    # discard frame data outside of time span
    for i in range(len(frame_data)):
      if frame_time[i] < start:  continue
      if frame_time[i] > end:  continue
      time.append(frame_time[i])
      data.append(frame_data[i])

  return time,data

# ==============================================================================
# Generate framecache in memory given lists of ifos and types
# ==============================================================================

def generate_cache(start,end,ifos,types,framecache=False):

  """
    If framecache==False, this function will return a glue.lal.Cache as
    found by ligo_data_find, otherwise if will return a
    pylal.dq.dqFrameUtils.FrameCache, given start and end time, and lists of
    ifos and types.
  """

  # construct span
  span    = segments.segment(start,end)

  # set options
  if framecache:
    cache = FrameCache()
    entry_class = FrameCacheEntry
    ldfopt = '--frame-cache'
  else:
    cache = LALCache()
    entry_class = LALCacheEntry
    ldfopt = '--lal-cache'

  # find ldf
  exe,err = make_external_call('which ligo_data_find')
  if err:
    print >>sys.stderr, err
  exe = os.path.abspath(exe.replace('\n',''))

  # if given strings, make single-element lists
  if isinstance(ifos,str):
    ifos=[ifos]
  if isinstance(types,str):
    types=[types]

  # loop over each ifo
  for ifo in ifos:
    # loop over each frame type
    for type in types:
      data_find_cmd = ' '.join([exe,'--gps-start-time',str(start),\
                                    '--gps-end-time',str(end),\
                                    '--observatory',ifo[0:1],\
                                    '--type',type,\
                                    '--url-type file',\
                                    '--gaps',\
                                    ldfopt])

      # run ligo_data_find and append each frame to the cache
      out,err = make_external_call(data_find_cmd,shell=True)
      if err:
        print >>sys.stderr, '\n'.join(err.splitlines())

      for line in out.splitlines():
        try:
          e = entry_class(line)
          if span.intersects(e.segment):
            cache.append(e)
        except ValueError:
          print >>sys.stderr, 'Could not convert %s to %s'\
                              % (filename,'.'.join([entry_class.__module__,\
                                                    entry_class.__name__]))

  cache.sort(key=lambda e: e.url)
  return cache

# =============================================================================
# Find types
# =============================================================================

def find_types(types,search='standard'):

  """
    This function will return a valid list of LIGO frame types given the list of
    type strings. The search option defines the breadth of the search, to speed
    up the search, the following search options are supported:
    'standard','short','full'. 

    The 'R', 'T', and 'M' (raw, raw second trends, and raw minute trends) are 
    treated as special cases, so as not to return all types containing those 
    letters. 

    Example:

    >>>find_types('H1_RDS')
    ['H1_RDS_C01_LX',
     'H1_RDS_C02_LX',
     'H1_RDS_C03_L1',
     'H1_RDS_C03_L2',
     'H1_RDS_C03_L2_ET',
     'H1_RDS_C03_L2_ET2',
     'H1_RDS_C03_L2_ET30',
     'H1_RDS_C04_LX',
     'H1_RDS_R_L1',
     'H1_RDS_R_L3',
     'H1_RDS_R_L4']

    >>>find_types(['H1_RDS','R'],search='short')
    ['H1_RDS_R_L1', 'H1_RDS_R_L3', 'H1_RDS_R_L4', 'R']
  """

  # make sure types is a list
  if types is None:  types = []
  if isinstance(types,str):  types = [types]

  # find ldf
  exe,err = make_external_call('which ligo_data_find')
  if err:
    print >>sys.stderr, err
    return []
  exe = os.path.abspath(exe.replace('\n',''))

  # set up search command
  find_cmd = ' '.join([exe,'-y'])

  # treat 'R','M' and 'T' as special cases,
  special_types = ['M','R','T']
  foundtypes = []

  # set list of ignored strings in `ligo_data_find -y`
  # there are thousands of GRBXXXXXX frame types, so ignore them
  if search!='full': 
    iglist = ['GRB']
  if search=='short':
    # all of these strings are part of frame types that can be ignored for a
    # short search
    short_iglist = ['CAL','BRST','Mon','SG','IMR','DuoTone','Concat',\
                    'BH','WNB','Lock','_M','_S5','Multi','Noise','_C0']
    iglist.extend(short_iglist)

  types_out = subprocess.Popen([exe,'-y'],stdout=subprocess.PIPE)
  for t in types_out.stdout.readlines():
    t = t.replace('\n','')
    ignore=False
    for ig in iglist:
      if re.search(ig,t):
        ignore=True
        break
    if ignore:  continue
    intypes=False
    if not types:
      # if no types have been specified, we find all types
      foundtypes.append(t)
    else:
      # else look for special types
      if t in types and t in special_types:
        foundtypes.append(t)
        continue
      # look for everything else
      for type in [tp for tp in types if tp not in special_types]:
        if re.search(type,t):
          foundtypes.append(t)

  types_out.stdout.close()

  return foundtypes

# =============================================================================
# Function to check ifos
# =============================================================================

#def find_ifos(channels,types,ifos):

#  """
#    Constructs an acceptable list of IFOs, parsing channel and frame data type
#    names, and any given ifos.
#
#    Example:
#
#    >>>find_ifos([H1:LSC-DARM_ERR,H0:PEM-ISCT1_ACCZ],
#                 [L1_RDS_R_L1],
#                 [V1])
#    [H0,H1,L1,V1]
#  """

#  accepted_ifos=['H','L','G','V','T',\
#                 'H0','H1','H2',\
#                 'L0','L1',\
#                 'G0','G1',\
#                 'V0','V1']
#
#  if isinstance(channels,str):
#    channels = [channels]
#  if isinstance(types,str):
#    types = [types]
#  if isinstance(ifos,str):
#    ifos = [ifos]
#
#  # if given no ifos, try to generate a list from the channels given  
#  if not ifos:
#    ifos=[]
#    if channels:
#      for channel in channels:
#        if channel.find(':')!=-1:
#          ifo = channel.split(':')[0]
#          if ifo in accepted_ifos and ifo not in ifos:
#            ifos.append(ifo)
#
#    if types:
#      for type in types:
#        ifo = type[0:2]
#        if ifo in accepted_ifos and ifo not in ifos:
#          ifos.append(ifo)
#
#  return ifos

# =============================================================================
# Function to find channels
# =============================================================================

def find_channels(channels=None,\
                  types=None,\
                  ifos=None,\
                  ex_channels=[],\
                  ignore=[],\
                  match=False,\
                  time=None,\
                  unique=False,\
                  verbose=False):

  """
    This function will use FrChannels to return all LIGO data channels matching
    the given list of 'channels' strings, whilst exluding the 'ex_channels'
    strings. Using and find_types() in the same module (if required),
    the search is performed over the given ifos for each given type.

    Use match=True to restrict the search to find channels that
    exactly match the given 'channels' list (i.e. not a partial match).
      Use time=True to search for channels in frame types defined at the given
    epoch.
      Use unique=True to return a unique list of channels, parsed using the
    parse_unique_channels() function, otherwise can return multiple instance of
    the same name string in different types.

    Returns a list of dqFrameUtils.Channel instances.

    Examples:

    >>>channels = find_channels(channels='DARM',types='H1_RDS_R_L1') 
    >>>for channel in channels:
         print channel.name,channel.type,channel.sampling
    LSC-DARM_CTRL H1_RDS_R_L1 16384.0
    LSC-DARM_ERR H1_RDS_R_L1 16384.0
    LSC-DARM_CTRL_EXC_DAQ H1_RDS_R_L1 16384.0
    LSC-DARM_GAIN H1_RDS_R_L1 16.0

    >>>channels = find_channels(channels='DARM_ERR',types=['H1_RDS_R_L1','H1_RDS_R_L3'])
    >>>for channel in channels:  print channel.name,channel.type,channel.sampling
    LSC-DARM_ERR H1_RDS_R_L1 16384.0
    LSC-DARM_ERR H1_RDS_R_L3 16384.0
  
    >>>channels = find_channels(channels='DARM_ERR',types=['H1_RDS_R_L1','H1_RDS_R_L3'],unique=True)
    >>>for channel in channels:  print channel.name,channel.type,channel.sampling
    LSC-DARM_ERR H1_RDS_R_L1 16384.0
  """
 
  # find ldf
  exe,err = make_external_call('which ligo_data_find')
  if err:
    print >>sys.stderr, err
    return []
  exe = os.path.abspath(exe.replace('\n',''))

  # find FrCheck
  frcheck,err = make_external_call('which FrCheck')
  if err or not frcheck:
    raise ValueError, "FrCheck not found."
  frcheck = frcheck.replace('\n','')

  # find types
  if not types:
    types = find_types(types)

  # check list status
  if isinstance(channels,str):
    channels = [channels]
  if isinstance(types,str):
    types = [types]
  if isinstance(ifos,str):
    ifos = [ifos]
  found_channels=[]

  # loop over each ifo
  for ifo in ifos:
    # set ligo_data_find frame search time
    if time is None:
      time = date.XLALUTCToGPS(datetime.datetime.now().timetuple())-(2*86400)
    if verbose:
      print_statement = \
          "Searching %s frame types for " % len(types)
      if channels:
        print_statement += ', '.join(channels)
      else:
        print_statement += "all channels"
      print_statement += ", in ifo %s" % ifo
      print >>sys.stdout, print_statement

    for type in types:
      count=0
      # skip empty frame types or those set for ignorance
      if type in ignore:  continue
      if type == '':  continue

      if verbose:
        print >>sys.stdout, "  Searching "+str(type)+"...",
      sys.stdout.flush()

      # find first frame file for type
      frame_cmd = ' '.join([exe,'--observatory',ifo[0:1],\
                                '--type',type,\
                                '--gps-start-time',str(time),\
                                '--gps-end-time',str(time),\
                                '--url-type','file',\
                                '--gaps'])
      out,err = make_external_call(frame_cmd,shell=True)

      if verbose and err: 
        print >>sys.stderr, '%s...' % ('\n'.join(err.splitlines())),

      frame=''
      for line in out.splitlines():
        if line.startswith('file://'):
          frame = line
          break

      # if frame is found:
      if frame:
        info = frame.split(' ')
        frame = urlparse.urlparse(info[-1])[2]

        # test frame for seg fault
        segtest = subprocess.Popen([frcheck,"-i",frame],stdout=subprocess.PIPE)
        if os.waitpid(segtest.pid,0)[1]==11:
          if verbose:
            print >>sys.stderr, "  Warning. Segmentation fault detected with "+\
                                "command:"
            print >>sys.stderr, "    FrCheck -i %s" % frame
          continue
        segtest.stdout.close()

        # get channels contained in frame, grepping for input channel string
        frchannels,err = make_external_call('FrChannels %s' % frame)
        if err:
          print >>sys.stderr, "  Failed to find channels for type %s," % (type),
          print >>sys.stderr, "using the following frame\n%s" % (frame)
          continue

        for line in frchannels.splitlines():
          # check match with ifo
          if not re.match(ifo,line):
            continue

          # split FrChannels output
          name,sampling = line.split(' ')

          # if channel matches any exlusion strings, skip it
          if ['match' for exchan in ex_channels if re.search(exchan,name)]:
            continue

          # if asked for exact match, check:
          if match and (name in channels):
            pass
          elif channels and\
              not ['match' for ch in channels if re.match('%s\Z' % ch,name)]:
            continue

          # generate structure and append to list  
          found_channel = Channel(name,type=type,sampling=sampling)
          found_channels.append(found_channel)
          count+=1
          sys.stdout.flush()

      # print channel count for data type
      if verbose:  print >>sys.stdout, count,"channels found"
    if verbose:  print >>sys.stdout

  if unique:
    found_channels = parse_unique_channels(found_channels)

  return found_channels

# ==============================================================================
# parse channels for uniqueness
# ==============================================================================

def parse_unique_channels(channels,type_order=None):

  """
    This function will parse a list of dqFrameUtils.Channel instances into a 
    unique list based on the given type_order, or the internal default (based on
    the S6 frame type convention), e.g for H1: 

    type_order=[H1_RDS_R_L1,H1_RDS_R_L3,R]
  
    Multiple instances of the same channel name will be tried for each type in 
    the order, if not matched the first instance will be chosen. This is an
    attempt to maximise performance by picking frame types that contain the
    least data.

    Example:
  
    If the list 'channels' includes H1:LSC-DARM_ERR from both the 'R' and 
    'H1_RDS_R_L3' frame types, the latter instance will be chosen.
  """

  # if given empty, return empty and hope user notices
  if channels==[]:
    return []

  # set up type preference order
  if type_order == None:
    type_order = ['H1_RDS_R_L1','L1_RDS_R_L1','H2_RDS_R_L1','R','RDS_R_L1']
  # sort channels by name
  channels.sort(key=lambda ch: ch.name)
  # set up loop variables
  uniq_channels = []
  channel_cluster = []
  previous_channel = Channel('tmp',type='tmp',sampling=0)
  # loop over channels
  for channel in channels:
    # if channel does not match previous channel, process previous cluster
    if channel.name != previous_channel.name:
      if channel_cluster!=[]:
        chosen = False
        # loop over types and channels to pick channel highest in type order
        for type in type_order:
          for element in channel_cluster:
            if element.type==type:
              chosen_channel = element
              chosen = True
              break
          if chosen == True:  break
        # if no channel found, take first in cluster
        if chosen == False:
          chosen_channel = channel_cluster[0]
          chosen = True
        # append chosen_channel to uniq channel list
        uniq_channels.append(chosen_channel)
      # reset the cluster to be the current working channel
      channel_cluster = [channel]

    # if channel name matches previous one, append to cluster
    else:
      channel_cluster.append(channel)
    # update previous_channel indicator
    previous_channel = channel

  # when loop is complete, analyse the final remaining cluster
  chosen = False
  for type in type_order:
    for element in channel_cluster:
      if element.type==type:
        chosen_channel = element
        chosen = True
        break
    if chosen == True:  break
  # if no channel found, take first in cluster
  if chosen == False:
    chosen_channel = channel_cluster[0]
    chosen = True
  # append final channel to uniq channel list
  uniq_channels.append(chosen_channel)

  return uniq_channels

# ==============================================================================
# Class to generate channel structure
# ==============================================================================

class Channel:
  """
  The Channel class defines objects to represent LIGO data channels. Each Channel
  has a 'name' attribute and can be assigned 'type' and 'sampling' attributes if
  relevant.

  Example:

  >>>GWChannel = Channel('H1:LSC-DARM_ERR,'H1_RDS_R_L3',4096)
  >>>GWChannel.name, GWChannel.type, GWChannel.sampling
  ('H1:LSC-DARM_ERR', 'H1_RDS_R_L3', 4096)
  """

  def __init__(self,name,type=None,sampling=None):
    """Initialise the dqFrameUtils.Channel object assuming the give name follows the standard convention:

IFO:SYSTEM-SUBSYSTEM_SIGNAL

For example:

c = dqDataUtils.Channel('H1:LSC-DARM_ERR')

will return

c.ifo = 'H1'
c.name = 'LSC-DARM_ERR'
c.site = 'H'
c.system = 'LSC'
c.subsystem = 'DARM'
c.signal = 'ERR'
    """

    attributes = ['ifo','site','name',\
                  'system','subsystem','signal',\
                  'type',\
                  'sampling']

    __slots__ = attributes

    # extract name attributes
    if ':' in name:
      self.ifo,self.name         = re.split(':',name,maxsplit=1)
      self.site                  = self.ifo[0]
    else:
      self.name = name

    tags = re.split('[-_]',self.name,maxsplit=3)

    self.system = tags[0]

    if len(tags)>1:
      self.subsystem = tags[1]

    if len(tags)>2:
      self.signal    = tags[2]

    if type:
      self.type = str(type)

    if sampling:
      self.sampling = float(sampling)

  def __getattribute__(self,name):

    return self.__dict__[name]

  def __str__( self ):

    return '%s:%s' % ( self.ifo, self.name )

# ==============================================================================
# Function to generate a framecache of /dmt types
# ==============================================================================
def dmt_cache(start,end,ifo,type,framecache=False):
  """
  This function will return a list of frame files in the given start and stop 
  time interval for the give IFO using the given DMT frame type. This is
  required if ligo_data_find will not return the dmt frames.

  Example:

  >>>dmt_cache(960000000,960010000,'H1','LockLoss_H1')
  ['/archive/frames/dmt/LHO/LockLoss_H1/H-M-960/H-LockLoss_H1_M-960001200-3600.gwf',
   '/archive/frames/dmt/LHO/LockLoss_H1/H-M-960/H-LockLoss_H1_M-960004800-3600.gwf']
  """

  # find dmt frames path
  host = getfqdn()
  if re.search('ligo-',host):
    dmt_dir = '/dmt'
  elif re.search('ligo.',host):
    site = {'H':'LHO','L':'LLO','V':'V1'}
    dmt_dir = os.path.join('/archive','frames','dmt',site[ifo[0]])


  span    = segments.segment(start,end)

  # set options
  if framecache:
    cache = FrameCache()
    entry_class = FrameCacheEntry
    ldfopt = '--frame-cache'
  else:
    cache = LALCache()
    entry_class = LALCacheEntry
    ldfopt = '--lal-cache'

  basedir = os.path.join(dmt_dir,type)

  # frames are 3600 seconds long, so round
  tmp = int(str(start)[0:3]+'000000')
  cache_start = tmp+3600*int((start-tmp)/3600)
  tmp = int(str(end)[0:3]+'000000')
  cache_end = tmp+3600*int((end-tmp)/3600)

  # directories are listed with three time digits
  start_three = int(str(cache_start)[0:3])
  end_three = int(str(cache_end)[0:3])
  first_three = numpy.arange(start_three,end_three+1,1)

  #loop over directories
  for t in first_three:
    querydir = os.path.join(basedir,'%s-%s-%s' % (ifo[0:1],'M',str(t)),'*')
    filenames = glob.glob(querydir)
    for filename in filenames:
      try:
        e = entry_class.from_T050017(filename)
        if span.intersects(e.segment):  cache.append(e)
      except ValueError:
        print >>sys.stderr, 'Could not convert %s to %s'\
                            % (filename,'.'.join([entry_class.__module__,\
                                                    entry_class.__name__]))

  return cache

# ==============================================================================
# Class for wCacheEntry
# ==============================================================================

class FrameCacheEntry(LALCacheEntry):
  """
    An object representing one line in a frame cache file.

    Each line in a frame cache identifies multiple files, and the line consists
    of six columns of white-space delimited text.

    The first column, "observatory", generally stores the name of an 
    observatory site or one or more instruments (preferably delimited by ",",
    but often there is no delimiter between instrument names).

    The second column, "description", stores a short string tag that is
    usually all capitals with "_" separating components, in the style of the
    description part of the LIGO-Virgo frame filename format.

    The third and fourth columns store the start time and stop time in GPS
    seconds of the interval spanned by the file identified by the cache line.

    The fifth column stored the duration of each frame identified in the cache
    line. 

    The sixth (last) column stores the file's URL.

    The values for these columns are stored in the .observatory,
    .description, .segment and .url attributes, respectively.  The
    .segment attribute stores a glue.segments.segment object describing
    the interval spanned by the file.  Any of these attributes except
    the URL is allowed to be None.
  """

  # How to parse a line in a frame cache file. Six white-space
  # delimited columns.
  _regex = re.compile(r"\A\s*(?P<obs>\S+)\s+(?P<dsc>\S+)\s+(?P<strt>\S+)\s+(?P<end>\S+)\s+(?P<dur>\S+)\s+(?P<url>\S+)\s*\Z")
  _url_regex = re.compile(r"\A((.*/)*(?P<obs>[^/]+)-(?P<dsc>[^/]+)-(?P<strt>[^/]+)-(?P<dur>[^/\.]+)\.[^/]+)\Z")

  def __init__(self, *args, **kwargs):
    """
    Intialize a FrameCacheEntry object. The arguments can take two
    forms:  a single string argument, which is interpreted and
    parsed as a line from a frame cache file, or four arguments
    used to explicitly initialize the observatory, description,
    segment and URL in that order.  When parsing a single line
    of text from a frame cache, an optional key-word argument
    "coltype" can be provided to set the type the start, end and
    durations are parsed as.  The default is glue.lal.LIGOTimeGPS.

    """
    if len(args) == 1:
      # parse line of text as an entry in a cache file
      match = self._regex.search(args[0])
      coltype = kwargs.pop("coltype", LIGOTimeGPS)

      try:
        match = match.groupdict()
      except AttributeError:
        raise ValueError, "could not convert %s to FrameCacheEntry"\
                          % repr(args[0])
      self.observatory = match["obs"]
      self.description = match["dsc"]
      start            = match["strt"]
      end              = match["end"]
      self.duration    = coltype(match["dur"])

      if start == "-" and end == "-":
        # no segment information
        self.segment = None
      else:
        self.segment = segments.segment(coltype(start),coltype(end))
      self.url = match["url"]

      if kwargs:
        raise TypeError, "unrecognized keyword arguments: %s" % ", ".join(kwargs)
    elif len(args) == 5:
      # parse arguments as observatory, description,
      # segment, duration, url
      if kwargs:
        raise TypeError, "invalid arguments: %s" % ", ".join(kwargs)
      self.observatory, self.description, self.segment, self.duration, self.url\
          = args
    else:
      raise TypeError, "invalid arguments: %s" % args

    # "-" indicates an empty column
    if self.observatory == "-":
      self.observatory = None
    if self.description == "-":
      self.description = None

  def __str__(self):
    """
    Convert the FrameCacheEntry to a string in the format of a line
    in a frame cache. Used to write the FrameCacheEntry to a file.

    """
    if self.segment is not None:
      start,end = [str(t) for t in self.segment]
    else:
      start    = "-"
      end      = "-"
      duration = "-"

    return "%s %s %s %s %s %s" % (self.observatory or "-", self.description or "-", start, end, self.duration, self.url)

  def __cmp__(self, other):
    """
    Compare two FrameCacheEntry objects by observatory, then
    description, then segment, then duration, then URL.
    """
    if type(other) != FrameCacheEntry:
      raise TypeError, "can only compare FrameCacheEntry to FrameCacheEntry"
    return cmp((self.observatory, self.description, self.segment, self.duration, self.url), (other.observatory, other.description, other.segment, other.duration, other.url))

  def get_files(self):
    """
    Return Find all files described by this FrameCacheEntry.
    """

    filenames = glob.glob(os.path.join(self.path(),'*'))
    cache = [e.path() for e in\
                 LALCache([LALCacheEntry.from_T050017(f) for f in filenames])\
             if e.observatory==self.observatory and\
                e.description==self.description and\
                self.segment.intersects(e.segment) and\
                abs(e.segment)==self.duration]

    return cache
    
  def from_T050017(cls, url, coltype = LIGOTimeGPS):      

    """      
    Parse a URL in the style of T050017-00 into a FrameCacheEntry.      
    The T050017-00 file name format is, essentially,  
  
    observatory-description-start-dur.ext  
  
    """      
    match = cls._url_regex.search(url)      
    if not match:      
            raise ValueError, "could not convert %s to CacheEntry" % repr(url)      
    observatory = match.group("obs")      
    description = match.group("dsc")      
    start = match.group("strt")      
    duration = match.group("dur")      
    if start == "-" and duration == "-":      
            # no segment information      
            segment = None      
    else:      
            segment = segments.segment(coltype(start), coltype(start) + coltype(duration))      
    return cls(observatory, description, segment, duration,\
               os.path.split(url)[0])    


  from_T050017 = classmethod(from_T050017)


# ==============================================================================
# Class for FrameCacheEntry 
# ==============================================================================

class FrameCache(LALCache):
  """
    An object representing a frame cache file. Currently it is possible to
    add anything to a FrameCache. This method should check that the thing you
    are adding is a FrameCacheEntry and throw and error if it is not.
  """
  entry_class = FrameCacheEntry
  
  def sieve(self, ifos=None, description=None, segment=None, duration=None,
    exact_match=False):
    """
    Return a FrameCache object with those FrameCacheEntries that contain the
    given patterns (or overlap, in the case of segment).  If
    exact_match is True, then non-None ifos, description, and
    segment patterns must match exactly.
    
    Bash-style wildcards (*?) are allowed for ifos and description.
    """
    if exact_match:
      segment_func = lambda e: e.segment == segment
    else:
      if ifos is not None: ifos = "*" + ifos + "*"
      if description is not None: description = "*" + description + "*"
      segment_func = lambda e: segment.intersects(e.segment)
    
    c = self
    
    if ifos is not None:
      ifos_regexp = re.compile(fnmatch.translate(ifos))
      c = [entry for entry in c if ifos_regexp.match(entry.observatory) is not None]
    
    if description is not None:
      descr_regexp = re.compile(fnmatch.translate(description))
      c = [entry for entry in c if descr_regexp.match(entry.description) is not None]
    
    if segment is not None:
      c = [entry for entry in c if segment_func(entry)]
    
    if duration is not None:
      c = [entry for entry in c if entry.duration==duration]

    return self.__class__(c)

def FrameCachetoLALCache( fcache ):

  lcache = FrameCache()

  for e in fcache:
    files = e.get_files()
    for f in files:
      lcache.append(LALCacheEntry.from_T050017( f ))
  

  return 0

def LALCachetoFrameCache( lcache ):

  lcache.sort( key=lambda e: (e.path(),e.segment[0]) )

  fcache = FrameCache()

  for e in lcache:

    matched = False

    dir = os.path.split(e.path())[0]

    # if path found in FrameCache try to coalesce with other entries
    dirs = [d.path() for d in fcache]

    if dir in dirs:

      pathentries = [fe for fe in fcache if fe.path()==dir]

      # test current entry against other entries for the same path in the
      # new cache
      for i,pe in enumerate( pathentries ):

        notdisjoint = e.segment[0] <= pe.segment[1]

        # if entries match in directory, duration, and are contiguous, append
        if pe.path()==os.path.split(e.path())[0]\
             and e.segment.__abs__() == pe.duration\
             and notdisjoint:
          seg = segments.segment( min(pe.segment[0], e.segment[0]),\
                                  max(pe.segment[1], e.segment[1]) )
          fcache[i].segment = seg

          matched = True
          break

      # if we haven't matched the entry to anything already in the cache add now
      if not matched:
        fe = FrameCacheEntry.from_T050017( e.path() )
        fcache.append(fe)

    # if from a new directory add
    else:
      fe = FrameCacheEntry.from_T050017( e.path() )
      fcache.append(fe)

  return fcache
