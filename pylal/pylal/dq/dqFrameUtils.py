#!/usr/bin/env python

from __future__ import division

from glue import git_version
import re
import os
import sys
from pylal import Fr
import numpy
import subprocess
__author__ = "Duncan Macleod <duncan.macleod@ligo.org>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date

"""
Module providing frame data utilities for DQ
"""

# =============================================================================
# Function to execute shell command and get output
# =============================================================================
def GetCommandOutput(command):
  # == function to execute bash commands and return the stdout and error status
  stdin, out, err = os.popen3(command)
  pid, status = os.wait()
  this_output = out.read()
  stdin.close()
  out.close()
  err.close()
  return this_output, status

# ==============================================================================
# Function to grab data from frames
# ==============================================================================
def grab_data(start,end,channel,type,\
              nds=False,verbose=False,dmt=False):
  """
  This function will return the frame data for the given channel of the given
  type in the given [start,end] time range and will construct a gps time vector
  to go with it. The nds option is not yet supported,
  and the dmt option will return data for dmt channels in frames not found by 
  ligo_data_find.

  >>>grab_data(960000000,960000001,'H1:IFO-SV_STATE_VECTOR','H1_RDS_R_L3')
  ([960000000.0,960000001.0,960000002.0,960000003.0,960000004.0,960000005.0],
   [15.0, 14.125, 13.0, 13.0, 13.0, 13.0])
  """

  time = []
  data = []

  #== find FrCheck
  p = subprocess.Popen(["which", "FrCheck"], stdout=subprocess.PIPE)
  frcheck = p.communicate()[0].replace('\n','')
  if p.returncode != 0:
    raise ValueError, "FrCheck not found."
  p.stdout.close()
  #== generate framecache
  if verbose:
    print >>sys.stdout, "Generating framecache..."
    sys.stdout.flush()
  if not dmt:
    cache = generate_cache(start,end,channel[0:1],type,return_files=True)
  else:
    cache = dmt_cache(start,end,channel[0:1],type)
  #== loop over frames in cache
  for frame in cache:
    #== check frame file exists
    if not os.path.isfile(frame):  continue
    #== check for Segmentation fault
    segtest = subprocess.Popen([frcheck,"-i",frame],stdout=subprocess.PIPE)
    if os.waitpid(segtest.pid,0)[1]==11:  
      print >>sys.stderr, "Warning. Segmentation fault detected with command:"
      print >>sys.stderr, "FrCheck -i "+frame
      continue
    segtest.stdout.close()
    #== try to extract data from frame
    try:
      frame_data,data_start,_,dt,_,_ = Fr.frgetvect1d(frame,channel)
      if frame_data==[]:
        print >>sys.stderr, "No data for "+channel+" in "+frame
        continue
      #== construct time array
      frame_length = float(dt)*len(frame_data)
      frame_time = data_start+dt*numpy.arange(len(frame_data))
      #== discard frame data outside of time span
      for i in range(len(frame_data)):
        if frame_time[i] < start:  continue
        if frame_time[i] > end:  continue
        time.append(frame_time[i])
        data.append(frame_data[i])
    except:
      print >>sys.stderr, "Failed to access frame:\n"+frame
      continue
  return time,data

# ==============================================================================
# Function to generate framecache in memory given lists of ifos and types
# ==============================================================================
def generate_cache(start_time,end_time,ifos,types,return_files=False):
  """
  This function will return a cache of files as found by ligo_data_find,
  given start and end time, and lists of ifos and types. If the return_files 
  option is given as 'True' the function will return a list of frames with 
  absolute paths, otherwise it will return a frame cache (as used in wpipline, 
  for example).

  Example:

  >>>generate_cache(961977615,962582415,R,H)
  ['H R 961977600 962000000 32 /archive/frames/S6/L0/LHO/H-R-9619'
   'H R 962000000 962064032 32 /archive/frames/S6/L0/LHO/H-R-9620']
  >>>generate_cache(961977615,962582415,R,H,return_files=True)
  [/archive/frames/S6/L0/LHO/H-R-9619/H-R-961977600-32.gwf,
   /archive/frames/S6/L0/LHO/H-R-9619/H-R-961977632-32.gwf,
   ...
   /archive/frames/S6/L0/LHO/H-R-9620/H-R-962064000-32.gwf]
  """

  cache = []

  #== find ldf
  p = subprocess.Popen(["which", "ligo_data_find"], stdout=subprocess.PIPE)
  ldf = p.communicate()[0].replace('\n','')
  if p.returncode != 0:
    raise ValueError, "ligo_data_find"
  p.stdout.close()
  ldf = os.path.realpath(ldf)
  #== if given strings, make single-element lists
  if isinstance(ifos,str):
    ifos=[ifos]
  if isinstance(types,str):
    types=[types]
  #== loop over each ifo
  for ifo in ifos:
    #== loop over each frame type
    for type in types:
      try:
        data_find_cmd =   ldf+\
                          ''' --gps-start-time '''+str(start_time)+\
                          ''' --gps-end-time '''+str(end_time)+\
                          ''' --observatory '''+ifo[0:1]+\
                          ''' --type '''+type+\
                          ''' --url-type file '''+\
                          ''' --frame-cache | sort'''
        #== run ligo_data_find and append each frame to the cache
        cache_out = subprocess.Popen(data_find_cmd,shell=True,\
                                     stdout=subprocess.PIPE)
        for line in cache_out.stdout.readlines():
          #== if line is not recognised in standard frame cache format, skip
          if len(line.split(' '))!=6:
            continue
          cache.append(line.replace('\n',''))
        cache_out.stdout.close()
      except:
        continue
  #== if no files:
  if cache==[]:
    print >>sys.stderr, "Warning: no frames found."
  #== if asked for the files, expand the cache
  if return_files:
    cache = expand_cache(cache)

  return cache

# =============================================================================
# Function to find types
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

  #== find ldf
  p = subprocess.Popen(["which", "ligo_data_find"], stdout=subprocess.PIPE)
  ldf = p.communicate()[0].replace('\n','')
  if p.returncode != 0:
    raise ValueError, "ligo_data_find"
  p.stdout.close()
  ldf = os.path.realpath(ldf)

  #== make sure types is a list
  if types is None:  types = []
  if isinstance(types,str):  types = [types]

  #== set up search command
  find_cmd = ldf+" -y | egrep "

  #== treat 'R','M' and 'T' as special cases,
  #== so not to grep for all types containing 'R'
  special_types = ['M','R','T']
  special_cases=[]

  #== set list of ignored strings in `ligo_data_find -y`
  #== there are thousands of GRBXXXXXX frame types, so ignore them
  if search!='full': 
    vgrep_list = ['GRB']
  if search=='short':
    #== all of these strings are part of frame types that can be ignored for a
    #== short search
    short_ignore_list = ['CAL','BRST','Mon','SG','IMR','DuoTone','Concat',\
                         'BH','WNB','Lock','_M','_S5','Multi','Noise','_C0']
    vgrep_list.extend(short_ignore_list)
  #== add each of those ignored strings to a vgrep command
  find_cmd+="-v '"
  for vstring in vgrep_list:
    find_cmd+=vstring+'|'
  #== take off last '|'
  find_cmd = find_cmd[0:-1] + "'"

  #== if given types
  if types:  find_cmd+=' | egrep '
  for type in types:
    #== if type is one of the special cases, save.
    if type in special_types:
      special_cases.append(type)
    #== otherwise add it to the grep command
    else:
      if not find_cmd.endswith('|'):
        find_cmd+= "'"
      find_cmd+=type+"|"
  #== take of the extra character
  if find_cmd[-1]=="|":
    find_cmd = find_cmd[0:-1] + "'"
  found_types = []
  #== if not searching only for special types, run the grep command
  if find_cmd != ldf+" -y | egrep '":
    found_types_out = subprocess.Popen(find_cmd,shell=True,\
                                       stdout=subprocess.PIPE)
    for line in found_types_out.stdout.readlines():
      if line=='\n':  continue
      found_type = line.replace('\n','')
      found_types.append(found_type)
    found_types_out.stdout.close()
  #== append all special cases to the list
  for type in special_cases:
    found_types.append(type)
  if found_types == ['']:
    print >>sys.stderr, "No data types found, exiting."
  return found_types

# =============================================================================
# Function to check ifos
# =============================================================================
def find_ifos(channels,types,ifos):
  """
  Constructs an acceptable list of IFOs, parsing channel and frame data type
  names, and any given ifos.

  Example:

  >>>find_ifos([H1:LSC-DARM_ERR,H0:PEM-ISCT1_ACCZ],
               [L1_RDS_R_L1],
               [V1])
  [H0,H1,L1,V1]
  """

  accepted_ifos=['H','L','G','V','T',\
                 'H0','H1','H2',\
                 'L0','L1',\
                 'G0','G1',\
                 'V0','V1']

  if isinstance(channels,str):
    channels = [channels]
  if isinstance(types,str):
    types = [types]
  if isinstance(ifos,str):
    ifos = [ifos]

  #== if given no ifos, try to generate a list from the channels given  
  if ifos is None:
    ifos=[]
    if channels is not None:
      for channel in channels:
        if channel.find(':')!=-1:
          ifo = channel.split(':')[0]
          if ifo in accepted_ifos and ifo not in ifos:
            ifos.append(ifo)

    if types is not None:
      for type in types:
        ifo = type[0:2]
        if ifo in accepted_ifos and ifo not in ifos:
          ifos.append(ifo)

    if ifos==[]:
      user_ask_string = \
          "Cannot determine IFO from channel names or frame types names.\n"+\
          "The list of accepted IFOs is:\n\n"
      for ifo in accepted_ifos:
        user_ask_string += ifo+','
      user_ask_string = user_ask_string[0:-1]+"\n\n"
      user_ask_string += "Please enter the releveant IFOs (comma separated):"
      user_input = raw_input(user_ask_string)
      ifos=user_input.split(',')
      if isinstance(ifos,str):
        ifos=[ifos]

  for ifo in ifos:
    if ifo not in accepted_ifos:
      print >>sys.stderr, 'Cannot recognise IFO: '+ifo+'. Please try again'
      ifos.remove(ifo)

  return ifos

# =============================================================================
# Function to find channels
# =============================================================================
def find_channels(channels=None,\
                  types=None,\
                  ifos=None,\
                  ex_channels=None,\
                  ignore=[],\
                  match=False,\
                  time=None,\
                  unique=False,\
                  verbose=False):

  """
  This function will use FrChannels to return all LIGO data channels matching
  the given list of 'channels' strings, whilst exluding the 'ex_channels'
  strings. Using find_ifos() and find_types() in the same module (if required),
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
  >>>for channel in channels:  print channel.name,channel.type,channel.sampling
  H1:LSC-DARM_CTRL H1_RDS_R_L1 16384.0
  H1:LSC-DARM_ERR H1_RDS_R_L1 16384.0
  H1:LSC-DARM_CTRL_EXC_DAQ H1_RDS_R_L1 16384.0
  H1:LSC-DARM_GAIN H1_RDS_R_L1 16.0

  >>>channels = find_channels(channels='DARM_ERR',types=['H1_RDS_R_L1','H1_RDS_R_L3'])
  >>>for channel in channels:  print channel.name,channel.type,channel.sampling
  H1:LSC-DARM_ERR H1_RDS_R_L1 16384.0
  H1:LSC-DARM_ERR H1_RDS_R_L3 16384.0
  
  >>>channels = find_channels(channels='DARM_ERR',types=['H1_RDS_R_L1','H1_RDS_R_L3'],unique=True)
  >>>for channel in channels:  print channel.name,channel.type,channel.sampling
  H1:LSC-DARM_ERR H1_RDS_R_L1 16384.0
  """
 
  #== find ldf
  p = subprocess.Popen(["which", "ligo_data_find"], stdout=subprocess.PIPE)
  ldf = p.communicate()[0].replace('\n','')
  if p.returncode != 0:
    raise ValueError, "ligo_data_find"
  p.stdout.close()
  ldf = os.path.realpath(ldf)

  #== cannot work with no ifos
  if ifos is None:
    ifos = find_ifos(channels,types,ifos)

  if types is None:
    types = find_types(types)

  #== check list status
  if isinstance(channels,str):
    channels = [channels]
  if isinstance(types,str):
    types = [types]
  if isinstance(ifos,str):
    ifos = [ifos]
  found_channels=[]
  #== loop over each ifo
  for ifo in ifos:
    #== set ligo_data_find frame search time
    if time is None:
      time = \
          str(GetCommandOutput('tconvert now -2 days')[0]).replace('\n','')
    if verbose:
      print_statement = \
          "Searching "+str(len(types))+" frame types for: "
      if channels is not None:
        for channel in channels:
          print_statement += channel+', '
        print_statement += " in ifo "+ifo
      else:
        print_statement+= "all channels, in ifo "+ifo
      print print_statement

    for type in types:
      count=0
      #== skip empty frame types or those set for ignorance
      if type in ignore:  continue
      if type == '':  continue

      if verbose:
        print >>sys.stdout, "  Searching "+str(type)+"...",
      sys.stdout.flush()

      #== find first frame file for type
      frame_cmd = ldf+''' --observatory '''+ifo[0:1]+\
                  ''' --type='''+type+\
                  ''' --gps-start-time '''+str(time)+\
                  ''' --gps-end-time '''+str(time)+\
                  ''' --url-type file'''
      frame_out = subprocess.Popen(frame_cmd,shell=True,stdout=subprocess.PIPE,\
                                   stderr=subprocess.PIPE)
      frame_status = 0
      frame=''
      for line in frame_out.stdout.readlines():
        if line[0:7]=='file://':  frame = line
        break
      frame_out.stdout.close()
      frame = frame.replace('\n','')
      #== if frame is found:
      if frame_status == 0 and frame != "":
        info = frame.split(' ')
        frame = info[-1].replace('file://localhost','')
        #== get channels contained in frame, grepping for input channel string
        channel_find_cmd = "FrChannels "+frame+" | grep "+ifo
        #== add grep options for each included channel
        if channels is not None:
          channel_find_cmd += " | egrep '"
          for channel in channels:
            channel_find_cmd += channel+"|"
          channel_find_cmd = channel_find_cmd[0:-1]+"'"
        #== add grep options for each excluded channel
        if ex_channels is not None:
          channel_find_cmd += " | egrep -v '"
          for ex_channel in ex_channels:
            channel_find_cmd += ex_channel+"|"
          channel_find_cmd = channel_find_cmd[0:-1]+"'"
        #== grab channels
        try:
          channel_list_out = subprocess.Popen(channel_find_cmd,shell=True,\
                                              stdout=subprocess.PIPE)
          for line in channel_list_out.stdout.readlines():
            data = line.replace('\n','')
            name,sampling = data.split(' ')
            #== if asked for exact match, check:
            if match:
              if name not in channels:  continue
            #== generate structure and append to list  
            found_channel = Channel(name,type=type,sampling=sampling)
            found_channels.append(found_channel)
            count+=1
            sys.stdout.flush()
          channel_list_out.stdout.close()
        except:
          print "  Failed to find channels for type "+type+", using the"+\
              " following frame\n"+frame
          continue
      #== print channel count for data type
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
  attempt to maximise performance by picking frame types that contain the least
  data.

  Example:
  
  If the list 'channels' includes H1:LSC-DARM_ERR from both the 'R' and 
  'H1_RDS_R_L3' frame types, the latter instance will be chosen.
  """
  #== if given empty, return empty and hope user notices
  if channels==[]:
    return []

  #== set up type preference order
  if type_order == None:
    type_order = ['H1_RDS_R_L1','L1_RDS_R_L1','H2_RDS_R_L1','R','RDS_R_L1']
  #== sort channels by name
  channels = sorted(channels,key=lambda ch: ch.name)
  #== set up loop variables
  uniq_channels = []
  channel_cluster = []
  previous_channel = Channel('tmp',type='tmp',sampling=0)
  #== loop over channels
  for channel in channels:
    #== if channel does not match previous channel, process previous cluster
    if channel.name != previous_channel.name:
      if channel_cluster!=[]:
        chosen = False
        #== loop over types and channels to pick channel highest in type order
        for type in type_order:
          for element in channel_cluster:
            if element.type==type:
              chosen_channel = element
              chosen = True
              break
          if chosen == True:  break
        #== if no channel found, take first in cluster
        if chosen == False:
          chosen_channel = channel_cluster[0]
          chosen = True
        #== append chosen_channel to uniq channel list
        uniq_channels.append(chosen_channel)
      #== reset the cluster to be the current working channel
      channel_cluster = [channel]

    #== if channel name matches previous one, append to cluster
    else:
      channel_cluster.append(channel)
    #== update previous_channel indicator
    previous_channel = channel

  #== when loop is complete, analyse the final remaining cluster
  chosen = False
  for type in type_order:
    for element in channel_cluster:
      if element.type==type:
        chosen_channel = element
        chosen = True
        break
    if chosen == True:  break
  #== if no channel found, take first in cluster
  if chosen == False:
    chosen_channel = channel_cluster[0]
    chosen = True
  #== append final channel to uniq channel list
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
    self.name = str(name)
    if type is not None:
      self.type = str(type)
    if sampling is not None:
      try:
        if sampling==int(sampling):
          self.sampling = int(sampling)
        else:
          self.sampling = float(sampling)
      except:
        self.sampling = float(sampling)

# ==============================================================================
# Function to generate a framecache of /dmt types
# ==============================================================================
def dmt_cache(start,end,ifo,type):
  """
  This function will return a list of frame files in the given start and stop 
  time interval for the give IFO using the given DMT frame type. This is
  required if ligo_data_find will not return the dmt frames.

  Example:

  >>>dmt_cache(960000000,960010000,'H1','LockLoss_H1')
  ['/archive/frames/dmt/LHO/LockLoss_H1/H-M-960/H-LockLoss_H1_M-960001200-3600.gwf',
   '/archive/frames/dmt/LHO/LockLoss_H1/H-M-960/H-LockLoss_H1_M-960004800-3600.gwf']
  """

  #== find dmt frames path
  host = GetCommandOutput('hostname -f')[0].replace('\n','')
  if host.find('ligo-')!=-1:
    dmt_dir = '/dmt'
  elif host.find('ligo.caltech')!=-1:
    site = {'H':'LHO','L':'LLO','V':'V1'}
    dmt_dir = os.path.join('/archive','frames','dmt',site[ifo[0]])

  cache=[]
  #== set base directory
  base_dir = os.path.join(dmt_dir,type)
  #== frames are 3600 seconds long, so round
  tmp = int(str(start)[0:3]+'000000')
  cache_start = tmp+3600*int((start-tmp)/3600)
  tmp = int(str(end)[0:3]+'000000')
  cache_end = tmp+3600*int((end-tmp)/3600)
  #== query with ls
  first_three = []
  start_three = int(str(cache_start)[0:3])
  end_three = int(str(cache_end)[0:3])
  first_three = numpy.arange(start_three,end_three+1,1)
  for t in first_three:
    ls_cmd = "ls "+os.path.join(base_dir,ifo[0:1]+'-M-'+str(t),'*')+" | "+\
             "awk -F - '($5>="+str(cache_start)+" && $5<="+str(cache_end)+")'"
    cache_out = os.popen(ls_cmd)
    for frame in cache_out.readlines():
      frame = frame.replace('\n','')
      cache.append(frame)
    cache_out.close()
  return cache

# ==============================================================================
# Function to expand an frame cache into frame files
# ==============================================================================
def expand_cache(cache):
  """
  Function to expand a frame cache (as given by ligo_data_find -W) into a list 
  of frame files with complete paths.

  Example:

  >>>expand_cache(['H R 963608000 963611296 32 /archive/frames/S6/L0/LHO/H-R-9636'])
  ['/archive/frames/S6/L0/LHO/H-R-9636/H-R-963608000-32.gwf',
   '/archive/frames/S6/L0/LHO/H-R-9636/H-R-963608032-32.gwf',
  ...
   '/archive/frames/S6/L0/LHO/H-R-9636/H-R-963611264-32.gwf']
  """
  
  frames=[]
  for line in cache:
    obs,type,start,end,duration,path = line.split(' ')
    start = int(start)
    end = int(end)
    duration = int(duration)
    cache_start = start
    while cache_start <=(end-duration):
      #== construct frame file name
      file = obs+'-'+type+'-'+str(cache_start)+'-'+str(duration)+'.gwf'
      #== construct full path and append to file
      frame = os.path.join(path,file)
      frames.append(frame)
      #== move to next frame
      cache_start+=duration

  return frames
