#!/usr/bin/env python

from __future__ import division

from glue import git_version

import os
import sys
from pylal import Fr
import numpy

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
  time = []
  data = []

  #== generate framecache
  if verbose:
    print >>sys.stdout, "Generating framecache..."
    sys.stdout.flush()
  if not dmt:
    cache = cache_gen(start,end,channel[0:1],type)
  else:
    cache = dmt_cache(start,end,channel[0:1],type)
  #== loop over frames in cache
  for frame in cache:
    try:
      frame_data,data_start,_,dt,_,_ = Fr.frgetvect1d(frame,channel)
      if frame_data==[]:
        print >>sys.stdout, "  No data for "+channel+" in "+frame
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
  [H R 962000000 962064032 32 /archive/frames/S6/L0/LHO/H-R-9620]
  [H R 961977600 962000000 32 /archive/frames/S6/L0/LHO/H-R-9619]
  >>>generate_cache(961977615,962582415,R,H,return_files=True)
  [/archive/frames/S6/L0/LHO/H-R-9619/H-R-961977600-32.gwf,
   /archive/frames/S6/L0/LHO/H-R-9619/H-R-961977632-32.gwf,
   ...
   /archive/frames/S6/L0/LHO/H-R-9620/H-R-962064000-32.gwf]
  """

  cache = []
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
        data_find_cmd = '''ligo_data_find '''+\
                        ''' --gps-start-time '''+str(start_time)+\
                        ''' --gps-end-time '''+str(end_time)+\
                        ''' --observatory '''+ifo[0:1]+\
                        ''' --type '''+type+\
                        ''' --url-type file '''+\
                        ''' --frame-cache '''
      #== run ligo_data_find and append each frame to the cache
        cache_out = Popen(data_find_cmd,stdout=PIPE)
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
def find_types(types,short=False,full=False):
  #== check for ldf
  ldf_exe='ligo_data_find'
  ldf_status = GetCommandOutput('which '+ldf_exe)[1]
  if ldf_status != 0:
    print >>sys.stderr, \
        "Error: ligo_data_find not found. Please ensure lscsoftrc is sourced"
    sys.exit()

  #== set up search command
  find_cmd = ldf_exe+" -y | egrep "

  #== treat 'R','M' and 'T' as special cases,
  #== so not to grep for all types containing 'R'
  special_types = ['M','R','T']
  special_cases=[]

  #== set list of ignored strings in `ligo_data_find -y`
  if types is None:
    #== there are thousands of GRBXXXXXX frame types, so ignore them ALWAYS
    vgrep_list = ['GRB']
    if short:
      #== all of these strings are part of frame types that can be ignored for a
      #== short search
      short_ignore_list = ['CAL','BRST','Mon','SG','IMR','DuoTone','Concat',\
                           'BH','WNB','Lock','_M','_S5','Multi','Noise']
      vgrep_list.extend(short_ignore_list)
    #== add each of those ignored strings to a vgrep command
    find_cmd+="-v '"
    for vstring in vgrep_list:
      find_cmd+=vstring+'|'
    #== take off last '|'
    find_cmd = find_cmd[0:-1] + "'"

  #== if given types
  else:
    for type in types:
      #== if type is one of the special cases, save for later
      if type in special_types:
        special_cases.append(type)
      #== otherwise add it to the grep command
      else:
        find_cmd+= "'"
        find_cmd+=type+"|"
    #== take of the extra character
    if find_cmd[-1]=="|":
      find_cmd = find_cmd[0:-1] + "'"

  found_types = []
  #== if not searching only for special types, run the grep command
  if find_cmd != ldf_exe+" -y | egrep '":
    found_types = GetCommandOutput(find_cmd)[0]
    found_types = found_types.split('\n')[0:-1]
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

  #== if given no ifos, try to generate a list from the channels given  
  if ifos is None:
    ifos=[]
    if channels is not None:
      for channel in channels:
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
                  ex_channels=None,\
                  types=None,\
                  ifos=None,\
                  ignore=[],\
                  match=False,\
                  time=None,\
                  unique=False,\
                  verbose=False):

  #== check for ldf
  ldf_exe='ligo_data_find'
  ldf_status = GetCommandOutput('which '+ldf_exe)[1]
  if ldf_status != 0:
    print "Error: ligo_data_find not found. Please ensure lscsoftrc is sourced"
    sys.exit()

  #== cannot work with no ifos
  if ifos is None:
    ifos = find_ifos(channels,types,ifos)

  #== check list status
  for input_set in [channels,types,ifos]:
    if isinstance(input_set,str):
      input_set = [input_set]
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
      frame_cmd = ldf_exe+''' --observatory '''+ifo[0:1]+\
                  ''' --type='''+type+\
                  ''' --gps-start-time '''+str(time)+\
                  ''' --gps-end-time '''+str(time)+\
                  ''' --url-type file --lal-cache | '''+\
                  ''' sort -g -k 3 -r | awk 'NR==1' '''
      frame,frame_status = GetCommandOutput(frame_cmd)
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
        channel_list,channel_status = \
                GetCommandOutput(channel_find_cmd)
        #== if grab was successful:
        if channel_status == 0:
          channel_list=channel_list.split('\n')
          for data in channel_list:
            if data=='':  continue
            data = data.replace('\n','')
            name,sampling = data.split(' ')
            #== if asked for exact match, check:
            if match:
              if name not in channels:  continue
            #== generate structure and append to list  
            found_channel = Channel(name,type=type,sampling=sampling)
            found_channels.append(found_channel)
            count+=1
            sys.stdout.flush()

      #== print channel count for data type
      if verbose:  print >>sys.stdout, count,"channels found"
    if verbose:  print >>sys.stdout

  if unique:
    found_channels = parse_unique_channels(found_channels)

  return found_channels

# ==============================================================================
# parse channels for uniqueness
# ==============================================================================
#== takes a list of DQutils.Channel structures and extracts a unique list of 
#== channels based on a given frame type preference order
#== NB: assumes channel list has been generated using Channel class below
def parse_unique_channels(channels,type_order=None):
  #== if given empty, return empty and hope user notices
  if channels==[]:
    return []

  #== set up type preference order
  if type_order == None:
    type_order = []
    #== set up list of ifos for type name construction
    ifos = []
    for channel in channels:
      if channel.name.find(':')==-1:
        ifos.append(channel.name[0:2])
    #== for each ifo add the *_RDS_R_L1 and *_RDS_R_L3 sets to the order
    for ifo in ifos:
      type_order.append(ifo+'_RDS_R_L1')
      type_order.append(ifo+'_RDS_R_L3')
    #== add the Raw data type
    type_order.append('R')
    #== add the S5 RDS type
    type_order.append('RDS_R_L1')

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
  cache=[]
  #== set base directory
  base_dir = os.path.join('/dmt',type)
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
  """
  
  frames=[]
  for line in cache:
    obs,type,start,end,duration,path = line.split(' ')
    cache_start = start
    while cache_start <=(end-duration):
      #== construct frame file name
      file = obs+'-'+type+'-'+cache_start+'-'+duration+'.gwf'
      #== construct full path and append to file
      frame = os.path.join(path,file)
      frames.append(frame)
      #== move to next frame
      cache_start+=duration

  return frames
