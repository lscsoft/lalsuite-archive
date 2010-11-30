#!/usr/bin/env python

from __future__ import division

from glue import git_version

import os
import sys
from pylal import Fr
sys.path.append('/archive/home/duncan.macleod/cvs/detchar/code/pyDQ')

__author__ = "Duncan Macleod <duncan.macleod@ligo.org>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date

"""
Module providing data find and manipulation utils for DQ
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
def data_grab(start_time,end_time,channel_name,frame_type,\
              nds=False,verbose=False):
  time = []
  data = []

  #== generate framecache
  if verbose:
    print >>sys.stdout, "Generating framecache..."
    sys.stdout.flush()
  cache = cache_gen(start_time,end_time,channel_name[0:1],frame_type)

  #== loop over frames in cache
  for frame in cache:
    frame_data,data_start,_,dt,_,_ = Fr.frgetvect1d(frame,channel_name)
    if frame_data==[]:
      print >>sys.stdout, "  No data for "+channel_name+" in "+frame
    #== construct time array
    frame_length = float(dt)*len(frame_data)
    frame_time = data_start+dt*numpy.arange(len(frame_data))
    #== discard frame data outside of time span
    for i in range(len(frame_data)):
      if frame_time[i] < start_time:  continue
      if frame_time[i] > end_time:  continue
      time.append(frame_time[i])
      data.append(frame_data[i])

  return time,data

# ==============================================================================
# Function to generate framecache in memory given lists of ifos and types
# ==============================================================================
def cache_gen(start_time,end_time,ifos,types,omega=False):
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
      data_find_cmd = '''ligo_data_find --observatory '''+ifo[0:1]+\
          ''' --gps-start-time '''+str(start_time)+\
          ''' --gps-end-time '''+str(end_time)+\
          ''' --url-type file --lal-cache --type '''+type
      #== if generating a framecache for use with wpipeline:
      if omega:
        data_find_cmd += '''| ~omega/opt/omega/bin/convertlalcache'''
      #== otherwise
      else:
        data_find_cmd += ''' | cut -d\  -f5 | cut -c17-'''
      #== run ligo_data_find and append each frame to the cache
      frame_cache_out = os.popen(data_find_cmd)
      for frame in frame_cache_out.readlines():
        cache.append(frame.replace('\n',''))
      frame_cache_out.close()

  if cache==[]:
    print >>sys.stderr, "Warning: no frames found."
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
  frame_type_find_cmd = ldf_exe+" -y | egrep "

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
    frame_type_find_cmd+="-v '"
    for vstring in vgrep_list:
      frame_type_find_cmd+=vstring+'|'
    #== take off last '|'
    frame_type_find_cmd = frame_type_find_cmd[0:-1] + "'"

  #== if given types
  else:
    for type in types:
      #== if type is one of the special cases, save for later
      if type in special_types:
        special_cases.append(type)
      #== otherwise add it to the grep command
      else:
        frame_type_find_cmd+= "'"
        frame_type_find_cmd+=type+"|"
    #== take of the extra character
    if frame_type_find_cmd[-1]=="|":
      frame_type_find_cmd = frame_type_find_cmd[0:-1] + "'"

  found_types = []
  #== if not searching only for special types, run the grep command
  if frame_type_find_cmd != ldf_exe+" -y | egrep '":
    found_types = GetCommandOutput(frame_type_find_cmd)[0]
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
def find_ifos(channel,types,ifos):
  """
  Constructs an acceptable list of IFOs given a list of channels, types and ifos
  """

  accepted_ifos=['H1','H0','H2','L0','L1','G0','G1','V0','V1']

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
      print >>sys.stderr, 'Unable to search for IFO '+ifo+'. Please try again'
      ifos.remove(ifo)

  return ifos

# =============================================================================
# Function to find channels
# =============================================================================
def find_channels(channels,types,ifos,ignore=[],match=False,time=None,\
                  verbose=False):
  #== check for ldf
  ldf_exe='ligo_data_find'
  ldf_status = GetCommandOutput('which '+ldf_exe)[1]
  if ldf_status != 0:
    print "Error: ligo_data_find not found. Please ensure lscsoftrc is sourced"
    sys.exit()

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
                  ''' --gps-start-time '''+time+\
                  ''' --gps-end-time '''+time+\
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
            found_channel = DQutils.Channel(name,type=type,sampling=sampling)
            found_channels.append(found_channel)
            count+=1
            sys.stdout.flush()

      #== print channel count for data type
      if verbose:  print >>sys.stdout, count,"channels found"
    if verbose:  print >>sys.stdout

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
      if sampling==int(sampling):
        self.sampling = int(sampling)
      else:
        self.sampling = float(sampling)

# =============================================================================
# Function to calculate blrms
# =============================================================================
def blrms(data,sampling,average=None,band=None,offset=0,w_data=None,\
          remove_mean=False):
  #== redefine None variables
  if average==None:
    average=len(data)/sampling
  if band==None:
    band=[0,sampling/2]
  #== calculate mean
  if remove_mean:
    mean = sum(data)/len(data)
    data = data-mean
  #== generate window
  window = pylab.hanning(len(data))
  data = numpy.multiply(data,window)
  #== Fourier transform
  fft_data = numpy.fft.rfft(data)
  #== PSD (homemade)
  psd_tmp = (8/3)/(pow(sampling,2)*average)*\
                numpy.multiply(fft_data,numpy.conj(fft_data))
  df = sampling/len(data)
  frequencies = list(numpy.arange(0,sampling/2,df))
  psd = {}
  #== set up psd as dictionary for ease
  for freq in frequencies:
    psd[freq] = psd_tmp[frequencies.index(freq)]

  #== define frequency band vector by removing psd frequencies outside of band
  for freq in frequencies:
    if freq < band[0]:
      del psd[freq]
    elif freq >= band[1]:
      del psd[freq]
  band_freq = sorted(psd.keys())
  print min(band_freq),max(band_freq)

  #== calculate banded weight function
  banded_weight = {}
  if w_data is not None:
    #== construct weight dictionary for ease
    w_frequencies = list(w_data[:,0])
    weight={}
    for freq in w_frequencies:
       weight[freq]=w_data[:,1][w_frequencies.index(freq)]
    #== calculate weight for each frequency in given band
    for freq in band_freq:
      w_index=-1
      #== if frequency is in the weighting function, use it
      if freq in w_frequencies:
        banded_weight[freq] = weight[freq]
      #== else, linearly extrapolate weight from weighting function 
      else:
        #== find weight frequency on either side using frequency list
        for w_freq in w_frequencies:
          #== find position of surrounding pair
          if w_freq>freq:
            w_index = w_frequencies.index(w_freq)-1
            if w_index==-1:  w_index-=1
            break
        #== if index not found, assign weight of one
        if w_index == -1:
          banded_weight[freq]=1
       #== unless not found because freq is below lowest weight freq, 
       #==   assign weight of lowest weight freq for consistency
        elif w_index ==-2:
          banded_weight[freq]=weight[w_frequencies[0]]
        else:
          wf_low,wf_high = w_frequencies[w_index],w_frequencies[w_index+1]
          #== calculate frequency weight linearly between weight on either side
          w_interval = weight[wf_high]-weight[wf_low]
          banded_weight[freq] = weight[wf_low] + \
              w_interval * (freq-wf_low)/(wf_high-wf_low)

  else:
    #== construct unity weight function
    for freq in band_freq:
      banded_weight[freq]=1

  #== restrict psd to band
  banded_psd=[]
  for freq in band_freq:
    banded_psd.append(psd[freq])

  #== calculate blrms
  blrms = sqrt(\
              (sum(\
                   numpy.multiply(banded_weight.values(),psd.values()))\
               + offset)\
              *df)

  return blrms

# =============================================================================
# Function to bandpass a time-series
# =============================================================================
def filter(data, f_low, f_high, sampling, padding=0, order=4):
  passband = [f_low*2/sampling,f_high*2/sampling]
  b,a = signal.butter(order,passband,btype='bandpass')
  data = signal.lfilter(b,a,data)
  #data = data[int(padding*sampling):] # removes padding from front end of data

  data = data[::-1] # reverses data

  b,a  = signal.butter(order,passband,btype='bandpass')
  data = signal.lfilter(b,a,data)
  #data = data[int(padding*sampling):] # removes padding from back end of data

  data = data[::-1]
  return data

# =============================================================================
# Function to generate an daily ihope cache 
# =============================================================================
def daily_ihope_cache(start,end,ifo,cluster=None):
  """
  Generates cache list of daily ihope INSPIRAL xml files for give ifo and clustering (None,'30ms','100ms', or '16s') between start and end time
  """

  #== daily path
  ihope_daily_path = '/archive/home/cbc/ihope_daily'

  #== set clustering tag
  if cluster==None:
    cluster_tag='UNCLUSTERED'
  elif cluster=='100ms':
    cluster_tag='100MILLISEC_CLUSTERED'
  elif cluster=='30ms':
    cluster_tag='30MILLISEC_CLUSTERED'
  elif cluster=='16s':
    cluster_tag='16SEC_CLUSTERED'

  #== work out days
  day_start = int(GetCommandOutput('tconvert `tconvert '+start+' -f %D`')[0])
  duration = end-start
  num_days = int(round((duration)/86400))
  #== generate array of days
  day_end = day_start+num_days*86400
  while day_end<end:
    day_end+=86400
  days = numpy.arange(day_start,day_end,86400)

  cache=[]
  #== loop over days gathering files
  for day in days:
    date = GetcommandOutput('tconvert '+day+' -f %Y%m%d')[0]
    day_path = os.path.join(ihope_daily_path,date[0:4],date)
    ls_cmd = 'ls '+day_path+'/'+ifo+'-INSPIRAL_UNCLUSTERED_'+cluster_tag+\
             '*.xml.gz'
    cache_out = os.popen(ls_cmd)
    for line in cache_out:
      trig_start = int(line.split('.xml')[0].split('-')[-2])
      duration = int(line.split('.xml')[0].split('-')[-2])
      if start<=trig_start<end or start<(trig_start+duration)<=end:
        cache.append(line.replace('\n',''))

  cache_out.close()

  return cache

# ==============================================================================
# Function to grab effective distance from ihope cache of trig files
# ==============================================================================
def grab_effective_distance(cache,time=False):
  distance=[]
  time=[]
  #== grab data
  for file in cache:
    data = GetCommandOutput('''ligolw_print -t summ_value '''+\
                            '''-c start_time -c end_time -c name -c value '''+\
                               file+''' | grep inspiral_effective_distance'''+\
                            ''' | 'NR==1 {print }' | cut -f1,2,4 -d, ''')
    start,end,dist = [float(num) for num in data.split(',')]
    time.append((start+end)/2)
    distance.append(dist)

  if time:  return distance,time
  else:  return distance
