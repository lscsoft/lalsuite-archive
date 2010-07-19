#!/usr/bin/env python

from __future__ import division
import types
import sys
import os

try:
  from glue.segments import segment, segmentlist
except:
  pass
from math import sqrt,log,fabs,exp,cos,pi

# Hey, scipy, shut up about your nose already.
import warnings
warnings.filterwarnings("ignore")
from scipy import signal as signal
from scipy.fftpack import fft, ifft, ifftshift, fft2, ifft2
from matplotlib import mlab, use
use('Agg')
import pylab
import numpy
from pylal import Fr

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

  #== grab data use nds2 server
  if nds==True:
    try:
      results = daq.get_data([channel_name],start_time,end_time)
    except:
      print >>sys.stdout, "Failed to get data for "+channel_name+" from nds2."
      sys.stdout.flush()
    #== add data      
    data = numpy.array(results[0].data)
    sampling = results[0].signal_rate
    dt = 1/sampling
    #== construct time array
    time = start_time + dt*numpy.arange(len(data))

  else:
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
  if isinstance(ifos,str):
    ifos=[ifos]
  if isinstance(types,str):
    types=[types]
  for ifo in ifos:
    for type in types:
      data_find_cmd = '''ligo_data_find --observatory '''+ifo[0:1]+\
          ''' --gps-start-time '''+str(start_time)+\
          ''' --gps-end-time '''+str(end_time)+\
          ''' --url-type file --lal-cache --type '''+type
      if omega:
        data_find_cmd += '''| ~omega/opt/omega/bin/convertlalcache'''
      else:
        data_find_cmd += ''' | cut -d\  -f5 | cut -c17-'''
      frame_cache_out = os.popen(data_find_cmd)
      for frame in frame_cache_out.readlines():
        cache.append(frame.replace('\n',''))
      frame_cache_out.close()
  if not cache:
    print >>sys.stdout, "No frames found."
  return cache

# =============================================================================
# Function to generate channel structure given names
# =============================================================================
def channel_grab(channels=None,excluded=None,types=None,\
                 ifos=None,unique=False,time=None,match=False):
  #== set up exe
  channel_data=[]
  channel_query_exe = 'ligo_channel_query'
  channel_query_args = ''

  #== set up --gps-time argument
  if time is not None:
    channel_query_args = ' --gps-time '+str(time)

  #== set up --channels argument
  if isinstance(channels,str):
    channels=[channels]
  if channels is not None:
    channel_query_args += ' --channels '
    for channel in channels:
      channel_query_args += channel+','
    channel_query_args = channel_query_args[0:-1]+' '

  #== set up --exclude-channels argument
  if isinstance(excluded,str):
    excluded=[excluded]
  if excluded is not None:
    channel_query_args += ' --exclude-channels '
    for channel in excluded:
      channel_query_args += channel+','
    channel_query_args = channel_query_args[0:-1]+' '

  #== set up --types argument
  if types is not None:
    if isinstance(types,str):
      types=[types]
    channel_query_args += '--types '
    for type in types:
      channel_query_args += type+','
      channel_query_args = channel_query_args[0:-1]+' '

  #== set up --ifos argument: generate ifos list if not given
  if ifos is None:
    accepted_ifos = ['H0','H1','H2','L0','L1','G1','V1']
    ifos = []
    for channel in channels:
      ifo = channel[0:2]
      if ifo not in ifos:
        if ifo not in accepted_ifos:  continue
        ifos.append(ifo[0:1])
  if isinstance(ifos,str):
    ifos=[ifos]
  channel_query_args += '--ifos '
  for ifo in ifos:
    channel_query_args += ifo+','
  channel_query_args = channel_query_args[0:-1]

  #== add match option
  if match:
    channel_query_args += ' --match'

  #== execute query
  channels=[]
  channel_query = os.popen(channel_query_exe+channel_query_args)
  try:
    for channel in channel_query.readlines():
      if channel=='':  continue
      channels.append(channel.replace('\n',''))
  except: 
    print >>sys.stdout, "Error in ligo_channel_query:\n\n"+\
        channel_query_exe+channel_query_args+"\n\nExiting."
    sys.exit()
  channel_query.close()

  #== if no channels, throw error, else, generate structures
  if channels==[]:
    print >>sys.stdout, "No channels found with cmd:\n\n"+\
        channel_query_exe+channel_query_args+"\n"
  else:
    for i in range(len(channels)):
      [name,type,sampling] = channels[i].split(' ')
      channels[i] = Channel(name,type=type,sampling=sampling)

  if unique == True:
    channels = parse_unique_channels(channels)
  
  return channels

# ==============================================================================
# parse channels for uniqueness
# ==============================================================================
#== takes an output from ligo_channel_query and extracts a unique list of 
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
# Class to generate channel structure from output of ligo_channel_query
# ==============================================================================
class Channel:
  def __init__(self,name,type=None,sampling=None):
    self.name = str(name)
    if type is not None:
      self.type = str(type)
    if sampling is not None:
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
  #band_freq = numpy.arange(band[0],band[1],1/average)

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

  #psd = psd[int(round(band[0]*average)):int(round(band[1]*average))]
  #== calculate blrms
  #blrms = numpy.multiply(banded_weight.values(),psd)
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
