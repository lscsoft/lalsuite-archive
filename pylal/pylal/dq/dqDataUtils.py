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

# =============================================================================
# Function to calculate blrms
# =============================================================================
def blrms(data,sampling,average=None,band=None,offset=0,w_data=None,\
          remove_mean=False):
  """
  This function will calculate the band-limited root-mean-square of the given
  data, using averages of the given length in the given [f_low,f_high) band.

  Options are included to offset the data, and weight frequencies given a 
  dict object of (frequency:weight) pairs.
  """
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
def filter(data, f_low, f_high, sampling, order=4):
  """
  This function will bandpass filter data in the given [f_low,f_high) band
  using the given order Butterworth filter.
  """

  #== construct passband
  passband = [f_low*2/sampling,f_high*2/sampling]
  #== construct filter
  b,a = signal.butter(order,passband,btype='bandpass')
  #== filter data forward then backward
  data = signal.lfilter(b,a,data)
  data = data[::-1]
  data = signal.lfilter(b,a,data)
  data = data[::-1]

  return data
