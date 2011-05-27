#!/usr/bin/env python

# =============================================================================
# Preamble
# =============================================================================

from __future__ import division
import re,numpy,math,subprocess

from glue.ligolw import ligolw,table,lsctables,utils
from glue.ligolw.utils import process as ligolw_process
from glue import segments

from pylal import llwapp
from pylal.xlal.datatypes.ligotimegps import LIGOTimeGPS
from pylal.dq import dqTriggerUtils

# Hey, scipy, shut up about your nose already.
import warnings
warnings.filterwarnings("ignore")
from scipy import signal as signal

from matplotlib import mlab
from glue import git_version

__author__  = "Duncan Macleod <duncan.macleod@astro.cf.ac.uk>"
__version__ = "git id %s" % git_version.id
__date__    = git_version.date

"""
This module provides a bank of useful functions for manipulating triggers and trigger files for data quality investigations.
"""

# =============================================================================
# Execute shell command and get output
# =============================================================================

def make_external_call(command):

  """
    Execute shell command and capture standard output and errors. 
    Returns tuple "(stdout,stderr)".
  """

  p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    shell=isinstance(command, str))
  out, err = p.communicate()

  return out,err

# =============================================================================
# Read injection files
# =============================================================================

def frominjectionfile( file, type, ifo=None, start=None, end=None):
  
  """
    Read generic injection file object file containing injections of the given
    type string. Returns an 'Sim' lsctable of the corresponding type.

    Arguments:
   
      file : file object
      type : [ "inspiral" | "burst" | "ringdown" ]

    Keyword arguments:

      ifo : [ "G1" | "H1" | "H2" | "L1" | "V1" ]
  """

  # read type
  type = type.lower()

  # read injection xml
  xml = re.compile('(xml$|xml.gz$)')
  if re.search(xml,file.name):
    xmldoc,digest = utils.load_fileobj(file)
    injtable = table.get_table(xmldoc,'sim_%s:table' % (type))

  # read injection txt
  else:
    cchar = re.compile('[#%<!()_\[\]{}:;\'\"]+')

    #== construct new Sim{Burst,Inspiral,Ringdown}Table
    injtable = lsctables.New(lsctables.__dict__['Sim%sTable' % (type.title())])
    if type=='inspiral':
      columns = ['geocent_end_time.geocent_end_time_ns',\
                 'h_end_time.h_end_time_ns',\
                 'l_end_time.l_end_time_ns',\
                 'v_end_time.v_end_time_ns',\
                 'distance'] 
      for line in file.readlines():
        if re.match(cchar,line):
          continue
        # set up siminspiral object
        inj = lsctables.SimInspiral()
        # split data
        sep = re.compile('[\s,]+')
        data = sep.split(line)
        # set attributes
        inj.geocent_end_time    = int(data[0].split('.')[0])
        inj.geocent_end_time_ns = int(data[0].split('.')[1])
        inj.h_end_time          = int(data[1].split('.')[0])
        inj.h_end_time_ns       = int(data[1].split('.')[1])
        inj.l_end_time          = int(data[2].split('.')[0])
        inj.l_end_time_ns       = int(data[2].split('.')[1])
        inj.v_end_time          = int(data[3].split('.')[0])
        inj.v_end_time_ns       = int(data[3].split('.')[1])
        inj.distance            = float(data[4])

        injtable.append(inj)

    if type=='burst':
      if file.readlines()[0].startswith('filestart'):
        # if given parsed burst file
        file.seek(0)

        snrcol = { 'G1':23, 'H1':19, 'L1':21, 'V1':25 }

        for line in file.readlines():
          inj = lsctables.SimBurst()
          # split data
          sep = re.compile('[\s,]+')
          data = sep.split(line)
          # set attributes

          # gps time
          if 'burstgps' in data:
            idx = data.index( 'burstgps' )+1
            geocent = LIGOTimeGPS(data[idx])

            inj.time_geocent_gps    = geocent.seconds
            inj.time_geocent_gps_ns = geocent.nanoseconds
          else:
            continue


          #inj.waveform            = data[4]
          #inj.waveform_number     = int(data[5])

          # frequency
          if 'freq' in data:
            idx = data.index( 'freq' )+1
            inj.frequency = float( data[idx] )
          else:
            continue

          # SNR a.k.a. amplitude
          if ifo and 'snr%s' % ifo in data:
            idx = data.index( 'snr%s' % ifo )+1
            inj.amplitude = float( data[idx] )
          elif 'rmsSNR' in data:
            idx = data.index( 'rmsSNR' )+1
            inj.amplitude = float( data[idx] )
          else:
            continue

          if 'phi' in data:
            idx = data.index( 'phi'  )+1
            inj.ra = float(data[idx])*24/(2*math.pi)       

          if 'theta' in data:
            idx = data.index( 'theta'  )+1 
            inj.ra = 90-(float(data[idx])*180/math.pi)

          if ifo and 'hrss%s' % ifo in data:
            idx = data.index( 'hrss%s' % ifo )+1
            inj.hrss = float( data[idx] )
          elif 'hrss' in data:
            idx = data.index( 'hrss' )+1
            inj.hrss = float( data[idx] )

          # extra columns to be added when I know how
          #inj.q = 0
          #inj.q                   = float(data[11])
          #h_delay = LIGOTimeGPS(data[41])
          #inj.h_peak_time         = inj.time_geocent_gps+h_delay.seconds
          #inj.h_peak_time_ns      = inj.time_geocent_gps_ns+h_delay.nanoseconds
          #l_delay = LIGOTimeGPS(data[43])
          #inj.l_peak_time         = inj.time_geocent_gps+l_delay.seconds
          #inj.l_peak_time_ns      = inj.time_geocent_gps_ns+l_delay.nanoseconds
          #v_delay = LIGOTimeGPS(data[43])
          #inj.v_peak_time         = inj.time_geocent_gps+v_delay.seconds
          #inj.v_peak_time_ns      = inj.time_geocent_gps_ns+v_delay.nanoseconds

          injtable.append(inj)

      else:
        # if given parsed burst file
        file.seek(0)
        for line in file.readlines():
          inj = lsctables.SimBurst()
          # split data
          sep = re.compile('[\s,]+')
          data = sep.split(line)
          # set attributes
          geocent = LIGOTimeGPS(data[0])
          inj.time_geocent_gps    = geocent.seconds
          inj.time_geocent_gps_ns = geocent.nanoseconds

          injtable.append(inj)

  injections = table.new_from_template( injtable )
  if not start:  start = 0
  if not end:    end   = 9999999999
  span = segments.segmentlist([ segments.segment(start, end) ])
  get_time = dqTriggerUtils.def_get_time( injections.tableName )
  injections.extend( inj for inj in injtable if get_time(inj) in span )

  return injections

# =============================================================================
# Calculate band-limited root-mean-square
# =============================================================================

def blrms(data,sampling,average=None,band=None,offset=0,w_data=None,\
          remove_mean=False):

  """
    This function will calculate the band-limited root-mean-square of the given
    data, using averages of the given length in the given [f_low,f_high) band.

    Options are included to offset the data, and weight frequencies given a 
    dict object of (frequency:weight) pairs.
  """

  # redefine None variables
  if average==None:
    average=len(data)/sampling
  if band==None:
    band=[0,sampling/2]
  # calculate mean
  if remove_mean:
    mean = sum(data)/len(data)
    data = data-mean
  # generate window
  window = pylab.hanning(len(data))
  data = numpy.multiply(data,window)
  # Fourier transform
  fft_data = numpy.fft.rfft(data)
  # PSD (homemade)
  psd_tmp = (8/3)/(pow(sampling,2)*average)*\
                numpy.multiply(fft_data,numpy.conj(fft_data))
  df = sampling/len(data)
  frequencies = list(numpy.arange(0,sampling/2,df))
  psd = {}
  # set up psd as dictionary for ease
  for freq in frequencies:
    psd[freq] = psd_tmp[frequencies.index(freq)]
  # define frequency band vector by removing psd frequencies outside of band
  for freq in frequencies:
    if freq < band[0]:
      del psd[freq]
    elif freq >= band[1]:
      del psd[freq]
  band_freq = sorted(psd.keys())
  #band_freq = numpy.arange(band[0],band[1],1/average)

  # calculate banded weight function
  banded_weight = {}
  if w_data is not None:
    # construct weight dictionary for ease
    w_frequencies = list(w_data[:,0])
    weight={}
    for freq in w_frequencies:
       weight[freq]=w_data[:,1][w_frequencies.index(freq)]
    # calculate weight for each frequency in given band
    for freq in band_freq:
      w_index=-1
      # if frequency is in the weighting function, use it
      if freq in w_frequencies:
        banded_weight[freq] = weight[freq]
      # else, linearly extrapolate weight from weighting function 
      else:
        # find weight frequency on either side using frequency list
        for w_freq in w_frequencies:
          # find position of surrounding pair
          if w_freq>freq:
            w_index = w_frequencies.index(w_freq)-1
            if w_index==-1:  w_index-=1
            break
        # if index not found, assign weight of one
        if w_index == -1:
          banded_weight[freq]=1
       # unless not found because freq is below lowest weight freq, 
       #   assign weight of lowest weight freq for consistency
        elif w_index ==-2:
          banded_weight[freq]=weight[w_frequencies[0]]
        else:
          wf_low,wf_high = w_frequencies[w_index],w_frequencies[w_index+1]
          # calculate frequency weight linearly between weight on either side
          w_interval = weight[wf_high]-weight[wf_low]
          banded_weight[freq] = weight[wf_low] + \
              w_interval * (freq-wf_low)/(wf_high-wf_low)

  else:
    # construct unity weight function
    for freq in band_freq:
      banded_weight[freq]=1

  # restrict psd to band
  banded_psd=[]
  for freq in band_freq:
    banded_psd.append(psd[freq])

  #psd = psd[int(round(band[0]*average)):int(round(band[1]*average))]
  # calculate blrms
  #blrms = numpy.multiply(banded_weight.values(),psd)
  blrms = math.sqrt(\
              (sum(\
                   numpy.multiply(banded_weight.values(),psd.values()))\
               + offset)\
              *df)

  return blrms

# =============================================================================
# Function to bandpass a time-series
# =============================================================================

def bandpass( data, f_low, f_high, sampling, order=4 ):

  """
    This function will bandpass filter data in the given [f_low,f_high) band
    using the given order Butterworth filter.
  """

  # construct passband
  passband = [f_low*2/sampling,f_high*2/sampling]
  # construct filter
  b,a = signal.butter(order,passband,btype='bandpass')
  # filter data forward then backward
  data = signal.lfilter(b,a,data)
  data = data[::-1]
  data = signal.lfilter(b,a,data)
  data = data[::-1]

  return data

# =============================================================================
# Calculate spectrum
# =============================================================================

def spectrum( data, sampling, NFFT=256, overlap=0.5,\
              window='hanning', detrender=mlab.detrend_linear,\
              sides='onesided', scale='PSD' ):

  numpoints  = len(data)
  numoverlap = int( sampling * (1.0 - overlap ))

  if isinstance(window,str):
    window=window.lower()

  win = signal.get_window(window, NFFT)

  # calculate PSD with given parameters
  spec,freq = mlab.psd( data, NFFT=NFFT, Fs=sampling, noverlap=numoverlap,\
                        window=win, sides=sides, detrend=detrender )

  # rescale data to meet user's request
  scale = scale.lower()
  if scale == 'asd':
    spec = numpy.sqrt(spec) * numpy.sqrt(2 / (sampling*sum(win**2)))
  elif scale == 'psd':
    spec = spec * 2 / (sampling*sum(win**2))
  elif scale == 'as':
    spec = nump.sqrt(spec) * numpy.sqrt(2) / sum(win)
  elif scale == 'ps':
    spec = spec * 2 / (sum(win)**2)

  return freq,spec
