#!/usr/bin/env python

# =============================================================================
# Preamble
# =============================================================================

from __future__ import division
import sys,os,re,numpy,math,shlex,subprocess,datetime,glob,tempfile
from socket import getfqdn
from glue.ligolw import ligolw,table,lsctables,utils
from glue.segments import segment, segmentlist
from glue.segmentdb import segmentdb_utils
from pylal import date
from pylal.xlal.datatypes.ligotimegps import LIGOTimeGPS

# Hey, scipy, shut up about your nose already.
import warnings
warnings.filterwarnings("ignore")
from scipy import signal as signal
from scipy.fftpack import fft, ifft, ifftshift, fft2, ifft2

from matplotlib import use
use('Agg')
import pylab

from glue import git_version

__author__  = "Duncan Macleod <duncan.macleod@astro.cf.ac.uk>"
__version__ = "git id %s" % git_version.id
__date__    = git_version.date

"""
This module provides a bank of useful functions for manipulating triggers and trigger files for data quality investigations.
"""

# =============================================================================
# Execute shell cmd and get output
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

# =============================================================================
# Convert list from text file into Sngl{Burst,Inspiral} object
# =============================================================================
def trigger(data,etg):

  """
    Reads the list object data and returns a Sngl{Burst,Inspiral} object
    with the correct columns seeded, based on the given string object etg.

    Arguments:

      data : [ string | list]
        string or list object containing the data to be parsed.

      etg : [ "ihope" | "kw" | "omega" | "omegaclustered" ]
        Defines how the data list is understood and parsed, based on the
        standard column layouts (comma-, space-, or tab-delimited) for each ETG.
        "ihope" is read as per the .csv files from ihope_daily. "omegaclustered"
        is read as per the .txt or .clst produced by the detchar scripts used to
        plot omega online events.

  """ 

  # read etg
  etg=etg.lower()

  etgcategories = {lsctables.SnglInspiral(): ['ihope'],\
                   lsctables.SnglBurst():    ['omega','omegaclustered','kw'],\
                   lsctables.SnglRingdown(): []}

  # set up trig object
  for obj,etgs in etgcategories.items():
    if etg in etgs:
      trig = obj

  # if given string, split on space, tab or comma
  if isinstance(data,str):
    sep = re.compile('[\t\s,]+')
    data = sep.split(data)

  # =====
  # ihope
  # =====
  if re.match('ihope',etg.lower()):
    # comma separated values are:
    # end_time,end_time_ns,ifo,snr,mass1,mass2,mtotal,eta,event_duration,
    # template_duration,eff_distance,chisq,chisq_dof,bank_chisq,bank_chisq_dof
    # cont_chisq,cont_chisq_dof

    trig.end_time              = int(data[0])
    trig.end_time_ns           = int(data[1])
    trig.ifo                   = str(data[2])
    trig.snr                   = float(data[3])
    trig.mass1                 = float(data[4])
    trig.mass2                 = float(data[5])
    trig.mtotal                = float(data[6])
    trig.eta                   = float(data[7])
    trig.event_duration        = float(data[8])
    trig.template_duration     = float(data[9])
    trig.eff_distance          = float(data[10])
    trig.chisq                 = float(data[11])
    trig.chisq_dof             = float(data[12])
    trig.bank_chisq            = float(data[13])
    try:
      trig.bank_chisq_dof      = float(data[14])
    except:
      trig.bank_chisq_dof = None
    try:
      trig.cont_chisq          = float(data[15])
      trig.cont_chisq_dof      = float(data[16])
    except ValueError:
      trig.cont_chisq = None
      trig.cont_chisq_dof = None

  # ============
  # kleine-welle
  # ============
  elif etg=='kw':
    # peak time
    peak  = LIGOTimeGPS(data[0])
    trig.peak_time             = peak.seconds
    trig.peak_time_ns          = peak.nanoseconds
    # start time
    start = LIGOTimeGPS(data[1])
    trig.start_time            = start.seconds
    trig.start_time_ns         = start.nanoseconds
    # end time
    end = LIGOTimeGPS(data[2])
    trig.end_time              = end.seconds
    trig.end_time_ns           = end.nanoseconds
    # others
    trig.central_freq          = float(data[3])
    trig.energy                = float(data[4])
    trig.norm_energy           = float(data[5])
    trig.n_pix                 = float(data[6])
    trig.significance          = float(data[7])
    trig.N                     = float(data[8])
    trig.snr                   = math.sqrt(trig.norm_energy-trig.n_pix)
 
  # =====
  # omega
  # =====
  if etg=='omega' or etg=='wpipe':
    # peak time
    peak = LIGOTimeGPS(data[0])
    trig.peak_time           = peak.seconds
    trig.peak_time_ns        = peak.nanoseconds
    # central frequency
    trig.central_freq        = float(data[1])
    trig.peak_frequency      = trig.central_freq
    # duration
    trig.duration            = LIGOTimeGPS(float(data[2]))
    # start time
    start = peak-trig.duration
    trig.start_time          = start.seconds
    trig.start_time_ns       = start.nanoseconds
    # end time
    stop = peak+trig.duration
    trig.stop_time           = stop.seconds
    trig.stop_time_ns        = stop.nanoseconds
    # bandwidth and flow,fhigh
    trig.bandwidth           = float(data[3])
    trig.flow                = trig.peak_frequency - 0.5*trig.bandwidth
    trig.fhigh               = trig.peak_frequency + 0.5*trig.bandwidth
    # energy
    trig.amplitude           = float(data[4])

    # cluster parameters
    #trig.cluster_size        = float(data[5])
    #trig.cluster_norm_energy = float(data[6])
    #trig.cluster_number      = float(data[7])

    # SNR
    trig.snr                 = math.sqrt(2*trig.amplitude)

  # ==============
  # omegaclustered
  # ==============
  # follow Cadonati's clustering output
  if etg=='omegaclustered' or etg=='wpipeclustered':
    # start time
    cstart = LIGOTimeGPS(data[0])
    trig.cluster_start_time              = cstart.seconds
    trig.cluster_start_time_ns           = cstart.nanoseconds
    # end time
    cstop  = LIGOTimeGPS(data[1])
    trig.cluster_end_time                = cstop.seconds
    trig.cluster_end_time_ns             = stop.nanoseconds
    # peak time
    peak  = LIGOTimeGPS(data[2])
    trig.peak_time               = peak.seconds
    trig.peak_time_ns            = peak.nanoseconds
    # bandwidth, and flow,fhigh,central_freq
    trig.flow                    = float(data[8])
    trig.fhigh                   = float(data[9])
    trig.bandwidth               = trig.fhigh - trig.flow
    trig.central_freq       = trig.flow + 0.5*trig.bandwidth
    # cluster params
    trig.cluster_flow            = float(data[3])
    trig.cluster_end_frequency   = float(data[4])
    trig.cluster_length          = float(data[5])
    trig.cluster_size            = float(data[10])
    trig.cluster_norm_energy     = float(data[11])
    trig.norm_energy             = float(data[12])
    trig.snr                     = math.sqrt(2*trig.norm_energy)

  return trig

# =============================================================================
# Function to write triggers to file in etg standard form
# =============================================================================
def totrigxml(file,triggers,etg):

  """
    Write the given table object triggers to the file object file, in the standard
    format for the given etg. INCOMPLETE (any suggestions are most welcome)
  """

  return 1

def totrigfile(file,table,etg,header=True):

  """
    Write the lsctables.table to the given file object file, in the standard
    format for the given etg
  """

  etg = etg.lower()

  # set columns
  if re.match('ihope',etg):
    # comma separated values are:
    columns = ['end_time','end_time_ns','ifo','snr','mass1','mass2','mtotal',\
               'eta','event_duration','template_duration','eff_distance',\
               'chisq','chisq_dof','bank_chisq','bank_chisq_dof','cont_chisq',\
               'cont_chisq_dof']

  elif re.match('omegaclustered',etg) or re.match('omega_clustered',etg):
    columns = ['cluster_start_time','cluster_end_time','peak_time',\
               'cluster_flow','cluster_fhigh',\
               'cluster_length','start_time','end_time','flow','fhigh'\
               'cluster_size','cluster_norm_energy','snr']

  elif re.match('omega',etg.lower()) or re.match('wpipe',etg.lower()):
    columns = ['peak_time','central_freq','duration',\
               'bandwidth','norm_energy',\
               'cluster_size','cluster_norm_energy','cluster_number']

  elif re.match('kw',etg.lower()):
    columns = ['peak_time','start_time','end_time','central_freq',\
               'energy','norm_energy','num_pixels','significance','N']

  # set delimiter
  if etg=='ihope':  d = ','
  else:  d=' '

  # print header
  if header:
    header = ' '.join(['#'].extend(columns))
    print >>file, header

  # print triggers
  for row in table:
    line = []
    for col in columns:
       if col not in table.columnnames:  continue
       entry = ''
       # if ihope, print column
       if re.match('ihope',etg.lower()):
         entry = str(row.__getattribute__(col))
       # if not ihope, check for time and print full GPS
       else:
         if col=='peak_time':
           entry = str(row.get_peak())
         elif col=='start_time':
           entry = str(row.get_start())
         elif col=='end_time':
           entry = str(row.get_end())
         else:
           entry = str(row.__getattribute__(col))

       line.append(entry)

    print >>file, d.join(line)
  
# =============================================================================
# Function to load triggers from cache
# =============================================================================
def fromtrigxml(file,tablename='sngl_inspiral:table',start=None,end=None):

  """
    Reads a trigger table from the given table from the xml
    file object file

    Arguments:

      file : file object
   
    Keyword arguments:

      start : [ float | int | LIGOTimeGPS ]
        GPS start time of requested period
      end : [ float | int | LIGOTimeGPS ]
        GPS end time of requested period
      tablename : string
        name of requested trigger table in xml file, defaults to
        'sngl_inspiral:table'
  """

  # set times
  if not start:
    start=0
  if not end:
    end=float('inf')

  # set tablename
  if not tablename.endswith(':table'):
    tablename = ':'.join([tablename,'table'])

  # crack open xml file
  xmldoc,digest = utils.load_fileobj(file,gz=file.name.endswith('gz'))
  alltriggers = table.get_table(xmldoc,tablename))

  triggers = lsctables.New(type(alltriggers))

  # parse triggers in time
  for line in alltriggers:
    if re.search('inspiral',alltriggers.tableName):  t = line.get_end()
    else:  t = line.get_peak()
    if start<=float(t)<=end:
      triggers.append(line)

  # sort table in time
  if re.search('inspiral',triggers.tableName):
    triggers.sort(key=lambda trig: trig.get_end())
  else:
    triggers.sort(key=lambda trig: trig.get_peak())

  return triggers

def fromtrigfile(file,etg,start=None,end=None,tabletype=None):

  """
    Reads the file object file containing standard columns for the given etg and
    returns either a corresponding lsctable.

    Arguments:

      file : file object
      etg : [ "ihope" | "kw" | "omega" | "omegaclustered" ]
        string defining how to read text file.
        "ihope" is read as per the .csv files from ihope_daily. "omegaclustered"
        is read as per the .txt or .clst produced by the detchar scripts used to
        plot omega online events.

    Keyword arguments:

      start : [ float | int | LIGOTimeGPS ]
        GPS start time of requested period
      end : [ float | int | LIGOTimeGPS ]
        GPS end time of requested period
      tabletype : type
        Specific ligolw table type for output. By default tables will be
        SnglInspiralTable or SnglBurstTable type depending on ETG
  """

  # define control character search
  cchar = re.compile('[#%<!()_\[\]{}:;\'\"]+')

  # set times
  if not start:
    start=0
  if not end:
    end=float('inf')

  if not tabletype:
    etgs = {'inspiral': ['ihope'],\
            'burst':    ['omega','omegaclustered','kw'],\
            'ringdown': []}
    # set up triggers table
    for search,etglist in etgs.items():
      if etg.lower() in etglist:
        tabletype = lsctables.__dict__['Sngl%sTable' % (search.title())]
        break

  triggers = lsctables.New(tabletype)

  # read table and append triggers
  for line in file.readlines():
    # if line starts with a control character:  continue
    if re.match(cchar,line):  continue 
    # read line as trigger
    trig = trigger(line,etg.lower())
    # append trig to table if within requested time
    if re.search('sngl_inspiral',triggers.tableName)\
        and start<=float(trig.get_end())<=end:
      triggers.append(trig)
    elif start<=float(trig.get_peak())<=end:
      triggers.append(trig)

  # sort triggers in time
  if re.search('sngl_inspiral',triggers.tableName):
    triggers.sort(key=lambda trig: trig.get_end())
  else:
    triggers.sort(key=lambda trig: trig.get_peak())

  return triggers

# =============================================================================
# Read injection files
# =============================================================================

def frominjectionfile(file,type,ifo=None):
  
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
        for line in file.readlines():
          inj = lsctables.SimBurst()
          # split data
          sep = re.compile('[\s,]+')
          data = sep.split(line)
          # set attributes
          geocent = LIGOTimeGPS(data[3])
          inj.time_geocent_gps    = geocent.seconds
          inj.time_geocent_gps_ns = geocent.nanoseconds
          inj.waveform            = data[4]
          inj.waveform_number     = int(data[5])
          inj.frequency           = float(data[9])

          # extra columns to be added when I know how

          #inj.q                   = float(data[11])
          #inj.hrss                = float(data[17])
          #inj.ra                  = float(data[19])*24/(2*math.pi)
          #inj.dec                 = 90-(float(data[21])*180/math.pi)
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

  return injtable

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
def bandpass(data, f_low, f_high, sampling, order=4):

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
# Function to write a list of files to a cache
# =============================================================================
def tocache(file,cachelist):

  """
    Reads a list of filepaths and writes each to the given file object cache.
    Filenames must be of the format IFO-USERTAG-STARTTIME-DURATION.ext. Printed
    columns are

    IFO USERTAG STARTTIME DURATION filepath

    for each file
  """

  for f in cachelist:
    lfn = os.path.basename(f)
    head,ext = os.path.splitext(lfn)
    if head.endswith('.xml'):
      head,ext = os.path.splitext(head)
    a,b,c,d = head.split('-')
    print >>file, "%s %s %s %s %s" % (a,b,c,d,f)

  file.flush() 

# =============================================================================
# Read a list of files from frame cache
# =============================================================================

def fromcache(file,ifo='',filetag='',start=0,end=float('inf')):

  """
    Reads a file object cache file containing the columns:

    IFO FILETAG GPSSTART DURATION FILEPATH

    and returns a list of filepaths.
  """

  cache = []
  cchar = re.compile('[#%<!()_\[\]{}:;\'\"]+')
  for line in file:
    if re.match(cchar,line):  continue
    a,b,c,d,f = line.split(' ')
    c = int(c)
    d = int(d)
    if re.match(ifo,a) and re.search(filetag,b)\
    and not segment(start,end).disjoint(segment(c,c+d)):
      cache.append(f.replace('\n',''))

  return cache

# =============================================================================
# Generate a daily ihope cache 
# =============================================================================

def daily_ihope_cache(start,end,ifo,cluster=None,filetype='xml',cat=0):

  """
    Generates cache list of daily ihope INSPIRAL files for give ifo and
    clustering between start and end time.

    Arguments:

      start : [ float | int | LIGOTimeGPS ]
        GPS start time of requested period
      end : [ float | int | LIGOTimeGPS ]
        GPS end time of requested period
      ifo : [ "H1" | "L1" | "V1" ]
        IFO

    Keyword arguments:
      cluster : [ "unclustered" | "100ms" | "30ms" | "16s" ]
        clustering time in human format
      filetype : [ "xml" | "csv" ]
        file format of desired cache contents
      cat : [ 0 | 1 | 2 | 3 | 4 ]
        int veto category of trigger files requested
  """

  # daily path
  ihope_daily_path = os.path.expanduser('~cbc/ihope_daily')

  # set clustering tag
  if cluster==None or cluster.upper()=='UNCLUSTERED':
    cluster_tag='UNCLUSTERED'
  elif cluster.upper()=='100MS':
    cluster_tag='100MILLISEC_CLUSTERED'
  elif cluster.upper()=='30MS':
    cluster_tag='30MILLISEC_CLUSTERED'
  elif cluster.upper()=='16S':
    cluster_tag='16SEC_CLUSTERED'

  # work out days
  days = gps_day_list(start,end) 
  cache=[]
  # loop over days gathering files
  for day in days:
    utc = datetime.datetime(*date.XLALGPSToUTC(day)[:6])
    day_path = os.path.join(ihope_daily_path,utc.strftime("%Y%m"),
                                             utc.strftime("%Y%m%d"))

    if filetype=='xml':
      files = glob.glob(os.path.join(day_path,
                                      ifo+'-INSPIRAL_'+cluster_tag+'*.xml.gz'))
      for line in files:
        if not line.startswith(ihope_daily_path):  continue
        a,b,c,d = re.split('.xml',line,maxsplit=1)[0].split('-')
        trig_start = int(c)
        duration = int(d)
        if start<=trig_start<end or start<(trig_start+duration)<=end:
          cache.append(line.replace('\n',''))

    elif filetype=='csv':
      csvfile = os.path.join(day_path,ifo+'-'+str(cat)+'-INSPIRAL_'+\
                                      cluster_tag+'.csv')
      if os.path.isfile(csvfile):
        cache.append(csvfile)

  return cache

# =============================================================================
# Function to generate an omega online cache
# =============================================================================

def omega_online_cache(start,end,ifo):

  """
    Returns a list of omega online trigger files between the given start and end
    time for the given ifo. For S6 triggers are only available for each IFO on
    it's own site cluster.

    Arguments:

      start : [ float | int | LIGOTimeGPS ]
        GPS start time of requested period
      end : [ float | int | LIGOTimeGPS ]
        GPS end time of requested period
      ifo : [ "H1" | "L1" | "V1" ]
        IFO
  """

  # verify host
  host = getfqdn()
  ifo_host = {'H1':'ligo-wa','H2':'ligo-wa','L1':'ligo-la'}
  if not re.search(ifo_host[ifo],host):
    print >>sys.stderr, "Error: Omega online files are not available for "+\
                        "IFO="+ifo+" on this host."
    return []

  cache = []
  basedir = os.path.expanduser('~omega/online/%s/archive/S6/segments'\
                               % (str(ifo)))

  # triggers are numbered from a given start time, which for S6 is:
  basetime = LIGOTimeGPS(931211808)
  triglength = 64
 
  start_time = int(start-math.fmod(start-basetime,triglength))
  t = start_time
  # loop over time segments constructing files and appending to the list
  while t < end:
    dirstart = "%.10d" % t
    dirend   = "%.10d" % (t+triglength)
    dirpath  = os.path.join(basedir,dirstart+'-'+dirend)
    trigfile = os.path.join(dirpath,'-'.join([ifo,\
                                              'OMEGA_TRIGGERS_CLUSTER',\
                                              dirstart,str(triglength)])+'.txt')
    if os.path.isfile(trigfile):
      cache.append(trigfile) 
    t+=triglength

  return cache 

# =============================================================================
# Function to generate a KW DARM_ERR cache
# =============================================================================
def kw_cache(start,end,ifo):
 
  """
    Returns a list of KW trigger files between the given start and end
    time for the given ifo. For S6 triggers are only available for each IFO on
    it's own site cluster.

    Arguments:

      start : [ float | int | LIGOTimeGPS ]
        GPS start time of requested period
      end : [ float | int | LIGOTimeGPS ]
        GPS end time of requested period
      ifo : [ "H1" | "L1" | "V1" ]
        IFO
  """

  # verify host
  host = getfqdn()
  ifo_host = {'H1':'ligo-wa','H2':'ligo-wa','L1':'ligo-la'}
  if not re.search(ifo_host[ifo],host):
    print >>sys.stderr, "Error: KW files are not available for "+\
                        "IFO="+ifo+" on this host."
    return []

  cache = []
  basedir = os.path.expanduser('~lindy/public_html/triggers/s6')

  # times are numbere from a given start, which for S6 is:
  basetime = LIGOTimeGPS(938736000)
  triglength = 86400
  
  start_time = int(start-math.fmod(start-base,triglength))
  t = start_time
  # loop over time segments constructing file paths and appending to the cache
  while t<end:
    dirstart = str(t)
    dirend   = str(t+triglength)
    dirpath  = os.path.join(basedir,dirstart+'_'+dirend)
    trigfile = os.path.join(dirpath,ifo+'_LSC-DARM_ERR_32_2048.trg')
    if os.path.isfile(trigfile):
      cache.append(trigfile)
    t+=triglength

  return cache

# ==============================================================================
# Function to construct list of days from start and end times
# ==============================================================================
def gps_day_list(start,end):

  """
    This script will construct a list of LIGOTimeGPS days encompassing the start and end GPS times.
  """

  start_d = date.XLALUTCToGPS(datetime.datetime(*date.XLALGPSToUTC(start)[:6])\
                                  .replace(hour=0,minute=0,second=0)\
                                  .timetuple())
  days = []
  day = start_d
  while day<=end:
    days.append(day)
    day+=86400

  return days 

# =============================================================================
# Function to cluster triggers
# =============================================================================

def cluster(triggers,clusterparam='time',width=1,rankparam='snr'):

  """
    Cluster the lsctable triggers in the given clusterparam(s) in order with the
    given clustering width.
  """

  if not isinstance(clusterparam,list) or isinstance(clusterparam,str):
    clusterparam = [clusterparam]
  if not isinstance(width,list):
    width = [width]

  # presort triggers in clusterparam
  if clusterparam[0]=='time':
    if re.search('sngl_inspiral',triggers.tableName):
      triggers.sort(key=lambda trig: trig.get_end())
    else:
      triggers.sort(key=lambda trig: trig.get_peak())
  else:
    triggers.sort(key=lambda trig: trig.__getattribute__(clusterparam[0]))

  cluster = []
  clustered = lsctables.New(type(triggers))
  prevtrig = None
  for trig in triggers:
    if not prevtrig:
      # add first trig to cluster and move on
      prevtrig = trig
      cluster = [trig]
      continue
    # for each clustering parameter, generate cluster and choose loudest
    # by rankparam
    for i in range(len(clusterparam)):
      if clusterparam[i]=='time':
        if re.search('sngl_inspiral',triggers.tableName):
          value = trig.get_end()
          prev = prevtrig.get_end()
        else:
          value = trig.get_peak()
          prev  = prevtrig.get_peak()
      else:
        value = trig.__getattribute__(clusterparam[i])
        prev  = prevtrig.__getattribute__(clusterparam[i])
      if math.fabs(value-prev)>=width[i]:
        cluster.sort(key=lambda item: item.__getattribute__(rankparam),\
                     reverse=True)
        event = cluster[0]
        # sort cluster attributes for SnglBurstTables
        #if triggers.tableName=='sngl_burst:table':
        #  event.cluster_flow   = min([t.flow for t in cluster])
        #  event.cluster_fhigh  = max([t.fhigh for t in cluster])
        #  event.cluster_length = len(cluster)
        clustered.append(event)
        cluster=[trig]
        break
      if i==(len(clusterparam)-1):
        cluster.append(trig)

    prevtrig = trig

  # process final cluster
  if cluster!=[]:
    cluster.sort(key=lambda item: item.__getattribute__(rankparam),\
                 reverse=True)
    event = cluster[0]
    # sort cluster attributes for SnglBurstTables
    #if triggers.tableName=='sngl_burst:table':
    #  event.cluster_flow   = min([t.flow for t in cluster])
    #  event.cluster_fhigh  = max([t.fhigh for t in cluster])
    #  event.cluster_length = len(cluster)
    clustered.append(event)

  return clustered

# ==============================================================================
# Calculate poisson significance of coincidences
# ==============================================================================
def coinc_significance(gwtriggers,auxtriggers,window=1,snrthresh=8,livetime=None):

  # get livetime
  if not livetime:
    if re.search('Inspiral',gwtriggers):
      start = min([t.get_end() for t in gwtriggers])
      start = max([t.get_end() for t in gwtriggers])
    else:
      start = min([t.get_peak() for t in gwtriggers])
      start = max([t.get_peak() for t in gwtriggers])

  # calculate probability of a GW trigger falling within the window
  gwprob = len(gwtriggers)*window/livetime

  # grab auxiliary triggers above threshold
  auxtriggers = lsctables.New(type(tableName))
  auxtriggers.extend([t for t in auxtriggers if t.snr>snrthresh])

  # calculate mean of Poisson distribution
  mu = gwprob * len(auxtriggers)

  # generate segments around auxiliary triggers
  if re.search('Inspiral',auxtriggers):
    coincsegs = segmentlist([segment(t.get_end()-window/2,t.get_end()+window/2)\
                             for t in triggers])

  else:
    coincsegs = segmentlist([segment(t.get_peak()-window/2,\
                                     t.get_peak()+window/2)\
                           for t in triggers])

  coinctriggers = lsctables.New(type(gwtriggers))
  if re.search('Inspiral',gwtriggers):
    coinctriggers.extend([g for g in gwtriggers if g.get_end() in coincsegs])
  else:
    coinctriggers.extend([g for g in gwtriggers if g.get_peak() in coincsegs])

  if len(coinctriggers)<1:
    significance = 0
  else:
    significance = -math.log(gammainc(len(coinctriggers),mu),10)

  return significance
