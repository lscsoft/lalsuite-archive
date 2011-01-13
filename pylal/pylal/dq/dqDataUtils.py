#!/usr/bin/env python

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

from matplotlib import mlab, use
use('Agg')
import pylab

# =============================================================================
# Execute shell cmd and get output
# =============================================================================
def make_external_call(cmd,shell=False):
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
  """ 

  # read etg
  etg=etg.lower()

  # if given string, split on space or comma
  if isinstance(data,str):
    sep = re.compile('[\s,]+')
    data = sep.split(data)

  # set up trig object
  if re.search('ihope',etg):
    trig = lsctables.SnglInspiral()
  else:
    trig = lsctables.SnglBurst()

  # ==========
  # IHOPEDAILY
  # ==========
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

  # ===========
  # KLIENE-WELLE
  # ===========
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
    trig.central_frequency     = float(data[3])
    trig.energy                = float(data[4])
    trig.norm_energy           = float(data[5])
    trig.n_pix                 = float(data[6])
    trig.significance          = float(data[7])
    trig.N                     = float(data[8])
    trig.snr                   = math.sqrt(trig.norm_energy-trig.n_pix)
 
  # ==========
  # OMEGA
  # ==========
  if etg=='omega' or etg=='wpipe':
    # peak time
    peak = LIGOTimeGPS(data[0])
    trig.peak_time           = peak.seconds
    trig.peak_time_ns        = peak.nanoseconds
    # central frequency
    trig.central_frequency   = float(data[1])
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
    trig.norm_energy         = float(data[4])
    # cluster parameters
    trig.cluster_size        = float(data[5])
    trig.cluster_norm_energy = float(data[6])
    trig.cluster_number      = float(data[7])
    # SNR
    trig.snr                 = math.sqrt(2*trig.norm_energy)
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
    # bandwidth, and flow,fhigh,central_frequency
    trig.flow                    = float(data[8])
    trig.fhigh                   = float(data[9])
    trig.bandwidth               = trig.fhigh - trig.flow
    trig.central_frequency       = trig.flow + 0.5*trig.bandwidth
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
#def totrigxml(file,trigs,options)

#  xmldoc = ligolw.Document()
 


def totrigfile(file,table,etg,header=True):
  """Write the lsctables.table to the given file object file, in the standard format for the given etg"""

  etg = etg.lower()

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
    columns = ['peak_time','central_frequency','duration',\
               'bandwidth','norm_energy',\
               'cluster_size','cluster_norm_energy','cluster_number']

  elif re.match('kw',etg.lower()):
    columns = ['peak_time','start_time','end_time','central_frequency',\
               'energy','norm_energy','num_pixels','significance','N']

  if etg=='ihope':  d = ','
  else:  d=' '

  # print header
  if header:
    header = '#'+' '+' '.join(columns)
    print >>file, header
  # print triggers
  for row in table:
    line = []
    for col in columns:
       if col not in table.columnnames:  continue
       entry = ''
       if re.match('ihope',etg.lower()):
         entry = str(row.__getattribute__(col))
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
def fromtrigxml(file,start=None,end=None,tablename='sngl_inspiral'):
  """Reads a list of Sngl{Burst,Inspiral}Table from the given table from the xml file object file"""

  if not start:
    start=0
  if not end:
    end=float('inf')

  xmldoc,digest = utils.load_fileobj(file,gz=file.name.endswith('gz'))
  trigs = table.get_table(xmldoc,':'.join([tablename,'table']))

  for line in trigs:
    if re.search('inspiral',tablename):  t = line.get_end()
    else:  t = line.get_peak()
    if not start<=float(t)<=end:
      trigs.remove(line)

  if re.search('inspiral',tablename):
    trigs.sort(key=lambda trig: trig.get_end())
  else:
    trigs.sort(key=lambda trig: trig.get_peak())

  return trigs

def fromtrigfile(file,start=None,end=None,etg="ihope"):
  """Reads the file object file containing standard columns for the given etg and returns either a SnglInspiralTable or SnglBurstTable."""

  # define control character search
  cchar = re.compile('[#%<!()_\[\]{}:;\'\"]+')

  # set times
  if not start:
    start=0
  if not end:
    end=float('inf')

  # set table
  if re.search('ihope',etg.lower()):
    t = 'SnglInspiralTable'
  else:
    t = 'SnglBurstTable'

  trigs = lsctables.New(lsctables.__dict__[t])
  for line in file.readlines():
    if re.match(cchar,line):  continue 
    trig = trigger(line,etg.lower())
    if re.search('sngl_inspiral',trigs.tableName)\
        and start<=float(trig.get_end())<=end:
      trigs.append(trig)
    elif start<=float(trig.get_peak())<=end:
      trigs.append(trig)

  if re.search('sngl_inspiral',trigs.tableName):
    trigs.sort(key=lambda trig: trig.get_end())
  else:
    trigs.sort(key=lambda trig: trig.get_peak())

  return trigs

# =============================================================================
# Read injection files
# =============================================================================
def frominjectionfile(file,type='inspiral',ifo=None):
  
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

    if injtype=='burst':
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
          inj.q                   = float(data[11])
          inj.hrss                = float(data[17])
          inj.ra                  = float(data[19])*24/(2*pi)
          inj.dec                 = 90-(float(data[21])*180/pi)
          h_delay = LIGOTimeGPS(data[41])
          inj.h_peak_time         = inj.time_geocent_gps+h_delay.seconds
          inj.h_peak_time_ns      = inj.time_geocent_gps_ns+h_delay.nanoseconds
          l_delay = LIGOTimeGPS(data[43])
          inj.l_peak_time         = inj.time_geocent_gps+l_delay.seconds
          inj.l_peak_time_ns      = inj.time_geocent_gps_ns+l_delay.nanoseconds
          v_delay = LIGOTimeGPS(data[43])
          inj.v_peak_time         = inj.time_geocent_gps+v_delay.seconds
          inj.v_peak_time_ns      = inj.time_geocent_gps_ns+v_delay.nanoseconds

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
  """Reads a list of filepaths and writes each to the given file object cache. Filenames must be of the format IFO-USERTAG-STARTTIME-DURATION.ext. Printed columns are

IFO USERTAG STARTTIME DURATION filepath

for each file"""

  for f in cachelist:
    lfn = os.path.basename(f)
    head,ext = os.path.splitext(lfn)
    if head.endswith('.xml'):
      head,ext = os.path.splitext(head)
    a,b,c,d = head.split('-')
    print >>file, "%s %s %s %s %s" % (a,b,c,d,f)

  file.flush() 

def fromcache(file,ifo='',filetag='',start=0,end=float('inf')):
  """Reads a file object cache file containing the columns:

IFO FILETAG GPSSTART DURATION FILEPATH

and returns a list of filepaths."""
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
# Function to generate an daily ihope cache 
# =============================================================================
def daily_ihope_cache(start,end,ifo,cluster=None,filetype='xml',cat=0):
  """
  Generates cache list of daily ihope INSPIRAL files for give ifo and clustering (None,'30ms','100ms', or '16s') between start and end time
  """

  # daily path
  ihope_daily_path = '/archive/home/cbc/ihope_daily'

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

  # verify host
  host = getfqdn()
  ifo_host = {'H1':'ligo-wa','H2':'ligo-wa','L1':'ligo-la'}
  if not re.search(ifo_host[ifo],host):
    print >>sys.stderr, "Error: Omega online files are not available for "+\
                        "IFO="+ifo+" on this host."
    sys.exit(1)

  # set variables
  cache = []
  basedir = os.path.expanduser('~omega/online/'+str(ifo)+'/archive/S6/segments')
  basetime = LIGOTimeGPS(931211808)
  triglength = 64
 
  start_time = int(start-math.fmod(start-basetime,triglength))
  t = start_time
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
 
  # verify host
  host = getfqdn()
  ifo_host = {'H1':'ligo-wa','H2':'ligo-wa','L1':'ligo-la'}
  if not re.search(ifo_host[ifo],host):
    print >>sys.stderr, "Error: KW files are not available for "+\
                        "IFO="+ifo+" on this host."
    sys.exit(1)

  # set variables
  cache = []
  basedir = os.path.expanduser('~lindy/public_html/triggers/s6')
  basetime = LIGOTimeGPS(938736000)
  triglength = 86400
  
  start_time = int(start-math.fmod(start-base,triglength))
  t = start_time
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
  """This script will construct a list of LIGOTimeGPS days encompassing the start and end GPS times."""

  start_d = date.XLALUTCToGPS(datetime.datetime(*date.XLALGPSToUTC(start)[:6])\
                                  .replace(hour=0,minute=0,second=0)\
                                  .timetuple())
  days = []
  day = start_d
  while day<=end:
    days.append(day)
    day+=86400

  return days 

# ==============================================================================
# Function to grab effective distance from ihope cache of trig files
# ==============================================================================
def grab_effective_distance(cache,time=False):
  distance=[]
  time=[]
  # grab data
  for file in cache:
    try:
      xmldoc = utils.load_filename(file,gz=file.endswith("gz"))
      value_table = table.get_table(xmldoc,lsctables.SummValueTable.tableName)
      for line in value_table:
        if line.name == 'inspiral_effective_distance':
          distance.append(line.value)
          time.append((line.start_time+line.end_time)/2)
          break
    except:
      continue
  if time:  return distance,time
  else:  return distance

# =============================================================================
# Function to cluster triggers
# =============================================================================
def cluster(trigs,clusterparam='time',width=1,rankparam='snr'):
  """Cluster given triggerlist in the given clusterparam(s) in order with the given clustering width."""

  if not isinstance(clusterparam,list) or isinstance(clusterparam,str):
    clusterparam = [clusterparam]
  if not isinstance(width,list):
    width = [width]

  # presort triggers in clusterparam
  if clusterparam[0]=='time':
    if re.search('sngl_inspiral',trigs.tableName):
      trigs.sort(key=lambda trig: trig.get_end())
    else:
      trigs.sort(key=lambda trig: trig.get_peak())
  else:
    trigs.sort(key=lambda trig: trig.__getattribute__(clusterparam[0]))

  cluster = []
  clustered = lsctables.New(type(trigs))
  prevtrig = None
  for trig in trigs:
    if not prevtrig:
      # add first trig to cluster and move on
      prevtrig = trig
      cluster = [trig]
      continue
    # for each clustering parameter, generate cluster and choose loudest
    # by rankparam
    for i in range(len(clusterparam)):
      if clusterparam[i]=='time':
        if re.search('sngl_inspiral',trigs.tableName):
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
        if trigs.tableName=='sngl_burst:table':
          event.cluster_flow   = min([t.flow for t in cluster])
          event.cluster_fhigh  = max([t.fhigh for t in cluster])
          event.cluster_length = len(cluster)
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
    if trigs.tableName=='sngl_burst:table':
      event.cluster_flow   = min([t.flow for t in cluster])
      event.cluster_fhigh  = max([t.fhigh for t in cluster])
      event.cluster_length = len(cluster)
    clustered.append(event)

  return clustered

# ==============================================================================
# Dump flags from segment database
# ==============================================================================
def dump_flags(squery="select ifos,name,version from segment_definer",\
               start=None,end=None,ifo=None,segment_url=None):

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
  seg_def_table = table.get_table(xmldoc,\
                                  lsctables.SegmentDefinerTable.tableName)
  
  flags = []
  for line in seg_def_table:
    if ifo and line.ifo!=ifo:  continue
    flag = ':'.join([line.ifos.rstrip(),line.name,str(line.version)])
    flags.append(flag)

  return flags
