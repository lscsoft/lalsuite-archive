#!/usr/bin/env python

from __future__ import division
import types,sys,os,re,numpy,copy,math,shlex,subprocess,datetime
from socket import getfqdn
from glue.ligolw import ligolw
from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw import utils
from glue.segments import segment, segmentlist
from pylal import Fr,date
from pylal.xlal.datatypes.ligotimegps import LIGOTimeGPS
from math import sqrt,log,fabs,exp,cos,pi

# Hey, scipy, shut up about your nose already.
import warnings
warnings.filterwarnings("ignore")
from scipy import signal as signal
from scipy.fftpack import fft, ifft, ifftshift, fft2, ifft2

from matplotlib import mlab, use
use('Agg')
import pylab

# =============================================================================
# Function to execute shell command and get output
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
# Class to define trigger object from data string
# =============================================================================
class trigger:
  def __init__(self,data,etg):
    etg=etg.lower()
    if isinstance(data,str):
      sep = re.compile('[\s,]+')
      data = sep.split(data)

    # ==========
    # IHOPE XML
    # ==========
    #== test if converting from inspiral xml trig to dqDataUtils.trigger
    if re.match('ihope',etg.lower()) and \
        'end_time' in self.__dict__.keys():
      pass
    # ==========
    # IHOPEDAILY
    # ==========
    if re.match('ihope',etg.lower()):
      #== comma separated values are:
      # end_time,end_time_ns,ifo,snr,mass1,mass2,mtotal,eta,event_duration,
      # template_duration,eff_distance,chisq,chisq_dof,bank_chisq,bank_chisq_dof
      # cont_chisq,cont_chisq_dof
      self.end_time              = int(data[0])
      self.end_time_ns           = int(data[1])
      self.peak_time             = self.end_time
      self.peak_time_ns          = self.end_time_ns
      self.ifo                   = str(data[2])
      self.snr                   = float(data[3])
      self.mass1                 = float(data[4])
      self.mass2                 = float(data[5])
      self.mtotal                = float(data[6])
      self.eta                   = float(data[7])
      self.event_duration        = float(data[8])
      self.template_duration     = float(data[9])
      self.eff_distance          = float(data[10])
      self.chisq                 = float(data[11])
      self.chisq_dof             = float(data[12])
      self.bank_chisq            = float(data[13])
      try:
        self.bank_chisq_dof      = float(data[14])
      except:
        self.bank_chisq_dof = None
      try:
        self.cont_chisq          = float(data[15])
        self.cont_chisq_dof      = float(data[16])
      except ValueError:
        self.cont_chisq = None
        self.cont_chisq_dof = None
    # ===========
    # KLIENE-WELLE
    # ===========
    elif etg=='kw':
      peak  = LIGOTimeGPS(data[0])
      self.peak_time             = peak.seconds
      self.peak_time_ns          = peak.nanoseconds
      start = LIGOTimeGPS(data[1])
      self.start_time            = start.seconds
      self.start_time_ns         = start.nanoseconds
      end = LIGOTimeGPS(data[2])
      self.end_time              = end.seconds
      self.end_time_ns           = end.nanoseconds
      self.frequency             = float(data[3])
      self.energy                = float(data[4])
      self.norm_energy           = float(data[5])
      self.n_pix                 = float(data[6])
      self.significance          = float(data[7])
      self.N                     = float(data[8])
      self.snr                   = sqrt(self.norm_energy-self.n_pix)
 
    # ==========
    # OMEGA
    # ==========
    if etg=='omega' or etg=='wpipe':
      peak = LIGOTimeGPS(data[0])
      self.peak_time           = peak.seconds
      self.peak_time_ns        = peak.nanoseconds
      self.frequency           = float(data[1])
      self.duration            = LIGOTimeGPS(float(data[2]))
      start = peak-self.duration
      self.start_time          = start.seconds
      self.start_time_ns       = start.nanoseconds
      stop = peak+self.duration
      self.end_time            = stop.seconds
      self.end_time_ns         = stop.nanoseconds
      self.bandwidth           = float(data[3])
      self.start_frequency     = self.frequency - 0.5*self.bandwidth
      self.end_frequency     = self.frequency + 0.5*self.bandwidth
      self.norm_energy         = float(data[4])
      self.cluster_size        = float(data[5])
      self.cluster_norm_energy = float(data[6])
      self.cluster_number      = float(data[7])
      self.snr                 = sqrt(2*self.norm_energy)
    #== follow Cadonati's clustering output
    if etg=='omegaclustered' or etg=='wpipeclustered':
      start = LIGOTimeGPS(data[0])
      self.start_time          = start.seconds
      self.start_time_ns       = start.nanoseconds
      stop  = LIGOTimeGPS(data[1])
      self.end_time            = stop.seconds
      self.end_time_ns         = stop.nanoseconds
      peak  = LIGOTimeGPS(data[2])
      self.peak_time           = peak.seconds
      self.peak_time_ns        = peak.nanoseconds
      self.cluster_N           = float(data[5])
      self.start_frequency     = float(data[8])
      self.end_frequency       = float(data[9])
      self.bandwidth           = self.end_frequency - self.start_frequency
      self.frequency           = self.start_frequency + 0.5*self.bandwidth
      self.cluster_size        = float(data[10])
      self.cluster_norm_energy = float(data[11])
      self.norm_energy         = float(data[12])
      self.snr                 = sqrt(2*self.norm_energy)

  def get_end(self):
    return LIGOTimeGPS(self.end_time,self.end_time_ns)

  def get_start(self):
    return LIGOTimeGPS(self.start_time,self.start_time_ns)

  def get_peak(self):
    return LIGOTimeGPS(self.peak_time,self.peak_time_ns)

  def __getattribute__(self,name):
    if name=='test':
      return 0.
    else:
      return self.__dict__[name]

# =============================================================================
# Function to write triggers to file in etg standard form
# =============================================================================
def totrigfile(file,trigs,etg,header=True):
  """Write the list of dqDataUtils.trigger objects to the given file object, in the standard format for the given etg"""

  if re.match('ihope',etg.lower()):
    #== comma separated values are:
    columns = ['end_time','end_time_ns','ifo','snr','mass1','mass2','mtotal',\
               'eta','event_duration','template_duration','eff_distance',\
               'chisq','chisq_dof','bank_chisq','bank_chisq_dof','cont_chisq',\
               'cont_chisq_dof']

  elif re.match('omega',etg.lower()) or re.match('wpipe',etg.lower()):
    columns = ['peak_time','frequency','duration','bandwidth','norm_energy',\
               'cluster_size','cluster_norm_energy','cluster_number']

  elif re.match('kw',etg.lower()):
    columns = ['peak_time','start_time','end_time','frequency','energy',\
               'norm_energy','num_pixels','significance','N']

  if etg=='ihope':  d = ','
  else:  d=' '

  #== print header
  if header:
    header = '#'+' '+' '.join(columns)
    print >>file, header
  #== print triggers
  for trig in trigs:
    line = []
    for col in columns:
       entry = ''
       if re.match('ihope',etg.lower()):
         if col in trig.__dict__.keys():
           entry = str(trig.__getattribute__(col))
       else:
         if col=='peak_time':
           entry = str(trig.get_peak())
         elif col=='start_time':
           entry = str(trig.get_start())
         elif col=='end_time':
           entry = str(trig.get_end())
         elif col in trig.__dict__.keys():
           entry = str(trig.__getattribute__(col))
       line.append(entry)

    print >>file, d.join(line)
  
# =============================================================================
# Function to load triggers from cache
# =============================================================================
def fromtrigxml(file,start=None,end=None,table='sngl_inspiral'):
  """Reads a list of dqDataUtils.trigger objects from the given table from the xml file object file"""

  if not start:
    start=0
  if not end:
    end=float('inf')

  trigs=[]
  xmldoc,digest = ligolw_utils.load_fileobj(file)
  trig_table = table.get_table(xmldoc,':'.join([table,'table']))

  for line in trig_table:
    trig = trigger(line,etg.lower(),filetype='xml')
    if start<=float(trig.get_peak())<=end:
      trigs.append(trig)

  trigs.sort(key=lambda trig: trig.get_peak())

  return trigs

def fromtrigfile(file,start=None,end=None,etg="ihope"):
  """Reads a list of dqDataUtils.trigger objects from the text file object file containing standard columns for the given etg."""
  #== define control character search
  cchar = re.compile('[#%<!()_\[\]{}:;\'\"]+')

  if not start:
    start=0
  if not end:
    end=float('inf')

  trigs = []
  for line in file.readlines():
    if re.match(cchar,line):  continue 
    trig = trigger(line,etg.lower())
    if start<=float(trig.get_peak())<=end:
      trigs.append(trig)
  
  trigs.sort(key=lambda trig: trig.get_peak())

  return trigs

# =============================================================================
# Function to load triggers from cache
# =============================================================================
def parse_trigs(trigs,segs,inclusive=True):
  """This function will parse a list of triggers using the given list of segments and will return: if inclusive==True, all triggers included in any segments; if inclusive==False, all triggers not found in any segments."""

  if inclusive:
    selected=[]
    #== work out vetoed triggers
    selected = [t for t in trigs if t.get_peak() in segs]
    #for seg in segs:
    #  selected.extend([t for t in trigs if seg[0]<=t.get_peak()<=seg[1]])
      
  else:
    #== work out unvetoed triggers
    selected = [t for t in trigs if not t.get_peak() in segs]
    #for seg in segs:
    #  selected = [t for t in trigs if t.get_peak()<seg[0]\
    #                               or t.get_peak()>seg[1]]
  return selected 

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
def bandpass(data, f_low, f_high, sampling, order=4):
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
    a,b,c,d = head.split('-')
    print >>file, "%s %s %s %s %s" % (a,b,c,d,f)

  file.flush() 

def fromcache(file):
  """Reads a file object cache file containing the columns:

IFO USERTAG GPSSTART DURATION FILEPATH

and returns a list of filepaths."""
  cache = []
  cchar = re.compile('[#%<!()_\[\]{}:;\'\"]+')
  for line in file:
    if re.match(cchar,line):  continue
    a,b,c,d,f = line.split(' ')
    cache.append(f.replace('\n',''))
  
  return cache

# =============================================================================
# Function to generate an daily ihope cache 
# =============================================================================
def daily_ihope_cache(start,end,ifo,cluster=None,filetype='xml',cat=0):
  """
  Generates cache list of daily ihope INSPIRAL files for give ifo and clustering (None,'30ms','100ms', or '16s') between start and end time
  """

  #== daily path
  ihope_daily_path = '/archive/home/cbc/ihope_daily'

  #== set clustering tag
  if cluster==None or cluster.upper()=='UNCLUSTERED':
    cluster_tag='UNCLUSTERED'
  elif cluster.upper()=='100MS':
    cluster_tag='100MILLISEC_CLUSTERED'
  elif cluster.upper()=='30MS':
    cluster_tag='30MILLISEC_CLUSTERED'
  elif cluster.upper()=='16S':
    cluster_tag='16SEC_CLUSTERED'

  #== work out days
  days = gps_day_list(start,end) 
  cache=[]
  #== loop over days gathering files
  for day in days:
    utc = datetime.datetime(*date.XLALGPSToUTC(day)[:6])
    day_path = os.path.join(ihope_daily_path,utc.strftime("%Y%m"),
                                             utc.strftime("%Y%m%d"))

    if filetype=='xml':
      ls_command = 'ls '+day_path+'/'+ifo+'-INSPIRAL_'+cluster_tag+\
             '*.xml.gz'
      cache_out = Popen(ls_command,shell=True,stdout=PIPE,stderr=PIPE)
      for line in cache_out.stdout.readlines():
        trig_start = int(line.split('.xml')[0].split('-')[-2])
        duration = int(line.split('.xml')[0].split('-')[-2])
        if start<=trig_start<end or start<(trig_start+duration)<=end:
          cache.append(line.replace('\n',''))

      cache_out.stdout.close()

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

  #== verify host
  host = getfqdn()
  ifo_host = {'H1':'ligo-wa','H2':'ligo-wa','L1':'ligo-la'}
  if not re.search(ifo_host[ifo],host):
    print >>sys.stderr, "Error: Omega online files are not available for "+\
                        "IFO="+ifo+" on this host."
    sys.exit(1)

  #== set variables
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
 
  #== verify host
  host = getfqdn()
  ifo_host = {'H1':'ligo-wa','H2':'ligo-wa','L1':'ligo-la'}
  if not re.search(ifo_host[ifo],host):
    print >>sys.stderr, "Error: KW files are not available for "+\
                        "IFO="+ifo+" on this host."
    sys.exit(1)

  #== set variables
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
  #== grab data
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

  #== presort triggers in clusterparam
  if clusterparam[0]=='time':
    trigs.sort(key=lambda trig: trig.get_peak())
  else:
    trigs.sort(key=lambda trig: trig.__getattribute__(clusterparam[0]))

  cluster = []
  clustered = []
  prevtrig = None
  for trig in trigs:
    if prevtrig == None:
      #== add first trig to cluster and move on
      prevtrig = trig
      cluster = [trig]
      continue
    for i in range(len(clusterparam)):
      if clusterparam[i]=='time':
        value = trig.get_peak()
        prev  = prevtrig.get_peak()
      else:
        value = trig.__getattribute__(clusterparam[i])
        prev  = prevtrig.__getattribute__(clusterparam[i])
      if math.fabs(value-prev)>=width[i]:
        cluster.sort(key=lambda item: item.__getattribute__(rankparam),\
                     reverse=True)
        clustered.append(cluster[0])
        cluster=[trig]
        break
      if i==(len(clusterparam)-1):
        cluster.append(trig)

    prevtrig = trig

  cluster.sort(key=lambda item: item.__getattribute__(rankparam),\
               reverse=True)
  clustered.append(cluster[0])
  return clustered    
