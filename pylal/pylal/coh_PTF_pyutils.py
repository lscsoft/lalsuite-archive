import os,sys,numpy,glob

from pylal import grbsummary,antenna,llwapp
from glue import segmentsUtils
from glue.ligolw import lsctables
from glue.ligolw.utils import process as ligolw_process

# define newsnr
def new_snr( snr, chisq, W, Q, N ):

  return snr / ( ( 1 + (chisq/W)**(Q/N))/2 ) ** (1/Q)

def get_signal_vetoes( trigger, bankq=0, bankn=0, autoq=0, auton=0, chiq=0,\
                       chin=0, fResp = None ):

  """
    Calculate and apply all signal based vetoes for the given trigger. This
    will generate eight different SNRs for the different stages of this
    process and return them as a dictionary.
    
    The signal based vetoes are applied in the following order:
    1) Chi-squared
    2) Coherent SNR < 6
    3) Bank chi-squared
    4) Auto veto (continuous chi-squared)
    5) Null SNR ( ( CoincSNR^2 - CohSNR^2 )^(1/2) )
    6) Single-detector SNR (from two most sensitive IFOs)
    7) None
    8) None
  """

  # set variables
  sbvs = {}
  q = bankq
  nhigh = bankn
  q2 = autoq
  nhigh2 = auton

  # get new SNR (from chisq)
  if trigger.chisq == 0:
    sbvs[1] = 0
  else:
    if trigger.chisq < 60:
      sbvs[1] = trigger.snr
    else:
      sbvs[1] = new_snr( trigger.snr, trigger.chisq, 60, chiq, chin )

  # get SNR2
  if trigger.snr < 6.:
    sbvs[2] = 0
  else:
    sbvs[2] = sbvs[1]

  # get new bank SNR
  if trigger.bank_chisq < 40:
    sbvs[3] = sbvs[2]
  else:
    bank_new_snr = new_snr( trigger.snr, trigger.bank_chisq, 40, q, nhigh )
    if bank_new_snr < 6.:
      sbvs[3] = 0
    else:
      sbvs[3] = sbvs[2]

  # get new auto SNR
  if trigger.cont_chisq < 40:
    sbvs[4] = sbvs[3]
  else:
    auto_new_snr = new_snr( trigger.snr, trigger.cont_chisq, 160, q2, nhigh2 )
    if auto_new_snr < 6.:
      sbvs[4] = 0
    else:
      sbvs[4] = sbvs[3]

  # get null affected SNR
  ifos = [ trigger.ifos[i*2:(i*2)+2]\
           for i in range(int(len(trigger.ifos)/2)) ]
  ifoAtt = { 'G1':'g', 'H1':'h1', 'H2':'h2', 'L1':'l', 'T1':'t', 'V1':'v' }
  # null SNR = sqrt( coincSNR^2 - cohSNR^2 )
  if len(ifos)<3:
    null_snr=0
  else:
    null_snr = ( sum([ trigger.__getattribute__('snr_%s' % ifoAtt[ifo])**2\
                       for ifo in ifos ]) - trigger.snr**2 )**0.5

  if null_snr > 3.5 and trigger.snr < 30:
    sbvs[5] = sbvs[4] * 1/(null_snr - 2.5)
  elif trigger.snr > 30:
    null_threshold = 3.5 + (trigger.snr -30)*50./700.
    if null_snr > null_threshold:
      sbvs[5]=sbvs[4]*1/(null_snr-(null_threshold-1))
    else:
      sbvs[5] = sbvs[4]
  else:
    sbvs[5] = sbvs[4]

  # apply IFO sensitivity cut
  ifoSens = []
  # get ifo sensitivity measure
  for ifo in ifos:
    ifoSens.append( (ifo, trigger.__getattribute__('sigmasq_%s' %ifoAtt[ifo])\
                              * fResp[ifo] ) )
  # rank and test two most sensitive IFOs
  ifoSens.sort( key=lambda (ifo,sens): sens, reverse=True )
  sbvs[8] = sbvs[5]
  for i in range(2):
    ifo = ifoSens[i][0]
    if trigger.__getattribute__('snr_%s' % ifoAtt[ifo]) < 4:
      sbvs[8] = 0
    
  # set old BestNRs
  sbvs[6] = sbvs[8]
  sbvs[7] = sbvs[6]

  return sbvs

def calculate_contours( bankq=0, bankn=0, autoq=0, auton=0, chiq=0, chin=0,\
                        nullt=0, nullo=0 ):

  cont_vals = [5.5,6,6.5,7,8,9,10,11]
  num_vals  = len(cont_vals)
  colors    = ['y-','k-','y-','y-','y-','y-','y-','y-']

  snr_vals      = numpy.arange(6,30,0.1)
  snr_high_vals = numpy.arange(30,500,1)
  snr_all       = []
  for snr in snr_vals:
    snr_all.append(snr)
  for snr in snr_high_vals:
    snr_all.append(snr)

  snr_vals = numpy.asarray(snr_all)

  bank_conts = [[],[],[],[],[],[],[],[]]
  auto_conts = [[],[],[],[],[],[],[],[]]
  chi_conts  = [[],[],[],[],[],[],[],[]]
  null_cont  = []

  for snr in snr_vals:
    for i in range(num_vals):
      bank_val = (snr/cont_vals[i])**bankq
      if (bank_val > 1):
        bank_val = (bank_val*2 - 1)**(bankn/bankq)
      else:
        bank_val = 1E-20
      bank_conts[i].append(bank_val*40)      

      auto_val = (snr/cont_vals[i])**autoq
      if (auto_val > 1):
        auto_val = (auto_val*2 - 1)**(auton/autoq)
      else:
        auto_val = 1E-20
      auto_conts[i].append(auto_val*160)

      chi_val = (snr/cont_vals[i])**chiq
      if (chi_val > 1):
        chi_val = (chi_val*2 - 1)**(chin/chiq)
      else:
        chi_val = 1E-20
      chi_conts[i].append(chi_val*60)

    if snr > nullo:
      null_cont.append(nullt+(snr-nullo)*50./700.)
    else:
      null_cont.append(nullt)

  return bank_conts,auto_conts,chi_conts,null_cont,snr_vals,colors

def plot_contours( axis, snr_vals, contours, colors ):

  for i in range(len(contours)):
    plot_vals_x = []
    plot_vals_y = []
    for j in range(len(snr_vals)):
      if contours[i][j] > 1E-15:
        plot_vals_x.append(snr_vals[j])
        plot_vals_y.append(contours[i][j])
    axis.plot(plot_vals_x,plot_vals_y,colors[i])

def readSegFiles(segdir):

  times = {}
  for name,fileName in\
      zip( ["buffer",       "off",             "on"],\
           ["bufferSeg.txt","offSourceSeg.txt","onSourceSeg.txt"] ):

    segs = segmentsUtils.fromsegwizard(open(os.path.join(segdir,fileName), 'r'))
    if len(segs)>1:
      raise AttributeError, 'More than one segment, an error has occured.'
    times[name] = segs[0]
  return times
       
def makePaperPlots():

  import pylab

  pylab.rcParams.update({
    "text.usetex": True,
    "text.verticalalignment": "center",
#    "lines.markersize": 12,
#    "lines.markeredgewidth": 2,
#    "lines.linewidth": 2.5,
    "figure.figsize": [8.0, 6.0],
    "font.size": 20,
    "axes.titlesize": 16,
    "axes.labelsize": 24,
    "xtick.labelsize": 18,
    "ytick.labelsize": 18,
    "legend.fontsize": 18,
    "font.family": "serif",
    "font.weight": "bold",
    })

def get_ra_dec(grbFile):

  """
    DEPRECATED
  """

  ext_trigs = grbsummary.load_external_triggers(grbFile)
  ra = ext_trigs[0].event_ra
  dec = ext_trigs[0].event_dec
  return ra,dec

def read_sigma_vals( sigmaFile ):

  """
    DEPRECATED
  """

  sigmaVals = {}
  file = open(sigmaFile,'r')
  for line in file:

    line = line.replace('\n','')
    ifo,min,max = line.split(' ')
    sigmaVals[ifo + 'min'] = float(min)
    sigmaVals[ifo + 'max'] = float(max)

  return sigmaVals
    
def get_det_response( ra, dec, triggerTime ):

  """
    Return detector response for complete set of IFOs for given sky location
    and time. Assumed inclination = 0, and polarization = 0. 
  """

  f_plus  = {}
  f_cross = {}
  inclination   = 0
  polarization  = 0
  for ifo in ['G1','H1','H2','L1','T1','V1']:
    f_plus[ifo],f_cross[ifo],_,_ = antenna.response( triggerTime, ra, dec,\
                                                     inclination, polarization,\
                                                     'degree', ifo )
  return f_plus,f_cross

def append_process_params( xmldoc, args, version, date ):

  """
    Construct and append process and process_params tables to ligolw.Document
    object xmldoc, using the given sys.argv variable args and other parameters.
  """

  xmldoc.childNodes[-1].appendChild(lsctables.New(lsctables.ProcessTable))
  xmldoc.childNodes[-1].appendChild(lsctables.New(lsctables.ProcessParamsTable))

  # build and seed process params
  progName = args[0]
  process = llwapp.append_process( xmldoc, program=progName,\
                                   version=version,\
                                   cvs_repository = 'lscsoft',\
                                   cvs_entry_time = date)
  params = []
  for i in range(len(args)):
    p = args[i]
    if not p.startswith('-'):
      continue
    v = ''
    if i < len(sys.argv)-1:
      v = sys.argv[i+1]
    params.append( map( unicode, (p,'string',v) ) )

  ligolw_process.append_process_params(xmldoc,process,params)

  return xmldoc

def identify_bad_injections(log_dir):

  files = glob.glob(os.path.join(log_dir,"*err"))

  badInjs = []

  for file in files:
    if os.stat(file)[6] != 0:
      fp = open(file,"r")
      conts = fp.read()
      if conts.find('terminated') != -1:
        conts=conts.split('\n')
        for line in conts:
          line = line.split(' ')
          line = [entry.replace(',','') for entry in line if entry]
          if 'terminated' in line:
            injDict = {}
            injDict['mass1'] = float(line[6])
            injDict['mass2'] = float(line[8])
            injDict['spin1x'] = float(line[10])
            injDict['spin1y'] = float(line[12])
            injDict['spin1z'] = float(line[14])
            injDict['spin2x'] = float(line[16])
            injDict['spin2y'] = float(line[18])
            injDict['spin2z'] = float(line[20])
            if not injDict in badInjs:
              badInjs.append(injDict)
  return badInjs

def remove_bad_injections(sims,badInjs):

  new_sims = []
  for sim in sims:
    for badInj in badInjs:
      if (abs(sim.mass1-badInj['mass1'])) < 0.001:
        if (abs(sim.mass2-badInj['mass2'])) < 0.001:
          if (abs(sim.spin1x-badInj['spin1x'])) < 0.001:
            if (abs(sim.spin1y-badInj['spin1y'])) < 0.001:
              if (abs(sim.spin1z-badInj['spin1z'])) < 0.001:
                if (abs(sim.spin2x-badInj['spin2x'])) < 0.001:
                  if (abs(sim.spin2y-badInj['spin2y'])) < 0.001:
                    if (abs(sim.spin2z-badInj['spin2z'])) < 0.001:
                      print "Removing injection:",sim.mass1,sim.mass2,sim.spin1x,sim.spin1y,sim.spin1z,sim.spin2x,sim.spin2y,sim.spin2z
                      break
    else:
      new_sims.append(sim)

  return new_sims
