import os,sys
import pylab
from pylal import grbsummary
from pylal import antenna
import scipy.stats

def get_signal_vetoes(trigger,bankq=0,bankn=0,autoq=0,auton=0,chiq=0,chin=0,sigmaVals = None,fResp = None):
  sbvs = {}
  q = bankq
  nhigh = bankn
  q2 = autoq
  nhigh2 = auton

  if trigger.chisq == 0:
    sbvs['BestNR1'] = 0
  else:
    if trigger.chisq < 60:
      sbvs['BestNR1'] = trigger.snr
    else:
      sbvs['BestNR1'] = trigger.snr/((1 + (trigger.chisq/60.)**(chiq/chin))/2.)**(1./chiq)

#  traceM = (27.5-6.5)/(30-6)
#  traceC = 6.5 - traceM * 6
#  if trigger.null_stat_degen > traceM * trigger.snr + traceC and trigger.snr < 30:
#    sbvs['BestNR2'] = 0
#  else:
#    sbvs['BestNR2'] = sbvs['BestNR1']
#  sbvs['BestNR2'] = sbvs['BestNR1'] * (trigger.snr/trigger.null_stat_degen)**2

  if trigger.snr < 6.:
    sbvs['BestNR2'] = 0
  else:
    sbvs['BestNR2'] = sbvs['BestNR1']

  if trigger.bank_chisq < 40:
    sbvs['BestNR3'] = sbvs['BestNR2']
  else:
    bank_new_snr = trigger.snr/((1 + (trigger.bank_chisq/40.)**(q/nhigh))/2.)**(1./q)
    if bank_new_snr < 6.:
      sbvs['BestNR3'] = 0
    else:
      sbvs['BestNR3'] = sbvs['BestNR2']

  if trigger.cont_chisq < 40:
    sbvs['BestNR4'] = sbvs['BestNR3']
  else:
    auto_new_snr = trigger.snr/((1 + (trigger.cont_chisq/160.)**(q2/nhigh2))/2.)**(1./q2)
    if auto_new_snr < 6.:
      sbvs['BestNR4'] = 0
    else:
      sbvs['BestNR4'] = sbvs['BestNR3']

  traceM = (21.5/24.)
  traceC = 27.5 - 30.*(traceM)

#  if trigger.null_stat_degen > traceM * trigger.snr + traceC and trigger.snr < 30:
#    sbvs['BestNR5'] = 0
#  else:
#    sbvs['BestNR5'] = sbvs['BestNR4'] 


#  compsList = [trigger.snr_h1,trigger.snr_l,trigger.snr_v]
#  compsList.sort()
#  if pylab.sqrt(compsList[1]) < (secondM*trigger.snr + secondC):
#    sbvs['BestNR6'] = 0
#  else:
#    sbvs['BestNR6'] = sbvs['BestNR5']

#  sbvs['BestNR6'] = sbvs['BestNR5']

#  if trigger.snr_h1 < 4 or trigger.snr_l < 4:
#    sbvs['BestNR6'] = 0
#  else:
#    sbvs['BestNR6'] = sbvs['BestNR5']

#  if pylab.sqrt(compsList[1]) < 4:
#    sbvs['BestNR6'] = 0
#  else:
#    sbvs['BestNR6'] = sbvs['BestNR5']

  if trigger.null_statistic > 3.5 and trigger.snr < 30:
    sbvs['BestNR5'] = sbvs['BestNR4'] * 1/(trigger.null_statistic - 2.5)
  elif trigger.snr > 30:
    null_threshold = 3.5 + (trigger.snr -30)*50./700.
    if trigger.null_statistic > null_threshold:
      sbvs['BestNR5']=sbvs['BestNR4']*1/(trigger.null_statistic-(null_threshold-1))
    else:
      sbvs['BestNR5'] = sbvs['BestNR4']
  else:
    sbvs['BestNR5'] = sbvs['BestNR4']

  if trigger.snr_h1 < 4 or trigger.snr_v < 4:
    sbvs['BestNR8'] = 0
  else:
    sbvs['BestNR8'] = sbvs['BestNR5']


  if sbvs['BestNR5'] == 0:
    sbvs['BestNR6'] = 0
  elif scipy.stats.ncx2.ppf(0.00135/2.,2,fResp['H1']*sigmaVals['H1min']*trigger.snr**2) > trigger.snr_h1**2:
    sbvs['BestNR6'] = 0
  elif scipy.stats.ncx2.ppf(1-0.00135/2.,2,fResp['H1']*sigmaVals['H1max']*trigger.snr**2) < trigger.snr_h1**2:
    sbvs['BestNR6'] = 0
  elif scipy.stats.ncx2.ppf(0.00135/2.,2,fResp['L1']*sigmaVals['L1min']*trigger.snr**2) > trigger.snr_l**2:
    sbvs['BestNR6'] = 0
  elif scipy.stats.ncx2.ppf(1-0.00135/2.,2,fResp['L1']*sigmaVals['L1max']*trigger.snr**2) < trigger.snr_l**2:
    sbvs['BestNR6'] = 0
  elif scipy.stats.ncx2.ppf(0.00135/2.,2,fResp['V1']*sigmaVals['V1min']*trigger.snr**2) > trigger.snr_v**2:
    sbvs['BestNR6'] = 0
  elif scipy.stats.ncx2.ppf(1-0.00135/2.,2,fResp['V1']*sigmaVals['V1max']*trigger.snr**2) < trigger.snr_v**2:
    sbvs['BestNR6'] = 0
  else: 
    sbvs['BestNR6'] = sbvs['BestNR5']

  if sbvs['BestNR6'] == 0:
    sbvs['BestNR7'] = 0
  elif scipy.stats.ncx2.ppf(0.0455/2.,2,fResp['H1']*sigmaVals['H1min']*trigger.snr**2) > trigger.snr_h1**2:
    sbvs['BestNR7'] = 0
  elif scipy.stats.ncx2.ppf(1-0.0455/2.,2,fResp['H1']*sigmaVals['H1max']*trigger.snr**2) < trigger.snr_h1**2:
    sbvs['BestNR7'] = 0
  elif scipy.stats.ncx2.ppf(0.0455/2.,2,fResp['L1']*sigmaVals['L1min']*trigger.snr**2) > trigger.snr_l**2:
    sbvs['BestNR7'] = 0
  elif scipy.stats.ncx2.ppf(1-0.0455/2.,2,fResp['L1']*sigmaVals['L1max']*trigger.snr**2) < trigger.snr_l**2:
    sbvs['BestNR7'] = 0
  elif scipy.stats.ncx2.ppf(0.0455/2.,2,fResp['V1']*sigmaVals['V1min']*trigger.snr**2) > trigger.snr_v**2:
    sbvs['BestNR7'] = 0
  elif scipy.stats.ncx2.ppf(1-0.0455/2.,2,fResp['V1']*sigmaVals['V1max']*trigger.snr**2) < trigger.snr_v**2:
    sbvs['BestNR7'] = 0
  else:
    sbvs['BestNR7'] = sbvs['BestNR6']

  return sbvs

def calculate_contours(bankq=0,bankn=0,autoq=0,auton=0,chiq=0,chin=0,nullt=0,nullo=0):
  cont_vals = [5.5,6,6.5,7,8,9,10,11]
  num_vals = len(cont_vals)
  colors = ['y-','k-','y-','y-','y-','y-','y-','y-']

  snr_vals = pylab.arange(6,30,0.1)
  snr_high_vals = pylab.arange(30,500,1)
  snr_all = []
  for snr in snr_vals:
    snr_all.append(snr)
  for snr in snr_high_vals:
    snr_all.append(snr)

  snr_vals = pylab.asarray(snr_all)

  bank_conts = [[],[],[],[],[],[],[],[]]
  auto_conts = [[],[],[],[],[],[],[],[]]
  chi_conts = [[],[],[],[],[],[],[],[]]
  null_cont = []

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

def plot_contours(snr_vals,contours,colors):
  for i in range(len(contours)):
    plot_vals_x = []
    plot_vals_y = []
    for j in range(len(snr_vals)):
      if contours[i][j] > 1E-15:
        plot_vals_x.append(snr_vals[j])
        plot_vals_y.append(contours[i][j])
    pylab.plot(plot_vals_x,plot_vals_y,colors[i])

def readSegFiles(segFileLocation):
  times = {}
  for name,fileName in zip(["buffer","off","on"],["bufferSeg.txt","offSourceSeg.txt","onSourceSeg.txt"]):
    file = open(os.path.join(segFileLocation,fileName),"r")
    for line in file:
      line = line.split('\t')
      if line[0] == '0':
        times[name+"_start"] = line[1]
        times[name+"_end"] = line[2]
  return times
       
def makePaperPlots():
  pylab.rcParams.update({
#    "text.usetex": True,
    "text.verticalalignment": "center",
#    "lines.markersize": 12,
#    "lines.markeredgewidth": 2,
#    "lines.linewidth": 2.5,
    "font.size": 218,
    "axes.titlesize": 24,
    "axes.labelsize": 24,
    "xtick.labelsize": 18,
    "ytick.labelsize": 18,
    "legend.fontsize": 18,
    })

def get_ra_dec(grbFile):
  ext_trigs = grbsummary.load_external_triggers(grbFile)
  ra = ext_trigs[0].event_ra
  dec = ext_trigs[0].event_dec
  return ra,dec

def read_sigma_vals(sigmaFile):
  sigmaVals = {}
  file = open(sigmaFile,'r')
  ifos = ['H1','H2','L1','V1']
  iter = 0
  for line in file:
    line = line.replace('\n','')
    line = line.split(' ')
    sigmaVals[ifos[iter] + 'min'] = float(line[0])
    sigmaVals[ifos[iter] + 'max'] = float(line[1])
    iter += 1
  return sigmaVals
    

def get_det_response(ra,dec,triggerTime):
  f_plus = {}
  f_cross = {}
  for det in ['H1','H2','L1','T1','V1','G1']:
    f_plus[det], f_cross[det],crap1,crap2 = antenna.response( triggerTime, ra,dec, 0,0, 'degree', det )
  return f_plus,f_cross
