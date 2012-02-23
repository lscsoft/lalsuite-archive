#!/usr/bin/python
# Written to load and manipulate output from different classifiers.
# Reed Essick (reed.essick@ligo.org), Young-Min Kim (young-min.kim@ligo.org)
# last modified: 12/07/2011

from optparse import *
import glob
import sys
import os
import matplotlib
matplotlib.use('Agg')
import pylab
import pdb
import numpy
from pylal import auxmvc_utils 
from pylal import git_version
import bisect
import pickle
from pylal import InspiralUtils


def CalculateFAPandEFF(total_data,classifier):
	"""
	Calculates fap, eff for the classifier and saves them into the corresponding columns of total_data.
	"""	
	if not classifier == 'cveto':
		glitch_ranks = total_data[numpy.nonzero(total_ranked_data['glitch'] == 1.0)[0],:][classifier+'_rank']
		clean_ranks = total_data[numpy.nonzero(total_ranked_data['glitch'] == 0.0)[0],:][classifier+'_rank']
		
		glitch_ranks.sort()
		clean_ranks.sort()
	else:
		print "Can not compute FAP and Efficiency for CVeto. Not an MVC classifier."
		sys.exit(1)
	  			  
	for trigger in total_data:
		trigger[classifier +'_fap'] = auxmvc_utils.compute_FAP(clean_ranks, glitch_ranks, trigger[classifier+'_rank'])
		trigger[classifier +'_eff'] = auxmvc_utils.compute_Eff(clean_ranks, glitch_ranks, trigger[classifier+'_rank'])
		
	return total_data


def compute_combined_rank(total_data):
	"""
	Computes combined rank defined as logarithm of highest EFF/FAP from the classifiers.
	"""

	for trigger in total_data:
		# We need to add here CVeto
		combined_mvc_rank = max(eff_over_fap(trigger,'mvsc'), eff_over_fap(trigger,'ann'), eff_over_fap(trigger,'svm'))
		# set highest combined rank to 1000.0
		if combined_mvc_rank == numpy.inf: 
			trigger['combined_rank'] = 1000.0
		else:
			trigger['combined_rank'] = combined_mvc_rank
	
	return total_data
	
def eff_over_fap(trigger, classifier):
  """
  Calculate ratio of eff over fap which approximates likelihood ratio.
  """	
  
  if trigger[classifier+'_fap'] == 0.0:
    return numpy.inf
  else:
    return numpy.log10(trigger[classifier+'_eff'] / trigger[classifier+'_fap'])
	
  
def PrateToRank(ranks,Prate):
	### convert a certain positive rate(Prate) to a corresponding rank in rank data
	ranks_sorted=numpy.sort(ranks)
	PositiveRates=[]
	for i,rank in enumerate(ranks_sorted):
		number_of_positive = len(ranks_sorted) - numpy.searchsorted(ranks_sorted,rank)
		PositiveRates.append(number_of_positive/float(len(ranks_sorted)))
	for i,ra in enumerate(PositiveRates):
		if ra <= Prate:
			break
	RankThr = ranks_sorted[i]

	return RankThr


def vetoGlitchesUnderFAP(glitch_data, rank_name, Rankthr, FAPthr):
	## veto triggers at Rankthr corresponding to FAPthr and return remained triggers
	glitch_data_sorted = numpy.sort(glitch_data,order=[rank_name])
	total_number_of_glitches = len(glitch_data_sorted[rank_name])
	number_of_vetoed_glitches = numpy.searchsorted(glitch_data_sorted[rank_name],Rankthr)
	number_of_true_alarms = total_number_of_glitches - number_of_vetoed_glitches
	efficiency = number_of_true_alarms/float(total_number_of_glitches)
	glitch_data_vetoed = glitch_data_sorted[:number_of_vetoed_glitches]

	return glitch_data_vetoed, efficiency


def ReadMVCRanks(list_of_files, classifier):
	"""
	Reads ranks from the classifier ouput files. Can read MVSC, ANN or SVM *.dat files.   
	"""
	data = auxmvc_utils.ReadMVSCTriggers(list_of_files)
	if classifier == 'mvsc':
		ranker = 'Bagger'
		rank_name = 'mvsc_rank'
	elif classifier == 'ann':
		ranker = 'glitch-rank'
		rank_name = 'ann_rank'
	else:
		ranker = 'SVMRank'
		rank_name = 'svm_rank'

	variables=['GPS','glitch',rank_name]
	formats=['g8','i','g8']
	gps_ranks = numpy.empty(len(data),dtype={'names':variables,'formats':formats})
	gps_ranks['glitch']=data['i']
	gps_ranks[rank_name]=data[ranker]
	gps_ranks['GPS'] = data['GPS_s'] + data['GPS_ms']*10**(-3)

	return gps_ranks


def BinToDec(binary):
	decimal=0
	power=len(str(binary)) - 1
	for i in str(binary):
		decimal += int(i)*2**power
		power -= 1
	return decimal


def ReadDataFromClassifiers(classifiers):
	"""
	This function reads MVCs(MVSC,ANN,SVM) ranked data and CVeto data and combines them into total_data array.
	When combining, it is assumed that input files contain same triggers. 
	In particular that sorting triggers in time is enough for cross-identification. 
	"""
	#reading all data from the first classifier
	data = auxmvc_utils.ReadMVSCTriggers(classifiers[0][1])
	n_triggers = len(data)

	# defining variables and formats for total_data array
	variables = ['GPS','glitch'] + list(data.dtype.names[5:-1])
	for var in ['mvsc','ann','svm','cveto']:
		variables += [var+j for j in ['_rank','_fap','_eff']]
	variables += ['cveto_chan','combined_rank', 'combined_eff', 'combined_fap']
	formats = ['g8','i']+['g8' for a in range(len(variables)-2)]

	# creating total_data array
	total_data = numpy.zeros(n_triggers, dtype={'names':variables, 'formats':formats})

	# populating columns using data from the first classifier
	total_data['glitch']=data['i']
	total_data[classifiers[0][0]+'_rank']=data[data.dtype.names[-1]]

	for name in data.dtype.names[5:-1]:
		total_data[name] = data[name]

	
	total_data['GPS'] = data['GPS_s'] + data['GPS_ms']*10**(-3)

	# sort total data first by glitch and then by GPS time.
	# numpy.lexsort() does inderct sort and return array's indices
	total_data=total_data[numpy.lexsort((total_data['GPS'], total_data['glitch']))]
	
	# looping over remaning classifiers
	if classifiers[1:]:
		for cls in classifiers[1:]:
			ranks_data=[]
			if not cls[0] == 'cveto':
				ranks_data=ReadMVCRanks(cls[1],cls[0])
				#sort trigers first in glitch then in GPS time to ensure correct merging
				ranks_data=ranks_data[numpy.lexsort((ranks_data['GPS'], ranks_data['glitch']))]
			
				# populating this classifier rank column in total_data array
				total_data[cls[0]+'_rank']=ranks_data[cls[0]+'_rank']
			
	
	# Reading in information from CVeto and populating corresponding columns in total_data
	if classifiers[-1][0] == 'cveto':
		# we first need to retrieve the CVeto data
		cveto_raw = auxmvc_utils.loadCV(classifiers[-1][1]) 
		# we want to iterate through all the glitches stored in total_data
		for glitch_index  in numpy.nonzero(total_data['glitch'] == 1.0)[0]:
			# we look for matching GW glitches in the CVeto data using the GPS time
			cveto_rank = GetCVetoRank(cveto_raw, total_data[glitch_index]['GPS'])
			if len(cveto_rank) > 0:
				total_data[glitch_index]['cveto_eff'] = cveto_rank[0]
				total_data[glitch_index]['cveto_fap'] = cveto_rank[1]
				total_data[glitch_index]['cveto_rank'] = cveto_rank[2]
				total_data[glitch_index]['cveto_chan'] = cveto_rank[3]
			else:
				pass
	return total_data


def GetCVetoRank(CVetoOutput, gwtcent, deltat = 0):
  '''
  This is meant to return the CVeto rank-data corresponding to a given central time (tcent) for a GW glitch. the CVeto data returned is a list of the form:
    [cveto_eff, cveto_fap, cveto_rank, cveto_chan]
  the output arguments correspond to internal CVeto data as follows:
    cveto_eff  : c_eff
    cveto_fap  : c_fct
    cveto_rank : 2-(rank/len(stats)) ... where rank is the row index in stats
    cveto_chan : vchan
  '''
  ###
  # the current incarnation will find all of the CVeto glitches that satisfy |tcent - gwtcent| < deltat, but will only return one of them. This function needs to be extended to decide between the multiple returns, perhaps using other informaion about the glitch.
  ###
  # define references for internal structure of CVetoOuptut
  Dic = {'tcent':0, 'vconfig':1, 'vstats':2}
  vconfigDic = {'vchan':0, 'vthr':1, 'vwin':2}
  vstatsDic = {'livetime':0, '#gwtrg':1, 'dsec':2, 'csec':3, 'vact':4, 'vsig':5, 'c_livetime':6, 'c_ngwtrg':7, 'c_dsec':8, 'c_csec':9, 'c_vact':10, 'rank':11}
  # define search variables
  begint = gwtcent - deltat
  endt = gwtcent + deltat
  cveto_dat = []
  max_rank = 0
  for index in range(len(CVetoOutput)):
    #find maximum rank (total number of vconfigs applied - 1)
    if max_rank < CVetoOutput[index][Dic['vstats']][vstatsDic['rank']]:
      max_rank = CVetoOutput[index][Dic['vstats']][vstatsDic['rank']]
    #check to see if the CVeto data matches gwtcent
    if CVetoOutput[index][Dic['tcent']] < begint:
      pass
    elif CVetoOutput[index][Dic['tcent']] > endt:
      pass
    else:
      # extract desired data from CVetoOutput
      h = CVetoOutput[index]
      cveto_eff = float(h[Dic['vstats']][vstatsDic['c_vact']]) / float(h[Dic['vstats']][vstatsDic['c_ngwtrg']])
      cveto_fap = float(h[Dic['vstats']][vstatsDic['c_csec']]) / float(h[Dic['vstats']][vstatsDic['c_livetime']])
      not_quite_rank = float(h[Dic['vstats']][vstatsDic['rank']])
      cveto_chan = h[Dic['vconfig']][vconfigDic['vchan']]
      cveto_dat += [ [ cveto_eff, cveto_fap, not_quite_rank, cveto_chan ] ]
#  print(cveto_dat)
  for index in range(len(cveto_dat)):
    # convert not_quite_rank to cveto_rank as defined above
    not_quite_rank = cveto_dat[index][3]
    cveto_dat[index][3] = 2 - float(rank)/max_rank
  # we now return the first entry in the list, BUT THIS SHOULD BE EXTENDED SO THAT WE PICK THE CORRECT ENTRY BASED ON OTHER INFORMATION ABOUT THE GLITCH
  if len(cveto_dat) > 1:
    print('found multiple glitches at:' + repr(gwtcent))
  return cveto_dat[0]



parser=OptionParser(usage="Generates comparison plots", version = "auxmvc")
parser.add_option("","--cveto-ranked-files", default=False, type="string", help="Provide the path for HVeto Results files which contain all trigger data and those ranks and globbing pattern")
parser.add_option("","--mvsc-ranked-files", default=False, type="string", help="Provide the path for MVSC *.dat files and globbing pattern")
parser.add_option("","--ann-ranked-files", default=False, type="string", help="Provide the path for ANN *.dat files and globbing pattern")
parser.add_option("","--svm-ranked-files", default=False, type="string", help="Provide the path for SVM *.dat files and globbing pattern")
parser.add_option("","--tag", default="auxmvc_results", help="filenames will be ROC_tag.png and efficiency_deadtime_tag.txt")
parser.add_option("","--fap-threshold", default=0.1,type="float", help="False Alarm Probability which is adapted to veto")
parser.add_option("","--output-dir", default=".", help="directory where output files will be written to")
(opts,files)=parser.parse_args()

try: os.mkdir(opts.output_dir)
except: pass


### Making ranked data to use glitch ranks and GW snr for plotting 
ranked_data={}

classifiers=[]
if opts.ann_ranked_files:
	#classifiers.append(['ann',glob.glob(opts.ann_ranked_files)])
	classifiers.append(['ann',opts.ann_ranked_files.split(',')])
if opts.mvsc_ranked_files:
	#classifiers.append(['mvsc',glob.glob(opts.mvsc_ranked_files)])
	classifiers.append(['mvsc',opts.mvsc_ranked_files.split(',')])
if opts.svm_ranked_files:
	#classifiers.append(['svm',glob.glob(opts.svm_ranked_files)])
	classifiers.append(['svm',opts.svm_ranked_files.split(',')])
if opts.cveto_ranked_files:
	#classifiers.append(['cveto',glob.glob(opts.cveto_ranked_files)])
	classifiers.append(['cveto',opts.cveto_ranked_files.split(',')])

#mvc_types=BinToDec(''.join(map(str,mvc_read)))

if not classifiers:
	print "Errors!! No Input Files(*.dat with MVCs' ranks and/or *.pickle with HVeto's ranks)"
	sys.exit()

# Reading and combining data from all classifers(MVSC,ANN,SVM,CVeto).
total_ranked_data = ReadDataFromClassifiers(classifiers)

# Computing FAP and Efficiency for MVCs
for cls in classifiers:
  if not cls[0] == 'cveto':
    total_ranked_data = CalculateFAPandEFF(total_ranked_data,cls[0])
	
# Computing combined rank		
total_ranked_data = compute_combined_rank(total_ranked_data)

# Computing FAP and Efficiency for combned rank
total_ranked_data = CalculateFAPandEFF(total_ranked_data,'combined')


#splitting data into glitch and clean samples
glitches = total_ranked_data[numpy.nonzero(total_ranked_data['glitch'] == 1.0)[0],:]
cleans = total_ranked_data[numpy.nonzero(total_ranked_data['glitch'] == 0.0)[0],:]


# write all data into test_file
test_file=open(opts.tag+'_data.dat','w')
test_file.write(' '.join(total_ranked_data.dtype.names)+'\n')

for da in total_ranked_data:
	test_file.write(' '.join(map(str,da))+'\n')

## Create scattered plots of mvc ranks vs. SNR and mvc ranks vs. Significance of GWTriggers
for cls in classifiers:
	## mvc ranks vs. SNR
	pylab.figure(1)
	pylab.plot(glitches[cls[0]+'_rank'],glitches['SNR'],'r*',label='glitches')
	pylab.plot(cleans[cls[0]+'_rank'],cleans['SNR'],'bo',label='clean samples')
	pylab.xlabel(cls[0]+' rank')
	pylab.ylabel('SNR of GW Triggers')
	pylab.yscale('log')
	pylab.legend()
	pylab.savefig(opts.tag+'_SNR_'+cls[0]+'rank.png')
	pylab.close()
	## mvc ranks vs. Significance
	pylab.figure(2)
	pylab.plot(glitches[cls[0]+'_rank'],glitches['signif'],'r*',label='glitches')
	pylab.plot(cleans[cls[0]+'_rank'],cleans['signif'],'bo',label='clean samples')
	pylab.xlabel(cls[0]+' rank')
	pylab.ylabel('Significance of GW Triggers')
	pylab.yscale('log')
	pylab.legend()
	pylab.savefig(opts.tag+'_Significance_'+cls[0]+'rank.png')
	pylab.close()

## Create Cumulative histograms of GW SNR for glitch triggers after vetoing at RankThr and FAPThr
FAPThr = opts.fap_threshold
eff_file=open(opts.tag+'_efficiency_at_fap'+str(FAPThr)+'.txt','w')
eff_file.write('Vetoed Results at FAP='+str(FAPThr)+'\n')
# histograms for SNR
pylab.figure(1)
pylab.hist(glitches['SNR'],400,histtype='step',cumulative=-1,label='before vetoing')
pylab.title("Cumulative histogram for SNR of glitches after vetoing at FAP "+str(FAPThr))
pylab.xlabel('SNR')
pylab.ylabel('Number of Glitches')
pylab.xscale('log')
pylab.yscale('log')
# histograms for Significance
pylab.figure(2)
pylab.hist(glitches['signif'],400,histtype='step',cumulative=-1,label='before vetoing')
pylab.title("Cumulative histogram for Significance of glitches after vetoing at FAP "+str(FAPThr))
pylab.xlabel('Significance')
pylab.ylabel('Number of Glitches')
pylab.xscale('log')
pylab.yscale('log')
glitches = total_ranked_data[numpy.nonzero(total_ranked_data['glitch'] == 1.0)[0],:]
for cls in classifiers:
	rank_name = cls[0]+'_rank'
	#RankThr = PrateToRank(cleans[rank_name],FAPThr)
	#glitches_vetoed, efficiency = vetoGlitchesUnderFAP(glitches,rank_name,RankThr,FAPThr)
	glitches_vetoed = glitches[numpy.nonzero(glitches[cls[0]+'_fap'] > FAPThr)[0],:]
	efficiency = min(glitches_vetoed[cls[0]+'_eff'])
	RankThr = max(glitches_vetoed[cls[0]+'_rank'])
	eff_file.write(cls[0]+' Efficiency : '+str(efficiency)+', threshold rank :'+str(RankThr)+'\n')
	# histograms for SNR
	pylab.figure(1)
	pylab.hist(glitches_vetoed['SNR'],400,histtype='step',cumulative=-1,label=cls[0])
	# histograms for Significance
	pylab.figure(2)
	pylab.hist(glitches_vetoed['signif'],400,histtype='step',cumulative=-1,label=cls[0])
pylab.figure(1)
pylab.legend()
pylab.savefig(opts.tag+'_cumul_hist_snr_fap'+str(FAPThr)+'.png')
pylab.close()
pylab.figure(2)
pylab.legend()
pylab.savefig(opts.tag+'_cumul_hist_signif_fap'+str(FAPThr)+'.png')
pylab.close()
