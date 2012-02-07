#!/usr/bin/python
# Written to load and manipulate output from different classifiers.
# Reed Essick (reed.essick@ligo.org), Young-Min Kim (young-min.kim@ligo.org)
# last modified: 12/07/2011

#from optparse import *
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
from scipy.interpolate import interp1d


def FAPToRank(ranks_data,FAP):
	ranks_data_sorted=numpy.sort(ranks_data)
	FAPs=[]
	for i,rank in enumerate(ranks_data_sorted):
		number_of_false_alarms = len(clean_ranks_sorted) - numpy.searchsorted(clean_ranks_sorted,rank)
		FAPs.append(number_of_false_alarms/float(len(clean_ranks_sorted)))
	ranks_func = interp1d(FAPs,ranks_data_sorted)

	return ranks_func(FAP)


def vetoGlitchesUnderFAP(glitch_data, rank_name, Rankthr, FAPthr):
	glitch_data_sorted = numpy.sort(glitch_data,order(rank_name))
	total_number_of_glitches = len(glitch_data_sorted[rank_name])
	number_of_vetoed_glitches = numpy.searchsorted(glitch_data_sorted[rank_name],Rankthr)
	number_of_true_alarms = total_number_of_glitches - number_of_vetoed_glitches
	efficiency = number_of_true_alarms/float(total_number_of_glitches)
	glitch_data_vetoed = glitch_data_sorted[number_of_vetoed_glitches+1:]

	return glitch_data_vetoed, efficiency


def ReadMVSCRanks(list_of_files, classifier):
	data = auxmvc_utils.ReadMVSCTriggers(list_of_files)
	if classifier == 'mvsc':
		ranker = 'Bagger'
		rank_name = 'mvsc_rank'
	elif classifier == 'ann':
		ranker = 'glitch-rank'
		rank_name = 'ann_rank'
	else:
		rank = 'SVMRank'
		rank_name = 'svm_rank'

	variables=['GPS','i',rank_name]
	formats=['g8','i','g8']
	gps_ranks = numpy.empty(len(data),dtype={'names':variables,'formats':formats})
	gps_ranks['i']=data['i']
	gps_ranks[rank_name]=data[ranker]
	for i in range(len(data)):
		#numpy.put(gps_ranks['GPS_s'],[i],[data['GPS_s'][i]])
		#numpy.put(gps_ranks['GPS_ms'],[i],[data['GPS_ms'][i]])
		numpy.put(gps_ranks['GPS'],[i],[str(int(data['GPS_s'][i]))+'.'+str(int(data['GPS_ms'][i]))])

	return gps_ranks


def BinToDec(binary):
	decimal=0
	power=len(str(binary)) - 1
	for i in str(binary):
		decimal += int(i)*2**power
		power -= 1
	return decimal


def generateTotalRankedTriggers(classifiers):
	## This function is to read MVCs(MVSC,ANN,SVM)' ranked data and CVeto data and combine them into total_data array.
	data = auxmvc_utils.ReadMVSCTriggers(classifiers[0][1])
	n_triggers = len(data)
	variables = ['GPS','glitch'] + list(data.dtype.names[5:-1])
	for var in ['mvsc','ann','svm','cveto']:
		variables += [var+j for j in ['_rank','_fap','_eff']]
	variables += ['cveto_chan']

	formats = ['g8','i']+['g8' for a in range(len(variables)-2)]
	total_data = numpy.zeros(n_triggers, dtype={'names':variables, 'formats':formats})
	total_data['glitch']=data['i']
	total_data[classifiers[0][0]+'_rank']=data[data.dtype.names[-1]]

	for name in data.dtype.names[5:-1]:
		total_data[name] = data[name]

	for n for range(n_triggers):
		numpy.put(total_data['GPS'],[n],[str(int(data['GPS_s'][n]))+'.'+str(int(data['GPS_ms'][n]))])

	if not classifiers[1:]:
		for cls in classifiers[1:]:
			ranks_data=[]
			if not cls[0] == 'cveto':
				ranks_data=ReadMVSCRanks(cls[1],cls[0])
				total_data[cls[0]+'_rank']=ranks_data[cls[0]+'_rank']
				#below 2 lines, GPS times are used to assign mvc ranks into right places. It needs to check !.
				#for i, gps in enumerate(ranks_data['GPS']):
				#	numpy.put(total_data[cls[0]+'_rank'],[numpy.searchsorted(total_data['GPS'],gps)],[ranks_data[cls[0]+'_rank'][i]])

	# this next bit of code was added by R. Essick. He doesn't really know what he's doing, so please check this
	if classifiers[-1][0] == 'cveto':
		# we first need to retrieve the CVeto data
		cveto_raw = auxmvc_utils.loadCV(classifiers[-1][1]) 
		# we want to iterate through all the glitches stored in total_data
		for n in range(len(total_data['GPS'][:])):
			# we look for matching GW glitches in the CVeto data using the GPS time
			cveto_stuff = giveMeCVeto(cveto_raw, total_data['GPS'][n])
			if len(cveto_stuff) > 0:
				total_data['cveto_eff'][n] = cveto_stuff[0]
				total_data['cveto_fap'][n] = cveto_stuff[1]
				total_data['cveto_rank'][n] = cveto_stuf[2]
				total_data['cveto_chan'][n] = cveto_stuff[3]
			else:
				pass
	return total_data
	

def giveMeCVeto(CVetoOutput, gwtcent, deltat = 0)
  '''
  this is meant to return the CVeto data corresponding to a given central time (tcent) for a GW glitch. the CVeto data returned is a list of the form:
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
    elseif CVetoOutput[index][Dic['tcent']] > endt:
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
parser.add_option("","--output-dir", default=".", help="directory where output files will be written to")
(opts,files)=parser.parse_args()

try: os.mkdir(opts.output_dir)
except: pass

### Making ranked data to use glitch ranks and GW snr for plotting 
ranked_data={}

classifiers=[]
if opts.mvsc_ranked_files:
	classifiers.append(['mvsc',glob.glob(opts.mvsc_ranked_files)])
if opts.ann_ranked_files:
	classifiers.append(['ann',glob.glob(opts.ann_ranked_files)])
if opts.svm_ranked_files:
	classifiers.append(['svm',glob.glob(opts.ann_ranked_files)])
if opts.cveto_ranked_files:
	classifiers.append(['cveto',glob.glob(opts.cveto_ranked_files)])

#mvc_types=BinToDec(''.join(map(str,mvc_read)))

if not classifiers:
	print "Errors!! No Input Files(*.dat with MVCs' ranks and/or *.pickle with HVeto's ranks)"
	sys.exit()

# Construction of total triggers with all ranks for all classifers(MVSC,ANN,SVM,CVeto).

total_ranked_data = generateTotalRankedTriggers(classifiers)

if classifiers[-1][0] == 'cveto':
	total_ranked_data['cveto_rank'] = 'please change this part for cveto_rank' 
	total_ranked_data['cveto_fap'] = 'plese change this part for cveto_fap'
	total_ranked_data['cveto_eff'] = 'please change this part for cfeto_eff'
	total_ranked_data['cveto_chan'] = 'please change this part for cveto_chan'


# create a cumulative histogram of snr 
#n, bins, patches = pylab.hist(snr,bins, normed=1, histtype='step', cumulative=True)
#y=pylab.normpdf(bins,mu,sigma).comsum()
#y /= y[-1]
#plot(bins,y,'k--',linewidth=1.5)
