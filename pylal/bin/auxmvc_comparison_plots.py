#!/usr/bin/python
# Reed Essick (reed.essick@ligo.org), Young-Min Kim (young-min.kim@ligo.org)



__author__ = "Young-Min Kim, Reed Essick, Ruslan Vaulin"
__version__ = "$Revision$"[11:-2]
__date__ = "$Date$"[7:-2]

__prog__="auxmvc_comparison_plots.py"
__Id__ = "$Id$"



from optparse import *
import glob
import sys
import os
import matplotlib
matplotlib.use('Agg')
from mpl_toolkits.axes_grid import ImageGrid
import pylab
import pdb
import numpy
from pylal import auxmvc_utils 
from pylal import git_version
import bisect
import pickle
from pylal import InspiralUtils
import math

def CalculateFAPandEFF(total_data,classifier):
	"""
	Calculates fap, eff for the classifier and saves them into the corresponding columns of total_data.
	"""	
#	if not classifier == 'ovl':
	glitch_ranks = total_data[numpy.nonzero(total_ranked_data['glitch'] == 1.0)[0],:][classifier+'_rank']
	clean_ranks = total_data[numpy.nonzero(total_ranked_data['glitch'] == 0.0)[0],:][classifier+'_rank']
		
	glitch_ranks.sort()
	clean_ranks.sort()
#	else:
#		print "Can not compute FAP and Efficiency for OVL. Not an MVC classifier."
#		sys.exit(1)
	  			  
	for trigger in total_data:
		trigger[classifier +'_fap'] = auxmvc_utils.compute_FAP(clean_ranks, glitch_ranks, trigger[classifier+'_rank'])
		trigger[classifier +'_eff'] = auxmvc_utils.compute_Eff(clean_ranks, glitch_ranks, trigger[classifier+'_rank'])
		
	return total_data


def compute_combined_rank(total_data):
	"""
	Computes combined rank defined as logarithm of highest EFF/FAP from the classifiers.
	"""

	for trigger in total_data:
		# We need to add here OVL
		combined_mvc_rank = max([eff_over_fap(trigger, cls[0]) for cls in classifiers if cls[0]!='ovl'])
#		combined_mvc_rank = max(eff_over_fap(trigger,'mvsc'), eff_over_fap(trigger,'ann'), eff_over_fap(trigger,'svm'))
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
	This function reads MVCs(MVSC,ANN,SVM) ranked data and OVL data and combines them into total_data array.
	When combining, it is assumed that input files contain same triggers. 
	In particular that sorting triggers in time is enough for cross-identification. 
	"""
	#reading all data from the first classifier
	data = auxmvc_utils.ReadMVSCTriggers(classifiers[0][1])
	n_triggers = len(data)

	# defining variables and formats for total_data array
	variables = ['GPS','glitch'] + list(data.dtype.names[5:-1])
	for var in ['mvsc','ann','svm','ovl']:
		variables += [var+j for j in ['_rank','_fap','_eff']]
	variables += ['ovl_chan','ovl_fdt','combined_rank', 'combined_eff', 'combined_fap']
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
			if not cls[0] == 'ovl':
				ranks_data=ReadMVCRanks(cls[1],cls[0])
				#sort trigers first in glitch then in GPS time to ensure correct merging
				ranks_data=ranks_data[numpy.lexsort((ranks_data['GPS'], ranks_data['glitch']))]
			
				# populating this classifier rank column in total_data array
				total_data[cls[0]+'_rank']=ranks_data[cls[0]+'_rank']
			
	
	# Reading in information from OVL and populating corresponding columns in total_data
	if classifiers[-1][0] == 'ovl':
		# we first need to retrieve the OVL data. Because we load both glitches and cleans into the same structure, we simply iterate over all the files listed in classifiers[-1][1]
		for file in classifiers[-1][1]:
			ovl_raw = auxmvc_utils.LoadOVL(file)
                        # figure out whether this is glitch or clean data
                        GorC = 1.0 
			for item in file.split('.'):
				if item == 'CLEAN':
					GorC = 0.0
			print 'ngwtrg_vtd from ' + file + ': ' + repr(len(ovl_raw))
			# we want to iterate through all the glitches stored in total_data that are glitches
			for glitch_index  in numpy.nonzero(total_data['glitch'] == GorC)[0]:
				# we look for matching GW glitches in the OVL data using the GPS time
				ovl_rank = GetOVLRank(ovl_raw, total_data[glitch_index]['GPS'], deltat = 0.0015)
				if len(ovl_rank) > 0:
					ovl_rank = ovl_rank[0]
					total_data[glitch_index]['ovl_eff'] = ovl_rank[0]
					total_data[glitch_index]['ovl_fdt'] = ovl_rank[1]
					total_data[glitch_index]['ovl_rank'] = ovl_rank[2]
					total_data[glitch_index]['ovl_chan'] = ovl_rank[3]
#					print '\nGPS:  ' + repr(total_data[glitch_index]['GPS'])
#	                                print 'ovl_eff: ' + repr(total_data[glitch_index]['ovl_eff'])
#					print 'ovl_fdt: ' + repr(total_data[glitch_index]['ovl_fdt'])
#					print 'ovl_rank: ' + repr(total_data[glitch_index]['ovl_rank'])
#					print 'ovl_chan: ' + repr(total_data[glitch_index]['ovl_chan'])
				else:
        	                        total_data[glitch_index]['ovl_eff'] = 1.0
                	                total_data[glitch_index]['ovl_fdt'] = 1.0
                        	        total_data[glitch_index]['ovl_chan'] = -1 

			print 'ngwtrg_vtd from total_data: ' + repr(len(total_data[numpy.nonzero(total_data['ovl_rank'] > 0)[0]]))
	return total_data


def GetOVLRank(OVLOutput, gwtcent, deltat = 0):
  '''
  This is meant to return the OVL rank-data corresponding to a given central time (tcent) for a GW glitch. the OVL data returned is a list of the form:
    [ovl_eff, ovl_fdt, ovl_rank, ovl_chan]
  the output arguments correspond to internal OVL data as follows:
    ovl_eff  : c_eff
    ovl_fdt  : c_fdt
    ovl_rank : 2-(rank/len(stats)) ... where rank is the row index in stats
    ovl_chan : vchan
  '''
  ###
  # the current incarnation will find all of the OVL glitches that satisfy |tcent - gwtcent| < deltat, but will only return one of them. This function needs to be extended to decide between the multiple returns, perhaps using other informaion about the glitch.
  ###
  # define references for internal structure of OVLOuptut
  Dic = {'tcent':0, 'vconfig':1, 'vstats':2}
  vconfigDic = {'vchan':0, 'vthr':1, 'vwin':2}
  vstatsDic = {'livetime':0, '#gwtrg':1, 'dsec':2, 'csec':3, 'vact':4, 'vsig':5, 'c_livetime':6, 'c_ngwtrg':7, 'c_dsec':8, 'c_csec':9, 'c_vact':10, 'rank':11}
  # define search variables
  begint = gwtcent - deltat
  endt = gwtcent + deltat
  ovl_dat = []
  max_rank = 0
  for index in range(len(OVLOutput)):
    #find maximum rank (total number of vconfigs applied - 1)
    if max_rank < OVLOutput[index][Dic['vstats']][vstatsDic['rank']]:
      max_rank = OVLOutput[index][Dic['vstats']][vstatsDic['rank']]
    #check to see if the OVL data matches gwtcent
    if (OVLOutput[index][Dic['tcent']] > begint) and (OVLOutput[index][Dic['tcent']] < endt):
      # extract desired data from OVLOutput
      h = OVLOutput[index]
      ovl_eff = float(h[Dic['vstats']][vstatsDic['c_vact']]) / float(h[Dic['vstats']][vstatsDic['c_ngwtrg']])
      ovl_fdt = float(h[Dic['vstats']][vstatsDic['c_dsec']]) / float(h[Dic['vstats']][vstatsDic['c_livetime']])
      not_quite_rank = float(h[Dic['vstats']][vstatsDic['rank']])
      ovl_chan = h[Dic['vconfig']][vconfigDic['vchan']]
      ovl_dat += [ [ ovl_eff, ovl_fdt, not_quite_rank, ovl_chan ] ]
    else:
      pass
  for index in range(len(ovl_dat)):
    # convert not_quite_rank to ovl_rank as defined above
    not_quite_rank = ovl_dat[index][2]
    ovl_dat[index][2] = 2 - float(not_quite_rank)/max_rank
  # we now return the first entry in the list, BUT THIS SHOULD BE EXTENDED SO THAT WE PICK THE CORRECT ENTRY BASED ON OTHER INFORMATION ABOUT THE GLITCH
  if len(ovl_dat) > 1:
    print('found multiple glitches at MVC GWtcent = ' + repr(gwtcent))
    print'\tcorresponding to OVL parameters: '
    for ovl in ovl_dat:
      print '\t\t' + repr(ovl) 
    return [ovl_dat[0]]
  else:
    return ovl_dat


def cluster(data, rank='signif', cluster_window=1.0):
	"""
	Clustering performed with the sliding window cluster_window keeping trigger with the highest rank;
	data is array of glitches, it is assumed to be sorted by GPS time in ascending order.
	Clustering algorithm is borrowed from pylal.CoincInspiralUtils cluster method.
	"""

	# initialize some indices (could work with just one)
	# but with two it is easier to read
	this_index = 0
	next_index = 1
	glitches = data[numpy.nonzero(data['glitch'] == 1.0)[0],:]
	while next_index < len(glitches):

		# get the time for both indices
		thisTime = glitches[this_index]['GPS']
		nextTime = glitches[next_index]['GPS']

		# are the two coincs within the time-window?
		if nextTime-thisTime < cluster_window:

			# get the ranks
			this_rank = glitches[this_index][rank]
			next_rank = glitches[next_index][rank]

			# and remove the trigger which has the lower rank
			if (next_rank > this_rank):
				glitches = numpy.delete(glitches, this_index, 0)
			else:
				glitches = numpy.delete(glitches, next_index, 0)

		else:
			# the two triggers are NOT in the time-window
			# so must increase index
			this_index+=1
			next_index+=1
			
	data = numpy.concatenate((glitches, data[numpy.nonzero(data['glitch'] == 0.0)[0],:]))
	return data
    

usage= """Written to load and manipulate output from different classifiers."""



parser=OptionParser(usage=usage, version = git_version.verbose_msg)
parser.add_option("","--ovl-ranked-files", default=False, type="string", help="Provide the path for HVeto Results files which contain all trigger data and those ranks and globbing pattern")
parser.add_option("","--mvsc-ranked-files", default=False, type="string", help="Provide the path for MVSC *.dat files and globbing pattern")
parser.add_option("","--ann-ranked-files", default=False, type="string", help="Provide the path for ANN *.dat files and globbing pattern")
parser.add_option("","--svm-ranked-files", default=False, type="string", help="Provide the path for SVM *.dat files and globbing pattern")
parser.add_option("","--fap-threshold", default=0.1,type="float", help="False Alarm Probability which is adapted to veto")
parser.add_option("","--cluster",action="store_true", default=False, help="cluster glitch samples")
parser.add_option("","--cluster-window", default=1.0, type="float", help="clustering window in seconds, default is 1 second.")
parser.add_option("","--cluster-rank", default='signif', type="string", help="rank used in clustering, default rank is significance of trigger in DARM.")
parser.add_option("-P","--output-path",action="store",type="string",default="",  metavar="PATH", help="path where the figures would be stored")	  
parser.add_option("-O","--enable-output",action="store_true", default="True",  metavar="OUTPUT", help="enable the generation of the html and cache documents")	 	  
parser.add_option("-u","--user-tag",action="store",type="string", default=None,metavar=" USERTAG", help="a user tag for the output filenames" )
parser.add_option("", "--figure-resolution",action="store",type="int", default=50, help="dpi of the thumbnails (50 by default)")
parser.add_option("", "--html-for-cbcweb",action="store", default=False, metavar = "CVS DIRECTORY", help="publish the html "\
      "output in a format that can be directly published on the cbc webpage "\
      "or in CVS. This only works IF --enable-output is also specified. The "\
      "argument should be the cvs directory where the html file will be placed "\
      "Example: --html-for-cbcweb protected/projects/s5/yourprojectdir")
parser.add_option("","--verbose", action="store_true", default=False, help="print information" )
parser.add_option("","--write-combined-data",action="store_true", default="True", help="write combined data to disk")
parser.add_option("","--DQ-ROC", type='string', default=False, help='plots DQ flags on ROC plots. DQ-ROC can be either S4 or S6') 

(opts,args)=parser.parse_args()

try: os.mkdir(opts.output_path)
except: pass


### Making ranked data to use glitch ranks and GW snr for plotting 
ranked_data={}

classifiers=[]
if opts.ann_ranked_files:
	classifiers.append(['ann',glob.glob(opts.ann_ranked_files)])
	#classifiers.append(['ann',opts.ann_ranked_files.split(',')])
if opts.mvsc_ranked_files:
	classifiers.append(['mvsc',glob.glob(opts.mvsc_ranked_files)])
	#classifiers.append(['mvsc',opts.mvsc_ranked_files.split(',')])
if opts.svm_ranked_files:
	classifiers.append(['svm',glob.glob(opts.svm_ranked_files)])
	#classifiers.append(['svm',opts.svm_ranked_files.split(',')])
if opts.ovl_ranked_files:
	classifiers.append(['ovl',glob.glob(opts.ovl_ranked_files)])
	#classifiers.append(['ovl',opts.ovl_ranked_files.split(',')])

#mvc_types=BinToDec(''.join(map(str,mvc_read)))

if not classifiers:
	print "Errors!! No Input Files(*.dat with MVCs' ranks and/or *.pickle with HVeto's ranks)"
	sys.exit()

if opts.verbose:
	print "Reading and combining data..."

# Reading and combining data from all classifers(MVSC,ANN,SVM,OVL).
total_ranked_data = ReadDataFromClassifiers(classifiers)

### TESTING OVL LOADING FUNCTIONALITY
'''
ovl_raw = auxmvc_utils.LoadOVL(classifiers[-1][1][0])
# we pull 10 random GPS times from ovl_raw[] and the corresponding data from total_ranked_data[]
# this is printed so we can check by eye whether they match
import random
deltaT = 0.0015
for ind in range(0,10):
  rind = random.randint(0,len(ovl_raw)-1)
  print '\n ovl_raw'
  print ovl_raw[rind] 
  print '\n total_ranked_data'
#  print total_ranked_data[numpy.nonzero( abs(ovl_raw[rind][0] - total_ranked_data['GPS']) <= deltaT )[0]]
  print 'GPS = ' + repr(total_ranked_data[numpy.nonzero( abs(ovl_raw[rind][0] - total_ranked_data['GPS']) <= deltaT )[0]]['GPS'])
  print 'ovl_eff = ' + repr(total_ranked_data[numpy.nonzero( abs(ovl_raw[rind][0] - total_ranked_data['GPS']) <= deltaT )[0]]['ovl_eff'])
  print 'ovl_fdt = ' + repr(total_ranked_data[numpy.nonzero( abs(ovl_raw[rind][0] - total_ranked_data['GPS']) <= deltaT )[0]]['ovl_fdt'])
  print 'ovl_chan = ' + repr(total_ranked_data[numpy.nonzero( abs(ovl_raw[rind][0] - total_ranked_data['GPS']) <= deltaT )[0]]['ovl_chan'])
'''
# cluster glitch samples if --cluster option is given
if opts.cluster:
	if opts.verbose:
		print "Number of glitch samples before clustering: ", len(total_ranked_data[numpy.nonzero(total_ranked_data['glitch'] == 1.0)[0],:])
		
	total_ranked_data = cluster(total_ranked_data, rank=opts.cluster_rank, cluster_window=opts.cluster_window)
	
	if opts.verbose:
		print "Number of glitch samples after clustering: ", len(total_ranked_data[numpy.nonzero(total_ranked_data['glitch'] == 1.0)[0],:])




# Computing FAP and Efficiency for MVCs
for cls in classifiers:
#  if not cls[0] == 'ovl':
  total_ranked_data = CalculateFAPandEFF(total_ranked_data,cls[0])
	
# Computing combined rank		
total_ranked_data = compute_combined_rank(total_ranked_data)

# Computing FAP and Efficiency for combned rank
total_ranked_data = CalculateFAPandEFF(total_ranked_data,'combined')

# add combined  to the list of classifiers
classifiers.append(['combined'])

#splitting data into glitch and clean samples
glitches = total_ranked_data[numpy.nonzero(total_ranked_data['glitch'] == 1.0)[0],:]
cleans = total_ranked_data[numpy.nonzero(total_ranked_data['glitch'] == 0.0)[0],:]

if opts.verbose:
	print "Done."

if opts.verbose:
	print "Writing combined data into a file..."

if opts.write_combined_data:
	# write glitch samples into a file
	glitch_file=open(opts.user_tag+'_glitch_data.dat','w')
	glitch_file.write(' '.join(glitches.dtype.names)+'\n')

	for da in glitches:
		glitch_file.write(' '.join(map(str,da))+'\n')	
	
	glitch_file.close()

	# write clean samples into a file
	clean_file=open(opts.user_tag+'_clean_data.dat','w')
	clean_file.write(' '.join(cleans.dtype.names)+'\n')

	for da in cleans:
		clean_file.write(' '.join(map(str,da))+'\n')
	
	clean_file.close()

if opts.verbose:
	print "Done."

################   PLOTS   #############################################################

# Initializing the html output
InspiralUtils.message(opts, "Initialisation...")
opts = InspiralUtils.initialise(opts, __prog__, __version__)
fnameList = []
tagList = []
fig_num = 0
comments = ""

###########################################################################################3
if opts.verbose:
  print 'Generating Plots...'

colorDIC = {'ann':'b', 'mvsc':'g', 'svm':'r', 'ovl':'c', 'combined':'m'}
labelDIC = {'ann':'ANN', 'mvsc':'MVSC', 'svm':'SVM', 'ovl':'OVL', 'combined':'MVC$_{\mathrm{max}}$'}

matplotlib.rc('text', usetex=True)

faircoin = numpy.linspace(10**-6,10**0.5,100)

# Scatter of ovl_fdt vs ovl_fap
if opts.verbose:
  print 'Scatter of ovl_fdt vs ovl_fap'
# this should help quantify how well fap estimates the fractional dead time with our random sampling
fig_num += 1
fig = pylab.figure(fig_num)
pylab.plot(total_ranked_data['ovl_fap'], total_ranked_data['ovl_fdt'], 'o', markerfacecolor = 'none', markeredgecolor = 'blue')
pylab.plot(faircoin, faircoin, 'k--')
pylab.grid(True)
pylab.xlabel('False Alarm Probability')
pylab.ylabel('Fractional Dead Time')
pylab.xscale('log')
pylab.yscale('log')
pylab.xlim(10**-4.1, 10**0.1)
pylab.ylim(10**-4.1, 10**0.1)
name = '_fdt_vs_fap'
fname = InspiralUtils.set_figure_name(opts, name)
fname_thumb = InspiralUtils.savefig_pylal(filename=fname, doThumb=True, dpi_thumb=opts.figure_resolution)
fnameList.append(fname)
tagList.append(name)
pylab.close()


# ROC curves from total_ranked_data[0]
if opts.verbose:
  print 'ROC curves from total_ranked_data[0]'
fig_num += 1
fig = pylab.figure(fig_num)
for cls in classifiers:
  total_data=total_ranked_data[numpy.lexsort((total_ranked_data[cls[0]+'_fap'], total_ranked_data[cls[0]+'_eff']))]
  if cls[0] == 'ovl':
    fap = total_data[numpy.nonzero(total_data['ovl_fap'] < 1)[0],:]
    pylab.plot(fap['ovl_fap'], fap['ovl_eff'], label = labelDIC[cls[0]], color=colorDIC[cls[0]])#, linewidth = 2)
  else:
    pylab.plot(total_data[cls[0]+'_fap'],total_data[cls[0]+'_eff'],label=labelDIC[cls[0]], color=colorDIC[cls[0]])#, linewidth = 2)
pylab.plot(faircoin, faircoin, '--k')
pylab.xlabel('False Alarm Probability')
pylab.ylabel('Efficiency')
pylab.xscale('log')
pylab.yscale('log')
pylab.xlim(10**-5, 10**0)
pylab.ylim(10**-3.5, 10**0)
pylab.title('Fig. '+str(fig_num)+': ROC Curves from total\_ranked\_data[]')
l1 = pylab.legend(loc = 'lower right')

name = '_ROC_log-log'
fname = InspiralUtils.set_figure_name(opts, name)
fname_thumb = InspiralUtils.savefig_pylal(filename=fname, doThumb=True, dpi_thumb=opts.figure_resolution)
fnameList.append(fname)
tagList.append(name)

fig_num +=1
pylab.title('Fig. '+str(fig_num)+': ROC Curves from total\_ranked\_data[]')
pylab.legend(loc = 'lower right')
pylab.yscale('linear')
pylab.xscale('linear')
pylab.xlim(0, 1)
pylab.ylim(0, 1)
name = '_ROC_lin-lin'
fname = InspiralUtils.set_figure_name(opts, name)
fname_thumb = InspiralUtils.savefig_pylal(filename=fname, doThumb=True, dpi_thumb=opts.figure_resolution)
fnameList.append(fname)
tagList.append(name)

fig_num += 1

if opts.DQ_ROC == 'S4':
  cbcS4DQfap = [2355./98147, 7672./98147, 19179./98147]
  cbcS4DQeff = [2357./16205, 6710./16205, 9445./16205]

  dq2,   = pylab.plot(cbcS4DQfap[0], cbcS4DQeff[0], linestyle = 'None', markerfacecolor = 'None', markeredgecolor = 'k', marker = '^', markersize = 6, label = 'CBC DQ\,2')
  dq23,  = pylab.plot(cbcS4DQfap[1], cbcS4DQeff[1], linestyle = 'None', markerfacecolor = 'None', markeredgecolor = 'k', marker = 's', markersize = 6, label = 'CBC DQ\,2+3')
  dq234, = pylab.plot(cbcS4DQfap[2], cbcS4DQeff[2], linestyle = 'None', markerfacecolor = 'None', markeredgecolor = 'k', marker = 'd', markersize = 6, label = 'CBC DQ\,2+3+4')
#  l2 = pylab.legend([dq2, dq23, dq234], ['CBC DQ\,2', 'CBC DQ\,2+3', 'CBC DQ\,2+3+4'], loc='upper right') # this removes l1 from the axes.
#  pylab.gca().add_artist(l1) # add l1 as a separate artist to the axes
  pylab.legend(loc = 'upper left')

if opts.DQ_ROC == 'S6':
  burstS6DQfap = [485./99869, 4028./99869, 4503./99869]
  burstS6DQeff = [547./2833, 721./2833, 748./2833]

  dq2,   = pylab.plot(burstS6DQfap[0], burstS6DQeff[0], linestyle = 'None', markerfacecolor = 'None', markeredgecolor = 'k', marker = '^', markersize = 6, label = 'Burst DQ\,2')
  dq23,  = pylab.plot(burstS6DQfap[1], burstS6DQeff[1], linestyle = 'None', markerfacecolor = 'None', markeredgecolor = 'k', marker = 's', markersize = 6, label = 'Burst DQ\,2+3')
  dq234, = pylab.plot(burstS6DQfap[2], burstS6DQeff[2], linestyle = 'None', markerfacecolor = 'None', markeredgecolor = 'k', marker = 'd', markersize = 6, label = 'Burst DQ\,2+3+4')
#  l2 = pylab.legend([dq2, dq23, dq234], ['Burst DQ\,2', 'Burst DQ\,2+3', 'Burst DQ\,2+3+4'], loc='upper right') # this removes l1 from the axes.
#  pylab.gca().add_artist(l1) # add l1 as a separate artist to the axes
  pylab.legend(loc = 'upper left')

pylab.title(opts.user_tag + '\nROC Curves')
pylab.legend(loc = 'upper left')
pylab.xscale('log')
pylab.xlim(10**-5, 10**0)
pylab.ylim(0,1)
	
name = '_ROC_lin-log'
fname = InspiralUtils.set_figure_name(opts, name)
fname_thumb = InspiralUtils.savefig_pylal(filename=fname, doThumb=True, dpi_thumb=opts.figure_resolution)
fnameList.append(fname)
tagList.append(name)
pylab.close()


# 5bit word plots
if opts.verbose:
  print '5bit word plots'
### we define a 5 bit word to be the concatenation of binary flags corresponding to whether or not a classifier removes a glitch at a set FAP
# 0: classifier did not remove glitch
# 1: classifier removed glitch

sigthr = 25

### this next bit is fragile, and only handles the special case of svm being absent. you should extend this to a more general case

if opts.svm_ranked_files:
  classify = [[['ann'], ['mvsc'], ['svm'], ['ovl'], ['combined']], [['ann'], ['mvsc'], ['svm'], ['ovl']], [['ann'], ['mvsc'], ['svm'], ['combined']], [['ann'], ['mvsc'], ['svm']], [['ovl'], ['combined']]]
else:
  classify = [[['ann'], ['mvsc'], ['ovl'], ['combined']], [['ann'], ['mvsc'], ['ovl']], [['ann'], ['mvsc'], ['combined']], [['ann'], ['mvsc']], [['ovl'], ['combined']]]

for FAPthr in [10**-3, 10**-2, 5*10**-2, 10**-1]:
  for clas in classify:
    #all glitches
    fbw = []
    for g in glitches:
      s = ""
      for cls in clas:
        if g[cls[0]+'_fap'] <= FAPthr:
          s += "1"
        else:
          s += "0"
      fbw += [BinToDec(s)]

    #glitches with signif <= sigthr
    g_rem = glitches[numpy.nonzero(glitches['signif'] <= sigthr)[0],:]
    fbwBelow = []
    for g in g_rem:
      s = ''
      for cls in clas:
        if g[cls[0]+'_fap'] <= FAPthr:
          s += '1'
        else:
          s += '0'
      fbwBelow += [BinToDec(s)]

    #glitches with signif >= sigthr
    g_rem = glitches[numpy.nonzero(glitches['signif'] >= sigthr)[0],:]
    fbwAbove = []
    for g in g_rem:
      s = ''
      for cls in clas:
        if g[cls[0]+'_fap'] <= FAPthr:
          s += '1'
        else:
          s += '0'
      fbwAbove += [BinToDec(s)]
    
    ### Add regular histograms 
    fig_num += 1
    fig = matplotlib.pyplot.figure()
    pylab.hist(fbw, bins = numpy.linspace(-0.5, 2**len(clas)-0.5, num = 2**len(clas) +1, endpoint=True), weights = numpy.ones(len(fbw))/len(fbw), histtype='step', log = True, label = 'all glitches')#, linewidth = 2)
    pylab.hist(fbwBelow, bins = numpy.linspace(-0.5, 2**len(clas)-0.5, num = 2**len(clas) +1, endpoint = True), weights = numpy.ones(len(fbwBelow))/len(fbwBelow), histtype = 'step', log = True, label = 'signif $\leq$ ' + str(sigthr))#,  linewidth = 2)
    (numA, bins, pathces) = pylab.hist(fbwAbove, bins = numpy.linspace(-0.5, 2**len(clas)-0.5, num = 2**len(clas) +1, endpoint = True), weights = numpy.ones(len(fbwAbove))/len(fbwAbove), histtype = 'step', log = True, label = 'signif $\geq$ ' + str(sigthr))#, linewidth = 2)
    pylab.ylabel('Fraction of Glitches')
    pylab.xlim(-0.5, 2**len(clas)-0.5)
    pylab.ylim(ymax = 1.0)
    pylab.title('Fig. '+str(fig_num)+': Histogram over Classifier Redundancy at FAP = '+str(FAPthr))
    pylab.legend(loc = 'upper center')
    pylab.grid(True, which = 'major')
    pylab.grid(True, which = 'minor')
  
    #convert decimals back to binaries and label the ticks accordingly
    tick_locs = numpy.linspace(0,2**len(clas)-1,num=2**len(clas))
    tick_labels = []
    for loc in tick_locs:
      s = ''
      l = range(len(clas) - 1)
      for ind in l[::-1]:
        if int(loc)/2**(ind+1) == 0:
          s += '0'
        else:
          s += '1'
          loc = int(loc)-2**(ind+1)
      if int(loc) == 0:
        s += '0'
      else:
        s += '1'
      tick_labels += [s]
    pylab.xticks(tick_locs, tick_labels)
    for label in fig.axes[0].xaxis.get_ticklabels():
      label.set_rotation(45)
  
    leg = "5 bit word is ordered in the following way:  "
    for cls in clas:
      leg += labelDIC[cls[0]]+' '
    pylab.figtext(0.5, 0.98, leg, ha = 'center', va = 'top')
  
    # print the fractions on top of the steps in the histograms
    for indn in range(len(numA)):
      loc = tick_locs[indn]
      if numA[indn] > 0:
        n = int(numA[indn]*10**3)/10.0
        pylab.text(loc,numA[indn],str(n)+'\,\%',ha='center',va='bottom')
 
    #adding to html page
    strclas = ''
    for cls in clas:
      strclas += cls[0] + '_'
    name = '_scatter_5bit-words_'+strclas+'_FAPthr_'+str(FAPthr)
    fname = InspiralUtils.set_figure_name(opts, name)
    fname_thumb = InspiralUtils.savefig_pylal(filename = fname, doThumb = True, dpi_thumb=opts.figure_resolution)
    fnameList.append(fname)
    tagList.append(name)
    pylab.close()

    ### add histograms over only those glitches that were removed
    fig_num += 1
    fig = matplotlib.pyplot.figure()
    fbwFound = []
    for fbwi in fbw:
      if fbwi > 0:
        fbwFound += [fbwi]
    fbwAboveFound = []
    for fbwA in fbwAbove:
      if fbwA > 0:
        fbwAboveFound += [fbwA]
    fbwBelowFound = []
    for fbwB in fbwBelow:
      if fbwB > 0:
        fbwBelowFound += [fbwB]
    pylab.hist(fbwFound, bins = numpy.linspace(0.5, 2**len(clas)-0.5, num = 2**len(clas), endpoint=True), weights = numpy.ones(len(fbwFound))/len(fbwFound), histtype='step', log = True, label = 'all glitches')#, linewidth = 2)
    pylab.hist(fbwBelowFound, bins = numpy.linspace(0.5, 2**len(clas)-0.5, num = 2**len(clas), endpoint = True), weights = numpy.ones(len(fbwBelowFound))/len(fbwBelowFound), histtype = 'step', log = True, label = 'signif $\leq$ ' + str(sigthr))#, linewidth = 2)
    (numA, bins, patches) = pylab.hist(fbwAboveFound, bins = numpy.linspace(0.5, 2**len(clas)-0.5, num = 2**len(clas), endpoint = True), weights = numpy.ones(len(fbwAboveFound))/len(fbwAboveFound), histtype = 'step', log = True, label = 'signif $\geq$ ' + str(sigthr))#, linewidth = 2)
    pylab.ylabel('Fraction of Glitches Found')
    pylab.xlim(0.5, 2**len(clas)-0.5)
    pylab.ylim(ymax = 1.0, ymin = 0.0)
    pylab.yscale('linear')
    pylab.title(opts.user_tag + '\nHistogram over Classifier Redundancy at FAP = '+str(FAPthr))
    pylab.legend(loc = 'upper right')
    pylab.grid(True, which = 'major')
    pylab.grid(True, which = 'minor')
  
    #convert decimals back to binaries and label the ticks accordingly
    tick_locs = numpy.linspace(1,2**len(clas)-1,num=2**len(clas)-1)
    tick_labels = []
    for loc in tick_locs:
      s = ''
      l = range(len(clas) - 1)
      for ind in l[::-1]:
        if int(loc)/2**(ind+1) == 0:
          s += '0'
        else:
          s += '1'
          loc = int(loc)-2**(ind+1)
      if int(loc) == 0:
        s += '0'
      else:
        s += '1'
      tick_labels += [s]
    pylab.xticks(tick_locs, tick_labels)
    for label in fig.axes[0].xaxis.get_ticklabels():
      label.set_rotation(45)
  
    # print the fractions on top of the steps in the histograms
    for indn in range(len(numA)):
      loc = tick_locs[indn]
      if numA[indn] > 0:
        n = int(numA[indn]*10**3)/10.0
        pylab.text(loc,numA[indn],str(n)+'\,\%',ha='center',va='bottom')
  
    leg = "Bit-word ordering:\n("
    for tmp in range(len(clas)-1):
      cls = clas[tmp]
      leg += labelDIC[cls[0]]+',\ '
    leg += labelDIC[clas[len(clas)-1][0]]+')'
    pylab.figtext(0.3, 0.8, leg, ha = 'center', va = 'center')
  
    #adding to html page
    strclas = ''
    for cls in clas:
      strclas += cls[0] + '_'
    name = '_scatter_5bit-words_FOUND_'+strclas+'_FAPthr_'+str(FAPthr)
    fname = InspiralUtils.set_figure_name(opts, name)
    fname_thumb = InspiralUtils.savefig_pylal(filename = fname, doThumb = True, dpi_thumb=opts.figure_resolution)
    fnameList.append(fname)
    tagList.append(name)
    pylab.close()


# Histogram over subsets of glitches from 5bit word plots
if opts.verbose:
  print 'Histogram over subsets of glitches from 5bit word plots'
# we pull out the glitches identified by certain classifiers but not by the others and histogram over their parameters:
#   GW signif
#   cls_rank (label the rank that corresponds to FAPthr)

sigthr = 25

### this next bit is fragile, and only handles the special case of svm being absent. you should extend this to a more general case

if opts.svm_ranked_files:
  classify = [[['ann'], ['mvsc'], ['svm']], [['ovl'], ['combined']]]
else:
  classify = [[['ann'], ['mvsc']], [['ovl'], ['combined']]]

for FAPthr in [10**-3, 10**-2, 5*10**-2, 10**-1]:
  for clas in classify:
    #all glitches
    fbw = []
    for g in glitches:
      s = ""
      for cls in clas:
        if g[cls[0]+'_fap'] <= FAPthr:
          s += "1"
        else:
          s += "0"
      fbw += [[g['GPS'],BinToDec(s)]]

    # separate and store GPS times for glitches that were identified by only one classifier
    only1 = [[] for i in range(len(clas))]
    for element in fbw:
      if (element[1] > 0) and (math.log(element[1],2)%1 == 0.0):
        only1[len(clas) - 1 - int(math.log(element[1],2))] += [element[0]]

    # histogram over parameters associated with GPS times in each list contained in only1[]
    for ind in range(len(only1)):
      fig_num += 1
      fig = matplotlib.pyplot.figure()
      # pull out all the information from glitches in only1[ind]
      GWsignifs = [g['signif'] for g in glitches for time in only1[ind] if g['GPS'] == time]
      if (len(GWsignifs) > 0):
        pylab.hist(GWsignifs, 100, histtype='step', label='all glitches')
        pylab.legend(loc='best')
      else:
        pylab.figtext(0.5, 0.5, 'no glitches identified by only '+labelDIC[clas[ind][0]], ha='center', va='center')
      pylab.ylabel('Number of Glitches')
      pylab.xlabel('Significance')
      pylab.title('Fig. '+str(fig_num)+': Histogram over Glitch Significance \nfor glitches found by ONLY '+labelDIC[clas[ind][0]]+' at FAP = '+str(FAPthr))

      #adding to html page
      name = '_hist_signif_glitches_removed_by_ONLY_'+clas[ind][0]+'_FAPthr_'+str(FAPthr)
      fname = InspiralUtils.set_figure_name(opts, name)
      fname_thumb = InspiralUtils.savefig_pylal(filename = fname, doThumb = True, dpi_thumb=opts.figure_resolution)
      fnameList.append(fname)
      tagList.append(name)
      pylab.close()

      for indd in range(len(only1)):
        fig_num += 1
        fig = matplotlib.pyplot.figure()
        rank = [g[clas[indd][0]+'_rank'] for g in glitches for t in only1[ind] if g['GPS'] == t]
        if (len(rank) > 0):
          pylab.hist(rank, 100, histtype='step', label='all glitches')
          hrank = [g[clas[indd][0]+'_rank'] for g in glitches for t in only1[ind] if g['GPS'] == t and g['signif'] >= sigthr]
          lrank = [g[clas[indd][0]+'_rank'] for g in glitches for t in only1[ind] if g['GPS'] == t and g['signif'] <= sigthr]
          if (len(lrank) > 0):
            pylab.hist(lrank, 100, histtype='step', label='signif $\leq$ '+str(sigthr))
          if (len(hrank) > 0):
            pylab.hist(hrank, 100, histtype='step', label='signif $\geq$ '+str(sigthr))
          pylab.legend(loc='best')
        else:
          pylab.figtext(0.5, 0.5, 'no glitces identified by only '+labelDIC[clas[ind][0]], ha='center', va='center')
        pylab.ylabel('Number of Glitches')
        pylab.xlabel(labelDIC[clas[indd][0]]+'\_rank')
        pylab.title('Fig. '+str(fig_num)+': Histogram over '+labelDIC[clas[indd][0]]+'\_rank\nfor glitches found by ONLY '+labelDIC[clas[ind][0]]+' at FAP = '+str(FAPthr))

        #adding to html page
        name = '_hist_'+clas[indd][0]+'_rank_glitches_removed_by_ONLY_'+clas[ind][0]+'_FAPthr_'+str(FAPthr)
        fname = InspiralUtils.set_figure_name(opts, name)
        fname_thumb = InspiralUtils.savefig_pylal(filename = fname, doThumb = True, dpi_thumb=opts.figure_resolution)
        fnameList.append(fname)
        tagList.append(name)
        pylab.close()


# statistics of glitches binned over cls_rank
if opts.verbose:
  print 'statistics of glitches binned over cls_rank'
### we bin glitches depending on cls_rank and then compute the statistics of glitches within those bins
### specifically, we calculate:
###  average significance
###  rms significance
###  the error in the estimate of the mean = (rms significance) / sqrt(num glitches in bin)

numBins = 100
maxminrank = {'mvsc':[0,1], 'ann':[0,1], 'svm':[0,1], 'ovl':[0,2],'combined':[0,1000]}

for cls in classifiers:
  for sigthr in [0, 25, 50]:
    #pull out the subset of glitches
    g_rem = glitches[numpy.nonzero(glitches['signif'] >= sigthr)[0],:]

    #define bin boundaries and list for signif values
    bins = numpy.linspace(maxminrank[cls[0]][0], maxminrank[cls[0]][1], num = numBins+1, endpoint = True)
    buckets = [[]]*numBins
    #iterate over glitches
    for g in g_rem:
      # iterate over bin boundaries. These are ordered, so we know that the first boundary that exceeds g[cls[0]+'_rank'] will determine to which bin the glitch belongs
      for ind in range(numBins):
        if (g[cls[0]+'_rank'] >= bins[ind]):
          if (g[cls[0]+'_rank'] <= bins[ind+1]):
            buckets[ind] = buckets[ind] + [g['signif']]
          else:
            pass
        else:
          pass

    # we now iterate over buckets and compute the statistics for glitches in each bucket
    aves = [0]*numBins
    stddevs = [0]*numBins
    ave_errs = [0]*numBins
    for ind in range(numBins):
      num_elements = float(len(buckets[ind]))
      if num_elements > 0:
        aves[ind] = sum(buckets[ind])/num_elements
        stddevs[ind] = (sum([(i - aves[ind])**2 for i in buckets[ind]])/(num_elements-1))**0.5
        ave_errs[ind] = stddevs[ind]/(num_elements)**0.5

    #generate figures
    fig_num+=1
    for ind in range(len(aves)):
      pylab.plot([(bins[ind]+bins[ind+1])/2., (bins[ind]+bins[ind+1])/2.], [aves[ind]+stddevs[ind], aves[ind]-stddevs[ind]],'-r')
      pylab.plot([bins[ind], bins[ind+1]], [aves[ind]+stddevs[ind], aves[ind]+stddevs[ind]], '-r')
      pylab.plot([bins[ind], bins[ind+1]], [aves[ind]-stddevs[ind], aves[ind]-stddevs[ind]], '-r')
      pylab.plot([(bins[ind]+bins[ind+1])/2., (bins[ind]+bins[ind+1])/2.], [aves[ind]+ave_errs[ind], aves[ind]-ave_errs[ind]],'-b')
      pylab.plot([bins[ind], bins[ind+1]], [aves[ind]+ave_errs[ind], aves[ind]+ave_errs[ind]], '-b')
      pylab.plot([bins[ind], bins[ind+1]], [aves[ind]-ave_errs[ind], aves[ind]-ave_errs[ind]], '-b')
    pylab.plot([(bins[i] + bins[i+1])/2. for i in range(len(bins)-1)], aves, 'ob')
    pylab.xlabel(labelDIC[cls[0]]+'\_rank')
    pylab.ylabel('average GW significance')
    pylab.title('Fig. '+str(fig_num)+': Statistics for glitches with signif $\geq$ '+str(sigthr)+' binned over '+labelDIC[cls[0]]+'\_rank')
    pylab.ylim(ymin = 0)
    #adding to html page
    name = '_statistics_binned_by_'+cls[0]+'_rank_sigthr_'+str(sigthr)
    fname = InspiralUtils.set_figure_name(opts, name)
    fname_thumb = InspiralUtils.savefig_pylal(filename = fname, doThumb = True, dpi_thumb=opts.figure_resolution)
    fnameList.append(fname)
    tagList.append(name)
    pylab.close()

# histograms of glitches over Classifier's rank
if opts.verbose:
  print 'histograms of glitches over Classifier rank'
for cls in classifiers:
  fig_num += 1
  pylab.figure(fig_num)
  pylab.hold(True)
  # histogram using all the glitches
  pylab.hist(glitches[cls[0] + '_rank'], 100, histtype='step', label = 'all glitches')
  # histogram using only a subset of glitches with sufficiently large 'signif'
  for sigthr in [10, 15, 25, 50]:
    glitches_removed = glitches[numpy.nonzero(glitches['signif'] >= sigthr)[0],:]
    if cls[0] == 'ovl':
      pylab.hist(glitches_removed[cls[0]+'_rank'], 100, histtype = 'step', label = 'signif $\geq$ ' + repr(sigthr), log=True)
    else:
      pylab.hist(glitches_removed[cls[0] + '_rank'], 100, histtype = 'step', label = 'signif $\geq$ ' + repr(sigthr), log=True)
  pylab.title('Fig. '+str(fig_num)+': Histogram for Glitches Based on ' + labelDIC[cls[0]] + '\_rank')
  pylab.xlabel(labelDIC[cls[0]]+ '\_rank')
  pylab.ylabel('Number of Glitches')
  pylab.legend(loc = 'upper center')
  
  #adding to the html page
  name = '_hist_glitches_' + cls[0] + '_rank_LOG'
  fname = InspiralUtils.set_figure_name(opts, name)
  fname_thumb = InspiralUtils.savefig_pylal(filename = fname, doThumb = True, dpi_thumb=opts.figure_resolution)
  fnameList.append(fname)
  tagList.append(name)

  #adding to the html page
  pylab.yscale('linear')
  name = '_hist_glitches_' + cls[0] + '_rank_LINEAR'
  fname = InspiralUtils.set_figure_name(opts, name)
  fname_thumb = InspiralUtils.savefig_pylal(filename = fname, doThumb = True, dpi_thumb=opts.figure_resolution)
  fnameList.append(fname)
  tagList.append(name)
  pylab.close()

# histograms of clean samples over Classifier's rank
if opts.verbose:
  print 'histograms of clean samples over Classifier rank'
for cls in classifiers:
  fig_num +=1
  pylab.figure(fig_num)
  pylab.hold(True)
  #histogram using all clean samples
  pylab.hist(cleans[cls[0]+'_rank'], 100, histtype='step',label = 'all', log = True)
  pylab.title('Fig. '+str(fig_num)+': Histogram for Clean Samples Based on '+labelDIC[cls[0]]+'\_rank')
  pylab.xlabel(labelDIC[cls[0]]+'\_rank')
  pylab.ylabel('Number of Samples')

  #adding to html page
  name = '_hist_cleans'+cls[0]+'_rank-LOG'
  fname = InspiralUtils.set_figure_name(opts, name)
  fname_thumb = InspiralUtils.savefig_pylal(filename = fname, doThumb = True, dpi_thumb=opts.figure_resolution)
  fnameList.append(fname)
  tagList.append(name)

  #adding to html page
  pylab.yscale('linear')
  name = '_hist_cleans'+cls[0]+'_rank-LINEAR'
  fname = InspiralUtils.set_figure_name(opts, name)
  fname_thumb = InspiralUtils.savefig_pylal(filename = fname, doThumb = True, dpi_thumb=opts.figure_resolution)
  fnameList.append(fname)
  tagList.append(name)
  pylab.close()


# Overlay of histograms of glitches and histograms of cleans over Classifier's rank
if opts.verbose:
  print 'Overlay of histogras of glitches and histograms of cleans over Classifier rank'
for cls in classifiers:
  fig_num += 1
  pylab.figure(fig_num)
  pylab.hold(True)
  # histogram using all the cleans
  pylab.hist(cleans[cls[0]+'_rank'], 100, histtype='step', weights = numpy.ones(len(cleans[cls[0]+'_rank']))/len(cleans[cls[0]+'_rank']), label = 'all cleans', log = True)
  # histogram using all the glitches
  pylab.hist(glitches[cls[0] + '_rank'], 100, histtype='step', weights = numpy.ones(len(glitches[cls[0]+'_rank']))/len(glitches[cls[0]+'_rank']), label = 'all glitches')
  # histogram using only a subset of glitches with sufficiently large 'signif'
  for sigthr in [10, 15, 25, 50]:
    glitches_removed = glitches[numpy.nonzero(glitches['signif'] >= sigthr)[0],:]
    pylab.hist(glitches_removed[cls[0]+'_rank'], 100, histtype = 'step', weights = numpy.ones(len(glitches_removed[cls[0]+'_rank']))/len(glitches_removed[cls[0]+'_rank']), label = 'signif $\geq$ ' + repr(sigthr), log=True)
  pylab.title('Fig. '+str(fig_num)+': Histogram for Glitches Based on ' + labelDIC[cls[0]] + '\_rank')
  pylab.xlabel(labelDIC[cls[0]]+ '\_rank')
  pylab.ylabel('Fraction of Glitches')
  pylab.legend(loc = 'upper center')

  #adding to the html page
  name = '_hist_glitches_and_cleans_' + cls[0] + '_rank_LOG'
  fname = InspiralUtils.set_figure_name(opts, name)
  fname_thumb = InspiralUtils.savefig_pylal(filename = fname, doThumb = True, dpi_thumb=opts.figure_resolution)
  fnameList.append(fname)
  tagList.append(name)

  #adding to the html page
  pylab.yscale('linear')
  name = '_hist_glitches_' + cls[0] + '_rank_LINEAR'
  fname = InspiralUtils.set_figure_name(opts, name)
  fname_thumb = InspiralUtils.savefig_pylal(filename = fname, doThumb = True, dpi_thumb=opts.figure_resolution)
  fnameList.append(fname)
  tagList.append(name)
  pylab.close()
	
## Create Cumulative histograms of GW Significance for glitch triggers after vetoing at RankThr and FAPThr
FAPThr = opts.fap_threshold
eff_file=open(opts.user_tag+'_efficiency_at_fap'+str(FAPThr)+'.txt','w')
eff_file.write('Vetoed Results at FAP='+str(FAPThr)+'\n')

# histograms for SNR
if opts.verbose:
  print 'histograms for SNR'
fig_num += 1
pylab.figure(fig_num)
pylab.hist(glitches['signif'],400,histtype='step',cumulative=-1,label='before vetoing', color='k')#, linewidth = 2)
for cls in classifiers:
	rank_name = cls[0]+'_rank'
	glitches_vetoed = glitches[numpy.nonzero(glitches[cls[0]+'_fap'] > FAPThr)[0],:]
	pylab.hist(glitches_vetoed['signif'],400,histtype='step',cumulative=-1,label=labelDIC[cls[0]], color=colorDIC[cls[0]])#, linewidth = 2)
	efficiency = min(glitches_vetoed[cls[0]+'_eff'])
	RankThr = max(glitches_vetoed[cls[0]+'_rank'])
	eff_file.write(cls[0]+' Efficiency : '+str(efficiency)+', threshold rank :'+str(RankThr)+'\n')
#pylab.title(r'Fig. '+str(fig_num)+'\nCumulative histogram for Significance  of glitches after vetoing at FAP '+str(fapthr))
pylab.title(opts.user_tag + '\nCumulative Histogram of Glitch DARM Significance at FAP = '+str(FAPThr))
pylab.xlabel('Significance')
pylab.ylabel('Number of Glitches')
pylab.xscale('log')
pylab.yscale('log')
pylab.xlim(xmin=min(glitches['signif']), xmax=4*10**2)
pylab.legend()
# adding to html page
name = '_cumul_hist_signif_fap'+str(FAPThr)
fname = InspiralUtils.set_figure_name(opts, name)
fname_thumb = InspiralUtils.savefig_pylal(filename=fname, doThumb=True, dpi_thumb=opts.figure_resolution)
fnameList.append(fname)
tagList.append(name)
pylab.close()

eff_file.close()

# histograms for SNR with pre-set FAPthr
for fapthr in [0.05, 0.1, 0.2, 0.4]:
  fig_num += 1
  pylab.figure(fig_num)
  pylab.hist(glitches['signif'],400,histtype='step',cumulative=-1,label='before vetoing', color='k')#, linewidth = 2)
  for cls in classifiers:
    rank_name = cls[0]+'_rank'
    glitches_vetoed = glitches[numpy.nonzero(glitches[cls[0]+'_fap'] > fapthr)[0],:]
    pylab.hist(glitches_vetoed['signif'],400,histtype='step',cumulative=-1,label=labelDIC[cls[0]], color = colorDIC[cls[0]])#, linewidth = 2)  
#  pylab.title(r'Fig. '+str(fig_num)+'\nCumulative histogram for Significance  of glitches after vetoing at FAP '+str(fapthr))
  pylab.title(opts.user_tag + '\nCumulative Histogram of Glitch DARM Significance at FAP = '+str(fapthr))
  pylab.xlabel('Significance')
  pylab.ylabel('Number of Glitches')
  pylab.xscale('log')
  pylab.yscale('log')
  pylab.xlim(xmin=min(glitches['signif']), xmax=4*10**2)
  pylab.legend()
  # adding to html page
  name = '_cumul_hist_signif_fap'+str(fapthr)
  fname = InspiralUtils.set_figure_name(opts, name)
  fname_thumb = InspiralUtils.savefig_pylal(filename=fname, doThumb=True, dpi_thumb=opts.figure_resolution)
  fnameList.append(fname)
  tagList.append(name)
  pylab.close()

# scatter plots of FAR from one classifier vs. FAR from another
if opts.verbose:
  print 'scatter plots of FAR from one classifier vs. FAR from another'
# we iterate through all pairs of classifiers.
start = 0
for cls in classifiers:
  start +=1
  for ind in range(start, len(classifiers)):
    cls2 = classifiers[ind]
    for sigthr in [10, 15, 25, 50]:
      g_rem = glitches[numpy.nonzero(glitches['signif'] >= sigthr)[0],:]
      fig_num += 1
      fig = matplotlib.pyplot.figure(fig_num)
      fig.hold(True)
      grid = ImageGrid(fig, 111, nrows_ncols = (2,2), axes_pad = 0.25, add_all = True)

      # histogram over cls_fap
      grid[0].hist(glitches[cls[0]+'_fap'],bins=numpy.logspace(-5,0,num=100,endpoint=True),histtype='step',weights = numpy.ones(len(glitches[cls[0]+'_fap']))/len(glitches[cls[0]+'_fap']), log=True)
      grid[0].hist(g_rem[cls[0]+'_fap'],bins=numpy.logspace(-5,0,num=100,endpoint=True),histtype='step',weights = numpy.ones(len(g_rem[cls[0]+'_fap']))/len(g_rem[cls[0]+'_fap']),color='red',log=True)

      grid[0].set_xscale('log')
      grid[0].set_yscale('log')
      grid[0].set_ylim(10**-4, 10**-0.5)
      grid[0].set_ylabel('Fraction of Glitches')

      # histogram over cls2_fap
      grid[3].hist(glitches[cls2[0]+'_fap'],bins=numpy.logspace(-5,0,num=100,endpoint=True),histtype='step',weights = numpy.ones(len(glitches[cls2[0]+'_fap']))/len(glitches[cls2[0]+'_fap']), orientation='horizontal', log=True)
      grid[3].hist(g_rem[cls2[0]+'_fap'],bins=numpy.logspace(-5,0,num=100,endpoint=True),histtype='step',weights = numpy.ones(len(g_rem[cls2[0]+'_fap']))/len(g_rem[cls2[0]+'_fap']), orientation='horizontal',color = 'red',log=True)
      grid[3].set_xscale('log')
      grid[3].set_yscale('log')
      grid[3].set_xlim(10**-4, 10**-0.5)
      grid[3].set_xlabel('Fraction of Glitches')

      # scatter of cls_fap vs cls2_fap
      grid[2].plot(glitches[cls[0]+'_fap'], glitches[cls2[0]+'_fap'], linestyle = 'none', marker = '.')
      grid[2].plot(g_rem[cls[0]+'_fap'], g_rem[cls2[0]+'_fap'], linestyle = 'none', marker = 'o', markerfacecolor = 'none', markeredgecolor = 'red')
      grid[2].plot([10**-5, 10**0], [10**-5, 10**0], '-k', label = '_nolegend')
      grid[2].set_xlabel(labelDIC[cls[0]]+'\_fap')
      grid[2].set_ylabel(labelDIC[cls2[0]]+'\_fap')
      grid[2].set_xscale('log')
      grid[2].set_yscale('log')
      
      matplotlib.pyplot.figtext(0.5, 0.98, 'Fig. '+str(fig_num)+': Scatter of '+labelDIC[cls[0]]+'\_fap vs. '+labelDIC[cls2[0]]+'\_fap for Glitches with signif $\geq$ '+str(sigthr), ha = 'center', va = 'top')
      matplotlib.pyplot.figtext(0.75, 0.75, 'blue : all glitches\nred : signif $\geq$ '+str(sigthr), ha = 'center', va = 'center')
      matplotlib.pyplot.setp(grid[1], visible=False)

      #adding to html page
      name = '_scatter_' + cls[0] + '_fap_vs_' + cls2[0] + '_fap_sigthr_'+str(sigthr)
      fname = InspiralUtils.set_figure_name(opts, name)
      fname_thumb = InspiralUtils.savefig_pylal(filename = fname, doThumb = True, dpi_thumb=opts.figure_resolution)
      fnameList.append(fname)
      tagList.append(name)
      pylab.close()

# scatter plots of Eff/FAR from one classifier vs. EFF/FAR from another
if opts.verbose:
  print 'scatter plots of Eff/FAR from one classifier vs. EFF/FAR from another'
# we iterate through all pairs of classifiers.
start = 0
for cls in classifiers:
  start +=1
  for ind in range(start, len(classifiers)):
    cls2 = classifiers[ind]
    for sigthr in [10, 15, 25, 50]:
      g_rem = glitches[numpy.nonzero(glitches['signif'] >= sigthr)[0],:]
      g_ebf = [g[cls[0]+'_eff']/g[cls[0]+'_fap'] for g in g_rem]
      g_ebf2 = [g[cls2[0]+'_eff']/g[cls2[0]+'_fap'] for g in g_rem]
      g_ebfALL = [glitch[cls[0]+'_eff']/glitch[cls[0]+'_fap'] for glitch in glitches]
      g_ebfALL2 = [glitch[cls2[0]+'_eff']/glitch[cls2[0]+'_fap'] for glitch in glitches]

      fig_num += 1
      fig = matplotlib.pyplot.figure(fig_num)
      fig.hold(True)
      grid = ImageGrid(fig, 111, nrows_ncols = (2,2), axes_pad = 0.25, add_all = True)

      # histogram over cls_eff/fap
      grid[0].hist(g_ebfALL,bins=numpy.logspace(0,4,num=100,endpoint=True),histtype='step', log=True)
      grid[0].hist(g_ebf,bins=numpy.logspace(0,4,num=100,endpoint=True),histtype='step',color='red',log=True)
      grid[0].set_xscale('log')
      grid[0].set_yscale('log')
      grid[0].set_ylim(ymin=10**0, ymax = 10**2.5)
      grid[0].set_ylabel('Number of Glitches')

      # histogram over cls2_eff/fap
      grid[3].hist(g_ebfALL2,bins=numpy.logspace(0,4,num=100,endpoint=True),histtype='step', orientation='horizontal', log=True)
      grid[3].hist(g_ebf2,bins=numpy.logspace(0,4,num=100,endpoint=True),histtype='step', orientation='horizontal',color = 'red',log=True)
      grid[3].set_xscale('log')
      grid[3].set_yscale('log')
      grid[3].set_xlim(xmin=10**0, xmax=10**2.5)
      grid[3].set_xlabel('Number of Glitches')

      # scatter of cls_eff/fap vs cls2_eff/fap
      grid[2].plot(g_ebfALL, g_ebfALL2, linestyle = 'none', marker = '.')
      grid[2].plot(g_ebf, g_ebf2, linestyle = 'none', marker = 'o', markerfacecolor = 'none', markeredgecolor = 'red')
      grid[2].plot([10**-5, 10**5], [10**-5, 10**5], '-k', label = '_nolegend')
      grid[2].set_xlabel(labelDIC[cls[0]]+'\_eff/fap')
      grid[2].set_ylabel(labelDIC[cls2[0]]+'\_eff/fap')
      grid[2].set_xscale('log')
      grid[2].set_yscale('log')
      grid[2].set_xlim(10**0, 10**3)
      grid[2].set_ylim(10**0, 10**3)

      matplotlib.pyplot.figtext(0.5, 0.98, 'Fig. '+str(fig_num)+': Scatter of '+labelDIC[cls[0]]+'\_eff/fap vs. '+labelDIC[cls2[0]]+'\_eff/fap for Glitches with signif $\geq$ '+str(sigthr), ha = 'center', va = 'top')
      matplotlib.pyplot.figtext(0.75, 0.75, 'blue : all glitches\nred : signif $\geq$ '+str(sigthr), ha = 'center', va = 'center')
      matplotlib.pyplot.setp(grid[1], visible=False)

      #adding to html page
      name = '_scatter_' + cls[0] + '_effbyfap_vs_' + cls2[0] + '_effbyfap_sigthr_'+str(sigthr)
      fname = InspiralUtils.set_figure_name(opts, name)
      fname_thumb = InspiralUtils.savefig_pylal(filename = fname, doThumb = True, dpi_thumb=opts.figure_resolution)
      fnameList.append(fname)
      tagList.append(name)
      pylab.close()


# Scatter of eff/fap for two classifiers using cleans
if opts.verbose:
  print 'Scatter of eff/fap for two classifiers using cleans'
# we iterate over all pairs of classifiers
start = 0
for cls in classifiers:
  start +=1
  for ind in range(start, len(classifiers)):
    cls2 = classifiers[ind]
    c_ebfALL = [g[cls[0]+'_eff']/g[cls[0]+'_fap'] for g in cleans]
    c_ebfALL2 = [g[cls2[0]+'_eff']/g[cls2[0]+'_fap'] for g in cleans]

    fig_num += 1
    fig = matplotlib.pyplot.figure(fig_num)
    fig.hold(True)
    grid = ImageGrid(fig, 111, nrows_ncols = (2,2), axes_pad = 0.25, add_all = True)

    # histogram over cls_eff/fap
    grid[0].hist(c_ebfALL,bins=numpy.logspace(0,4,num=100,endpoint=True),weights = numpy.ones(len(c_ebfALL))/len(c_ebfALL), histtype='step', log=True)
    grid[0].set_xscale('log')
    grid[0].set_yscale('log')
    grid[0].set_ylim(ymin=10**-4, ymax = 10**0.0)
    grid[0].set_ylabel('Fraction of Glitches')

    # histogram over cls2_eff/fap
    grid[3].hist(c_ebfALL2,bins=numpy.logspace(0,4,num=100,endpoint=True), weights = numpy.ones(len(c_ebfALL2))/len(c_ebfALL2), histtype='step', orientation='horizontal', log=True)
    grid[3].set_xscale('log')
    grid[3].set_yscale('log')
    grid[3].set_xlim(xmin=10**-4, xmax=10**0)
    grid[3].set_xlabel('Fraction of Glitches')

    # scatter of cls_eff/fap vs cls2_eff/fap
    grid[2].plot(c_ebfALL, c_ebfALL2, linestyle = 'none', marker = 'o')
    grid[2].plot([10**-5, 10**5], [10**-5, 10**5], '-k', label = '_nolegend')
    grid[2].set_xlabel(labelDIC[cls[0]]+'\_eff/fap')
    grid[2].set_ylabel(labelDIC[cls2[0]]+'\_eff/fap')
    grid[2].set_xscale('log')
    grid[2].set_yscale('log')
    grid[2].set_xlim(10**0, 10**3)
    grid[2].set_ylim(10**0, 10**3)

    matplotlib.pyplot.figtext(0.5, 0.98, 'Fig. '+str(fig_num)+': Scatter of '+labelDIC[cls[0]]+'\_eff/fap vs. '+labelDIC[cls2[0]]+'\_eff/fap for Clean Samples', ha = 'center', va = 'top')
    matplotlib.pyplot.setp(grid[1], visible=False)

    #adding to html page
    name = '_scatter_' + cls[0] + '_effbyfap_vs_' + cls2[0] + '_effbyfap_cleans'
    fname = InspiralUtils.set_figure_name(opts, name)
    fname_thumb = InspiralUtils.savefig_pylal(filename = fname, doThumb = True, dpi_thumb=opts.figure_resolution)
    fnameList.append(fname)
    tagList.append(name)
    pylab.close()


## Create scattered plots of Significance vs  SNR for GW Triggers
if opts.verbose:
  print 'scatter plots of signif vs snr for GW trigs'
fig_num += 1
pylab.figure(fig_num)
pylab.plot(glitches['signif'],glitches['SNR'],'rx',label='glitches')
pylab.xlabel('Significance of GW Triggers')
pylab.ylabel('SNR of GW Triggers')
pylab.yscale('log')
pylab.xscale('log')
pylab.title(r"Fig. "+str(fig_num)+": significance vs SNR of GW trigger")
pylab.legend()
# adding to html page
name = 'Sig_vs_SNR'
fname = InspiralUtils.set_figure_name(opts, name)
fname_thumb = InspiralUtils.savefig_pylal(filename=fname, doThumb=True, dpi_thumb=opts.figure_resolution)
fnameList.append(fname)
tagList.append(name)
pylab.close()

## Create scattered plots of mvc ranks vs. SNR and mvc ranks vs. Significance of GWTriggers
if opts.verbose:
  print 'scatter plots of rank vs. SNR and rank vs Signif'
for cls in classifiers:
        ## mvc ranks vs. SNR
        #fig_num += 1
        #pylab.figure(fig_num)
        #pylab.plot(cleans[cls[0]+'_rank'],cleans['SNR'],'k+',label='clean samples')
        #pylab.plot(glitches[cls[0]+'_rank'],glitches['SNR'],'rx',label='glitches')
        #pylab.xlabel(cls[0]+' rank')
        #pylab.ylabel('SNR of GW Triggers')
        #pylab.yscale('log')
        #pylab.xscale('log')
        #pylab.title(r"Fig. "+str(fig_num)+": " + cls[0]+" rank vs SNR of GW trigger")
        #pylab.legend()
        # adding to html page
        #name = '_SNR_'+cls[0]+'rank'
        #fname = InspiralUtils.set_figure_name(opts, name)
        #fname_thumb = InspiralUtils.savefig_pylal(filename=fname, doThumb=True, dpi_thumb=opts.figure_resolution)
        #fnameList.append(fname)
        #tagList.append(name)
        #pylab.close()


        ## mvc ranks vs. Significance
        fig_num += 1
        pylab.figure(fig_num)
        if cls[0] == 'combined':
                off_scale_glitches = glitches[numpy.nonzero(glitches[cls[0]+'_rank'] == 1000.0)]
                other_glitches = glitches[numpy.nonzero(glitches[cls[0]+'_rank'] != 1000.0)]
                pylab.plot(other_glitches[cls[0]+'_rank'],other_glitches['signif'],'rx',label='glitches')
                pylab.plot(numpy.ones(len(off_scale_glitches)) + max(other_glitches[cls[0]+'_rank']),off_scale_glitches['signif'],'g+',label='off-scale glitches')
        else:
                pylab.plot(glitches[cls[0]+'_rank'],glitches['signif'],'rx',label='glitches')
        pylab.xlabel(cls[0]+' rank')
        pylab.ylabel('Significance of GW Triggers')
        pylab.yscale('log')
        #pylab.xscale('log')
        pylab.title(r"Fig. "+str(fig_num)+": " + cls[0]+' rank vs significance of GW trigger')
        pylab.legend()
        # adding to html page
        name = '_Sig_'+cls[0]+'rank'
        fname = InspiralUtils.set_figure_name(opts, name)
        fname_thumb = InspiralUtils.savefig_pylal(filename=fname, doThumb=True, dpi_thumb=opts.figure_resolution)
        fnameList.append(fname)
        tagList.append(name)
        pylab.close()

##############################################################################################################
### this next big will help to print the options nicely on the html page
opts_dict = vars(opts)

html_filename = InspiralUtils.write_html_output(opts, [" --"+key.replace('_','-')+"="+str(opts_dict[key]) for key in opts_dict.keys()], fnameList, tagList, comment=comments)
InspiralUtils.write_cache_output(opts, html_filename, fnameList)

if opts.html_for_cbcweb:
  html_filename_publish = InspiralUtils.wrifacte_html_output(opts, args, fnameList, tagList, cbcweb=True)

##############################################################################################################
if opts.verbose:
  print 'Done.'

