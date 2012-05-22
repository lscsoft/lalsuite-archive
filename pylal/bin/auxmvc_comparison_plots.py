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

# Histogram over subsets of glitches from 5bit word plots
# we pull out the glitches identified by certain classifiers but not by the others and histogram over their parameters: 
#   GW signif
#   cls_rank (label the rank that corresponds to FAPthr)
#   else??

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
        pylab.figtext(0.5, 0.5, 'NO GLITCHES IDENTIFIED BY '+labelDIC[clas[ind][0]]+' ALONE', ha='center', va='center')
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
          pylab.figtext(0.5, 0.5, 'NO GLITCHES IDENTIFIED BY '+labelDIC[clas[ind][0]]+' ALONE', ha='center', va='center')
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

  
"""

# 3D surface
# we iterate through all pairs of classifiers.
from mpl_toolkits.mplot3d import axes3d, Axes3D
istart = 0
for cls in classifiers:
  start +=1
  for ind in range(start, len(classifiers)):
    cls2 = classifiers[ind]
    for sigthr in [10, 15, 25, 50]:
      g_rem = glitches[numpy.nonzero(glitches['signif'] >= sigthr)[0],:]
      bins=numpy.logspace(-5,0,num=10**2,endpoint=True)
      mdpts = [(bins[i]+bins[i+1])/2.0 for i in range(len(bins)-1)]
      (cls_fap, cls2_fap) = numpy.meshgrid(mdpts, mdpts)
      N = numpy.zeros((len(cls_fap), len(cls2_fap)))
      for g in g_rem:
        for nx in range(len(bins)-1):
          if (g[cls[0]+'_fap'] > bins[nx]) and (g[cls[0]+'_fap'] < bins[nx+1]):
            for ny in range(len(bins)-1):
              if (g[cls2[0]+'_fap'] > bins[ny]) and (g[cls2[0]+'_fap'] < bins[ny+1]):
                N[ny][nx] += 1

      fig_num += 1
      fig = pylab.figure(fig_num)

      ax = Axes3D(fig)
      ax.plot_surface(cls_fap, cls2_fap, N, rstride=20, cstride=20, alpha=0.3)
      cset = ax.contour(cls_fap, cls2_fap, N, zdir='z', offset=-100)
      cset = ax.contour(cls_fap, cls2_fap, N, zdir='x', offset=-40)
      cset = ax.contour(cls_fap, cls2_fap, N, zdir='y', offset=40)

      ax.set_xlabel(cls[0]+'\_fap')
      ax.set_xlim(10**-5, 10**0)
      ax.set_ylabel(cls2[0]+'\_fap')
      ax.set_ylim(10**-5, 10**0)
      ax.set_zlabel('Count')
#      ax.set_zlim(-10, 100)

      ax.set_xscale('log')
      ax.set_yscale('log')

      matplotlib.pyplot.figtext(0.5, 0.98, 'Fig. '+str(fig_num)+': Contour Plot of '+labelDIC[cls[0]]+'\_fap vs. '+labelDIC[cls2[0]]+'\_fap for Glitches with signif $\geq$ '+str(sigthr), ha = 'center', va = 'top')

      #adding to html page
      name = '_surface_' + cls[0] + '_fap_vs_' + cls2[0] + '_fap_sigthr_'+str(sigthr)
      fname = InspiralUtils.set_figure_name(opts, name)
      fname_thumb = InspiralUtils.savefig_pylal(filename = fname, doThumb = True, dpi_thumb=opts.figure_resolution)
      fnameList.append(fname)
      tagList.append(name)
      pylab.close()

# 2D Contour Plots of fap vs fap
# we iterate through all pairs of classifiers.
start = 0
for cls in classifiers:
  start +=1
  for ind in range(start, len(classifiers)):
    cls2 = classifiers[ind]
    for sigthr in [10, 15, 25, 50]:
      g_rem = glitches[numpy.nonzero(glitches['signif'] >= sigthr)[0],:]
      bins=numpy.logspace(-5,0,num=10**2,endpoint=True)
      mdpts = [(bins[i]+bins[i+1])/2.0 for i in range(len(bins)-1)]
      (cls_fap, cls2_fap) = numpy.meshgrid(mdpts, mdpts)
      N = numpy.zeros((len(cls_fap), len(cls2_fap)))
      for g in g_rem:
        for nx in range(len(bins)-1):
          if (g[cls[0]+'_fap'] > bins[nx]) and (g[cls[0]+'_fap'] < bins[nx+1]):
            for ny in range(len(bins)-1):
              if (g[cls2[0]+'_fap'] > bins[ny]) and (g[cls2[0]+'_fap'] < bins[ny+1]):
                N[ny][nx] += 1

      fig_num += 1
      pylab.figure(fig_num)
      CS = matplotlib.pyplot.contour(cls_fap, cls2_fap, N, 15)
#      matplotlib.pyplot.clabel(CS, fontsize=9, inline=1)
      pylab.xlabel(cls[0]+'\_fap')
      pylab.ylabel(cls2[0]+'\_fap')
      pylab.xscale('log')
      pylab.yscale('log')
      pylab.xlim(10**-5, 10**0)
      pylab.ylim(10**-5, 10**0)

      matplotlib.pyplot.figtext(0.5, 0.98, 'Fig. '+str(fig_num)+': Contour Plot of '+labelDIC[cls[0]]+'\_fap vs. '+labelDIC[cls2[0]]+'\_fap for Glitches with signif $\geq$ '+str(sigthr), ha = 'center', va = 'top')

      #adding to html page
      name = '_contour_' + cls[0] + '_fap_vs_' + cls2[0] + '_fap_sigthr_'+str(sigthr)
      fname = InspiralUtils.set_figure_name(opts, name)
      fname_thumb = InspiralUtils.savefig_pylal(filename = fname, doThumb = True, dpi_thumb=opts.figure_resolution)
      fnameList.append(fname)
      tagList.append(name)
      pylab.close()

# 2D Contour Plots of fap vs fap
# we iterate through all pairs of classifiers.
start = 0
for cls in classifiers:
  start +=1
  for ind in range(start, len(classifiers)):
    cls2 = classifiers[ind]
    for sigthr in [10, 15, 25, 50]:
      g_rem = glitches[numpy.nonzero(glitches['signif'] >= sigthr)[0],:]
      bins=numpy.logspace(-5,0,num=10**2,endpoint=True)
      mdpts = [(bins[i]+bins[i+1])/2.0 for i in range(len(bins)-1)]
      (cls_fap, cls2_fap) = numpy.meshgrid(mdpts, mdpts)
      N = numpy.zeros((len(cls_fap), len(cls2_fap)))
      for g in g_rem:
        for nx in range(len(bins)-1):
          if (g[cls[0]+'_fap'] > bins[nx]) and (g[cls[0]+'_fap'] < bins[nx+1]):
            for ny in range(len(bins)-1):
              if (g[cls2[0]+'_fap'] > bins[ny]) and (g[cls2[0]+'_fap'] < bins[ny+1]):
                N[ny][nx] += 1

      fig_num += 1
      pylab.figure(fig_num)
      pylab.pcolor(cls_fap, cls2_fap, N)
      pylab.colorbar()
      pylab.xlabel(cls[0]+'\_fap')
      pylab.ylabel(cls2[0]+'\_fap')
      pylab.xscale('log')
      pylab.yscale('log')
      pylab.xlim(10**-5, 10**0)
      pylab.ylim(10**-5, 10**0)

      matplotlib.pyplot.figtext(0.5, 0.98, 'Fig. '+str(fig_num)+': Contour Plot of '+labelDIC[cls[0]]+'\_fap vs. '+labelDIC[cls2[0]]+'\_fap for Glitches with signif $\geq$ '+str(sigthr), ha = 'center', va = 'top')

      #adding to html page
      name = '_color_' + cls[0] + '_fap_vs_' + cls2[0] + '_fap_sigthr_'+str(sigthr)
      fname = InspiralUtils.set_figure_name(opts, name)
      fname_thumb = InspiralUtils.savefig_pylal(filename = fname, doThumb = True, dpi_thumb=opts.figure_resolution)
      fnameList.append(fname)
      tagList.append(name)
      pylab.close()

# scatter plots of FAR from one classifier vs. FAR from another
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

"""

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

