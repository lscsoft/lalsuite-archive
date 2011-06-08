#!/usr/bin/python
from optparse import *
import glob
import sys
import matplotlib
matplotlib.use('Agg')
import pylab
import pdb
import numpy
from pylal import auxmvc_utils 
import bisect

parser=OptionParser(usage="""
Generates summary plots for events processed by StatPatternRecognition
A single SprOutputWriterApp job produces a single .dat file
This code accepts multiple .dat files so that you can aggregate results for multiple sets of events (i.e. round robin evaluation sets)
instructions:
mvsc_ROC_and_histograms.py --tag test [--histograms] *.dat
""", version = "Kari Hodge")
parser.add_option("","--histograms", action="store_true", default=False, help="use if you want to produce histograms for each dimension")
parser.add_option("","--tag", help="filenames will be ROC_tag.png and efficiency_tag.txt")
(opts,files)=parser.parse_args()

tmp=None
print files
for f in files:
	flines = open(f).readlines()
	variables = flines[0].split()
	formats = ['i','i']+['g8' for a in range(len(variables)-2)]
	if tmp != None:
		data = numpy.concatenate((data,numpy.loadtxt(f,skiprows=1, dtype={'names': variables,'formats':formats})),axis=0)
	else:
		tmp = numpy.loadtxt(f,skiprows=1, dtype={'names': variables,'formats':formats})
		data = tmp
print variables
background_data = data[numpy.nonzero(data['i']==0)[0],:]
signal_data = data[numpy.nonzero(data['i']==1)[0],:]

all_ranks = numpy.concatenate((background_data['Bagger'],signal_data['Bagger']))
all_ranks_sorted = numpy.sort(all_ranks)
background_ranks_sorted = numpy.sort(background_data['Bagger'])
signal_ranks_sorted = numpy.sort(signal_data['Bagger'])
number_of_false_alarms=[]
number_of_true_alarms=[]
deadtime=[] #fraction
efficiency=[] #fraction

print len(all_ranks_sorted), "all ranks"
print len(signal_ranks_sorted), "signal ranks"
print len(background_ranks_sorted), "background ranks"


FAP = []
DP = []
# classic ROC curve
for i,rank in enumerate(background_ranks_sorted):
	# get the number of background triggers with rank greater than or equal to a given rank
	number_of_false_alarms = len(background_ranks_sorted) -	numpy.searchsorted(background_ranks_sorted,rank)
	# get the number of signals with rank greater than or equal to a given rank
	number_of_true_alarms = len(signal_ranks_sorted) -	numpy.searchsorted(signal_ranks_sorted,rank)
	# calculate the total deadime if this given rank is used as the threshold
	FAP.append( number_of_false_alarms / float(len(background_ranks_sorted)))
	# calculate the fraction of correctly flagged signals
	DP.append(number_of_true_alarms / float(len(signal_ranks_sorted)))

pylab.figure(1)
pylab.plot(FAP,DP, linewidth = 2.0)
pylab.hold(True)
x = numpy.arange(min(DP), max(DP) + (max(DP) - min(DP))/1000.0, (max(DP) - min(DP))/1000.0)
pylab.plot(x,x, linestyle="dashed", linewidth = 2.0)
pylab.xlabel('False Alarm Fraction')
pylab.ylabel('Efficiency')
pylab.xlim([0,1])
pylab.ylim([0,1])
pylab.savefig('ROC_'+opts.tag+'.png')	
pylab.close()

#FAP is a list that is naturally sorted in reverse order (highest to lowest),
#we need to turn it into a regularly sorted list so that we can find the DP for
#fiducial FAPs
FAP.sort()
edfile = open('efficiency_'+opts.tag+'.txt','w')
for threshold in [.01,.05,.1]:
	tmpindex=bisect.bisect_left(FAP,threshold)
	edfile.write("deadtime: "+str(FAP[tmpindex])+" efficiency: "+str(DP[len(FAP)-tmpindex-1])+"\n")

if opts.histograms:
	for i,var in enumerate(variables):
		pylab.figure(i)
		print var
		pylab.hist(background_data[var],100)
		pylab.savefig("hist_twosided"+var+"_background"+opts.tag) 
		pylab.figure(100+i)
		pylab.hist(signal_data[var],bins=100)
		pylab.savefig("hist_twosided"+var+"_signals"+opts.tag) 
