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

parser=OptionParser(usage="Generates summary plots for KW aux-triggers processed by StatPatternRecognition", version = "Kari Hodge")
parser.add_option("","--histograms", action="store_true", default=False, help="use if you want to produce histograms for each dimension")
parser.add_option("","--tag", help="filenames will be ROC_tag.png and efficiency_deadtime_tag.txt")
(opts,files)=parser.parse_args()

print files
data = auxmvc_utils.ReadMVSCTriggers(files)        
print data.shape
clean_data = data[numpy.nonzero(data['i']==0)[0],:]
glitch_data = data[numpy.nonzero(data['i']==1)[0],:]

all_ranks = numpy.concatenate((clean_data['Bagger'],glitch_data['Bagger']))
all_ranks_sorted = numpy.sort(all_ranks)
clean_ranks_sorted = numpy.sort(clean_data['Bagger'])
glitch_ranks_sorted = numpy.sort(glitch_data['Bagger'])
number_of_false_alarms=[]
number_of_true_alarms=[]
deadtime=[] #fraction
efficiency=[] #fraction

print len(all_ranks_sorted), "all ranks"
print len(glitch_ranks_sorted), "glitch ranks"
print len(clean_ranks_sorted), "clean ranks"


FAP = []
DP = []
# classic ROC curve
for i,rank in enumerate(clean_ranks_sorted):
	# get the number of clean triggers with rank greater than or equal to a given rank
	number_of_false_alarms = len(clean_ranks_sorted) -	numpy.searchsorted(clean_ranks_sorted,rank)
	# get the number of glitches with rank greater than or equal to a given rank
	number_of_true_alarms = len(glitch_ranks_sorted) -	numpy.searchsorted(glitch_ranks_sorted,rank)
	# calculate the total deadime if this given rank is used as the threshold
	FAP.append( number_of_false_alarms / float(len(clean_ranks_sorted)))
	# calculate the fraction of correctly flagged glitches
	DP.append(number_of_true_alarms / float(len(glitch_ranks_sorted)))

pylab.figure(1)
pylab.plot(FAP,DP, linewidth = 2.0)
pylab.hold(True)
x = numpy.arange(min(DP), max(DP) + (max(DP) - min(DP))/1000.0, (max(DP) - min(DP))/1000.0)
pylab.plot(x,x, linestyle="dashed", linewidth = 2.0)
pylab.xlabel('False Alarm Probability')
pylab.ylabel('Efficiency')
pylab.xlim([0,1])
pylab.ylim([0,1])
pylab.savefig('ROC_'+opts.tag+'.png')	
pylab.close()

print FAP 
#FAP is a list that is naturally sorted in reverse order (highest to lowest),
#we need to turn it into a regularly sorted list so that we can find the DP for
#fiducial FAPs
print len(FAP)
print FAP[len(FAP)-36-1]
print DP[len(FAP)-36-1]
FAP.sort()
edfile = open('efficiency_deadtime_'+opts.tag+'.txt','w')
for threshold in [.01,.05,.1]:
	tmpindex=bisect.bisect_left(FAP,threshold)
	edfile.write("deadtime: "+str(FAP[tmpindex])+" efficiency: "+str(DP[len(FAP)-tmpindex-1])+"\n")

if opts.histograms:
	for i,var in enumerate(variables):
		pylab.figure(i)
		print var
		pylab.hist(clean_data[var],100)
		pylab.savefig("hist_twosided"+var+"_cleantimes") 
		pylab.figure(100+i)
		pylab.hist(glitch_data[var],bins=100)
		pylab.savefig("hist_twosided"+var+"_glitches") 
