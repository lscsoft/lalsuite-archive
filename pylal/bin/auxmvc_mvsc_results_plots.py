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

parser=OptionParser(usage="Generates summary plots for KW aux-triggers processed by StatPatternRecognition", version = "Kari Hodge")
parser.add_option("","--histograms", action="store_true", default=False, help="use if you want to produce histograms for each dimension")
(opts,files)=parser.parse_args()

print files
data = auxmvc_utils.ReadMVSCTriggers(files)        
print data.shape
clean_data = data[numpy.nonzero(data['i']==0)[0],:]
glitch_data = data[numpy.nonzero(data['i']==1)[0],:]
print clean_data.shape, 'background (clean events)'
print glitch_data.shape, 'signal (glitches)'


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
	# get the number of background triggers with rank greater than or equal to a given rank
	number_of_false_alarms = len(clean_ranks_sorted) -	numpy.searchsorted(clean_ranks_sorted,rank)
	# get the number of signal (glitch) triggers with rank greater than or equal to a given rank
	number_of_true_alarms = len(glitch_ranks_sorted) -	numpy.searchsorted(glitch_ranks_sorted,rank)
	# calculate the total deadime if this given rank is used as the threshold
	FAP.append( number_of_false_alarms / float(len(clean_ranks_sorted)))
	# calculate the fraction of correctly flagged signals (glitches)
	DP.append(number_of_true_alarms / float(len(glitch_ranks_sorted)))

pylab.figure(1)
pylab.plot(FAP,DP, linewidth = 2.0)
pylab.hold(True)
x = numpy.arange(min(DP), max(DP) + (max(DP) - min(DP))/1000.0, (max(DP) - min(DP))/1000.0)
pylab.plot(x,x, linestyle="dashed", linewidth = 2.0)
pylab.xlabel('False Alarm Probability')
pylab.ylabel('Efficiency')
pylab.savefig('classic_ROC.png')	
pylab.close()

	
#print clean_data
if opts.histograms:
	for i,var in enumerate(variables):
		pylab.figure(i)
		#print clean_data[var]
		print var
		pylab.hist(clean_data[var],100)
		pylab.savefig("hist_twosided"+var+"_background") 
		pylab.figure(100+i)
		pylab.hist(glitch_data[var],bins=100)
		pylab.savefig("hist_twosided"+var+"_signal") 
