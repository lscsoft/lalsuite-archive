#!/usr/bin/python
from optparse import *
import glob
import sys
import matplotlib
matplotlib.use('Agg')
import pylab
import math
import numpy

parser=OptionParser(usage="usage", version = "kari")
parser.add_option("","--files", default="H1L1*evaluation*.dat", help="your evaluation/testing files, warning: do not use training files here")
parser.add_option("","--zerolag", default="H1L1*zerolag*.dat")
parser.add_option("","--tag", default="histogram", help="label with ifo combo, run tag")
parser.add_option("","--title", default="preliminary", help="this will go on the top of all the plots")
parser.add_option("","--histograms",action="store_true", default=False, help="use if you want to turn on histogram production")
parser.add_option("","--comparison",help="name of variable (look at choices in top row of .dat files) to compare MVSC to in ROC. optional.")
(opts,args)=parser.parse_args()

files = glob.glob(opts.files)
print "These are the .dat files you are plotting from:", files
print "Is this all of your evaluation files? It should be..."
zerolag = glob.glob(opts.zerolag)
tag=opts.tag

variables = []
tmp = None

for f in files:
	flines = open(f).readlines()
	variables = flines[0].split()
	formats = ['i','i']+['g8' for a in range(len(variables)-2)]
	if tmp != None:
		data = numpy.concatenate((data,numpy.loadtxt(f,skiprows=1, dtype={'names': variables,'formats':formats})),axis=0)
	else:
		tmp = numpy.loadtxt(f,skiprows=1, dtype={'names': variables,'formats':formats})
		data = tmp
print "variables:", variables
print "total events", data.shape
bk_data = data[numpy.nonzero(data['i']==0)[0],:]
sg_data = data[numpy.nonzero(data['i']==1)[0],:]
print "total background events (Class 0)", bk_data.shape
print "total signal events (Class 1)", sg_data.shape

for f in zerolag:
	flines = open(f).readlines()
	variables = flines[0].split()
	formats = ['i','i']+['g8' for a in range(len(variables)-2)]
	zl_data = numpy.loadtxt(f,skiprows=1, dtype={'names': variables,'formats':formats})
zl_data = zl_data[numpy.nonzero(zl_data['i']==0)[0],:]

if opts.histograms:
	# make loglog histograms
	for i,var in enumerate(variables):
		pylab.figure(i)
		#print var
		min_x=min(numpy.min(numpy.log10(bk_data[var])),numpy.min(numpy.log10(sg_data[var])),numpy.min(numpy.log10(zl_data[var])))
		max_x=max(numpy.max(numpy.log10(bk_data[var])),numpy.max(numpy.log10(sg_data[var])),numpy.max(numpy.log10(zl_data[var])))
		#print min_x, max_x
		try:
			n,bin,patches = pylab.hist(numpy.log10(bk_data[var]),bins=100,range=(min_x,max_x),normed=True,log=True,label="background")
		except ValueError:
			continue
		pylab.setp(patches,'facecolor','k','alpha',.5)
		pylab.hold(1)
		n,bin,patches = pylab.hist(numpy.log10(sg_data[var]),bins=100,range=(min_x,max_x),normed=True,log=True,label="signal")
		pylab.setp(patches,'facecolor','r','alpha',.5)
		n,bin,patches = pylab.hist(numpy.log10(zl_data[var]),bins=100,range=(min_x,max_x),normed=True,log=True,label="zerolag")
		pylab.setp(patches,'facecolor','b','alpha',.5)
		pylab.xlabel("log10 "+var,fontsize=18)
		pylab.title(opts.title,fontsize=18)
		pylab.ylabel("density",fontsize=18)
		pylab.legend(loc ="upper right")
		pylab.savefig("hist_"+var+"_"+tag+"_logx")

	# make histograms
	for i,var in enumerate(variables):
		pylab.figure(i+100)
		print var
		min_x=min(numpy.min(bk_data[var]),numpy.min(sg_data[var]),numpy.min(zl_data[var]))
		max_x=max(numpy.max(bk_data[var]),numpy.max(sg_data[var]),numpy.max(zl_data[var]))
		print "min: ", min_x, " max: ", max_x
		try:
			n,bin,patches = pylab.hist(bk_data[var],bins=100,range=(min_x,max_x),normed=True,log=True,label="background")
		except ValueError:
			continue
		pylab.setp(patches,'facecolor','k','alpha',.5)
		pylab.hold(1)
		n,bin,patches = pylab.hist(sg_data[var],bins=100,range=(min_x,max_x),normed=True,log=True,label="signal")
		pylab.setp(patches,'facecolor','r','alpha',.5)
		n,bin,patches = pylab.hist(zl_data[var],bins=100,range=(min_x,max_x),normed=True,log=True,label="zerolag")
		pylab.setp(patches,'facecolor','b','alpha',.5)
		pylab.xlabel(var,fontsize=18)
		pylab.title(opts.title,fontsize=18)
		pylab.ylabel("density",fontsize=18)
		pylab.legend(loc ="upper right")
		pylab.savefig("hist_"+var+"_"+tag)

#prepare to make ROC plot
all_ranks = numpy.concatenate((bk_data['Bagger'],sg_data['Bagger']))
all_ranks_sorted = numpy.sort(all_ranks)
background_ranks_sorted = numpy.sort(bk_data['Bagger'])
signal_ranks_sorted = numpy.sort(sg_data['Bagger'])
number_of_false_alarms=[]
number_of_true_alarms=[]
deadtime=[] #fraction
efficiency=[] #fraction

print len(all_ranks_sorted), "all ranks"
print len(signal_ranks_sorted), "signal ranks"
print len(background_ranks_sorted), "background ranks"

FAP = [] # false alarm percentage
TAP = [] # true alarm percentage
# classic ROC curve
for i,rank in enumerate(background_ranks_sorted):
	# get the number of background triggers with rank greater than or equal to a given rank
	number_of_false_alarms = len(background_ranks_sorted) - numpy.searchsorted(background_ranks_sorted,rank)
	# get the number of signals with rank greater than or equal to a given rank
	number_of_true_alarms = len(signal_ranks_sorted) - numpy.searchsorted(signal_ranks_sorted,rank)
	# calculate the total deadime if this given rank is used as the threshold
	FAP.append( number_of_false_alarms / float(len(background_ranks_sorted)))
	# calculate the fraction of correctly flagged signals
	TAP.append(number_of_true_alarms / float(len(signal_ranks_sorted)))

pylab.figure(1000)
pylab.semilogx(FAP,TAP, linewidth = 2.0,label='MVSC')
pylab.hold(True)
#x = numpy.arange(0,1,.00001)
#pylab.semilogy(x,x, linestyle="dashed", linewidth = 2.0)
pylab.xlabel('False Alarm Fraction',fontsize=18)
pylab.ylabel('Efficiency',fontsize=18)
pylab.xlim([0,1])
pylab.ylim([0,1])

if opts.comparison:
	all_ranks = numpy.concatenate((bk_data[opts.comparison],sg_data[opts.comparison]))
	all_ranks_sorted = numpy.sort(all_ranks)
	background_ranks_sorted = numpy.sort(bk_data[opts.comparison])
	signal_ranks_sorted = numpy.sort(sg_data[opts.comparison])
	number_of_false_alarms=[]
	number_of_true_alarms=[]
	deadtime=[] #fraction
	efficiency=[] #fraction
	FAP = [] # false alarm percentage
	TAP = [] # true alarm percentage
	# classic ROC curve
	for i,rank in enumerate(background_ranks_sorted):
		# get the number of background triggers with rank greater than or equal to a given rank
		number_of_false_alarms = len(background_ranks_sorted) - numpy.searchsorted(background_ranks_sorted,rank)
		# get the number of signals with rank greater than or equal to a given rank
		number_of_true_alarms = len(signal_ranks_sorted) - numpy.searchsorted(signal_ranks_sorted,rank)
		# calculate the total deadime if this given rank is used as the threshold
		FAP.append( number_of_false_alarms / float(len(background_ranks_sorted)))
		# calculate the fraction of correctly flagged signals
		TAP.append(number_of_true_alarms / float(len(signal_ranks_sorted)))
	pylab.semilogx(FAP,TAP, linewidth = 2.0,label=opts.comparison)
pylab.legend(loc='lower right')
pylab.title(opts.title,fontsize=18)
pylab.savefig('ROC_'+opts.tag+'.png')
pylab.close()
