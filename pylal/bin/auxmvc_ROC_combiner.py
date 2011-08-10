#!/bin/python

from optparse import *
import matplotlib
matplotlib.use('Agg')
import pylab
import pickle

parser=OptionParser(usage="supply several pickle files with false alarm percentages v. true alarm percentages, this code will plot them all on the same image", version="Kari Hodge")
(opts,picklefiles)=parser.parse_args()

labels=[]
for file in picklefiles:
	data = pickle.load(open(file))
	pylab.figure(1)
	thislabel=file.split('.')[0]
	labels.append(thislabel)
	pylab.loglog(data[0],data[1],label=thislabel)
	pylab.xlabel('False Alarm Probability')
	pylab.ylabel('Efficiency')
	pylab.xlim([0,1])
	pylab.ylim([0,1])
	pylab.hold(True)
pylab.legend(tuple(labels))
pylab.savefig('combined_ROC.png')
