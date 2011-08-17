#!/bin/python

from optparse import *
import matplotlib
matplotlib.use('Agg')
import pylab
import pickle

parser=OptionParser(usage="supply several pickle files with false alarm percentages v. true alarm percentages, this code will plot them all on the same image", version="Kari Hodge")
parser.add_option("","--tag", help="filenames will be combined_ROC_tag.png")
parser.add_option("","--output-dir", help="directory where output files will be written to")

(opts,picklefiles)=parser.parse_args()

labels=[]

matplotlib.rcParams.update({
       'axes.color_cycle':['b','g','y','r','m','c',[(0.7,0.5,0)],[(0.7,0.,0.5)],'k'],
       #'figure.subplot.right': 0.6,
       #'figure.subplot.top': 0.6,
       'legend.fontsize': 8
})

for file in picklefiles:
	data = pickle.load(open(file))
	fig = pylab.figure(1)
	ax = pylab.axes()
	print type(ax)
	thislabel=file.split('.')[0].split('/')[1]
	print thislabel
	dqlabel=thislabel.split('_S4')[0]
	typelabel='S4'+thislabel.split('_S4')[1]
	print dqlabel, typelabel
	labels.append(dqlabel)
	pylab.loglog(data[0],data[1])
	pylab.xlabel('False Alarm Probability')
	pylab.ylabel('Efficiency')
	pylab.xlim([0,2])
	pylab.ylim([0,2])
	pylab.hold(True)
	pylab.legend(labels,loc=4)

pylab.text(0.005,1.0,typelabel,horizontalalignment='center')
pylab.savefig(opts.output_dir+'/combined_ROC'+opts.tag+'.png')

