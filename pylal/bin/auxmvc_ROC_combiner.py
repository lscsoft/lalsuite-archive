#!/usr/bin/python

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

fig = pylab.figure(1)
ax = pylab.axes()
matplotlib.rcParams.update({'legend.fontsize': 8})
ax.set_color_cycle(['b','c','g','y',(0.9,0.5,0.2),'r','m',(0.5,0.,0.8),'k'])

for file in picklefiles:
	data = pickle.load(open(file))
	thislabel=file.split('.')[0].split('/')[1]
	dqlabel=thislabel.split('_')[1]
	typelabel=opts.tag
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

