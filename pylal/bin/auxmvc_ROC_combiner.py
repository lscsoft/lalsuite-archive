#!/usr/bin/python

from optparse import *
import matplotlib
matplotlib.use('Agg')
import pylab
import pickle

parser=OptionParser(usage="""
Supply several pickle files with false alarm percentages v. true alarm percentages, this code will plot them all on the same image
Usage: 
auxmvc_ROC_combiner.py --tag S6_all_data --output-dir auxmvc_results_plots file1:label1 file2:label2 file3:label3 ...
""", version="Kari Hodge")
parser.add_option("","--tag", help="filenames will be combined_ROC_tag.png")
parser.add_option("","--title", help="make an informative title here")
parser.add_option("","--output-dir", help="directory where output files will be written to")

(opts,args)=parser.parse_args()

picklefiles=[]
labels=[]
for arg in args:
	picklefiles.append(arg.split(':')[0])
	labels.append(arg.split(':')[1])

fig = pylab.figure(1)
ax = pylab.axes()
matplotlib.rcParams.update({'legend.fontsize': 8})
ax.set_color_cycle(['b','c','g','y',(0.9,0.5,0.2),'r','m',(0.5,0.,0.8),'k'])

for file in picklefiles:
	data = pickle.load(open(file))
	pylab.loglog(data[0],data[1])
	pylab.xlabel('False Alarm Probability')
	pylab.ylabel('Efficiency')
	pylab.xlim([0,2])
	pylab.ylim([0,2])
	pylab.hold(True)

pylab.legend(labels,loc=4)
pylab.title(opts.title)
#pylab.text(0.005,1.0,typelabel,horizontalalignment='center')
pylab.savefig(opts.output_dir+'/combined_ROC_'+opts.tag+'.png')

