#!/usr/bin/python
# Written to load and manipulate output from different classifiers.
# Reed Essick (reed.essick@ligo.org), Young-Min Kim (young-min.kim@ligo.org)
# last modified: 12/07/2011

from optparse import *
import glob
import sys
import os
import matplotlib
matplotlib.use('Agg')
import pylab
import pdb
import numpy
from pylal import auxmvc_utils 
import bisect
import pickle

def vetoGlitchesUnderFAP(MVSCTriggers,rank_name, FAPthr):
	clean_data = data[numpy.nonzero(data['i']==0)[0],:]
	glitch_data = data[numpy.nonzero(data['i']==1)[0],:]
	
	clean_rank_sorted = numpy.sort(clean_data[rank_name])
	glitches_sorted = numpy.sort(glitches,order(rank_name))
	for i,rank in enumerate(clean_rank_sorted):
		number_of_false_alarms = len(clean_ranks_sorted) -	numpy.searchsorted(clean_ranks_sorted,rank)
		number_of_true_alarms = len(glitch_ranks_sorted) -	numpy.searchsorted(glitches_sorted[rank_name],rank)
		FAP = number_of_false_alarms / float(len(clean_ranks_sorted))
		DP = number_of_true_alarms / float(len(glitch_ranks_sorted))
		if FAP < FAPthr:
			break

	return glitches_sorted[:i-1]

def ReadMVSCRanks(list_of_files, rank_name):
	data = auxmvc_utils.ReadMVSCTriggers(list_of_files)
	variables=['GPS','i',rank_name]
	formats=['g8','i','g8']
	gps_ranks = numpy.empty(len(data),dtype={'names':variable,'formats':formats})
	gps_ranks['i']=data['i']
	gps_ranks[rank_name]=data[rank_name]
	for i in range(len(data)):
		numpy.put(gps_ranks['GPS'],[i],[data['GPS_s'][i]+'.'+data['GPS_ms'][i]])
	
	return gps_ranks

	

parser=OptionParser(usage="Generates comparison plots", version = "auxmvc")
parser.add_option("","--hveto-ranked-files", default=False, help="Provide the path for HVeto Results files which contain all trigger data and those ranks and globbing pattern")
parser.add_option("","--mvsc-ranked-files", default=False, help="Provide the path for MVSC *.dat files and globbing pattern")
parser.add_option("","--ann-ranked-files", default=False, help="Provide the path for ANN *.dat files and globbing pattern")
parser.add_option("","--tag", help="filenames will be ROC_tag.png and efficiency_deadtime_tag.txt")
parser.add_option("","--output-dir", help="directory where output files will be written to")
(opts,files)=parser.parse_args()

try: os.mkdir(opts.output_dir)
except: pass

### Making ranked data to use glitch ranks and GW snr for plotting 
ranked_data={}

if opts.mvsc_ranked_files:
	list_of_mvsc=glob.glob(opts.mvsc_ranked_files)
	ranked_data["MVSC"] = auxmvc_utils.ReadMVSCTriggers(list_of_mvsc)

if opts.ann_ranked_files:
	list_of_ann=glob.glob(opts.ann_ranked_files)
	ranked_data["ANN"] = auxmvc_utils.ReadMVSCTriggers(list_of_ann)

if opts.hveto_ranked_files:
	ranked_data["HVeto"] = auxmvc_utils.loadHVefbydt(opts.hveto_ranked_files)

classifiers=ranked_data.keys()

