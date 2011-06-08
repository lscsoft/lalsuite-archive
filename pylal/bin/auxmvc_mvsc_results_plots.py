#!/usr/bin/python
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

parser=OptionParser(usage="Generates summary plots for KW aux-triggers processed by StatPatternRecognition", version = "Kari Hodge")
parser.add_option("","--histograms", action="store_true", default=False, help="use if you want to produce histograms for each dimension")
parser.add_option("","--plot-rank-vs-significance", action="store_true", default=False, help="Make scatter plots of ranke vs Aux KW significance")
parser.add_option("","--tag", help="filenames will be ROC_tag.png and efficiency_deadtime_tag.txt")
parser.add_option("","--output-dir", help="directory where output files will be written to")
(opts,files)=parser.parse_args()

try: os.mkdir(opts.output_dir)
except: pass

print files
data = auxmvc_utils.ReadMVSCTriggers(files)        

clean_data = data[numpy.nonzero(data['i']==0)[0],:]
glitch_data = data[numpy.nonzero(data['i']==1)[0],:]

all_ranks = numpy.concatenate((clean_data['Bagger'],glitch_data['Bagger']))
all_ranks_sorted = numpy.sort(all_ranks)


FAP, DP = auxmvc_utils.ROC(clean_data['Bagger'], glitch_data['Bagger'])


# Plot ROC curve
pylab.figure(1)
pylab.loglog(FAP,DP, linewidth = 2.0)
pylab.hold(True)
x = numpy.arange(min(DP), max(DP) + (max(DP) - min(DP))/1000.0, (max(DP) - min(DP))/1000.0)
pylab.loglog(x,x, linestyle="dashed", linewidth = 2.0)
pylab.xlabel('False Alarm Probability')
pylab.ylabel('Efficiency')
pylab.xlim([0,1])
pylab.ylim([0,1])
pylab.savefig(opts.output_dir + "/"+'ROC_'+opts.tag+'.png')	
pylab.close()

# save ROC curve in a file
roc_file = open(opts.output_dir + "/"+"ROC_" + opts.tag + ".pickle", "w")
pickle.dump(FAP, roc_file)
pickle.dump(DP, roc_file)
roc_file.close()

#print FAP 
#FAP is a list that is naturally sorted in reverse order (highest to lowest),
#we need to turn it into a regularly sorted list so that we can find the DP for
#fiducial FAPs
#print len(FAP)
FAP.sort()
edfile = open(opts.output_dir + "/"+'efficiency_deadtime_'+opts.tag+'.txt','w')
for threshold in [.01,.05,.1]:
	tmpindex=bisect.bisect_left(FAP,threshold)
	edfile.write("deadtime: "+str(FAP[tmpindex])+" efficiency: "+str(DP[len(FAP)-tmpindex-1])+"\n")



if opts.plot_rank_vs_significance:

  # plot rank vs significance
  variables = data.dtype.names
  sig_variables = [s for s in variables if s.endswith("_sig")]

  #form arrays with significances for clean data sets
  clean_sig_rss = []
  clean_sig_max = []
  for i in range(len(clean_data)):
    sig_rss = 0
    sig_max = 0
    for var in sig_variables:
      sig_rss += clean_data[var][i]**2
      if clean_data[var][i] > sig_max:
        sig_max = clean_data[var][i]
    clean_sig_rss.append(sig_rss**0.5)
    clean_sig_max.append(sig_max)
    

  #form arrays with significances for glitch data sets
  glitch_sig_rss = []
  glitch_sig_max = []
  for i in range(len(glitch_data)):
    sig_rss = 0
    sig_max = 0
    for var in sig_variables:
      sig_rss += glitch_data[var][i]**2
      if glitch_data[var][i] > sig_max:
        sig_max = glitch_data[var][i]
    glitch_sig_rss.append(sig_rss**0.5)
    glitch_sig_max.append(sig_max)
    

  pylab.figure(1)
  pylab.semilogx(clean_sig_max,clean_data['Bagger'], "kx", label="clean")
  pylab.hold(True)
  pylab.semilogx(glitch_sig_max,glitch_data['Bagger'], "r+", label="glitch")
  pylab.hold(True)
  pylab.xlabel('Max Aux KW significance')
  pylab.ylabel('MVSC rank')
  pylab.savefig(opts.output_dir + "/"+'rank_vs_sig_max_'+opts.tag+'.png')	
  pylab.close()


  pylab.figure(1)
  pylab.semilogx(clean_sig_rss,clean_data['Bagger'], "kx", label="clean")
  pylab.hold(True)
  pylab.semilogx(glitch_sig_rss,glitch_data['Bagger'], "r+", label="glitch")
  pylab.hold(True)
  pylab.xlabel('rss Aux KW significance')
  pylab.ylabel('MVSC rank')
  pylab.savefig(opts.output_dir + "/"+'rank_vs_sig_rss_'+opts.tag+'.png')	
  pylab.close()




if opts.histograms:
  variables = data.dtype.names
  for i,var in enumerate(variables):
    pylab.figure(i)
    print var
    pylab.hist(clean_data[var],100)
    pylab.savefig(opts.output_dir + "/"+"hist_twosided"+var+"_cleantimes") 
    pylab.figure(100+i)
    pylab.hist(glitch_data[var],bins=100)
    pylab.savefig(opts.output_dir + "/"+"hist_twosided"+var+"_glitches") 
