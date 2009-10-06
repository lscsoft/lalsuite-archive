#!/usr/bin/python
__author__ = "Ruslan Vaulin <vaulin@gravity.phys.uwm.edu>"
__version__ = "$Revision: 1.7 $"[11:-2]
__date__ = "$Date: 2009/03/18 22:18:59 $"[7:-2]
__prog__="lv_final_clustering.py"
__Id__ = "$Id: lv_final_clustering.py,v 1.7 2009/03/18 22:18:59 jclayton Exp $"


#loading standard modules
from optparse import *
import glob
import sys
import os
import cPickle
import shutil

#loading modules used for input/output of data 
from glue import lal
from glue.ligolw import lsctables
from glue.ligolw import table
from glue.ligolw import utils
from glue.ligolw import ligolw
from glue.ligolw.utils import ligolw_add
from glue import segments
from glue import segmentsUtils
from  glue import iterutils

from pylal import CoincInspiralUtils
from pylal import SnglInspiralUtils
from pylal import SimInspiralUtils
from pylal import InspiralUtils
from pylal import SearchSummaryUtils
import numpy


def getSnglInspiralTableFromCoincTable(coinctable):
  """
  Return the all sngl inspiral triggers from the coinc inspiral table
  """
   
  sngl_table = lsctables.New(lsctables.SnglInspiralTable)
  for coinc in coinctable:
    ifos, ifolist = coinc.get_ifos()
    for ifo in ifolist:
      sngl_table.append(getattr(coinc, ifo))
  return sngl_table

# Functions for un-sliding triggers taken form SnglInspiral Utils with very minor modifications
# we need to get rings from coire files, so we added them top the function as input

def ReadSnglInspiralSlidesFromFiles(fileList, coire_file_list, shiftVector, vetoFile=None,
  mangleEventId=False, verbose=False, non_lsc_tables_ok=False):
  """
  Function for reading time-slided single inspiral triggers
  with automatic resliding the times, given a list of input files.

  @param fileList: List of files containing single inspiral time-slided
                   triggers
  @param shiftVector: Dictionary of time shifts to apply to triggers
                      keyed by IFO  
  @param vetoFile: segwizard formatted file used to veto all triggers
  @param mangleEventId: ID remapping is necessary in cases where multiple
    files might have event_id collisions (ex: exttrig injections)
  @param verbose: print ligolw_add progress
  """

  # read raw triggers
  inspTriggers = SnglInspiralUtils.ReadSnglInspiralFromFiles(\
    fileList, verbose=verbose, mangle_event_id=mangleEventId, non_lsc_tables_ok=non_lsc_tables_ok )
  if inspTriggers:
    # get the rings
    segDict = SearchSummaryUtils.GetSegListFromSearchSummaries(coire_file_list)
    rings = segments.segmentlist(iterutils.flatten(segDict.values()))
    rings.sort()
    # FIXME:  remove with thinca's ring boundary bug is fixed
    rings = segments.segmentlist(segments.segment(ring[0], ring[1] + 1e-9) for ring in rings)

    # perform the veto
    if vetoFile is not None:
      segList = segmentsUtils.fromsegwizard(open(vetoFile))
      inspTriggers = inspTriggers.veto(segList)

    # now slide all the triggers within their appropriate ring
    SnglInspiralUtils.slideTriggersOnRingWithVector(inspTriggers, shiftVector, rings)

  # return the re-slided triggers
  return inspTriggers





################################################################################
# Main program
################################################################################
usage= """
usage: %prog [options]

This code calculates false alarm rate (FAR) in untis [number of events per year] for  coincident events in the input file and saves its inverse in alpha column of sngl_inspiral table.
"""
###############################################################################
# Options to read in Input
###############################################################################
def parse_command_line():

  """
  Parser function dedicated
  """

  parser = OptionParser( usage=usage, version="%prog CVS $Id: inspiral_search_summary,v 1.7 2009/03/18 22:18:59 jclayton Exp $ " )

  parser.add_option("","--slides-glob",action="store",type="string",\
      default=None, metavar=" GLOB",help="GLOB time slides files to read" )
  
  parser.add_option("","--zerolag-glob",action="store",type="string",\
      default=None, metavar=" GLOB",help="GLOB zerolag files to read" )
	  
  parser.add_option("","--statistic",action="store",default='snr',\
      type="string",\
      help="choice of statistic used in building coinc table, valid arguments are: snr (DEFAULT), snr_over_chi, s3_snr_chi_stat, effective_snr, bitten_l, bitten_lsq, ifar")
	  
  parser.add_option("","--num-slides", action="store",type="int",\
      default = 0, metavar="numslides", help="number of time slides performed, must match the corresponding parameter from the .ini file of the search" )

  parser.add_option("","--cluster-window",action="store",type="int",default=10,\
      metavar="SEC",help="clustering time" )
	  
	  
  parser.add_option("","--gps-start-time", action="store",type="int",\
      default = 0, help="start time of the analysis, used in the names of the output files" )

  parser.add_option("","--gps-end-time", action="store",type="int",\
      default = 0, help="end time of the analysis, used in the names of the output files" )


  parser.add_option("-u","--user-tag",action="store",type="string",\
      default=None, metavar=" USERTAG",\
      help="The user tag used in the name of the figures" )

  parser.add_option("","--ignore-IFO-times",action="store",type="string",\
      default='', metavar=" USERTAG",\
      help="comma separated list of IFOTYPE_IFOTIME that should not be included in efficiency calculation e.g. H1H2_H1H2,H2L1_H1H2L1. This option will work only with s5 LV search files that are named IFOTIME_IFOTYPE-*xml" )
  
  parser.add_option("","--output-dir",action="store",\
      type="string",default=None,  metavar="PATH",\
      help="directory for reclustered files")
	  
  parser.add_option("", "--h1-triggers",action="store_true", default=False,\
      help="input files contain triggers from H1")

  parser.add_option("", "--h2-triggers",action="store_true", default=False,\
      help="input files contain triggers from H2")

  parser.add_option("", "--l1-triggers",action="store_true", default=False,\
      help="input files contain triggers from L1")

  parser.add_option("", "--g1-triggers",action="store_true", default=False,\
      help="input files contain triggers from G1")

  parser.add_option("", "--v1-triggers",action="store_true", default=False,\
      help="input files contain triggers from V1")
	  
  parser.add_option("","--output-path",action="store",\
      type="string",default=None,  metavar="PATH",\
      help="path where the figures would be stored")

  parser.add_option("", "--figure-resolution",action="store",type="int",\
      default=50, metavar="FIGURERESOLUTION", \
      help="resolution of the thumbnails (50 by default)" )
	  
  parser.add_option("","--verbose", action="store_true",\
      default=False, help="print information" )
	  
  parser.add_option("-O","--enable-output",action="store_true",\
      default="True",  metavar="OUTPUT",\
      help="enable the generation of the html and cache documents")

  parser.add_option("", "--html-for-cbcweb",action="store",\
      default=False, metavar = "CVS DIRECTORY", help="publish the html "\
      "output in a format that can be directly published on the cbc webpage "\
      "or in CVS. This only works IF --enable-output is also specified. The "\
      "argument should be the cvs directory where the html file will be placed "\
      "Example: --html-for-cbcweb protected/projects/s5/yourprojectdir")


  (opts,args) = parser.parse_args()


  return opts, sys.argv[1:]
#####################################################################
opts, args = parse_command_line()

# Initializing the html output
InspiralUtils.message(opts, "Initialisation...")
opts = InspiralUtils.initialise(opts, __prog__, __version__)
fnameList = []
tagList = []
fig_num = 0
comments = ""


#Calculating statistic for coincidences
statistic = CoincInspiralUtils.coincStatistic(opts.statistic) 

# constructing the list of the IFO's
ifo_list = [ifo for ifo in ("G1", "H1", "H2", "L1", "V1") \
            if getattr(opts, "%s_triggers" % ifo.lower())]

#constructing the list of all possible IFO combinations (doubles, triples etc) 
ifo_combos = CoincInspiralUtils.get_ifo_combos(ifo_list)

ifo_times = []
for ifo_combo in ifo_combos:
  #ifo_combo = ifo_combo.sort()
  ifo_times.append("".join(ifo_combo))
  
# create output directory for clustered likelihood files
try: os.mkdir(opts.output_dir)
except: pass

output_dir = opts.output_dir


# constructing list of zerolag files
zero_lag_files = glob.glob(opts.zerolag_glob)


# contsructing list of  time slides files
slidesfiles = glob.glob(opts.slides_glob)

if not len(slidesfiles) > 0:
  print >>sys.stderr, "List of time slides files is empty: your sieve (glob) pattern may be wrong or files do not exist in the location given by the cache file"
  #sys.exit(1)

# get rid of certain IFO times if necessary
if opts.ignore_IFO_times:

  ifo_times_to_ignore = opts.ignore_IFO_times.split(",")
  
  
  
  
  # zero lag
  # Assumes that ifo time is the first thing in the file name
  new_zero_lag_files = []
  for file in zero_lag_files:
    tmpifotime=file.split("/")[-1].split("_")[0]
    tmpifotype=file.split("/")[-1].split("_")[1].split("-")[0]
    category="_".join([tmpifotype,tmpifotime])
    if not (category in ifo_times_to_ignore):
      new_zero_lag_files.append(file)

  zero_lag_files = []
  zero_lag_files = new_zero_lag_files   
  
  # time slides
  # Assumes that ifo time is the first thing in the file name
  new_slidesfiles = []
  for file in slidesfiles:
	tmpifotime=file.split("/")[-1].split("_")[0]
	tmpifotype=file.split("/")[-1].split("_")[1].split("-")[0]
	category="_".join([tmpifotype,tmpifotime])
	if not (category in ifo_times_to_ignore):
	  new_slidesfiles.append(file)

  slidesfiles = []
  slidesfiles = new_slidesfiles


############# CLUSTERING ##############################################################

# Cluster zero lag

if opts.verbose:
  print "Reading zero-lag files ..."
# read in sngl inspiral table
Triggers = SnglInspiralUtils.ReadSnglInspiralFromFiles(zero_lag_files, mangle_event_id=False, non_lsc_tables_ok=True)
# construct coincidence 
CoincTriggers = CoincInspiralUtils.coincInspiralTable(Triggers, statistic)

if opts.verbose:
  print "There were " + str(len(CoincTriggers)) + " coincident triggers in the input files."
# cluster them 

# save clustered triggers
CoincTriggers.cluster(opts.cluster_window)

if opts.verbose:
  print " There are " + str(len(CoincTriggers)) + " coincident triggers after clustering." 

# open the output file for writing
output_file_name = output_dir + "/" + "H1H2L1V1_H1H2L1V1-LIKELIHOOD_CLUSTERED-" + str(opts.gps_start_time) + "-" + str(opts.gps_end_time - opts.gps_start_time) + ".xml"
output_file =  open(output_file_name, "w")

# create xml doc
output_doc = ligolw.Document()

#create LIGO_LW element
output_ligo_lw = ligolw.LIGO_LW()

#append it to xml doc
output_doc.appendChild(output_ligo_lw)

# get clustered coincs as sngl_inspiral table
sngl_table = getSnglInspiralTableFromCoincTable(CoincTriggers)
if len(sngl_table) > 0:
  output_doc.childNodes[0].appendChild(sngl_table)

#writing xml doc to the output file
output_doc.write(output_file)
output_file.close()


# Cluster time slides
    
if opts.verbose:
  print "reading time slides  triggers ..."
slidesTriggers = SnglInspiralUtils.ReadSnglInspiralFromFiles(slidesfiles, non_lsc_tables_ok=True)
 

# construct the time slides coincs
slidesCoincTriggers = CoincInspiralUtils.coincInspiralTable(slidesTriggers, statistic)

if opts.verbose:
  print "There were " + str(len(slidesCoincTriggers)) + " coincident triggers in the input files."
# cluster  triggers in time slides
slidesCoincTriggers.cluster(cluster_window = opts.cluster_window, numSlides = 50)
  
if opts.verbose:
  print " There are " + str(len(slidesCoincTriggers)) + " coincident triggers after clustering." 

  
# open the output file for writing
output_file_name = output_dir + "/" + "H1H2L1V1_H1H2L1V1-LIKELIHOOD_SLIDES_CLUSTERED-"+ str(opts.gps_start_time) + "-" + str(opts.gps_end_time - opts.gps_start_time) + ".xml"
output_file =  open(output_file_name, "w")

# create xml doc
output_doc = ligolw.Document()

#create LIGO_LW element
output_ligo_lw = ligolw.LIGO_LW()

#append it to xml doc
output_doc.appendChild(output_ligo_lw)

# get clustered coincs as sngl_inspiral table
sngl_table = getSnglInspiralTableFromCoincTable(slidesCoincTriggers)
if len(sngl_table) > 0:
  output_doc.childNodes[0].appendChild(sngl_table)

#writing xml doc to the output file
output_doc.write(output_file)
output_file.close()



# Make plots
 
slides_numbers = numpy.asarray(range(-opts.num_slides, opts.num_slides +1))


# Change to Agg back-end if show() will not be called 
# thus avoiding display problem
import matplotlib
matplotlib.use('Agg')
import pylab

fig_num = 0

# make histograms for underclustered triggers
ifo_combos = CoincInspiralUtils.get_ifo_combos(ifo_list)
for ifo_combo in ifo_combos:
  slides_coincs_from_combo = slidesCoincTriggers.coinctype(ifo_combo)
  zerolag_coincs_from_combo = CoincTriggers.coinctype(ifo_combo)
  
  num_triggers_in_slides = numpy.zeros(len(slides_numbers))
  for (i, slide) in enumerate(slides_numbers) :
	if slide != 0:
	  num_triggers_in_slides[i] = len(slides_coincs_from_combo.getslide(slide))
	else:
	  num_triggers_in_slides[i] = len(zerolag_coincs_from_combo.getslide(slide))
  
  fig_num += 1
  pylab.figure(fig_num)
  legend_text = []

  pylab.bar(slides_numbers, num_triggers_in_slides, width=(0.8), color='b', edgecolor="k")
  pylab.hold(True)
  pylab.bar(slides_numbers[(len(num_triggers_in_slides)-1)/2], num_triggers_in_slides[(len(num_triggers_in_slides)-1)/2], width=(0.8), color='r', edgecolor="k")
  
  pylab.title(" Number of triggers in time slides: " + "".join(ifo_combo))
  name = "_hist_slides_" + "".join(ifo_combo) + "_in_"
  fname = InspiralUtils.set_figure_name(opts, name)
  fname_thumb = InspiralUtils.savefig_pylal(filename=fname, doThumb=True, dpi_thumb=opts.figure_resolution)
  fnameList.append(fname)
  tagList.append("Histogram of triggers in time slides: " + "".join(ifo_combo))
  pylab.close()
	
##############################################################################################################


html_filename = InspiralUtils.write_html_output(opts, args, fnameList, tagList, comment=comments)
InspiralUtils.write_cache_output(opts, html_filename, fnameList)

if opts.html_for_cbcweb:
  html_filename_publish = InspiralUtils.write_html_output(opts, args, fnameList, tagList, cbcweb=True)

##############################################################################################################
  













