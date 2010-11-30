#!/usr/bin/python

__author__ = "Jacob Slutsky <jslutsky@phys.lsu.edu>, Chad Hanna <channa@phys.lsu.edu>, Romain Gouaty <romain@phys.lsu.edu>"

import sys
import os
import copy
import re
import exceptions
import glob
import fileinput
import linecache
import string
import random
from optparse import *
from types import *
import matplotlib
matplotlib.use('Agg')

from pylab import *
from glue import segments
from glue import segmentsUtils
from glue.ligolw import ligolw
from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw import utils
from pylal import git_version
from pylal import CoincInspiralUtils

# this imports the follow up utilities necessary to run this code
from pylal.fu_utils import *
from pylal.fu_dosomething import *
from pylal.fu_writeXMLparams import *
##############################################################################
# help message
##############################################################################
usage = """\
followup.py [options]

  Do some awesome stuff!
  This program assumes that you have lalapps_coire in the path
  It also assumes that you are using some linux os with sensible
  things like grep etc...

  ----------------------------------------------------------------------------
  --verbose:			Make a lot of text everywhere!	
  --statistic: 			Coinc stat to use [effective_snrsq]
  --bitten-l-a:                 bitten l a [3]
  --bitten-l-b:                 bitten l b [3]
  --xml-glob:  			A glob of xml files 
  --injection-glob:  		A glob of the astrophysical injections
  --cluster-time:		The cluster time for triggers [10]
  --injection-coinc-time:  	How close to consider a found injection [0.1]
  --data-type		        playground_only | all_data | exclude_play 
  --num-trigs			The number of triggers to follow up [1]
"""



##############################################################################
# This Section parses the options given on the command line.
# The goal is to keep the options to a minimum keep this in mind
# when developing all modules!
##############################################################################
parser = OptionParser( usage=usage, version=git_version.verbose_msg)

parser.add_option("-V","--verbose",action="store_true",default=False,\
    help="print additional information when running" )

# input files - Let's not try to input every type of xml file differently
# We'll write the code to be smart and figure it out by itself.  It should
# know how to handle each XML it is given.  That is the point of using
# XML in the first place.
parser.add_option("-I","--injection-glob",action="store",type="string",\
    default=None,metavar=" INJ_GLOB",\
    help="GLOB of files containing astrophysically distributed injections")

parser.add_option("-g","--xml-glob",action="store",type="string",\
    default=None, metavar=" XML_GLOB", \
    help="GLOB of xml files to read" )

# Statistic - We don't need to know the number of slides - we can get
# that from the xml - we do need to know what statistic should be
# used for clustering, sorting, etc...
# let's simplify the handling of bitten l statistics by including it
# with statistic option as "--statistic bitten-l <a> <b>"

parser.add_option("-K","--statistic",action="store",type="string",\
    default="effective_snrsq",metavar=" STAT",\
    help="coincident statistic (default = effective_snr)")

parser.add_option("-a","--bitten-l-a",action="store",type="float",\
    default=3,metavar=" BLA",\
    help="bitten l a parameter")

parser.add_option("-b","--bitten-l-b",action="store",type="float",\
    default=3,metavar=" BLB",\
    help="bitten l b parameter")


parser.add_option("-c","--cluster-time",action="store",type="float",\
    default=10,metavar=" CLUST_TIME",\
    help="cluster time in seconds (default = 10)")

parser.add_option("-i","--injection-coinc-time",action="store",type="float",\
    default=10,metavar=" INJ_TIME",\
    help="injection coincidence time in seconds (default = .10)")

parser.add_option("-d","--data-type",action="store",type="string",\
    default="all_data",metavar=" DATA_TYPE",\
    help="data type (default = all_data)")

parser.add_option("-n","--num-trigs",action="store",type="float",\
    default=1,metavar=" NUM_TRIGS",\
    help="number of triggers to follow up (default = 1)")

parser.add_option("-p","--page",action="store",type="string",\
    default="investigations/s5/people/followups/",metavar=" PAGE",\
    help="web page path (default 'investigations/s5/people/followups/'")


# Now that things are set up we will actually parse the options
(opts,args) = parser.parse_args()

##############################################################################
# Although we don't have many options we need to make sure they are sane
##############################################################################

if not opts.xml_glob:
  print >> sys.stderr, "Must specify a GLOB of xmls to read"
  sys.exit(1)

if not opts.statistic:
  print >> sys.stderr, "Must specify a statistic to use"
  sys.exit(1)



##############################################################################
# redefine the SimInspiral columns of interest
##############################################################################
lsctables.SimInspiralTable.loadcolumns = [
    "waveform",
    "geocent_end_time",
    "geocent_end_time_ns",
    "h_end_time",
    "h_end_time_ns",
    "l_end_time",
    "l_end_time_ns",
    "source",
    "mass1",
    "mass2",
    "mchirp",
    "eta",
    "distance",
    "spin1x",
    "spin1y",
    "spin1z",
    "spin2x",
    "spin2y",
    "spin2z",
    "eff_dist_h",
    "eff_dist_l",
    "eff_dist_g",
    "eff_dist_t",
    "eff_dist_v"]

##############################################################################
# redefine the SnglInspiral columns of interest
##############################################################################
lsctables.SnglInspiralTable.loadcolumns = [
    "ifo",
    "end_time",
    "end_time_ns",
    "eff_distance",
    "mass1",
    "mass2",
    "mchirp",
    "eta",
    "snr",
    "chisq",
    "chisq_dof",
    "sigmasq",
    "event_id"]

##############################################################################
# function to read in a list of files and extract the simInspiral tables
# and sngl_inspiral tables
##############################################################################
def readFiles(fileGlob,statistic=None):
  """
  read in the Sngl and SimInspiralTables from a list of files
  if Sngls are found, construct coincs, add injections (if any)
  also return Sims (if any)
  @param fileGlob: glob of input files
  @param statistic: statistic to use in creating coincs
  """
  #if fileGlob is empty return empty structures
  if not fileGlob:
    if opts.verbose:
      print "Warning: No glob specified, returning empty structures..."
    return None, CoincInspiralUtils.coincInspiralTable()

  # if there aren't any files globbed exit
  fList = glob.glob(fileGlob)
  if not fList:
    print >>sys.stderr, "The glob for " + fileGlob + " returned no files"
    sys.exit(1)

  sims = None
  coincs = None
  for thisFile in fList:
    doc = utils.load_filename(thisFile)
    # extract the sim inspiral table
    try:
      simInspiralTable = \
          table.get_table(doc, lsctables.SimInspiralTable.tableName)
      if sims: sims.extend(simInspiralTable)
      else: sims = simInspiralTable
    except: simInspiralTable = None

    # extract the sngl inspiral table, construct coincs
    try: snglInspiralTable = \
      table.get_table(doc, lsctables.SnglInspiralTable.tableName)
    except: snglInspiralTable = None
    if snglInspiralTable:
      coincInspiralTable = \
        CoincInspiralUtils.coincInspiralTable(snglInspiralTable,statistic)
      if simInspiralTable:
        coincInspiralTable.add_sim_inspirals(simInspiralTable)
      if coincs: coincs.extend(coincInspiralTable)
      else: coincs = coincInspiralTable
  return sims,coincs

##############################################################################
# MAIN PROGRAM - I like this heading.
##############################################################################

#Set up the directories to output files to
setupdirs()

#Glob for both injections and trigger files
fList,iList = globxmlandinj(opts)

#Coire the files
coire(opts,fList,iList)

#Read in the coired files
#In general we want different behavior for injections and triggers
#Namely we will follow up loud triggers and loud missed injections
#But quiet found injections...
(found, missed, coincs) = (None,None,None)
if opts.injection_glob:
  found, coincs = readFiles("found/FOUND"+opts.xml_glob,getstatistic(opts))
  missed, junk = readFiles("missed/MISSED"+opts.xml_glob,getstatistic(opts))
else:
  junk,coincs = readFiles("clustered/CLUSTER" + str(opts.cluster_time) + 
                          "s_" + opts.xml_glob, getstatistic(opts))

# this will get either the found injections (if injections were specified)
# and the missed injections 
# or just the coincidence triggers 

#print type(coincs)
#print type(coincs[0].sim)
#print type(missed)
#print coincs[0].sim.__class__.__dict__
#print coincs[0].__slots__

followuptrigs = getfollowuptrigs(opts,coincs,missed)

summaryHTMLlist = []


for trig in followuptrigs:
  trig.is_trigs()
  # Set up the summary table
  summaryHTML = summaryHTMLTable(trig)
  # Call the followup functions they should all return
  # an HTMLcontainer class
  summaryHTML.containers.append(dosomething(trig))
  summaryHTML.containers.append(writeXMLparams(trig))
  # ....................................................
  # Add your functions Here in th following way        :
  # summaryHTML.containers.append(<yourfunction>(trig)): 
  # ....................................................
  # Append the summary table to the list of tables 
  summaryHTMLlist.append(summaryHTML)
  
writeHTMLTables(summaryHTMLlist)

publishToIULGroup(opts.page)
sys.exit(1)
