#!/usr/bin/python

# $Id$
__author__ = "Jacob Slutsky <jslutsky@phys.lsu.edu>, Chad Hanna <channa@phys.lsu.edu>, Romain Gouaty <romain@phys.lsu.edu>"
__version__ = "$Revision$"[11:-2]
__date__ = "$Date$"[7:-2]

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
from pylal import CoincInspiralUtils

# this imports the follow up utilities necessary to run this code
from pylal.fu_utils import *

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
  --statistic: 			Coinc stat to use
  --bitten-l-a:                 bitten l a
  --bitten-l-b:                 bitten l b
  --xml-glob:  			A glob of xml files 
  --injection-glob:  		A glob of the astrophysical injections
  --cluster-time:		The cluster time for triggers
  --injection-coinc-time:  	How close to consider a found injection
  --data-type		        playground_only | all_data | exclude_play 
  --num-trigs			The number of triggers to follow up
"""



##############################################################################
# This Section parses the options given on the command line.
# The goal is to keep the options to a minimum keep this in mind
# when developing all modules!
##############################################################################
parser = OptionParser( usage=usage, version= "%prog CVS $Id$")

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
    help="injection cluster time in seconds (default = .10)")

parser.add_option("-d","--data-type",action="store",type="string",\
    default="all_data",metavar=" DATA_TYPE",\
    help="data type (default = all_data)")

parser.add_option("-n","--num-trigs",action="store",type="float",\
    default=1,metavar=" NUM_TRIGS",\
    help="number of triggers to follow up (default = 1)")


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

#############################################################################
# function to set up directories
#############################################################################
#def setupdirs():

#  try:
#    os.chdir("clustered")
#    os.chdir('..')
#  except: os.mkdir("clustered")

#  try:
#    os.chdir("found")
#    os.chdir('..')
#  except: os.mkdir("found")

#  try:
#    os.chdir("missed")
#    os.chdir('..')
#  except: os.mkdir("missed")

#  try:
#    os.chdir("followuptrigs")
#    os.chdir('..')
#  except: os.mkdir("followuptrigs")

#  try:
#    os.chdir("followupfound")
#    os.chdir('..')
#  except: os.mkdir("followupfound")

#  try:
#    os.chdir("followupmissed")
#    os.chdir('..')
#  except: os.mkdir("followupmissed")

#############################################################################
# function to glob files
#############################################################################
#def globxmlandinj():

#  #Glob for the files both injections and trigger files
#  fList = glob.glob(opts.xml_glob)
#  if not fList:
#    print >>sys.stderr, "The glob for " + opts.xml_glob + " returned no files"
#    sys.exit(1)

#  if opts.injection_glob:
#    iList = glob.glob(opts.injection_glob)
#    if not iList:
#      print >> sys.stderr, "The glob for " + opts.injection_glob + " returned no files"
#      sys.exit(1)
#  else:
#    iList = None
#
#  return fList,iList


#############################################################################
# function to return the number of slides in a file (as a string)
#############################################################################
#def getslidenum(fName):
#  command = "grep -m 1 -e '--num-slides' " + fName + " | sed -e 's@\"@\'@g"
#  fStr = os.popen(command).readlines()
#  if fStr:
#    fStrlist = fStr[0].rsplit(",")
#    return fStrlist[-2]
#  else:
#   return "0"


##############################################################################
# function to coire the input files 
##############################################################################
#def coire(fList,iList=None):
#  """
#  do coires
#  """
#  sims = None
#  coincs = None
#
#  # we need to get slide num and data type from the files
#  # this is just hard coded now!!!
#  for xmls in fList:
#    if opts.verbose: 
#      print "running lalapps_coire on " + xmls 
#    numslides = getslidenum(xmls)
#    command = "lalapps_coire " + \
#              "--glob '" + xmls + "'" + \
#              " --data-type " + str(opts.data_type) + \
#              " --num-slides " + numslides + " --coinc-stat " + \
#              opts.statistic + " --cluster-time " + \
#              str(1000*opts.cluster_time)
#            
#    if iList:
#      command += " --injection-file " + iList[0] 
#      command += " --injection-window " + str(opts.injection_coinc_time)
#      command += " --missed-injections missed/MISSED" + xmls 
#      command += " --output found/FOUND" + xmls
#    else:
#      command += " --output clustered/CLUSTER" 
#      command += str(opts.cluster_time) + "s_" + xmls 
#    if opts.bitten_l_a and opts.bitten_l_b:
#      command += " --h1-bittenl-a " + str(opts.bitten_l_a)
#      command += " --h2-bittenl-a " + str(opts.bitten_l_a)
#      command += " --l1-bittenl-a " + str(opts.bitten_l_a)
#      command += " --h1-bittenl-b " + str(opts.bitten_l_b)
#      command += " --h2-bittenl-b " + str(opts.bitten_l_b)
#      command += " --l1-bittenl-b " + str(opts.bitten_l_b)
#
#    os.system(command) 	

##############################################################################
# function to extract the statistic information
##############################################################################
#def getstatistic():
#
#  if opts.statistic == "effective_snrsq":
#    newstat = "effective_snr"
#  else: 
#    newstat = opts.statistic 
#
#  statistic=CoincInspiralUtils.coincStatistic( newstat, opts.bitten_l_a, 
#                                               opts.bitten_l_b)
#  return statistic  

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

#############################################################################
# Follow up list class definition
#############################################################################
class followUpList:
  """
  Class to hold gps and ifo pairs to send to subsequent functions
  It also holds an instance of the coinc class which contains 
  All of the relevant xml information
  """
  def __init__(self,Coincs = None, Missed = None ):
    self.gpsTime = {"H1" : None, "H2" : None, "L1" : None,
                  "G1" : None, "V1" : None, "T1" : None}
    self.coincs = Coincs
    self.missed = Missed
  def add_coincs(self,Coincs):
    setattr(self,"coincs",Coincs)
  def add_missed(self,Missed):
    setattr(self,"missed",Missed)
  def is_trigs(self):
    if isinstance(self.coincs,CoincInspiralUtils.coincInspiralTable.row):
      return 1
  #def is_found(self):
  #def is_missed(self):
  
#############################################################################
# Function to return the follow up list of coinc triggers
#############################################################################
def getfollowuptrigs(coincs=None,missed=None):
  
  followups = []

  if coincs:
    coincs.sort()
    numTrigs = 0
    for ckey in coincs:
      fuList = followUpList()
      fuList.add_coincs(ckey)
      try:
        getattr(ckey,'H1')
        fuList.gpsTime["H1"] = (float(getattr(ckey,'H1').end_time_ns)/1000000000)+float(getattr(ckey,'H1').end_time)  
      except: fuList.gpsTime["H1"] = None
      try:
        getattr(ckey,'H2')
        fuList.gpsTime["H2"] = (float(getattr(ckey,'H2').end_time_ns)/1000000000)+float(getattr(ckey,'H2').end_time)
      except: fuList.gpsTime["H2"] = None
      try:
        getattr(ckey,'L1')
        fuList.gpsTime["L1"] = (float(getattr(ckey,'L1').end_time_ns)/1000000000)+float(getattr(ckey,'L1').end_time)
      except: fuList.gpsTime["L1"] = None
      try:
        getattr(ckey,'G1')
        fuList.gpsTime["G1"] = (float(getattr(ckey,'G1').end_time_ns)/1000000000)+float(getattr(ckey,'G1').end_time)
      except: fuList.gpsTime["G1"] = None
      try:
        getattr(ckey,'V1')
        fuList.gpsTime["V1"] = (float(getattr(ckey,'V1').end_time_ns)/1000000000)+float(getattr(ckey,'V1').end_time)
      except: fuList.gpsTime["V1"] = None
      try:
        getattr(ckey,'T1')
        fuList.gpsTime["T1"] = (float(getattr(ckey,'T1').end_time_ns)/1000000000)+float(getattr(ckey,'T1').end_time)
      except: fuList.gpsTime["T1"] = None
      followups.append(fuList) 
      numTrigs += 1
      if numTrigs >= opts.num_trigs:
        break
  
  # the missed stuff doesnt work yet!!!
  if missed:
    followups = None

  return followups

  
#############################################################################
# Class to hold summary HTML information for all of the functions
#############################################################################
class summaryHTMLTable:

  def __init__(self,trig):
    if trig.is_trigs():
      self.eventID = trig.coincs.event_id
      self.statValue = trig.coincs.stat
    else:
      self.eventID = None
      self.statValue = None
    self.H1time = trig.gpsTime["H1"]
    self.H2time = trig.gpsTime["H2"]
    self.L1time = trig.gpsTime["L1"]
    self.G1time = trig.gpsTime["G1"]
    self.V1time = trig.gpsTime["V1"]
    self.T1time = trig.gpsTime["T1"]
    self.containers = []
    
class HTMLcontainer:
  
  def __init__(self,trig):
    # The injections dont work yet!!!
    if trig.is_trigs():
      self.path = "followuptrigs/"
    else: self.path = ""
    self.name = __name__
    self.image = self.path + __name__ + ".png"
    self.text = "click here"
    self.link = self.path + self.name + ".html"

##############################################################################
# Function to write the HTML tables to pages
##############################################################################
def writeHTMLTables(summaryHTMLlist):
  for table in summaryHTMLlist:
    print table.eventID
    for container in table.containers:
      print container.name
      # this is the wrong value for function name - it give __main__


##############################################################################
# Examples of follow up functions
##############################################################################
def dosomething(trig):
  container = HTMLcontainer(trig)
  return container

def dosomethingelse(trig):
  container = HTMLcontainer(trig)
  # lets say I don't have an image for this function or a link
  # I need to set their value to None -  that will prevent
  # the HTML from displaying an empty image or a bogus link!
  # Otherwise I should make sure and produce and image with the
  # default name and create a file with the default link name!
  container.image = None
  container.link = None
  return container

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
  found, coincs = readFiles("found/FOUND"+opts.injection_glob,getstatistic(opts))
  missed, junk = readFiles("missed/MISSED"+opts.injection_glob,getstatistic(opts))
else:
  junk,coincs = readFiles("clustered/CLUSTER" + str(opts.cluster_time) + 
                          "s_" + opts.xml_glob, getstatistic(opts))

# this will get either the found injections (if injections were specified)
# and the missed injections 
# or just the coincidence triggers 
followuptrigs = getfollowuptrigs(coincs,missed)

summaryHTMLlist = []


for trig in followuptrigs:
  trig.is_trigs()
  # Set up the summary table
  summaryHTML = summaryHTMLTable(trig)
  # Call the followup functions they should all return
  # an HTMLcontainer class
  summaryHTML.containers.append(dosomething(trig))
  summaryHTML.containers.append(dosomethingelse(trig))
  # ....................................................
  # Add your functions Here in th following way        :
  # summaryHTML.containers.append(<yourfunction>(trig)): 
  # ....................................................
  # Append the summary table to the list of tables 
  summaryHTMLlist.append(summaryHTML)
  
writeHTMLTables(summaryHTMLlist)

sys.exit(1)
