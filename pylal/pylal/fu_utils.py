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

#############################################################################
# function to set up directories
#############################################################################
def setupdirs():

  try:
    os.chdir("clustered")
    os.chdir('..')
  except: os.mkdir("clustered")

  try:
    os.chdir("found")
    os.chdir('..')
  except: os.mkdir("found")

  try:
    os.chdir("missed")
    os.chdir('..')
  except: os.mkdir("missed")

  try:
    os.chdir("followuptrigs")
    os.chdir('..')
  except: os.mkdir("followuptrigs")

  try:
    os.chdir("followupfound")
    os.chdir('..')
  except: os.mkdir("followupfound")

  try:
    os.chdir("followupmissed")
    os.chdir('..')
  except: os.mkdir("followupmissed")


#############################################################################
# function to return the number of slides in a file (as a string)
#############################################################################
def getslidenum(fName):
  command = "grep -m 1 -e '--num-slides' " + fName + " | sed -e 's@\"@\'@g"
  fStr = os.popen(command).readlines()
  if fStr:
    fStrlist = fStr[0].rsplit(",")
    return fStrlist[-2]
  else:
   return "0"

#############################################################################
# function to glob files
#############################################################################
def globxmlandinj(opts):

  #Glob for the files both injections and trigger files
  fList = glob.glob(opts.xml_glob)
  if not fList:
    print >>sys.stderr, "The glob for " + opts.xml_glob + " returned no files"
    sys.exit(1)

  if opts.injection_glob:
    iList = glob.glob(opts.injection_glob)
    if not iList:
      print >> sys.stderr, "The glob for " + opts.injection_glob + " returned no files"
      sys.exit(1)
  else:
    iList = None

  return fList,iList

##############################################################################
# function to coire the input files
##############################################################################
def coire(opts,fList,iList=None):
  """
  do coires
  """
  sims = None
  coincs = None

  # we need to get slide num and data type from the files
  # this is just hard coded now!!!
  for xmls in fList:
    if opts.verbose:
      print "running lalapps_coire on " + xmls
    numslides = getslidenum(xmls)
    command = "lalapps_coire " + \
              "--glob '" + xmls + "'" + \
              " --data-type " + str(opts.data_type) + \
              " --num-slides " + numslides + " --coinc-stat " + \
              opts.statistic + " --cluster-time " + \
              str(1000*opts.cluster_time)

    if iList:
      command += " --injection-file " + iList[0]
      command += " --injection-window " + str(opts.injection_coinc_time)
      command += " --missed-injections missed/MISSED" + xmls
      command += " --output found/FOUND" + xmls
    else:
      command += " --output clustered/CLUSTER"
      command += str(opts.cluster_time) + "s_" + xmls
    if opts.bitten_l_a and opts.bitten_l_b:
      command += " --h1-bittenl-a " + str(opts.bitten_l_a)
      command += " --h2-bittenl-a " + str(opts.bitten_l_a)
      command += " --l1-bittenl-a " + str(opts.bitten_l_a)
      command += " --h1-bittenl-b " + str(opts.bitten_l_b)
      command += " --h2-bittenl-b " + str(opts.bitten_l_b)
      command += " --l1-bittenl-b " + str(opts.bitten_l_b)

    os.system(command)

##############################################################################
# function to extract the statistic information
##############################################################################
def getstatistic(opts):

  if opts.statistic == "effective_snrsq":
    newstat = "effective_snr"
  else:
    newstat = opts.statistic

  statistic=CoincInspiralUtils.coincStatistic( newstat, opts.bitten_l_a,
                                               opts.bitten_l_b)
  return statistic


