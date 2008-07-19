#!/usr/bin/env python
"""
Something

$Id$

This program reads the coire files containing the list of coincident triggers, build separated lists of single ifo triggers from these coincidences, and sort the triggers by snr or effective snr.
"""

__author__ = 'Romain Gouaty <romain@phys.lsu.edu>'
__date__ = '$Date$'
__version__ = '$Revision$'[11:-2]

##############################################################################
# import standard modules and append the lalapps prefix to the python path
import sys, os, copy, math, random
import socket, time
import re, string
import commands
from optparse import *
import tempfile
import ConfigParser
import urlparse
import urllib
from UserDict import UserDict
sys.path.append('@PYTHONLIBDIR@')

##############################################################################
# import the modules we need to build the pipeline
from glue import lal
from glue import segments
from glue import segmentsUtils
from glue.ligolw import ligolw
from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw import utils
from pylal import CoincInspiralUtils
from pylal.fu_utils import *
from pylal.fu_writeXMLparams import *
from pylal.webUtils import *
from pylal import Fr
from lalapps import inspiral

##############################################################################
#
#  INPUTS
#
##############################################################################

inputDir = "/archive/home/dkeppel/post_processing/lowcbc/20051104-20061114/upperlimits_v99/corse_dag/corse/"
#inputDir = "/archive/home/romain/Projects/LowMassCBC/20051104-20061114/triggers/triggers_2_8_17_35/"
mass_bin_list = ["mchirp_2_8","mchirp_8_17","mchirp_17_35"]
cat_dir = "full_data_cat123"
#ifo_coincs = ["H1H2L1","H1L1","H2L1","H1H2"]
ifo_coincs = ["H1H2L1","H1L1","H2L1"]
triggerFileString = "815160323-32395247"
nb_loudest = 100


outputFile = "/archive/home/romain/Projects/LowMassCBC/20051104-20061114/triggers/upperlimits_v99_july08/triggerPerIfos/loudest_triggers_upperlimits_v99_070808-"

##############################################################################
#
#  MAIN PROGRAM
#
##############################################################################

trigger_list_H1 = []
trigger_list_H2 = []
trigger_list_L1 = []

for ifo_coinc in ifo_coincs:
  if ifo_coinc == "H1H2L1":
    ifo_times = ["H1H2L1"]
  if ifo_coinc == "H1L1":
    ifo_times = ["H1H2L1","H1L1"]
  if ifo_coinc == "H2L1":
    ifo_times = ["H1H2L1","H2L1"]
  if ifo_coinc == "H1H2":
    ifo_times = ["H1H2L1","H1H2"]

  for bin in mass_bin_list:
    for time in ifo_times:
      inputFile = inputDir + bin + "/" + cat_dir + "/" + ifo_coinc + "/" + ifo_coinc + "-CORSE_" + time + "-" + triggerFileString + ".xml.gz"
      #inputFile = inputDir + bin + "/hipecoire/" + cat_dir + "/" + ifo_coinc + "/" + ifo_coinc + "-COIRE_in_*_CLUST_10s.xml.gz"
      found, coincs, search = readFiles(inputFile,getstatistic("effective_snrsq",None,None))


      if coincs:
        for coinc_trigger in coincs:
          for j in range(0,len(ifo_coinc)-1,2):
            ifo = ifo_coinc[j:j+2]
            eval("trigger_list_" + ifo).append(getattr(coinc_trigger,ifo))

for IFO in ["H1","H2","L1"]:
  exec("sorted_trigger_list_" + IFO + " = [ (trigger.get_effective_snr(), trigger) for trigger in trigger_list_" + IFO + "]")
  eval("sorted_trigger_list_" + IFO + ".sort()")
  eval("sorted_trigger_list_" + IFO + ".reverse()")
  #print eval("sorted_trigger_list_" + IFO)

for IFO in ["H1","H2","L1"]:
  file = open(outputFile + IFO + ".txt",'w')
  file.write("rank" + "\t" + "end_time" + "\t" + "end_time_ns" + "\t" + "effective_snr" + "\t" + "snr" + "\t" + "FAR" + "\t" + "mchirp" + "\t" + "eta" + "\t" + "mass1" + "\t" + "mass2" + "\t" + "eff_distance" + "\n")
  numTrigs = 0
  for (effsnr,trig) in eval("sorted_trigger_list_" + IFO):
    numTrigs +=1
    if numTrigs > nb_loudest:
      break
    file.write(str(numTrigs) + "\t" + str(trig.end_time) + "\t" + str(trig.end_time_ns) + "\t" + str(trig.get_effective_snr()) + "\t" + str(trig.snr) + "\t" + str(trig.alpha) + "\t" + str(trig.mchirp) + "\t" + str(trig.eta) + "\t" + str(trig.mass1) + "\t" + str(trig.mass2) + "\t" + str(trig.eff_distance) + "\n")
  file.close()

