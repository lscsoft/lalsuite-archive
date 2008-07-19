#!/usr/bin/env python
"""
Something

$Id$

This program generates a table of the loudest candidates in wiki format.
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
#  INPUTS TO BE UPDATED APPROPRIATELY
#
##############################################################################

inputDir = "/archive/home/dkeppel/post_processing/lowcbc/20051104-20061114/upperlimits_v99/corse_dag/corse/"
mass_bin_list = ["mchirp_2_8","mchirp_8_17","mchirp_17_35"]
cat_veto_list = ["Cat1","Cat12","Cat123","Cat1234"]
cat_dir = ["full_data_cat1","full_data_cat12","full_data_cat123","full_data_cat1234"]
ifo_times = ["H1H2L1","H1L1","H2L1","H1H2"]
nbCoincCat = 7 #number of categories of coincidences (counting ifo times and coincidence types)
outputFile = "upperlimits_v99_July2008_candidate_table_v1.txt"
ifar_url = "http://ldas-jobs.ligo.caltech.edu/~dkeppel/results/s5_1yr_lowcbc/upperlimits_v99/corse_dag/Images/"
listOfRedTimes = "./candidate_redlist.txt"
far_threshold = 0.10
tableHeader1 = "|||| Time:  ||||||||<:> !H1H2L1 ||<:> !H1L1 ||<:> !H2L1 ||<:> !H1H2 ||"
tableHeader2 = "|||| Coincidence:  ||<:> !H1H2L1 ||<:> !H1L1 ||<:> !H2L1 ||<:> !H1H2||<:> !H1L1 ||<:> !H2L1 ||<:> !H1H2 ||"


##############################################################################
#
#  MAIN PROGRAM
#
##############################################################################

# call the method "listFromFile" from the pylal module "fu_utils.py" to read the text file "listOfRedTimes"
timeList = listFromFile(listOfRedTimes)

file = open(outputFile,'w')
file.write(tableHeader1 + "\n")
file.write("||||<:> all mass bins||<-" + str(nbCoincCat) + ">||\n")

for cat_veto in cat_veto_list:

  file.write("||||" + cat_veto + " ||||||||")
  for ifo_time in ifo_times:
    if not ifo_time == "H1H2":
      ifar_c = ifar_url + ifo_time + "-plotifar_ALL_DATA" + "_" + cat_veto.upper() + "_cumhist_ifar_combined-815160323-32395247.png"
      file.write("<:> [[ImageLink(" + ifar_c + ",width=150)]]||")
  file.write(" ||" + "\n")


file.write(tableHeader2 + "\n")

for mass_bin in mass_bin_list:
  file.write("||||<:> " + mass_bin + "||<-7>||\n")
  
  for k,cat_veto in enumerate(cat_veto_list):
    file.write("|| ||<:> " + cat_veto + "||")
    
    for ifo_time in ifo_times:
      if ifo_time == "H1H2L1":
        coincs = ["H1H2L1","H1L1","H2L1","H1H2"]
        nbCat = (len(coincs) - 1) * len(mass_bin_list)
      if ifo_time == "H1L1":
        coincs = ["H1L1"]
        nbCat = len(mass_bin_list)
      if ifo_time == "H2L1":
        coincs = ["H2L1"]
        nbCat = len(mass_bin_list)
      if ifo_time == "H1H2":
        coincs = ["H1H2"]
        nbCat = len(mass_bin_list)

      for coinc in coincs:
        inputFile = inputDir + mass_bin + "/" + cat_dir[k] + "/" + coinc + "/" + coinc + "-CORSE_" + ifo_time + "-815160323-32395247.xml.gz"
        found, coincs, search = readFiles(inputFile,getstatistic("effective_snrsq",None,None))
        followuptrigs = getfollowuptrigs("1",None,coincs,None,None,None)
        if len(followuptrigs) > 0:
          gps_time = followuptrigs[0].gpsTime[followuptrigs[0].ifolist_in_coinc[0]] 
        #stat = followuptrigs[0].statValue**2
        #stat_string = "%0.2f"%stat
          FAR_ul = getattr(followuptrigs[0].coincs,followuptrigs[0].ifolist_in_coinc[0]).alpha

          if not coinc == "H1H2":
            FAR_ul_c = getattr(followuptrigs[0].coincs,followuptrigs[0].ifolist_in_coinc[0]).alpha * float(nbCat)

          candidate_checklist = "http://www.lsc-group.phys.uwm.edu/ligovirgo/cbc/protected/projects/s5/followup/official_candidates/LowMassCBC/20051104-20061114/full_zerolag/" + "followup_" + str(int(gps_time)) + ".html"
          #cumulative_histo = "http://ldas-jobs.ligo.caltech.edu/~romain/LowMassCBC/20051104-20061114/candidates_ethinca_rerun/" + cat_veto.lower() + "_ifar/" + mass_bin + "/" + coinc + "_in_" + ifo_time + "/" + "plotthinca_" + coinc + "_in_" + ifo_time + "_" + coinc + "_cum_hist_effective_snr--.png"
          #histo = "http://ldas-jobs.ligo.caltech.edu/~romain/LowMassCBC/20051104-20061114/candidates_ethinca_rerun/plots/histograms/" + cat_veto.lower() + "_ifar/" + mass_bin + "/" + coinc + "_in_" + ifo_time + "/" + ifo_time + "-plotthinca_hist_effective_snr-unspecified-gpstime.png"
          #ifar = "http://www.lsc-group.phys.uwm.edu/ligovirgo/cbc/protected/projects/s5/lowcbc/20051104-20061114/ethinca_rerun/s5_1yr_run_" + cat_veto.lower() + "_ifar/ifar_upperlimits/" + mass_bin + "/corse/full_data/" + ifo_time + "-plotifar_" + mass_bin.upper() + "-unspecified-gpstime.html"
          ifar = ifar_url + ifo_time + "-plotifar_" + mass_bin.upper() + "_" + cat_veto.upper() + "_cumhist_ifar-815160323-32395247.png"
          if not coinc == "H1H2":
            ifar_c = ifar_url + ifo_time + "-plotifar_ALL_DATA" + "_" + cat_veto.upper() + "_cumhist_ifar_combined-815160323-32395247.png" 

          # decide whether this candidate is red by checking if its GPS matches an element of the list of red times
          try:
            timeList.index(str(int(gps_time)))
            color = "<#FF7744>"
          except:
            color = ""

          # decide whether the candidate is significant and needs to be printed in large font
          if FAR_ul <= far_threshold:
            tuneFont = "~+"
          else:
            tuneFont = ""

          file.write(color + tuneFont + "[" + candidate_checklist + " " + str(int(gps_time)) + "]")
          #file.write(" [" + cumulative_histo + " " + stat_string + "]")
          #file.write(";[" + histo + " histo]")
          file.write(";[" + ifar + " FAR=" + "]" + "%0.3f"%FAR_ul)
          if not coinc == "H1H2":
            file.write(" [" + ifar_c + " FARc=" + "]" + "%0.3f"%FAR_ul_c)
          file.write(tuneFont[::-1])
        else:
          file.write("None")
        file.write("||")
    file.write("\n")
file.close()
