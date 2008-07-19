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

# "inputDir" indicates the main part of the path where the trigger files are written. These trigger files (xml) contain the list of candidates to analyze.
inputDir = "/archive/home/dkeppel/post_processing/lowcbc/20051104-20061114/upperlimits_v99/corse_dag/"

# lists the different veto categories
cat_veto_list = ["Cat123"]

# list of ifo times being analyzed
ifo_times = ["H1H2L1","H1L1","H2L1"]

# Name of the outputFile (containing the candidates' table in wiki format)
outputFile = "upperlimits_v99_July2008_candidate_table_threeloudest.txt"

# page containing the IFAR plots
ifar_url = "http://ldas-jobs.ligo.caltech.edu/~dkeppel/results/s5_1yr_lowcbc/upperlimits_v99/corse_dag/Images/"

# input text file which contains a list of GPS times (integers). This must be the list of trigger times which are ruled out by the detection checklist. The program use this list of times to decide which cells of the candidates' table need to be put in red.
listOfRedTimes = "./candidate_redlist.txt"

# threshold on the FAR value used to determine whether the candidate is significant and needs to be printed in large font
far_threshold = 0.10

# number of candidates to followup in each ifo_time
nb_loudest = 3

# first row of the wiki table (lists all the ifo times)
tableHeader1 = "|||| Time:  || H1H2L1 || H1L1 || H2L1 ||"


##############################################################################
#
#  MAIN PROGRAM
#
##############################################################################

##########################################################################
# First we read the trigger files and store the sorted list of triggers in
# variables
##########################################################################

# loop over the list of veto categories.
for cat_veto in cat_veto_list:
    
  # loop over the ifo_times
  for ifo_time in ifo_times:

    # define path and name of trigger file to read
    inputFile = inputDir + "/" + ifo_time + "-CORSE_" + cat_veto.upper() + "-815160323-32395247.xml.gz"
    # read the trigger file with the "readFiles" method of fu_utils.py (PYLAL)
    found, coincs, search = readFiles(inputFile,getstatistic("effective_snrsq",None,None))

    # sort the list of triggers ("coincs") by FAR
    # "coinc_trigger.get_ifos()[1][0]" is a string corresponding to the name of the first interferometer of the coincidence
    exec("coincs_list_" + ifo_time + "_" + cat_veto.upper() + " = [ (getattr(coinc_trigger,coinc_trigger.get_ifos()[1][0]).alpha, coinc_trigger) for coinc_trigger in coincs ]")
    eval("coincs_list_" + ifo_time + "_" + cat_veto.upper() + ".sort()")


##########################################################################
# We will now build the wiki table of candidates
##########################################################################

# call the method "listFromFile" from the pylal module "fu_utils.py" to read the text file "listOfRedTimes" (this file contains the list of ruled out candidate times)
timeList = listFromFile(listOfRedTimes)

# prepare the file in which the table of candidates will be written
file = open(outputFile,'w')
file.write(tableHeader1 + "\n")

# loop over the list of veto categories. The table of candidates contain one row per category.
for cat_veto in cat_veto_list:
  file.write("||||<:> " + cat_veto + "||<-" + str(len(ifo_times)) + ">||\n")

  # loop over the "nb_loudest" candidates
  for k in range(nb_loudest):

    file.write("||||<:>||")

    # Loop over the ifo_times. For each ifo_times we need to specify how to normalize the FAR by giving the number of category in each ifo_time ("nbCat" variable)
    for ifo_time in ifo_times:
      if ifo_time == "H1H2L1":
        nbCat = 9
      if ifo_time == "H1L1":
        nbCat = 3
      if ifo_time == "H2L1":
        nbCat = 3

      coincs_list = eval("coincs_list_" + ifo_time + "_" + cat_veto.upper())

      # if the element "coincs_list[k]" does not exist, we should just report "None" in this cell of the table (there are no candidates of this type). This is the purpose of the following try
      try: 
        coincs_list[k]

        # get the time of the candidate (integer part) in the first ifo
        gps_time = getattr(coincs_list[k][1],coincs_list[k][1].get_ifos()[1][0]).end_time

        # get the FAR of the candidate
        FAR_ul   = getattr(coincs_list[k][1],coincs_list[k][1].get_ifos()[1][0]).alpha
        FAR_ul_c = FAR_ul*nbCat

        # define weblink to the checklist of this candidate (assuming it exists)
        candidate_checklist = "http://www.lsc-group.phys.uwm.edu/ligovirgo/cbc/protected/projects/s5/followup/official_candidates/LowMassCBC/20051104-20061114/full_zerolag/" + "followup_" + str(int(gps_time)) + ".html"
        # defibe weblink to the corresponding IFAR plot 
        ifar_c = ifar_url + ifo_time + "-plotifar_ALL_DATA" + "_" + cat_veto.upper() + "_cumhist_ifar_combined-815160323-32395247.png" 

        # decide whether this candidate is red by checking if its GPS matches an element of the list of red times (red is supposed to mean ruled out)
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

        # Now fill one row of the table of candidates using the information obtained previously
        file.write(color + tuneFont + "[" + candidate_checklist + " " + str(int(gps_time)) + "]")

        file.write(" [" + ifar_c + " FARc=" + "]" + "%0.3f"%FAR_ul_c)
        file.write(tuneFont[::-1])

      except:
        file.write("None")
      file.write("||")
    file.write("\n")
    k+=1
file.close()
