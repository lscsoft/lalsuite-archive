#!/usr/bin/env python

# This script launches the python program "buildCandidateList_CORSE.py" to generate an html page containing the list of loudest triggers for a given mass bin, a given coincidence type, and a given ifo time. The html page contains information about the inspiral parameters of the triggers as well as their FAR in each of the DQ vetoes categories.

# import required modules
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
#INPUT VARIABLES
##############################################################################

# give the main input path to the xml trigger files
#input = "/archive/home/dkeppel/post_processing/lowcbc/20051104-20061114/upperlimits_v99/corse_dag/corse/"
input = "/home/romain/S5_followup/1stCalendarYear/fulldata/LowMassCBC_triggers/upperlimits_v99_april04/"

# give the list of ifo times
timeList = ["H1H2L1","H1L1","H2L1","H1H2"]

# give the list of mass bins (the strings must match the real directory names in the path to the xml trigger files)
massbins = ["mchirp_2_8","mchirp_8_17","mchirp_17_35"]

# indicate the normalization constant to be used to compute the combined FAR from the FAR: "combined FAR = FAR * far_normalization_???" 
# "far_normalization_???" should be the product of len(massbins) by the number of types of coincidences in a given ifo time category.
# For example for a trigger found in triple times, there are 3 mass bins and 3 coincidence types (H1H2L1, H1L1, H2L1), which gives "far_normalization_triples=9". For a trigger found in double times: "far_normalization_doubles=3".
# if the trigger is a H1H2 coincidence the argument "far_normalization_???" will be ignored, as we don't compute the combined FAR for H1H2 coincidences.
far_normalization_triples = 9
far_normalization_doubles = 3

# Give the number of candidates you want to be displayed in the html page
nb_candidates = "30"

# This is a string to be used in the naming of the output file
#userString = "upperlimits_v99_july08"
userString = "upperlimits_v99_april04"

# give the output path where the html page will be generated
#outputPath = "/archive/home/romain/Projects/LowMassCBC/20051104-20061114/triggers/upperlimits_v99_july08/triggerFiles/"
outputPath = "/home/romain/python_scripts/"

# give the list of DQ veto categories: use strings that match the directories containing the xml trigger files
categories_dir = ["full_data_cat1","full_data_cat12","full_data_cat123","full_data_cat1234"]
# list of DQ veto categories. These strings are used to set up option names for the program "buildCandidateList_CORSE.py"
categories = ["cat1","cat12","cat123","cat1234"]


##############################################################################
#MAIN PROGRAM
##############################################################################

for ifo_time in timeList:
  if ifo_time == "H1H2L1":
    coincList = ["H1H2L1","H1L1","H2L1","H1H2"]
    far_normalization = far_normalization_triples
  if ifo_time == "H1L1":
    coincList = ["H1L1"]
    far_normalization = far_normalization_doubles
  if ifo_time == "H2L1":
    coincList = ["H2L1"]
    far_normalization = far_normalization_doubles
  if ifo_time == "H1H2":
    coincList = ["H1H2"]
    far_normalization = far_normalization_doubles

  for ifo_coinc in coincList:

    #zerolagFile = ifo_coinc + "-CORSE_" + ifo_time + "-815160323-32395247.xml.gz"
    zerolagFile = ifo_coinc + "-CORSE_" + ifo_time + "-loudest.xml"

    for massbin in massbins:
      outputFile = "candidates_summary_" + ifo_coinc + "_" + ifo_time + "_" + userString + "_" + massbin + ".html"
      title = "List_of_triggers_in_" + massbin + "," + ifo_coinc + "_coincidence," + ifo_time + "_times"

      zerolag_files = ""
      for j,cat in enumerate(categories):
        zerolag_files = zerolag_files + "--zerolag-" + cat + " " + input + massbin + "/" + categories_dir[j] + "/" + ifo_coinc + "/" + zerolagFile + " "

      command = "python buildCandidateList_CORSE.py " + \
          "-l ./. " + \
          "-o " + outputPath + " " + \
          "--header " + title + " " + \
          "--page http://www.lsc-group.phys.uwm.edu/ligovirgo/cbc/protected/projects/s5/followup/official_candidates/LowMassCBC/20051104-20061114/full_zerolag/ " + \
          "--web-output-file " + outputFile + " " + \
          "-s effective_snrsq " + \
          "-n " + nb_candidates + " " + \
          zerolag_files + " " + \
          "--far-normalization " + str(far_normalization)

      print "we are running the following command: "
      print command
      commands.getoutput(command)
      print "done"

