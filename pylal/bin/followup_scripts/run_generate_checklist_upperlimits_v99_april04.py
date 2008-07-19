#!/usr/bin/env python

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

###########################################################################
# What is needed to run the script ?
###########################################################################

# You will need the executable generate_checklist.py which is called by this routine.
# You need an installation of GLUE, LAL, LALAPPS, PYLAL

###########################################################################
# Inputs to be specified
###########################################################################

# provide the list of mass bins (example: mass_bins = ['mchirp_2_8','mchirp_8_17','mchirp_17_35'])
mass_bins = ['mchirp_2_8','mchirp_8_17','mchirp_17_35']
# provide the list of DQ flags (example: cat_vetoes = ['cat12'])
cat_vetoes = ['cat1','cat12','cat123','cat1234']
# coincs = ['H1L1','H2L1']
coincs = ['H2L1']

# Specify the path to the executable
executable = "./generate_checklist.py"

# path to the COIRE files containing the list of the loudest candidates
input_path = "/home/romain/S5_followup/1stCalendarYear/fulldata/LowMassCBC_triggers/upperlimits_v99_april04/"
# link the automated followup web pages
automated_page = "http://www.lsc-group.phys.uwm.edu/ligovirgo/cbc/protected/projects/s5/followup/LowMassCBC/20051104-20061114/full_zerolag/upperlimits_v99_april04/"
# link to the web page containing calculated False Alarm Probabilities 
stat_page = "http://www.lsc-group.phys.uwm.edu/ligovirgo/cbc/protected/projects/s5/followup/official_candidates/LowMassCBC/20051104-20061114/full_zerolag/"
statpage_string = "upperlimits_v99_april04"

# link to cumulative histograms
#cumul_histo_page = "http://ldas-jobs.ligo.caltech.edu/~romain/LowMassCBC/20051104-20061114/candidates_ethinca_rerun/"
cumul_histo_page = "http://www.lsc-group.phys.uwm.edu/ligovirgo/cbc/protected/people/dkeppel/s5_1yr_lowcbc_20051104-20061114/box/"

# link to IFAR plots
#ifar_page = "http://www.lsc-group.phys.uwm.edu/ligovirgo/cbc/protected/projects/s5/lowcbc/20051104-20061114/ethinca_rerun/"
ifar_page = "http://ldas-jobs.ligo.caltech.edu/~dkeppel/results/s5_1yr_lowcbc/upperlimits_v99/corse_all_data/"

# link to qscan webpages
qscan_page = "http://ldas-jobs.ligo.caltech.edu/~romain/S5_followup/LowMassCBC/20051104-20061114/full_zerolag/"

# Number of coincident triggers to follow-up in each COIRE file
nb_triggers = 1

############################################################################
# Main program
############################################################################

for cat_veto in cat_vetoes:
  for mass_bin in mass_bins:
    for coinc in coincs:
      if coinc == "H1H2L1":
        ifo_times = ["H1H2L1"]
      if coinc == "H1L1":
        ifo_times = [["H1H2L1","triples"],["H1L1","doubles"]]
      if coinc == "H2L1":
        ifo_times = [["H1H2L1","triples"],["H2L1","doubles"]]

      for ifo_time in ifo_times:

        command = "python " + executable + \
        " --xml-file " + input_path + mass_bin + "/full_data_" + cat_veto + "/" + coinc + "/" + coinc + "-CORSE_" + ifo_time[0] + "-loudest.xml" \
        + " --num-triggers " + str(nb_triggers) \
        + " --statistic effective_snrsq " + \
        "--automated-page " + automated_page + coinc + "/" + mass_bin + "/" + ifo_time[1] + "_" + cat_veto + "/" + \
        " --statistic-page " + stat_page +  "candidates_summary_" + coinc + "_" + ifo_time[0] + "_" + statpage_string + "_" + mass_bin + "_" + cat_veto + ".html" + \
        " --cumulhisto-page " + cumul_histo_page + mass_bin + "_" + cat_veto + "_final_v99/plots/" + coinc + "/" + coinc + "_in_" + ifo_time[0] + "_cum_hist_effective_snr.png" + \
        " --ifar-page " + ifar_page + mass_bin + "_" + cat_veto + "/Images/" + ifo_time[0] + "-plotifar_cumhist_ifar-unspecified-gpstime.png" + \
        " --qscan-page " + qscan_page + \
        " --string-id " + cat_veto + "_" + mass_bin.split("_")[1] + "-" + mass_bin.split("_")[2]
        #" --cumulhisto-page " + cumul_histo_page + cat_veto + "_ifar/" + mass_bin  + "/" + coinc + "_in_" + ifo_time[0] + "/plotthinca_" + coinc + "_in_" + ifo_time[0] + "_" + coinc + "_cum_hist_effective_snr--.png" + \
        #" --histo-page " + cumul_histo_page + "plots/histograms/" + cat_veto + "_ifar/" + mass_bin  + "/" + coinc + "_in_" + ifo_time[0] + "/" + ifo_time[0] + "-plotthinca_hist_effective_snr-unspecified-gpstime.png" + \
        #" --ifar-page " + ifar_page + "/s5_1yr_run_" + cat_veto + "_ifar/ifar_upperlimits/" + mass_bin + "/corse/full_data/" + ifo_time[0] + "-plotifar_" + mass_bin.upper() + "-unspecified-gpstime.html" + \

        result = commands.getoutput(command)
        print result

