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
# provide the list of coincidence types to analyze (example: coincs = ['H1L1','H2L1'])
coincs = ['H1H2L1','H1L1','H2L1','H1H2']

# Specify the path to the executable
executable = "/archive/home/romain/cvs/pylal/bin/followup_scripts/generate_checklist.py"

# path to the COIRE or CORSE files containing the list of the loudest candidates
input_path = "/archive/home/dkeppel/post_processing/lowcbc/20051104-20061114/upperlimits_v99/corse_dag/corse/"

# link the automated followup web pages
automated_page = "http://www.lsc-group.phys.uwm.edu/ligovirgo/cbc/protected/projects/s5/followup/LowMassCBC/20051104-20061114/full_zerolag/upperlimits_v99_july08/" 

# link to the web page containing calculated False Alarm Probabilities 
stat_page = "http://www.lsc-group.phys.uwm.edu/ligovirgo/cbc/protected/projects/s5/followup/official_candidates/LowMassCBC/20051104-20061114/full_zerolag/"
statpage_string = "upperlimits_v99_july08"

# link to cumulative histograms
#cumul_histo_page = "http://ldas-jobs.ligo.caltech.edu/~romain/LowMassCBC/20051104-20061114/candidates_ethinca_rerun/"
cumul_histo_page = "http://ldas-jobs.ligo.caltech.edu/~dkeppel/results/s5_1yr_lowcbc/upperlimits_v99/corse_dag/"

# link to IFAR plots
#ifar_page = "http://www.lsc-group.phys.uwm.edu/ligovirgo/cbc/protected/projects/s5/lowcbc/20051104-20061114/ethinca_rerun/"
ifar_page = "http://ldas-jobs.ligo.caltech.edu/~dkeppel/results/s5_1yr_lowcbc/upperlimits_v99/corse_dag/Images/"

# link to qscan webpages
qscan_page = "http://ldas-jobs.ligo.caltech.edu/~romain/S5_followup/LowMassCBC/20051104-20061114/full_zerolag/"

# Link to qscan output directory (used to get the table of DQ flags and Science segments)
qscan_dir = "/archive/home/romain/public_html/S5_followup/LowMassCBC/20051104-20061114/full_zerolag/"

# Number of coincident triggers to follow-up in each COIRE file
nb_triggers = 3

############################################################################
# Main program
############################################################################

#loop over the DQ veto categories
for cat_veto in cat_vetoes:
  # loop over the mass bins
  for mass_bin in mass_bins:
    # loop over the coincidence types
    for coinc in coincs:
      if coinc == "H1H2L1":
        ifo_times = [["H1H2L1","triples"]]
      if coinc == "H1L1":
        ifo_times = [["H1H2L1","triples"],["H1L1","doubles"]]
      if coinc == "H2L1":
        ifo_times = [["H1H2L1","triples"],["H2L1","doubles"]]

      for ifo_time in ifo_times:

        command = "python " + executable + \
        " --xml-file " + input_path + mass_bin + "/full_data_" + cat_veto + "/" + coinc + "/" + coinc + "-CORSE_" + ifo_time[0] + "-815160323-32395247.xml.gz" \
        + " --num-triggers " + str(nb_triggers) \
        + " --statistic effective_snrsq " + \
        "--automated-page " + automated_page + coinc + "/" + mass_bin + "/" + ifo_time[1] + "_" + cat_veto + "/" + \
        " --statistic-page " + stat_page +  "candidates_summary_" + coinc + "_" + ifo_time[0] + "_" + statpage_string + "_" + mass_bin + ".html" + \
        " --cumulhisto-page " + cumul_histo_page + mass_bin + "/plots_" + cat_veto + "/Images/" + ifo_time[0] + "-plotthinca_" + coinc + "_" + coinc + "_cum_hist_effective_snr-unspecified-gpstime_thumb.png" + \
        " --histo-page " + cumul_histo_page + mass_bin + "/plots_" + cat_veto + "/Images/" + ifo_time[0] + "-plotthinca_" + coinc + "_" + "hist_effective_snr-unspecified-gpstime_thumb.png" + \
        " --ifar-page " + ifar_page + ifo_time[0] + "-plotifar_" + mass_bin.upper() + "_" + cat_veto.upper() + "_cumhist_ifar-815160323-32395247.png" + \
        " --ifar-combined-page " + ifar_page + ifo_time[0] + "-plotifar_ALL_DATA_" + cat_veto.upper() + "_cumhist_ifar_combined-815160323-32395247.png" + \
        " --qscan-page " + qscan_page + \
        " --string-id " + cat_veto + "_" + mass_bin.split("_")[1] + "-" + mass_bin.split("_")[2] + \
        " --qscan-dir " + qscan_dir

        result = commands.getoutput(command)
        print result

