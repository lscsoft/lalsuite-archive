#!/usr/bin/env python
"""
$Id$

This script is used to set up a correct followup_pipe.ini file in each of the subdirectories where the followup dags will be submitted. Since there might be 84 categories of candidates, the followup is organised in a tree of directories where each of these categories of candidates are followed up. THIS SCRIPT MUST BE RUN IN THE DIRECTORY WHERE YOU WANT TO SET UP THE FOLLOWUP TREE.
"""

__author__ = 'Romain Gouaty <romain@phys.lsu.edu>'
__date__ = '$Date$'
__version__ = '$Revision$'[11:-2]

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

######################## OPTION PARSING  #####################################
usage = """usage:
This script is used to set up a correct followup_pipe.ini file in each of the subdirectories where the followup dags will be submitted. Since there might be 84 categories of candidates, the followup is organised in a tree of directories where each of these categories of candidates are followed up. THIS SCRIPT MUST BE RUN IN THE DIRECTORY WHERE YOU WANT TO SET UP THE FOLLOWUP TREE.
"""

parser = OptionParser( usage )

parser.add_option("-v", "--version",action="store_true",default=False,\
    help="print version information and exit")

parser.add_option("-f","--input-file",action="store",type="string",\
    metavar=" FILE",help="The input file should be an instance of the followup configuration file (followup_pipe.ini). The .ini file provided as input should already be properly set up for one example of category of candidates. The goal of this program is only to propagate the .ini file through the followup directory tree and to update appropriately the [hipe-cache], [triggers] and [ouput] sections")

parser.add_option("-x","--xml-glob-template",action="store",type="string",\
    metavar=" FILE",help="This argument is used to specify the general format of the field \"xml-glob\" of the [triggers] section. For example this could be: \"/archive/home/dkeppel/post_processing/lowcbc/20051104-20061114/upperlimits_v99/corse_dag/corse/MASSBIN/full_data_CAT/COINC/COINC-CORSE_IFOTIME-815160323-32395247.xml.gz\". Note that the script will interprete the keywords as follows: MASSBIN is one of [\"mchirp_2_8\",\"mchirp_8_17\",\"mchirp_17_35\"], CAT is one of [\"cat1\",\"cat12\",\"cat123\",\"cat1234\"], COINC and IFOTIME can be one of [\"H1H2L1\",\"H1L1\",\"H2L1\",\"H1H2\"]")

parser.add_option("-i","--hipe-path-template",action="store",type="string",\
    metavar=" FILE",help="This argument specifies where the hipe intermediate xml files can be found. This will be used to fill the field \"second-inspiral-path\" of the section [hipe-cache] (the other fields \"tmpltbank-path\", \"trigbank-path\", ... will be filled also but they are not currently relevant). An example could be \"/archive/home/dkeppel/analysis/lowcbc/20051104-20061114/full_data/slides_CAT_final_v99/\" where the key word CAT will be replaced within the script by one of [\'cat1\',\'cat12\',\'cat123\',\'cat1234\']")

parser.add_option("-o","--output-template",action="store",type="string",\
    metavar=" FILE",help="This argument specifies where the followup results will be written (used to fill the field \"page\" of the section [output]). The default that we will use is something like: \"/ligovirgo/cbc/protected/projects/s5/followup/LowMassCBC/20051104-20061114/full_zerolag/upperlimits_v99_july08/COINC/MASSBIN/TIME_CAT/\", where the key words are: COINC (the code will replace it by one of [\'H1H2L1\',\'H1L1\',\'H2L1\',\'H1H2\'], MASSBIN is one of [\'mchirp_2_8\',\'mchirp_8_17\',\'mchirp_17_35\'], TIME is one of [\'doubles\',\'triples\'], CAT is one of [\'cat1\',\'cat12\',\'cat123\',\'cat1234\'])")

command_line = sys.argv[1:]
(opts,args) = parser.parse_args()

#################################
# if --version flagged
if opts.version:
  print "$Id: generate_checklist.py, v 1.0 2008/05/20 07:00:00 romain Exp"
  sys.exit(0)

##############################################################################
# main program

iniOrig = ConfigParser.ConfigParser()
iniOrig.read(opts.input_file)

for coinc in ["H1H2L1","H1L1","H2L1","H1H2"]:

  if coinc == "H1H2L1":
    timeList = [["H1H2L1","triples"]]
  if coinc == "H1L1":
    timeList = [["H1L1","triples"],["H1L1","doubles"]]
  if coinc == "H2L1":
    timeList = [["H2L1","triples"],["H2L1","doubles"]]
  if coinc == "H1H2":
    timeList = [["H1H2","triples"],["H1H2","doubles"]]

  xmlfile = opts.xml_glob_template.replace("COINC",coinc)
  outputdir = opts.output_template.replace("COINC",coinc)

  for massBin in ["mchirp_2_8","mchirp_8_17","mchirp_17_35"]:
    xmlfile1 = xmlfile.replace("MASSBIN",massBin)
    outputdir1 = outputdir.replace("MASSBIN",massBin)

    for time in timeList:
      xmlfile2 = xmlfile1.replace("IFOTIME",time[0])
      outputdir2 = outputdir1.replace("TIME",time[1])

      for cat in ["cat1","cat12","cat123","cat1234"]:
        xmlfile3 = xmlfile2.replace("CAT",cat)
        outputdir3 = outputdir2.replace("CAT",cat)
        hipefile = opts.hipe_path_template.replace("CAT",cat)

        iniOrig.set('triggers','xml-glob',xmlfile3)
        iniOrig.set('output','page',outputdir3)
        iniOrig.set('hipe-cache','tmpltbank-path',hipefile)
        iniOrig.set('hipe-cache','trigbank-path',hipefile)
        iniOrig.set('hipe-cache','first-inspiral-path',hipefile)
        iniOrig.set('hipe-cache','second-inspiral-path',hipefile)
        iniOrig.set('hipe-cache','first-coinc-path',hipefile)
        iniOrig.set('hipe-cache','second-coinc-path',hipefile)
        
        relativePath = coinc + "/" + massBin + "/" + time[1] + "_" + cat + "/"

        fp = file(relativePath + "followup_pipe.ini","w")
        iniOrig.write(fp)
        fp.close()

