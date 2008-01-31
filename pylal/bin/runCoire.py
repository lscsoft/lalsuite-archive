#!/usr/bin/python

import sys
import os
import copy
import shutil
from optparse import *
import exceptions
import glob
import ConfigParser
import string

#######################################################################
usage = """
usage: %prog [options] 

Script for running the coire steps on the GRB injections.

"""

parser = OptionParser( usage )
#
# a c i p n s u 
#

# path
parser.add_option("-n","--number",action="store",type="int",\
    default=0, metavar=" NUMBER", help="number of injections")

parser.add_option("-s","--start",action="store",type="int",\
    default=1, metavar=" START", help="start number of injections")

parser.add_option("-S","--second",action="store_true",\
    default=False, help="running coire over the second THINCA files")

parser.add_option("-I","--inspiral",action="store_true",\
    default=False, help="running sire on the single IFO inspiral files")

parser.add_option("-d","--dataset",action="store",type="string",\
    default="H1H2L1", metavar=" DATASET", help="what dataset to analyze")

parser.add_option("-T","--test",action="store_true",\
    default=False, help="only print the commands to be run")

# mass range options
parser.add_option("-t","--data-type",action="store",type="string",\
    default="all_data", metavar=" DATA_TYPE",help="data type") 

parser.add_option("-w","--injection-window",action="store",type="float",\
    default=15,metavar=" INJWINDOW", help="injection window [10]" )

parser.add_option("-c","--cluster-time",action="store",type="float",\
    default=30, metavar=" CLUSTERTIME", help=" cluster-time [20]. If set to zero, no clustering will be used." )

parser.add_option("-u","--user-tag",action="store",type="string",\
    default="", metavar=" USERTAG",help=" user tag" )

parser.add_option("-o","--output-tag",action="store",type="string",\
    default="", metavar=" OUTPUTTAG", help=" output tag" )

parser.add_option("-C","--condor",action="store_true",\
    default=False, help=" to run it as condor job" )

parser.add_option("-A","--all",action="store_true",\
    default=False, help=" use all combinations" )

#######################################################################
# check options and initialize
#######################################################################
(opts,args)=parser.parse_args()

# two datasets for now
dataSets=['H1H2L1','H1L1']

number=opts.number
if number==0:
  number=opts.start

# create name for the output-files
cfile="COIRE-"
if opts.second:
  cfile="COIRE2-"

for inj in range(opts.start,number+1):
  
  if opts.inspiral:

     cfile="-SIRE-"
     for ifo in ['H1','H2','L1']:
       command = "./lalapps_sire --glob '" + ifo + "-INSPIRAL"

       # use thinca from the second coinc step?
       if opts.second:
          command+="_"+set
       if opts.output_tag:
         outputPart=ifo +cfile+opts.output_tag+'-'+str(inj)
       else:
         outputPart=ifo +cfile+str(inj)

       command+="_" +opts.user_tag+str(inj) + "-*.xml'"+\
              " --output  "+ outputPart+"_FOUND.xml"+\
              " --missed-injections  "+ outputPart+"_MISSED.xml"+\
              " --summary-file  "+ outputPart+".txt"+\
              " --injection-file HL-INJECTION."+str(inj)+".xml"+\
              " --injection-window "+ str(opts.injection_window) +\
              " --data-type " +opts.data_type
       if opts.cluster_time>0:
         command+=" --cluster-time "+str(opts.cluster_time)+" --cluster-algorithm snr"

       # should not use this, condor creates some problems...
       if opts.condor:
          command="condor_run \""+command+"\"&"
       if opts.test:
          print command
       else:
          os.system(command)


  else:
    dataSets=[opts.dataset]
   
    for set in dataSets:
       if set=="H1H2L1":
         coincs=['H1H2L1','H1H2','H1L1','H2L1']
       else:
         coincs=[set]

       if opts.all:
         coincs=['ALL']

       for coinc in coincs:
          command = "./lalapps_coire --glob '" + set + "-THINCA"

          # use thinca from the second coinc step?
          if opts.second:
            command+="_"+set
          if opts.output_tag:
            outputPart=coinc +"-H1H2L1-"+cfile+opts.output_tag+'-'+str(inj)
          else:
            outputPart=coinc +"-H1H2L1-"+cfile+str(inj)

          command+="_" +opts.user_tag+str(inj) + "-*.xml'"+\
              " --output  "+outputPart+"_FOUND.xml"+\
              " --missed-injections  "+ outputPart+"_MISSED.xml"+\
              " --summary-file  "+ outputPart+".txt"+\
              " --injection-file HL-INJECTION."+str(inj)+".xml"+\
              " --injection-window "+ str(opts.injection_window) +\
              " --data-type " +opts.data_type
          if opts.cluster_time>0:
            command+=" --cluster-time "+str(opts.cluster_time)+" --coinc-stat snrsq"
          if not opts.all:
              command+=" --coinc-cut "+coinc

          # should not use this, condor creates some problems...
          if opts.condor:
            command="condor_run \""+command+"\"&"

          if opts.test:
            print command
          else:
            os.system(command)

