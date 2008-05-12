#!/usr/bin/python
__author__ = "Romain Gouaty"
__version__ = "$Revision$"[11:-2]
__date__ = "$Date$"[7:-2]

import os
import sys
import commands

from optparse import *

usage = """ %prog [options]
"""

parser = OptionParser( usage )

parser.add_option("-v","--version",action="store_true",default=False,\
    help="display version information and exit")

parser.add_option("-f","--mcmc-file",action="store",type="string",\
    metavar=" FILE",help="text file containing the mcmc data used as input")

parser.add_option("-r","--plot-routine",action="store",type="string",\
    metavar=" FILE",help="path to the \"mcmcsummary.R\" script")

parser.add_option("-e","--executable",action="store",type="string",\
    metavar=" FILE",help="path to the R executable")

parser.add_option("-b","--burnin",action="store",type="string",\
    metavar=" VALUE",help="number of mcmc iterations to disregard")

parser.add_option("-t","--reference-time",action="store",type="string",\
    metavar=" GPS",help="GPS time to be used as reference")

parser.add_option("-C","--reference-mchirp",action="store",type="string",\
    metavar=" VALUE",help="VALUE chirp mass to be used as reference")

parser.add_option("-E","--reference-eta",action="store",type="string",\
    metavar=" VALUE",help="VALUE eta to be used as reference")

parser.add_option("-d","--reference-distance",action="store",type="string",\
    metavar=" VALUE",help="VALUE distance to be used as reference")

parser.add_option("-P","--reference-phi",action="store",type="string",\
    metavar=" VALUE",help="VALUE phi to be used as reference")

parser.add_option("-o","--output-path",action="store",type="string",\
    metavar=" PATH",help="use output path PATH")

parser.add_option("-i","--identity",action="store",type="string",\
    metavar=" STRING",help="job identity STRING")

command_line = sys.argv[1:]
(opts,args) = parser.parse_args()

#################################
# if --version flagged
if opts.version:
  print "$Id$"
  sys.exit(0)

#################################
# Sanity check of input arguments

if not opts.mcmc_file:
  print >> sys.stderr, "No mcmc file specified."
  print >> sys.stderr, "Use --mcmc-file FILE to specify location."
  sys.exit(1)

if not opts.plot_routine:
  print >> sys.stderr, "No R routine specified."
  print >> sys.stderr, "Use --plot-routine FILE to specify location."
  sys.exit(1)

if not opts.executable:
  print >> sys.stderr, "No R executable specified."
  print >> sys.stderr, "Use --executable FILE to specify location."
  sys.exit(1)

if not opts.burnin:
  print >> sys.stderr, "No burnin specified."
  print >> sys.stderr, "Use --burnin VALUE to specify it."
  sys.exit(1)

if not opts.reference_time:
  print >> sys.stderr, "No reference time specified."
  print >> sys.stderr, "Use --reference-time GPS to specify location."
  sys.exit(1)

if not opts.reference_mchirp:
  print >> sys.stderr, "No reference mchirp specified."
  print >> sys.stderr, "Use --reference-mchirp VALUE to specify location."
  sys.exit(1)

if not opts.reference_eta:
  print >> sys.stderr, "No reference eta specified."
  print >> sys.stderr, "Use --reference-eta VALUE to specify location."
  sys.exit(1)

if not opts.reference_distance:
  print >> sys.stderr, "No reference distance specified."
  print >> sys.stderr, "Use --reference-distance VALUE to specify location."
  sys.exit(1)

if not opts.reference_phi:
  print >> sys.stderr, "No reference phi specified."
  print >> sys.stderr, "Use --reference-phi VALUE to specify location."
  sys.exit(1)

if not opts.output_path:
  print >> sys.stderr, "No output path specified."
  print >> sys.stderr, "Use --output-path PATH to specify location."
  sys.exit(1)

if not opts.identity:
  print >> sys.stderr, "No identity specified."
  print >> sys.stderr, "Use --identity STRING to specify it."
  sys.exit(1)

#################################
# Main program

file = open(opts.output_path + "/" + opts.identity + ".R",'w')

file.write("# load the \"mcmcsummary\" code:\n")
file.write("source(\"" + opts.plot_routine + "\")\n\n")

file.write("# load the data written in the MCMC txt file\n")
file.write("input <- read.table(\"" + opts.mcmc_file + "\",header=TRUE)\n\n")

file.write("# prepare the input data for the plotting routine\n")
file.write("dataset = cbind(input[,1],input[,4:8],input[,12:16],input[,20:24],input[,28:32],input[,36:40],input[,44:48])\n\n")

file.write("# enter injected or inspiral parameters\n")
file.write("injpar <- c(\"mc\"=" + opts.reference_mchirp + ",\"eta\"=" + opts.reference_eta + ",\"tc\"=" + opts.reference_time + ",\"phi\"=" + opts.reference_phi + ",\"dl\"=" + opts.reference_distance + ")\n")
file.write("trueparList = c(injpar,injpar,injpar,injpar,injpar,injpar)\n\n")

file.write("# execute the \"mcmcsummary\" code:\n")
file.write("mcmcsummary(data=dataset[,-1],targetdirectory=\"" + opts.output_path + "/" + opts.identity + "\",iteration=dataset[,1],burnin=" + opts.burnin + ",varnames=colnames(dataset[,-1]),truevalue=trueparList,graphicsformats = c(\"png\"), overwrite = T)")

file.close()

command = opts.executable + " --slave --vanilla --file=" + opts.output_path + "/" + opts.identity + ".R"

result = commands.getstatusoutput(command)
#result = commands.getoutput(command)

if not result[0] == 0:
  print sys.stderr, "status=" + str(result[0])
  print sys.stderr, result[1]
  sys.exit(1)

