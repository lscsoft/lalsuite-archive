#!/usr/bin/python

import sys
import os
import copy
import shutil
from optparse import *
import exceptions
import glob
import ConfigParser

#######################################################################
usage = """
usage: %prog [options] 

Script for getting an upper limit from the end of the pipeline.  It
uses hipecoire, plotnumgalaxies, lalapps_compute_posterior and some
other stuff to properly do the upper limit.  

There are options to skip each step if it has already been done. This
is very useful when dealing with something like the BBH search where
the number of triggers is very large.

"""

parser = OptionParser( usage )
#
# a c i p n s u 
#

# path
parser.add_option("-a","--trig-path",action="store",type="string",\
    default='', metavar=" DIR",help="directory containing the triggers" )

parser.add_option("-C","--config-file",action="store",type="string",\
    default='', metavar=" INI File",\
    help="ini file with information about run directories" )

parser.add_option("-R","--results-dir",action="store",type="string",\
    default='', metavar=" RESULTS_DIR", help="directory for the results" )

# general options
parser.add_option("-s","--skip-setup",action="store_true",default=False,\
    help="skip setup stage, i.e. mkdirs, copy files, run inspinj")

parser.add_option("-c","--skip-coiredata",action="store_true",default=False,\
    help="skip coiring of data")

parser.add_option("-p","--skip-population",action="store_true",default=False,\
    help="skip population generation")

parser.add_option("-i","--skip-coireinj",action="store_true",default=False,\
    help="skip coiring of injections")

parser.add_option("-n","--skip-png",action="store_true",default=False,\
    help="skip plotnumgalaxies")

parser.add_option("-u","--skip-upper-limit",action="store_true",default=False,\
    help="skip upper limit")

parser.add_option("-t","--usertag",action="store",type="string",\
    default=None, metavar=" USERTAG", help="the user tag used for injections")

parser.add_option("-T","--test",action="store_true",\
    default=False, help="only print the commands to be run")

# mass range options
parser.add_option("-d","--min-mass",action="store",type="float",\
    default=None, metavar=" MIN_MASS",help="minimum mass covered by the search" )

parser.add_option("-e","--max-mass",action="store",type="float",\
    default=None, metavar=" MAX_MASS",help="maximum mass covered by the search" )

parser.add_option("-f","--mean-mass",action="store",type="string",\
    default=None, metavar=" MEAN_MASS",help="mean mass for Gaussian distribution" )

parser.add_option("-g","--stdev-mass",action="store",type="string",\
    default=None, metavar=" STD_MASS",help="std mass for Gaussian distribution" )

parser.add_option("-m","--m-dm",action="store",type="string",\
    default=None, metavar=" DM",help="mass interval for rate v mass" )

# hierarchy pipeline
parser.add_option("-j","--second-coinc",action="store_true",default=False,\
    help="use the second stage thinca files as input")

# rsqr parameters
parser.add_option("-k","--rsq-threshold",action="store",type="string",\
    default=None, metavar=" RSQ_THRESH", help="rsq threshold")

parser.add_option("-l","--rsq-max-snr",action="store",type="string",\
    default=None, metavar=" RSQ_MAX_SNR", help="rsq max snr")

(opts,args) = parser.parse_args()

#######################################################################
# check options and initialize
#######################################################################

if not (opts.min_mass and opts.max_mass and opts.mean_mass and opts.stdev_mass):
  print >> sys.stderr, "Must specify all mass options"
  sys.exit(1)

if (opts.rsq_threshold or opts.rsq_max_snr):
  if not (opts.rsq_threshold and opts.rsq_max_snr):
    print >> sys.stderr, "Must specify both rsq options"
    sys.exit(1)

if not opts.config_file:
  print >> sys.stderr, "Must specify --config-file"
  sys.exit(1)

if not opts.results_dir:
  print >> sys.stderr, "Must specify --results-dir"
  sys.exit(1)

##############################################################################
# create the config parser object and read in the ini file
cp = ConfigParser.ConfigParser()
cp.read(opts.config_file)
fulldata = cp.items("fulldata")
injdata = cp.items("injdata")
vetodata = cp.items("vetodata")
sourcedata = cp.items("sourcedata")

#S4DIR="/home/htdocs/uwmlsc/root/iulgroup/investigations/s4"
#SEARCHDIR="bns"
#LALAPPS_LOCATION=getenv("LALAPPS_LOCATION")
#USERTAG=""

#######################################################################
# do the set up
#######################################################################
MYRESULTSDIR=opts.results_dir
minmass = str(opts.min_mass)
maxmass = str(opts.max_mass)
minmtotal = str(2.0*opts.min_mass)
maxmtotal = str(2.0*opts.max_mass)

if not opts.skip_setup:
  os.mkdir(MYRESULTSDIR)
  os.chdir(MYRESULTSDIR)

  os.symlink(fulldata[0][1], "./full_data")
  for injdir in injdata:
    os.symlink(injdir[1], "./" + injdir[0])
  for myfile in vetodata:
    shutil.copy(myfile[1], "./")
  for myfile in sourcedata:
    shutil.copy(myfile[1], "./inspsrcs.new")

  os.chdir(MYRESULTSDIR)

  # write out the options that were supplied
  myfile = open("upperlimit.log","w")
  myfile.write( str(opts) )
  myfile.close()

#######################################################################
# generate population
#######################################################################
if not opts.skip_population:
  print "** Generating the population with Gaussian mass distribution **"
  command = "lalapps_inspinj \
    --source-file inspsrcs.new \
    --gps-start-time 793130413 --gps-end-time 795679213 \
    --time-step 8.000000e+00 \
    --m-distr 2 --min-mass " + minmass + " --max-mass " + maxmass + " \
    --mean-mass " + opts.mean_mass + " --stdev-mass " + opts.stdev_mass + "  \
    --enable-milkyway 1.700000e+00"
  if opts.test:
    print command + "\n"
  else:
    os.system( command )
  
  print "** Generating the population with uniform total mass distribution **"
  command = "lalapps_inspinj \
    --user-tag UNIFORM \
    --source-file inspsrcs.new \
    --gps-start-time 793130413 --gps-end-time 795679213 \
    --time-step 8.000000e+00 \
    --m-distr 0 --min-mass " + minmass + " --max-mass " + maxmass + " \
    --enable-milkyway 1.700000e+00"
  if opts.test:
    print command + "\n"
  else:
    os.system( command )

#######################################################################
# sire and coire full data
#######################################################################
if not opts.skip_coiredata:
  print "** Processing full data set"
  command = "hipecoire --trig-path " + MYRESULTSDIR + "/full_data/ --ifo H1 --ifo H2 --ifo L1 \
    --num-slides 50 --zero-data all_data --cluster-time 10 \
    --veto-file " + MYRESULTSDIR + "/combinedVetoesH1-23.list \
    --veto-file " + MYRESULTSDIR + "/combinedVetoesH2-23.list \
    --veto-file " + MYRESULTSDIR + "/combinedVetoesL1-23.list "
  if opts.second_coinc:
    command += " --second-coinc "
  if opts.rsq_threshold and opts.rsq_max_snr:
    command += " --rsq-threshold " + opts.rsq_threshold + " --rsq-max-snr " +\
        opts.rsq_max_snr
  if opts.test:
    print command + "\n"
  else:
    os.makedirs( MYRESULTSDIR + "/hipecoire/full_data" )
    os.chdir( MYRESULTSDIR + "/hipecoire/full_data" )
    os.system( command )
    # link to the files needed for the upper limit
    os.chdir( MYRESULTSDIR + "/hipecoire" )
    for file in glob.glob("full_data/H*SLIDE*.xml"):
      tmpdest = os.path.splitext( os.path.basename(file) )
      os.symlink( file, tmpdest[0] + "_slides.xml" )
    for file in glob.glob("full_data/H*THINCA_CLUST*.xml"):
      tmpdest = os.path.splitext( os.path.basename(file) )
      os.symlink( file, tmpdest[0] + "_zero.xml" )

os.chdir(MYRESULTSDIR)

#######################################################################
# sire and coire injections
#######################################################################
if not opts.skip_coireinj:
# loop over the injection directories running hipecoire and linking
# the appropriate files. 
  print "** Processing the injections"
  for mydir in glob.glob( "injections*" ):
    print "** Processing " + mydir
    injectionfile=glob.glob( MYRESULTSDIR + "/" + mydir + "/HL-INJECTIONS*.xml")
    print injectionfile
    command = "hipecoire --trig-path " + MYRESULTSDIR + "/" + mydir +\
        " --ifo H1 --ifo H2 --ifo L1 --injection-file " + injectionfile[0] +\
        " --injection-window 10 --cluster-time 10 \
        --veto-file " + MYRESULTSDIR + "/combinedVetoesH1-23.list \
        --veto-file " + MYRESULTSDIR + "/combinedVetoesH2-23.list \
        --veto-file " + MYRESULTSDIR + "/combinedVetoesL1-23.list "
    if opts.second_coinc:
      command += " --second-coinc "
    if opts.rsq_threshold and opts.rsq_max_snr:
      command += " --rsq-threshold " + opts.rsq_threshold + " --rsq-max-snr " +\
        opts.rsq_max_snr
    if opts.usertag:
      command += " --usertag " + opts.usertag

    if opts.test:
      print command + "\n"
    else:
      os.makedirs( MYRESULTSDIR + "/hipecoire/" + mydir )
      os.chdir( MYRESULTSDIR + "/hipecoire/" + mydir ) 
      os.system( command )
      os.chdir( MYRESULTSDIR + "/hipecoire/" )
      for file in glob.glob( mydir + "/H*FOUND.xml"):
        tmpdest = os.path.splitext( os.path.basename(file) )
        os.symlink( file, tmpdest[0] + "_" + mydir + ".xml" )
      for file in glob.glob( mydir + "/H*MISSED.xml"):
        tmpdest = os.path.splitext( os.path.basename(file) )
        os.symlink( file, tmpdest[0] + "_" + mydir + ".xml" )
      os.chdir(MYRESULTSDIR)

# generate the background/foreground number of events plot
command = "plotthinca --glob 'H1*-THINCA_*.xml' --plot-slides --num-slides 50 \
  --figure-name 'summary' --add-zero-lag --snr-dist --min-snr 8.5 \
  --max-snr 11.5 --statistic effective_snr"
if opts.test:
  print command + "\n"
else:
  os.chdir( MYRESULTSDIR + "/hipecoire/full_data" )
  os.system( command )

#######################################################################
# skip plotnumgalaxies
#######################################################################
if not opts.skip_png:
  timetypes = ["H1H2L1", "H1H2", "H1L1", "H2L1"]
  for times in timetypes:
    print "running plotnumgalaxies for " + times 
    command = "plotnumgalaxies \
      --slide-glob '" + MYRESULTSDIR + "/hipecoire/H*slides.xml' \
      --zero-glob '" + MYRESULTSDIR + "/hipecoire/H*zero.xml' \
      --found-glob '" + MYRESULTSDIR + "/hipecoire/" + times + "-THINCA*FOUND*.xml' \
      --missed-glob '" + MYRESULTSDIR + "/hipecoire/" + times + "-THINCA*MISSED*.xml' \
      --source-file '" + MYRESULTSDIR + "/inspsrcs.new' \
      --injection-glob '" + MYRESULTSDIR + "/HL-INJECTIONS_1-793130413-2548800.xml'" 
    if times == "H1H2":
      command += " --num-slides 50 --plot-cum-loudest --plot-pdf-loudest \
        --x-value chirp_dist_h --x-max 40.0 --plot-ng --plot-efficiency \
        --cum-search-ng --mc-errors --figure-name H1H2 --verbose --nbins 20 \
        --distance-error positive --magnitude-error positive \
        --waveform-systematic 0.1 --h-calibration 0.08"
    else:
      command += " --num-slides 50 --plot-cum-loudest --plot-pdf-loudest \
        --x-value chirp_dist_h --x-max 40.0 --y-value chirp_dist_l --axes-square \
        --plot-2d-ng --plot-effcontour --cum-search-2d-ng --mc-errors \
        --figure-name " + times + " --verbose --nbins 20 --distance-error positive \
        --magnitude-error positive --waveform-systematic 0.1 --h-calibration 0.08 \
        --l-calibration 0.05"
    if opts.test:
      print command + "\n"
    else:
      os.makedirs( MYRESULTSDIR + "/plotnumgalaxies/" + times )
      os.chdir( MYRESULTSDIR + "/plotnumgalaxies/" + times )
      os.system( command )


#######################################################################
# skip plotnumgalaxies
#######################################################################
if not opts.skip_upper_limit:
  print "** computing the upper limit"
  command = "lalapps_compute_posterior \
      -f H1H2/png-output.ini -t 0.0143802 \
      -f H1H2L1/png-output.ini -t 0.0416896 \
      -f H1L1/png-output.ini -t 0.00524945 \
      -f H2L1/png-output.ini -t 0.0044068 \
      -p uniform -m 500 -d 0.01 -F s4-gauss -x 5 \
      --calibration-error  --magnitude-error  --montecarlo-error \
      --waveform-error  --distance-error -n 16000"
  if opts.test:
    print command + "\n"
  else:
    os.chdir( MYRESULTSDIR + "/plotnumgalaxies/" )
    os.system( command )


#######################################################################
# skip plotnumgalaxies;  this is the second run for the uniform mass
# distribution in order give an upper limit as a function of mass
#######################################################################
if not opts.skip_png:
  timetypes = ["H1H2L1", "H1H2", "H1L1", "H2L1"]
  for times in timetypes:
    print "running plotnumgalaxies for " + times 
    os.makedirs( MYRESULTSDIR + "/plotnumgalaxies/" + times + "_uniform")
    os.chdir( MYRESULTSDIR + "/plotnumgalaxies/" + times + "_uniform")
    command = "plotnumgalaxies \
      --slide-glob '" + MYRESULTSDIR + "/hipecoire/H*slides.xml' \
      --zero-glob '" + MYRESULTSDIR + "/hipecoire/H*zero.xml' \
      --found-glob '" + MYRESULTSDIR + "/hipecoire/" + times + "-THINCA*FOUND*.xml' \
      --missed-glob '" + MYRESULTSDIR + "/hipecoire/" + times + "-THINCA*MISSED*.xml' \
      --source-file '" + MYRESULTSDIR + "/inspsrcs.new' \
      --injection-glob '" + MYRESULTSDIR + "/HL-INJECTIONS_1_UNIFORM-793130413-2548800.xml'" 
    if times == "H1H2":
      command += " --num-slides 50 --plot-cum-loudest --plot-pdf-loudest \
        --x-value chirp_dist_h --x-max 40.0 --plot-ng --plot-efficiency \
        --cum-search-ng --mc-errors --figure-name H1H2 --verbose --nbins 20 \
        --distance-error positive --magnitude-error positive \
        --waveform-systematic 0.1 --h-calibration 0.08 \
        --m-low " + minmtotal + " --m-high " + maxmtotal + " --m-dm " +\
        opts.m_dm
    else:
      command += " --num-slides 50 --plot-cum-loudest --plot-pdf-loudest \
        --x-value chirp_dist_h --x-max 40.0 --y-value chirp_dist_l --axes-square \
        --plot-2d-ng --plot-effcontour --cum-search-2d-ng --mc-errors \
        --figure-name " + times + " --verbose --nbins 20 --distance-error positive \
        --magnitude-error positive --waveform-systematic 0.1 --h-calibration 0.08 \
        --l-calibration 0.05 \
        --m-low " + minmtotal + " --m-high " + maxmtotal + " --m-dm " +\
        opts.m_dm
    if opts.test:
      print command + "\n"
    else:
      os.system( command )


#######################################################################
# skip plotnumgalaxies
#######################################################################
if not opts.skip_upper_limit:
  print "** computing the upper limit"
  command = "lalapps_compute_posterior \
      -f H1H2_uniform/png-output.ini -t 0.0143802 \
      -f H1H2L1_uniform/png-output.ini -t 0.0416896 \
      -f H1L1_uniform/png-output.ini -t 0.00524945 \
      -f H2L1_uniform/png-output.ini -t 0.0044068 \
      -p uniform -m 500 -d 0.01 -F 's4-uniform' -x 5 \
      --calibration-error  --magnitude-error  --montecarlo-error \
      --waveform-error  --distance-error -n 16000"
  if opts.test:
    print command + "\n"
  else:
    os.chdir( MYRESULTSDIR + "/plotnumgalaxies/" )
    os.system( command )

