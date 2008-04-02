#!/usr/bin/python

# $Id$
__author__ = "Patrick Brady <patrick@gravity.phys.uwm.edu>"
__version__ = "$Revision$"[11:-2]
__date__ = "$Date$"[7:-2]

import sys
import os
import copy
import shutil
from optparse import *
import exceptions
import glob
import ConfigParser

# just testing

def mkdirsafe( directory ):
  """
     Creates a directory, does not nag about when it already exists
  """
  try:
     os.makedirs(directory)
  except OSError, (errno, strerror):
    if errno==17:
      print "WARNING: directory %s already exist" % (directory) 
    else:
      raise
  

def symlinksafe( target, linkname ):
  """
     Creates a link, does not nag about when it already exists
  """
  try:
     os.symlink( target, linkname )
  except OSError, (errno, strerror):
    if errno==17:
      print "WARNING: link %s already exist" % (linkname)
    else:
      raise

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

parser = OptionParser( usage=usage, version="%prog CVS $Id$" )
# used letters:
# a c d e f g i j k l m n p r s t u v x y z
# A B C H L R S T
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

parser.add_option("","--skip-injcut",action="store_true",default=False,\
    help="skip the injection file pruning by mass")

parser.add_option("","--skip-coiremasscut",action="store_true",default=False,\
    help="skip the mass cut")

parser.add_option("-p","--skip-population",action="store_true",default=False,\
    help="skip population generation")

parser.add_option("-c","--skip-coiredata",action="store_true",default=False,\
    help="skip coiring of data")

parser.add_option("-i","--skip-coireinj",action="store_true",default=False,\
    help="skip coiring of injections")

parser.add_option("-L","--skip-plotthinca",action="store_true",default=False,\
    help="skip plotthinca step")

parser.add_option("-n","--skip-png",action="store_true",default=False,\
    help="skip plotnumgalaxies")

parser.add_option("-u","--skip-upper-limit",action="store_true",default=False,\
    help="skip upper limit")

parser.add_option("-t","--usertag",action="store",type="string",\
    default=None, metavar=" USERTAG", help="the user tag used for injections")

parser.add_option("-T","--test",action="store_true",\
    default=False, help="only print the commands to be run")

parser.add_option("-v","--no-veto",action="store_true",\
    default=False, \
    help="do not apply any vetoes (in the case they already have been applied)")

parser.add_option("-r","--remove-h2l1",action="store_true",\
    default=False, help="remove H2-L1 triggers from three ifo times")

parser.add_option("-z","--output-to-file",action="store_true",\
    default=False, 
    help="pipe output from plotnumgalaxies and compute_posterior to files")

# hierarchy pipeline
parser.add_option("-j","--second-coinc",action="store_true",default=False,\
    help="specify that the input comes from second stage thinca files")

parser.add_option("","--clustered-files",action="store_true",default=False,\
    help="specify that the input comes from clustered files")

(opts,args) = parser.parse_args()

#######################################################################
# check options and initialize
#######################################################################

if not opts.config_file:
  print >> sys.stderr, "Must specify --config-file"
  sys.exit(1)

if not opts.results_dir:
  print >> sys.stderr, "Must specify --results-dir"
  sys.exit(1)

#######################################################################
# create the config parser object and read in the ini file
cp = ConfigParser.ConfigParser()
cp.read(opts.config_file)
fulldata = cp.items("fulldata")
injdata = cp.items("injdata")
vetodata = cp.items("vetodata")
sourcedata = cp.items("sourcedata")
inspinj_options = cp.items("inspinjdata")
distributions = cp.items("distdata")
gaussianmass_options = cp.items("gaussianmassdata")
totalmass_options = cp.items("totalmassdata")
componentmass_options = cp.items("componentmassdata")
masses = cp.items("massdata")
stat_options = cp.items("statistic")
analyzedtimes = cp.items("analyzedtimes")
coire_options = cp.items("coire")
coireinj_options = cp.items("coireinj")
coiredata_options = cp.items("coiredata")
injcut_options = cp.items("injcut")
coiremasscut_options = cp.items("coiremasscut")
plotthinca_options=cp.items("plotthinca")
png_options = cp.items("plotnumgalaxies")
png_y_options = cp.items("plotnumgalaxies-y")
posterior_options = cp.items("upperlimit")

# convert the statistic options into a dictionary
stat_dict=dict()
for opt in stat_options:
    stat_dict[opt[0]]=opt[1]

# convert the distributon options into a dictionary
dist_dict=dict()
for opt in distributions:
  dist_dict[opt[0]]=int(opt[1])

# convert the gaussian mass options into a dictionary
gmass_dict=dict()
for opt in gaussianmass_options:
    gmass_dict[opt[0]]=opt[1]

# convert the total mass options into a dictionary
tmass_dict=dict()
for opt in totalmass_options:
    tmass_dict[opt[0]]=opt[1]

# convert the component mass options into a dictionary
cmass_dict=dict()
for opt in componentmass_options:
    cmass_dict[opt[0]]=opt[1]

# convert the mass options into a dictionary
mass_dict=dict()
for opt in masses:
    mass_dict[opt[0]]=opt[1]

# convert the coiredata options into a dictionary
coiredata_dict=dict()
for opt in coiredata_options:
    coiredata_dict[opt[0]]=opt[1]

# convert the injcut options into a dictionary
injcut_dict=dict()
for opt in injcut_options:
    injcut_dict[opt[0]]=opt[1]

#######################################################################
# do the set up
#######################################################################
MYRESULTSDIR=opts.results_dir

if not opts.skip_setup:
  mkdirsafe(MYRESULTSDIR)
  os.chdir(MYRESULTSDIR)

  symlinksafe(fulldata[0][1], "./full_data")
  for injdir in injdata:
    symlinksafe(injdir[1], "./" + injdir[0])
  vetodict=dict()  
  for myfile in vetodata:
    vetodict[myfile[0]]=myfile[1]
  shutil.copy( vetodict["h1vetofile"], "./h1veto.list" )
  shutil.copy( vetodict["h2vetofile"], "./h2veto.list" )
  shutil.copy( vetodict["l1vetofile"], "./l1veto.list" )  
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
  os.chdir(MYRESULTSDIR)
  if dist_dict["gaussian"]:
    print "** Generating the population with Gaussian mass distribution **"
    command = "lalapps_inspinj" + \
      " --user-tag GAUSSIANMASS" + \
      " --source-file inspsrcs.new" + \
      " --gps-start-time 793130413 --gps-end-time 795679213" + \
      " --time-step 8.000000e+00" + \
      " --write-compress"
    for opt in gaussianmass_options:
      command += " --" + opt[0] + " " + opt[1]
    for opt in inspinj_options:
      command += " --" + opt[0] + " " + opt[1]
    if opts.test:
      print command + "\n"
    else:
      os.system( command )

  if dist_dict["totalmass"]:
    print "** Generating the population with uniform total mass distribution **"
    command = "lalapps_inspinj" + \
      " --user-tag TOTALMASS" + \
      " --source-file inspsrcs.new" + \
      " --gps-start-time 793130413 --gps-end-time 795679213" + \
      " --time-step 8.000000e+00" + \
      " --write-compress"
    for opt in totalmass_options:
      command += " --" + opt[0] + " " + opt[1]
    for opt in inspinj_options:
      command += " --" + opt[0] + " " + opt[1]
    if opts.test:
      print command + "\n"
    else:
      os.system( command )

  if dist_dict["componentmass"]:
    print "** Generating the population with uniform component mass distribution **"
    command = "lalapps_inspinj" + \
      " --user-tag COMPONENTMASS" + \
      " --source-file inspsrcs.new" + \
      " --gps-start-time 793130413 --gps-end-time 795679213" + \
      " --time-step 8.000000e+00" + \
      " --write-compress"
    for opt in componentmass_options:
      command += " --" + opt[0] + " " + opt[1]
    for opt in inspinj_options:
      command += " --" + opt[0] + " " + opt[1]
    if opts.test:
      print command + "\n"
    else:
      os.system( command )

#######################################################################
# create mass cut injection files
#######################################################################
if not (opts.skip_injcut):
  print "** Creating mass cut injection files"
  os.chdir(MYRESULTSDIR)
  mkdirsafe( MYRESULTSDIR + "/injcut" )
  for mydir in glob.glob( "injections*" ):
    print "** Processing " + mydir
    mkdirsafe( MYRESULTSDIR + "/injcut/" + mydir )
    injectionfile=glob.glob(MYRESULTSDIR + "/" + mydir + \
        "/HL-INJECTIONS*.xml.gz")
    injFileName = os.path.basename( injectionfile[0] )
    command = "lalapps_injcut --injection-file " + injectionfile[0]
    for key in injcut_dict.keys():
      if key in ['mass-cut', 'mass-range-low', 'mass-range-high',
          'mass2-range-low', 'mass2-range-high']:
        command += " --" + key + " " + injcut_dict[key]
    command +=  " --output " + MYRESULTSDIR + "/injcut/" + mydir + \
        "/" + injFileName
    if opts.test:
      print command + "\n"
    else:
      os.system( command )
    if not injectionfile:
      print "ERROR in injcut: No injection-file (HL*.xml.gz) \
          in the injections-directory. Exiting..."
      sys.exit(1)

#######################################################################
# sire and coire full data
#######################################################################
if not opts.skip_coiredata:
  print "** Processing full data set"
  mkdirsafe( MYRESULTSDIR + "/hipecoire/full_data" )
  command = "hipecoire --trig-path " + MYRESULTSDIR + \
      "/full_data/ --ifo H1 --ifo H2 --ifo L1" + \
      " --cluster-infinity --coinc-stat " + stat_dict["statistic"]
  if stat_dict["statistic"][-3:] == "snr":
    # coire takes "snrsq"/"effective_snrsq"
    command += "sq"
  if "bitten_l" in stat_dict["statistic"]:
    command += " --h1-bittenl-a "+stat_dict["bittenl_a"] + \
        " --h1-bittenl-b "+stat_dict["bittenl_b"] + \
        " --h2-bittenl-a "+stat_dict["bittenl_a"] + \
        " --h2-bittenl-b "+stat_dict["bittenl_b"] + \
        " --l1-bittenl-a "+stat_dict["bittenl_a"] + \
        " --l1-bittenl-b "+stat_dict["bittenl_b"] 
  if not opts.no_veto:
    command+=" --veto-file " + MYRESULTSDIR + "/h1veto.list" + \
        " --veto-file " + MYRESULTSDIR + "/h2veto.list" + \
        " --veto-file " + MYRESULTSDIR + "/l1veto.list"
  if opts.second_coinc:
    command += " --second-coinc"
  if opts.clustered_files:
    command += " --clustered-files"
  for opt in coiredata_options:
    command += " --" + opt[0] + " " + opt[1]
  for opt in coire_options:
    command += " --" + opt[0] + " " + opt[1]
  if not opts.skip_coiremasscut:
    for opt in coiremasscut_options:
      command += " --" + opt[0] + " " + opt[1]
  if opts.test:
    print command + "\n"
  else:
    os.chdir( MYRESULTSDIR + "/hipecoire/full_data" )
    os.system( command )
    
    # link to the files needed for the upper limit
    os.chdir( MYRESULTSDIR + "/hipecoire" )
    if not opts.remove_h2l1:
      for file in glob.glob("full_data/H*SLIDE*.xml.gz"):
        tmpdest = os.path.splitext( os.path.splitext( \
            os.path.basename(file) )[0] )
        symlinksafe( file, tmpdest[0] + "_slides.xml.gz" )
      for file in ( glob.glob("full_data/H*COIRE_CLUST*.xml.gz") + \
          glob.glob("full_data/H*COIRE_LOUDEST*.xml.gz") ):
        tmpdest = os.path.splitext( os.path.splitext( \
            os.path.basename(file) )[0] )
        symlinksafe( file, tmpdest[0] + "_zero.xml.gz" )
    else:
      for file in ( glob.glob("full_data/H???-COIRE_SLIDE*.xml.gz") + \
          glob.glob("full_data/H1*/H1*SLIDE_in_H1H2L1*.xml.gz") ):
        tmpdest = os.path.splitext( os.path.splitext( \
            os.path.basename(file) )[0] )
        symlinksafe( file, tmpdest[0] + "_slides.xml.gz" )
      for file in ( glob.glob("full_data/H???-COIRE_CLUST*.xml.gz") + \
          glob.glob("full_data/H???-COIRE_LOUDEST*.xml.gz") + \
          glob.glob("full_data/H1*/H1*-*COIRE_in_H1H2L1*.xml.gz") ):
        tmpdest = os.path.splitext( os.path.splitext( \
            os.path.basename(file) )[0] )
        symlinksafe( file, tmpdest[0] + "_zero.xml.gz" )

#######################################################################
# sire and coire injections
#######################################################################
if not opts.skip_coireinj:
  # loop over the injection directories running hipecoire and linking
  # the appropriate files. 
  print "** Processing the injections"
  os.chdir(MYRESULTSDIR)
  for mydir in glob.glob( "injections*" ):
    print "** Processing " + mydir
    mkdirsafe( MYRESULTSDIR + "/hipecoire/" + mydir )
    if not opts.skip_injcut:
      injectionfile=glob.glob(MYRESULTSDIR + "/injcut/" + mydir + \
          "/HL-INJECTIONS*.xml.gz")
    else:
      injectionfile=glob.glob(MYRESULTSDIR + "/" + mydir + \
          "/HL-INJECTIONS*.xml.gz")
    injFileName = os.path.basename( injectionfile[0] )
    symlinksafe( injectionfile[0], MYRESULTSDIR + "/hipecoire/" + mydir + \
        "/" + injFileName )
    if not injectionfile:
      print "ERROR in coireinj: No injection-file (HL*.xml.gz) \
          in the injections-directory. Exiting..."
      sys.exit(1)
    command = "hipecoire --trig-path " + MYRESULTSDIR + "/" + mydir + \
        " --ifo H1 --ifo H2 --ifo L1 --injection-file " + injFileName + \
        " --coinc-stat " + stat_dict["statistic"]
    if stat_dict["statistic"][-3:] == "snr":
      # coire takes "snrsq"/"effective_snrsq"
      command += "sq"
    if "bitten_l" in stat_dict["statistic"]:
      command += " --h1-bittenl-a " + stat_dict["bittenl_a"] + \
          " --h1-bittenl-b " + stat_dict["bittenl_b"] + \
          " --h2-bittenl-a " + stat_dict["bittenl_a"] + \
          " --h2-bittenl-b " + stat_dict["bittenl_b"] + \
          " --l1-bittenl-a " + stat_dict["bittenl_a"] + \
          " --l1-bittenl-b " + stat_dict["bittenl_b"] 
    for opt in coireinj_options:
      command += " --" + opt[0] + " " + opt[1]
    if not opts.no_veto:
      command += " --veto-file " + MYRESULTSDIR + "/h1veto.list" + \
          " --veto-file " + MYRESULTSDIR + "/h2veto.list" + \
          " --veto-file " + MYRESULTSDIR + "/l1veto.list"
    if opts.second_coinc:
      command += " --second-coinc"
    if opts.clustered_files:
      command += " --clustered-files"
    for opt in coire_options:
      command += " --" + opt[0] + " " + opt[1]
    if not opts.skip_coiremasscut:
      for opt in coiremasscut_options:
        command += " --" + opt[0] + " " + opt[1]
    if opts.test:
      print command + "\n"
    else:
      os.chdir( MYRESULTSDIR + "/hipecoire/" + mydir )
      os.system( command )
      os.chdir( MYRESULTSDIR + "/hipecoire/" )
      for file in glob.glob( mydir + "/H*FOUND.xml.gz"):
        tmpdest = os.path.splitext( os.path.splitext( \
            os.path.basename(file) )[0] )
        symlinksafe( file, tmpdest[0] + "_" + mydir + ".xml.gz" )
      for file in glob.glob( mydir + "/H*MISSED.xml.gz"):
        tmpdest = os.path.splitext( os.path.splitext( \
            os.path.basename(file) )[0] )
        symlinksafe( file, tmpdest[0] + "_" + mydir + ".xml.gz" )
      os.chdir(MYRESULTSDIR)

#######################################################################
# plotthinca step
#######################################################################
if not opts.skip_plotthinca:
  print "running plotthinca to generate foreground/background plots"
  # generate the background/foreground number of events plot
  command = "plotthinca --glob '" +\
      MYRESULTSDIR + "/hipecoire/*CLUST*zero* " + \
      MYRESULTSDIR + "/hipecoire/*CLUST*slides*'" + \
      " --plot-slides --num-slides " + coiredata_dict["num-slides"] + \
      " --figure-name 'summary' --add-zero-lag --snr-dist" + \
      " --statistic " + stat_dict["statistic"]
  for opt in plotthinca_options:
    command+=" --"+opt[0]+" "+opt[1]
  if "bitten_l" in stat_dict["statistic"]:
    command += " --bittenl_a " + stat_dict["bittenl_a"] + \
     " --bittenl_b "+stat_dict["bittenl_b"]
  if opts.test:
    print command + "\n"
  else:
    mkdirsafe( MYRESULTSDIR + "/plotthinca" )
    os.chdir( MYRESULTSDIR + "/plotthinca" )
    os.system( command )

#######################################################################
# plotnumgalaxies step
#######################################################################
popdistr={}
if dist_dict["gaussian"]:
  popdistr["gaussian"] = {}
  popdistr["gaussian"]["file"] = \
      "/HL-INJECTIONS_1_GAUSSIANMASS-793130413-2548800.xml.gz"
  popdistr["gaussian"]["popType"] = "gaussian"
  popdistr["gaussian"]["useMassInfo"] = 0
if dist_dict["totalmass"]:
  popdistr["totalmass"] = {}
  popdistr["totalmass"]["file"] = \
      "/HL-INJECTIONS_1_TOTALMASS-793130413-2548800.xml.gz"
  popdistr["totalmass"]["popType"] = "totalmass"
  popdistr["totalmass"]["useMassInfo"] = 1
  popdistr["totalmass"]["min-mass"] = tmass_dict["min-mtotal"]
  popdistr["totalmass"]["max-mass"] = tmass_dict["max-mtotal"]
if dist_dict["componentmass"]:
  popdistr["componentmass"] = {}
  popdistr["componentmass"]["file"] = \
      "/HL-INJECTIONS_1_COMPONENTMASS-793130413-2548800.xml.gz"
  popdistr["componentmass"]["popType"] = "componentmass"
  popdistr["componentmass"]["useMassInfo"] = 1
  popdistr["componentmass"]["min-mass"] = cmass_dict["min-mass1"]
  popdistr["componentmass"]["max-mass"] = cmass_dict["max-mass1"]

if not opts.skip_png:
  for key in popdistr:
    pop=popdistr[key]["file"]
    for times in analyzedtimes:
      print "running plotnumgalaxies for " + times[0].upper() + \
          " using population from " + pop
      command = "plotnumgalaxies" + \
          " --slide-glob '" + MYRESULTSDIR + \
          "/hipecoire/H*LOUDEST*slides.xml.gz'" + \
          " --zero-glob '" + MYRESULTSDIR + \
          "/hipecoire/H*LOUDEST*zero.xml.gz'" + \
          " --found-glob '" + MYRESULTSDIR + "/hipecoire/" + \
          times[0].upper() + "-COIRE*FOUND*.xml.gz'" + \
          " --missed-glob '" + MYRESULTSDIR + "/hipecoire/" + \
          times[0].upper() + "-COIRE*MISSED*.xml.gz'" + \
          " --source-file '" + MYRESULTSDIR + "/inspsrcs.new'" + \
          " --population-glob '" + MYRESULTSDIR + pop + \
          "' --figure-name " + times[0].upper() + \
          " --plot-cum-loudest --plot-pdf-loudest" + \
          " --num-slides " + coiredata_dict["num-slides"] + \
          " --statistic " + stat_dict["statistic"] + \
          " --plot-efficiency"
      if  "bitten_l" in stat_dict["statistic"]:
        command += " --bittenl_a " + stat_dict["bittenl_a"] + \
            " --bittenl_b "+stat_dict["bittenl_b"]
      for opt in png_options:
        command += " --" + opt[0] + " " + opt[1]
      if times[0].upper() == "H1H2":
        command +=  " --plot-ng --plot-efficiency --cum-search-ng" + \
            " --plot-ng-vs-stat"
      else:
        command += " --axes-square --plot-2d-ng" + \
            " --plot-effcontour --cum-search-2d-ng" + \
            " --plot-2d-ng-vs-stat"
        for opt in png_y_options:
          command += " --" + opt[0] + " " + opt[1]
      command += " --population-type " + popdistr[key]["popType"]
      if not popdistr[key]["popType"] == "gaussian":
        command += " --m-low " + popdistr[key]["min-mass"] + \
            " --m-high " + popdistr[key]["max-mass"] + \
            " --m-dm " + mass_dict["m-dm"] + " --cut-inj-by-mass " + \
            mass_dict["cut-inj-by-mass"]
      if opts.output_to_file:
        command += " > png-output-" + times[0].upper() + ".log"
      if opts.test:
        print command + "\n"
      else:
        dir = MYRESULTSDIR + "/plotnumgalaxies/" + times[0].upper()
        dir = dir + "_" + popdistr[key]["popType"]
        mkdirsafe( dir )
        os.chdir( dir )
        os.system( command )

#######################################################################
# calculate upper limit
#######################################################################
if not opts.skip_upper_limit:
  secondsInAYear = 31556736.0
  for key in popdistr:
    pop = popdistr[key]["popType"]
    print
    print "** Computing the upper limit for " + pop + " population"
    command = "lalapps_compute_posterior"
    for times in analyzedtimes:
      command += " -f " + times[0].upper()
      command += "_" + pop
      command += "/png-output.ini -t " + str(float(times[1])/secondsInAYear)
    command += " --figure-name " + pop
    for opt in posterior_options:
       command += " --" + opt[0] + " " + opt[1]
    command += " --mass-region-type " + pop
    if opts.output_to_file:
      command += " > ul-output-" + pop.lower() + ".log"
    if opts.test:
      print command + "\n"
    else:
      os.chdir( MYRESULTSDIR + "/plotnumgalaxies/" )
      os.system( command )
