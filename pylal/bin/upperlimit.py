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
    help="use the second stage thinca files as input")

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
masses = cp.items("massdata")
stat_options = cp.items("statistic")
analyzedtimes = cp.items("analyzedtimes")
coire_options = cp.items("coire")
plotthinca_options=cp.items("plotthinca")
png_options = cp.items("plotnumgalaxies")
png_y_options = cp.items("plotnumgalaxies-y")
posterior_options=cp.items("upperlimit")

# convert the statistic options into a dictionary
stat_dict=dict()
for opt in stat_options:
    stat_dict[opt[0]]=opt[1]

#######################################################################
# do the set up
#######################################################################
MYRESULTSDIR=opts.results_dir

# dealing with the masses
mass_dict=dict()
for opt in masses:
    mass_dict[opt[0]]=opt[1]

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
  print "** Generating the population with Gaussian mass distribution **"
  command = "lalapps_inspinj " + \
    "--source-file inspsrcs.new " + \
    "--gps-start-time 793130413 --gps-end-time 795679213 " + \
    "--time-step 8.000000e+00 " + \
    "--m-distr 2 --min-mass " + mass_dict["min-mass"] + \
    " --max-mass " +  mass_dict["max-mass"] + \
    " --mean-mass " +  mass_dict["mean-mass"] + \
    " --stdev-mass " +  mass_dict["stdev-mass"] + \
    " --enable-milkyway 1.700000e+00"
  if opts.test:
    print command + "\n"
  else:
    os.system( command )
  
  print "** Generating the population with uniform total mass distribution **"
  command = "lalapps_inspinj " + \
    " --user-tag UNIFORM " +\
    " --source-file inspsrcs.new" + \
    " --gps-start-time 793130413 --gps-end-time 795679213" + \
    " --time-step 8.000000e+00" + \
    " --m-distr 0 --min-mass " +  mass_dict["min-mass"] + \
    " --max-mass " + mass_dict["max-mass"] + \
    " --enable-milkyway 1.700000e+00"
  if opts.test:
    print command + "\n"
  else:
    os.system( command )

#######################################################################
# sire and coire full data
#######################################################################
if not opts.skip_coiredata:
  print "** Processing full data set"
  command = "hipecoire --trig-path " + MYRESULTSDIR + \
      "/full_data/ --ifo H1 --ifo H2 --ifo L1 " + \
      "--num-slides 50 --zero-data exclude_play  " + \
      "--cluster-infinity --coinc-stat " + stat_dict["statistic"]
  if stat_dict["statistic"][-3:] == "snr":
    # coire takes "snrsq"/"effective_snrsq"
    command += "sq"
  if stat_dict["statistic"] == "bitten_l":
    command += " --h1-bittenl-a "+stat_dict["bittenl_a"] + \
        " --h1-bittenl-b "+stat_dict["bittenl_b"] + \
        " --h2-bittenl-a "+stat_dict["bittenl_a"] + \
        " --h2-bittenl-b "+stat_dict["bittenl_b"] + \
        " --l1-bittenl-a "+stat_dict["bittenl_a"] + \
        " --l1-bittenl-b "+stat_dict["bittenl_b"] 
  if not opts.no_veto:
    command+=" --veto-file " + MYRESULTSDIR + "/h1veto.list " + \
        " --veto-file " + MYRESULTSDIR + "/h2veto.list " + \
        " --veto-file " + MYRESULTSDIR + "/l1veto.list"
  if opts.second_coinc:
    command += " --second-coinc "
  for opt in coire_options:
       command+=" --"+opt[0]+" "+opt[1]
  if opts.test:
    print command + "\n"
  else:
    mkdirsafe( MYRESULTSDIR + "/hipecoire/full_data" )
    os.chdir( MYRESULTSDIR + "/hipecoire/full_data" )
    os.system( command )
    
    # link to the files needed for the upper limit
    os.chdir( MYRESULTSDIR + "/hipecoire" )
    if not opts.remove_h2l1:
      for file in glob.glob("full_data/H*SLIDE*.xml"):
        tmpdest = os.path.splitext( os.path.basename(file) )
        symlinksafe( file, tmpdest[0] + "_slides.xml" )
      for file in ( glob.glob("full_data/H*THINCA_CLUST*.xml") + \
          glob.glob("full_data/H*THINCA_LOUDEST*.xml") ):
        tmpdest = os.path.splitext( os.path.basename(file) )
        symlinksafe( file, tmpdest[0] + "_zero.xml" )
    else:
      for file in ( glob.glob("full_data/H???-THINCA_SLIDE*.xml") + \
          glob.glob("full_data/H1*/H1*SLIDE_in_H1H2L1*.xml") ):
        tmpdest = os.path.splitext( os.path.basename(file) )
        symlinksafe( file, tmpdest[0] + "_slides.xml" )
      for file in ( glob.glob("full_data/H???-THINCA_CLUST*.xml") + \
          glob.glob("full_data/H???-THINCA_LOUDEST*.xml") + \
          glob.glob("full_data/H1*/H1*-*THINCA_in_H1H2L1*.xml") ):
        tmpdest = os.path.splitext( os.path.basename(file) )
        symlinksafe( file, tmpdest[0] + "_zero.xml" )


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
    injectionfile=glob.glob( MYRESULTSDIR + "/" + mydir + "/HL-INJECTIONS*.xml")
    if not injectionfile:
      print "ERROR in coireinj: No injection-file (HL*.xml) \
          in the injections-directory. Exiting..."
      sys.exit(1)
    command = "hipecoire --trig-path " + MYRESULTSDIR + "/" + mydir +\
        " --ifo H1 --ifo H2 --ifo L1 --injection-file " + injectionfile[0] +\
        " --injection-window 10 --coinc-stat " + stat_dict["statistic"]
    if stat_dict["statistic"][-3:] == "snr":
      # coire takes "snrsq"/"effective_snrsq"
      command += "sq"
    if stat_dict["statistic"] == "bitten_l":
      command += " --h1-bittenl-a " + stat_dict["bittenl_a"] + \
          " --h1-bittenl-b " + stat_dict["bittenl_b"] + \
          " --h2-bittenl-a " + stat_dict["bittenl_a"] + \
          " --h2-bittenl-b " + stat_dict["bittenl_b"] + \
          " --l1-bittenl-a " + stat_dict["bittenl_a"] + \
          " --l1-bittenl-b " + stat_dict["bittenl_b"] 
    if not opts.no_veto:
      command+=" --veto-file " + MYRESULTSDIR + "/h1veto.list " + \
          " --veto-file " + MYRESULTSDIR + "/h2veto.list " + \
          " --veto-file " + MYRESULTSDIR + "/l1veto.list"
    if opts.second_coinc:
      command += " --second-coinc "
    for opt in coire_options:
       command+=" --"+opt[0]+" "+opt[1]
    if opts.test:
      print command + "\n"
    else:
      mkdirsafe( MYRESULTSDIR + "/hipecoire/" + mydir )
      os.chdir( MYRESULTSDIR + "/hipecoire/" + mydir ) 
      os.system( command )
      os.chdir( MYRESULTSDIR + "/hipecoire/" )
      for file in glob.glob( mydir + "/H*FOUND.xml"):
        tmpdest = os.path.splitext( os.path.basename(file) )
        symlinksafe( file, tmpdest[0] + "_" + mydir + ".xml" )
      for file in glob.glob( mydir + "/H*MISSED.xml"):
        tmpdest = os.path.splitext( os.path.basename(file) )
        symlinksafe( file, tmpdest[0] + "_" + mydir + ".xml" )
      os.chdir(MYRESULTSDIR)




#######################################################################
# plotthinca step
#######################################################################
if not opts.skip_plotthinca:
  print "running plotthinca to generate foreground/background plots"
  # generate the background/foreground number of events plot
  command = "plotthinca --glob '" + MYRESULTSDIR + "/hipecoire/*CLUST*zero* " \
      + MYRESULTSDIR + "/hipecoire/*CLUST*slides*' " + \
      "--plot-slides --num-slides 50 " \
      "--figure-name 'summary' --add-zero-lag --snr-dist"+\
      " --statistic " + stat_dict["statistic"]
  for opt in plotthinca_options:
    command+=" --"+opt[0]+" "+opt[1]
  if stat_dict["statistic"] == "bitten_l":
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
if not opts.skip_png:
  popfiles = ["/HL-INJECTIONS_1-793130413-2548800.xml",
      "/HL-INJECTIONS_1_UNIFORM-793130413-2548800.xml"]
  for nr in range(0,2):
    pop=popfiles[nr]
    for times in analyzedtimes:
      print "running plotnumgalaxies for " + times[0].upper() + \
          " using population from " + pop
      command = "plotnumgalaxies " + \
          "--slide-glob '" + MYRESULTSDIR + "/hipecoire/H*LOUDEST*slides.xml' "\
          + "--zero-glob '" + MYRESULTSDIR + "/hipecoire/H*LOUDEST*zero.xml' "\
          + "--found-glob '" + MYRESULTSDIR + "/hipecoire/" \
          + times[0].upper() + "-THINCA*FOUND*.xml' " \
          + "--missed-glob '" + MYRESULTSDIR + "/hipecoire/" \
          + times[0].upper() + "-THINCA*MISSED*.xml' " \
          + "--source-file '" + MYRESULTSDIR + "/inspsrcs.new' " \
          + "--injection-glob '" + MYRESULTSDIR + pop + \
          "' --figure-name " + times[0].upper() + \
          " --plot-cum-loudest --plot-pdf-loudest" + \
          " --num-slides 50 --statistic " + stat_dict["statistic"]
      if  stat_dict["statistic"] == "bitten_l":
        command += " --bittenl_a " + stat_dict["bittenl_a"] + \
            " --bittenl_b "+stat_dict["bittenl_b"]
      for opt in png_options:
        command+=" --"+opt[0]+" "+opt[1]
      if times[0].upper() == "H1H2":
        command +=  " --plot-ng --plot-efficiency --cum-search-ng " \
            + "--plot-ng-vs-stat"
      else:
        command += " --axes-square --plot-2d-ng "\
            + "--plot-effcontour --cum-search-2d-ng "\
            + "--plot-2d-ng-vs-stat"
        for opt in png_y_options:
          command+=" --" + opt[0] + " " + opt[1]
      if (nr==1):
        command+=" --m-low " + str(2.0*float(mass_dict["min-mass"])) + \
            " --m-high " + str(2.0*float(mass_dict["max-mass"])) + \
            " --m-dm " + mass_dict["m-dm"] + " --cut-inj-by-mass " + \
            mass_dict["cut-inj-by-mass"]
      if opts.output_to_file:
        command+=" > png-output-" + times[0].upper() + ".log"
      if opts.test:
        print command + "\n"
      else:
        dir=MYRESULTSDIR + "/plotnumgalaxies/" + times[0].upper()
        if  (nr==1):
          dir=dir+"_uniform"
        mkdirsafe( dir )
        os.chdir( dir )
        os.system( command )


#######################################################################
# calculate upper limit
#######################################################################
if not opts.skip_upper_limit:
   secondsInAYear=31556736.0
   popfiles = ["Gaussian","Uniform"]
   for pop in popfiles:
     print
     print "** Computing the upper limit for " + pop + " population"
     command = "lalapps_compute_posterior"
     for times in analyzedtimes:
       command+=" -f "+times[0].upper()
       if (pop=="Uniform"):
         command+="_uniform"
       command+="/png-output.ini -t "+ str(float(times[1])/secondsInAYear)
     command+=" --figure-name "+pop
     for opt in posterior_options:
        command+=" --"+opt[0]+" "+opt[1]
     if opts.output_to_file:
       command += " > ul-output-" + pop.lower() + ".log"
     if opts.test:
       print command + "\n"
     else:
       os.chdir( MYRESULTSDIR + "/plotnumgalaxies/" )
       os.system( command )
