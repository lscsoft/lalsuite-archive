#!/usr/bin/env @PYTHONPROG@
"""
this code contains some functions used by ihope (and followup_pipe)

$Id$
"""
__author__ = 'Stephen Fairhurst <sfairhur@gravity.phys.uwm.edu>'
__date__ = '$Date$'
__version__ = '$Revision$'[11:-2]

##############################################################################
# import standard modules and append the lalapps prefix to the python path
import os, sys, copy
import ConfigParser
import optparse
import tempfile
import urllib
sys.path.append('@PYTHONLIBDIR@')

##############################################################################
# import the modules we need to build the pipeline
from glue import segments
from glue import segmentsUtils
from glue import pipeline

dq_url_pattern = "http://ldas-cit.ligo.caltech.edu/segments/S5/%s/dq_segments.txt"

##############################################################################
# Functions used in setting up the dag:
def make_external_call(command, show_stdout=False, show_command=False):
  """
  Run a program on the shell and print informative messages on failure.
  """
  if show_command: print command

  stdin, out, err = os.popen3(command)
  pid, status = os.wait()

  if status != 0:
      print >>sys.stderr, "External call failed."
      print >>sys.stderr, "  status: %d" % status
      print >>sys.stderr, "  stdout: %s" % out.read()
      print >>sys.stderr, "  stderr: %s" % err.read()
      print >>sys.stderr, "  command: %s" % command
      sys.exit(status)
  if show_stdout:
      print out.read()
  stdin.close()
  out.close()
  err.close()

##############################################################################
# Function to set up the segments for the analysis
def science_segments(ifo, config, opts):
  """
  generate the segments for the specified ifo
  """
  segFindFile = ifo + "-SCIENCE_SEGMENTS-" + str(opts.gps_start_time) + "-" + \
      str(opts.gps_end_time - opts.gps_start_time) + ".txt"

  # if not generating segments, all we need is the name of the segment file
  if not opts.generate_segments: return segFindFile

  executable = config.get("condor", "segfind")
  if executable[0] != "/": executable = "../" + executable

  # run segFind to determine science segments
  segFindCall = executable + " --interferometer=" + ifo + \
      " --type=\"" + config.get("segments", "analyze") + "\""\
      " --gps-start-time=" + str(opts.gps_start_time) + \
      " --gps-end-time=" + str(opts.gps_end_time) + " > " + segFindFile
  make_external_call(segFindCall)
  return segFindFile

##############################################################################
# Function to set up the segments for the analysis
def veto_segments(ifo, config, segmentList, dqSegFile, categories, opts):
  """
  generate veto segments for the given ifo
  
  ifo         = name of the ifo
  segmentList = list of science mode segments
  dqSegfile   = the file containing dq flags
  categories  = list of veto categories 
  """
  executable = config.get("condor", "query_dq")
  if executable[0] != "/": executable = "../" + executable

  vetoFiles = {}

  for category in categories:    
    dqFile = config.get("segments", ifo.lower() + "-cat-" + str(category) + \
        "-veto-file")
    if dqFile[0] != "/": dqFile = "../" + dqFile

    vetoFile = ifo + "-CATEGORY_" + str(category) + "_VETO_SEGS-" + \
        str(opts.gps_start_time) + "-" + \
        str(opts.gps_end_time - opts.gps_start_time) + ".txt"

    dqCall = executable + " --ifo " + ifo + " --dq-segfile " + dqSegFile + \
        " --segfile " + segmentList + " --flagfile " + dqFile + \
        " --outfile " + vetoFile

    # generate the segments
    make_external_call(dqCall)

    # if there are previous vetoes, generate combined
    try: previousSegs = \
        segmentsUtils.fromsegwizard(open(vetoFiles[category-1]))
    except: previousSegs = None

    if previousSegs:
      combinedFile = ifo + "-COMBINED_CAT_" + str(category) + "_VETO_SEGS-" + \
          str(opts.gps_start_time) + "-" + \
          str(opts.gps_end_time - opts.gps_start_time) + ".txt"

      vetoSegs = segmentsUtils.fromsegwizard(open(vetoFile)).coalesce()
      vetoSegs |= previousSegs
      segmentsUtils.tosegwizard(file(combinedFile,"w"), vetoSegs)
      vetoFiles[category] = combinedFile

    else: vetoFiles[category] = vetoFile

  return vetoFiles

