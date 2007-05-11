#!/usr/bin/python

import copy
import math
import sys
import os
import socket, time
from optparse import *
import re, string
import exceptions
import glob
import tempfile
import ConfigParser
import urlparse
from types import *
from pylab import *
from pylal import frgetvect
from pylal import viz
sys.path.append('@PYTHONLIBDIR@')


##############################################################################
# function to plot the SNR and CHISQ time series

def plotsnrchisq(gpsTime,frameFile,ligoName,chisqBins,chisqDelta,outputPath):
  
    # BE CAREFUL ! The length of the frame file itself (without taking into account the path) is supposed to be 30 characters !
    ind1 = len(frameFile) - 30  # this variable is also used in the definition of the plot title
    ind2 = ind1 + 2
    ifoName = frameFile[ind1:ind2]
    print ifoName

    chanString = ifoName + ':' + ligoName + '_SNRSQ_0'
    print chanString

    # find the start time of the first channel
    # BE CAREFUL ! it is assumed that the sampling frequency is higher than 200 Hz
    testOnFirstChannel = frgetvect.frgetvect(frameFile,chanString,-1,0.01,0)
    gpsStart = testOnFirstChannel[3]
    
    # This actually prints only one digit after the .
    print gpsStart

    # find the channel which contains the data we want to look at
    # BE CAREFUL ! it is assumed that the segment length is 128 s
    segmentLength = 128.
    position = (float(gpsTime) - float(gpsStart) - segmentLength/2.)/segmentLength
    position = int(position)
    chanNumber = str(position)
    chanNameSnr = ifoName + ':' + ligoName + '_SNRSQ_' + chanNumber
    print chanNameSnr
    chanNameChisq = ifoName + ':' + ligoName + '_CHISQ_' + chanNumber

    # now, read the data !!
    # The window width should be an input argument maybe ?
    duration = 2.0
    startWindow = float(gpsTime) - duration/2.
    squareSnr_tuple = frgetvect.frgetvect(frameFile,chanNameSnr,startWindow,duration,0)
    print squareSnr_tuple[0]

    squareChisq_tuple = frgetvect.frgetvect(frameFile,chanNameChisq,startWindow,duration,0)

    # compute the snr vector
    snr_vector = sqrt(squareSnr_tuple[0])
    print snr_vector
    snr_time = squareSnr_tuple[1]
    print snr_time

    # compute the r^2
    rsq_vector = squareChisq_tuple[0]
    chisq_time = squareChisq_tuple[1]
    print rsq_vector

    # compute the normalized chisq
    chisqNorm_vector = rsq_vector/(1 + chisqDelta/chisqBins*squareSnr_tuple[0])
    print chisqNorm_vector
    
    
    # Now plot the snr time serie !!
    figure(1)
    plot(snr_time - duration/2.,snr_vector)
    xlabel('time (s)',size='x-large')
    ylabel(r'$\rho$',size='x-large')
    grid(1)
    title(ifoName + ' trigger: ' + gpsTime)   
    savefig(ifoName + '_' + str(int(float(gpsTime))) + '_snr.png')

    # Now plot the r^2 time serie !!
    figure(2)
    plot(chisq_time - duration/2.,rsq_vector)
    xlabel('time (s)',size='x-large')
    ylabel(r'$r^2$',size='x-large')
    grid(1)
    title(ifoName + ' trigger: ' + gpsTime)   
    savefig(ifoName + '_' + str(int(float(gpsTime))) + '_rsq.png')

    # Now plot the normalized chisq time serie !!
    figure(3)
    plot(chisq_time - duration/2.,chisqNorm_vector)
    xlabel('time (s)',size='x-large')
    ylabel(r'$\chi^2 / (p + \delta^2\rho^2)$',size='x-large')
    grid(1)
    title(ifoName + ' trigger: ' + gpsTime)   
    savefig(ifoName + '_' + str(int(float(gpsTime))) + '_chisq.png')

    
##############################################################################
#
#  MAIN PROGRAM
#
##############################################################################
usage = """ %prog [options]
"""

parser = OptionParser( usage )

parser.add_option("-v","--version",action="store_true",default=False,\
    help="display version information and exit")

parser.add_option("-f","--frame-file",action="store",type="string",\
    metavar=" FILE",help="use frame files list FILE")

parser.add_option("-t","--gps",action="store",type="string",\
    metavar=" GPS",help="use gps time GPS")

parser.add_option("-o","--output-path",action="store",type="string",\
    metavar=" PATH",help="use output path PATH")

parser.add_option("-l","--ligo-channel",action="store",type="string",\
    metavar=" CHANNEL",help="use ligo channel CHANNEL")

parser.add_option("-b","--chisq-dof",action="store",type="float",\
    metavar=" BINS",help="use chisq bins BINS")

parser.add_option("-d","--chisq-delta",action="store",type="float",\
    metavar=" DELTA",help="use chisq delta DELTA")


command_line = sys.argv[1:]
(opts,args) = parser.parse_args()

#################################
# if --version flagged
if opts.version:
  print "$Id: plotsnrchisq.py,v 1.1 2007/04/28 22:23:19 channa Exp $"
  sys.exit(0)

#################################
# Sanity check of input arguments

if not opts.frame_file:
  print >> sys.stderr, "No frame file specified."
  print >> sys.stderr, "Use --frame-file FILE to specify location."
  sys.exit(1)

if not opts.gps:
  print >> sys.stderr, "No gps time specified."
  print >> sys.stderr, "Use --gps GPS to specify location."
  sys.exit(1)

if not opts.output_path:
  print >> sys.stderr, "No output path time specified."
  print >> sys.stderr, "Use --output-path PATH to specify location."
  sys.exit(1)

if not opts.ligo_channel:
  print >> sys.stderr, "No ligo_channel specified."
  print >> sys.stderr, "Use --ligo-channel CHANNEL to specify location."
  sys.exit(1)  

if not opts.chisq_dof:
  print >> sys.stderr, "No chisq bins specified."
  print >> sys.stderr, "Use --chisq-dof BINS to specify location."
  sys.exit(1)
  
if not opts.chisq_delta:
  print >> sys.stderr, "No chisq delta specified."
  print >> sys.stderr, "Use --chisq-delta DELTA to specify location."
  sys.exit(1)

#################################

# if we want to take as insput a list of gwf files

# frameFileList = open(opts.frame_files_list,"r")
# frameFiles = []
# frameFiles = frameFileList.readlines()

# if not len(frameFiles):
  # print >> sys.stderr, "No frame file found in FILE"
  # print >> sys.stderr, "Is the first line blank?"
  # sys.exit(1)

# frameFile = []

# for file in frameFiles:
  # file = file.replace ( "\n", "" )
  # print file
  # frameFile.append(file)

plotsnrchisq(opts.gps,opts.frame_file,opts.ligo_channel,opts.chisq_dof,opts.chisq_delta,opts.output_path)

#################################
