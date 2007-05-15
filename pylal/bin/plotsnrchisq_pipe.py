#!/usr/bin/python

import matplotlib
matplotlib.use('Agg')
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
sys.path.append('/archive/home/channa/opt/pylal/lib64/python2.4/site-packages')
from pylal import viz
from pylal import Fr
sys.path.append('/archive/home/channa/opt/glue/lib64/python2.4/site-packages')
from glue import segments
from glue import segmentsUtils
from glue.ligolw import ligolw
from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw import utils
from pylal import CoincInspiralUtils
from glue import pipeline
from glue import lal


#import matplotlib
#matplotlib.use('Agg')
#from pylab import *
#from pylab import *
#from pylal import viz


sys.path.append('@PYTHONLIBDIR@')

rc('text', usetex=False)

def plotsnrchisq(gpsTime,frameFile,outputPath,inspProcParams,tableFile):
 
    rsqThreshold = 0;

    for row in inspProcParams:
      if row.param == "--channel-name":
        chanStringBase = row.value
        ifoName = str(row.value).split(":",1)
      if row.param == "--segment-length":
        segLen = eval(row.value)
      if row.param == "--sample-rate":
        sampleRate = eval(row.value)
      if row.param == "--segment-overlap":
        segOverlap = eval(row.value)
      if (row.param == "--chisq-delta") or (row.param == "--minimal-match"):
        chisqDelta = eval(row.value)
      if row.param =="--chisq-bins":
        chisqBins = eval(row.value)
      if row.param == "--snr-threshold":
        snrThreshold = eval(row.value)
      if row.param == "--chisq-threshold":
        chisqThreshold = eval(row.value)
      if row.param == "--rsq-veto-threshold":
        rsqThreshold = eval(row.value)
      if row.param == "--trig-start-time":
        trigStart = eval(row.value)
      if row.param == "--gps-start-time":
        gpsStart = eval(row.value)

    segLenSec = segLen / sampleRate
    segOverlapSec = segOverlap / sampleRate
     
    if (trigStart):
      trigPosition = int((trigStart - gpsStart - segOverlapSec ) / (segLenSec -segOverlapSec))
    else: trigPosition = 0

    gpsPosition = int((eval(gpsTime) - gpsStart - segOverlapSec/2.) / (segLenSec -segOverlapSec))
    position = gpsPosition - trigPosition
    chanNumber = str(position)

    chanNameSnr = chanStringBase + "_SNRSQ_" + chanNumber
    chanNameChisq = chanStringBase + "_CHISQ_" + chanNumber

    # now, read the data !!
    # The window width should be an input argument maybe ?
    duration = 2.0
   
    squareSnr_tuple = Fr.frgetvect(frameFile,chanNameSnr,-1,segLenSec,0)
    
    squareChisq_tuple = Fr.frgetvect(frameFile,chanNameChisq,-1,segLenSec,0)

    snr_position = eval(gpsTime) - (gpsStart + gpsPosition* (segLenSec - segOverlapSec) )
    chisq_position = snr_position
    
    # compute the snr vector
    snr_vector = sqrt(squareSnr_tuple[0])
    # print squareSnr_tuple
    snr_time = array(range(0, segLen)) * squareSnr_tuple[3][0] - snr_position
      
    # compute the r^2
    rsq_vector = squareChisq_tuple[0]

    chisq_time = array(range(0, segLen)) * squareChisq_tuple[3][0] - chisq_position

    # compute the normalized chisq
    if(chisqBins > 0):
        chisqNorm_vector = rsq_vector/(1 + chisqDelta/chisqBins*squareSnr_tuple[0])
    else:
        chisqNorm_vector = rsq_vector
        
    # define the time boundaries of the plot
    snr_start =  (snr_position - duration/2.) * 1/squareSnr_tuple[3][0]
    snr_stop = (snr_position + duration/2.) * 1/squareSnr_tuple[3][0]
    lengthSnr = int(snr_stop) - int(snr_start)
	
    chisq_start =  (chisq_position - duration/2.) * 1/squareChisq_tuple[3][0]
    chisq_stop = (chisq_position + duration/2.) * 1/squareChisq_tuple[3][0]	
    lengthChisq = int(chisq_stop) - int(chisq_start) 

    #prepare the thresholds to be plotted 
    
    snrThreshVect = lengthSnr * [snrThreshold]
    if(chisqThreshold < 100.):	
    	chisqThreshVect = lengthChisq * [chisqThreshold]
    if((rsqThreshold > 0) & (rsqThreshold < 100.)):
	rsqThreshVect = lengthChisq * [rsqThreshold]
   

    # Now plot the snr time serie !!
    
    figure(1)
    plot(snr_time[int(snr_start):int(snr_stop)],snr_vector[int(snr_start):int(snr_stop)])    
    hold(1)
    plot(snr_time[int(snr_start):int(snr_stop)],snrThreshVect)
    hold(0)
    xlim(-duration/2., duration/2.)
    xlabel('time (s)',size='x-large')
    ylabel(r'snr',size='x-large')
    grid()
    title(ifoName[0] + ' trigger: ' + gpsTime)
    figName = ifoName[0] + '_' + str(gpsTime).replace(".","_") + '_snr.png'
    savefig(outputPath +"/" + figName)
    tableFile = open(container.locallink,'a')
    table = HTMLTable()
    table.add_column('<img width=400 src="' + figName +'">','SNR')

    figure(11)
    plot(snr_time[int(snr_start):int(snr_stop)],snr_vector[int(snr_start):int(snr_stop)])
    hold(1)
    plot(snr_time[int(snr_start):int(snr_stop)],snrThreshVect)
    hold(0)
    xlim(-duration/20., duration/20.)
    xlabel('time (s)',size='x-large')
    ylabel(r'snr',size='x-large')
    grid()
    title(ifoName[0] + ' trigger: ' + gpsTime + ' Zoom')
    figName = ifoName[0] + '_' + str(gpsTime).replace(".","_") + '_snr_zoom.png'
    savefig(outputPath +"/" + figName)
    table.add_column('<img width=400 src="' + figName +'">','SNR ZOOM')
    table.write(tableFile)
    tableFile.close()

    ### END CHADS CHANGES ###

    
    # Now plot the r^2 time serie !!
    figure(2)
    plot(chisq_time[int(chisq_start):int(chisq_stop)],rsq_vector[int(chisq_start):int(chisq_stop)])
    if(chisqThreshold < 100.):    
	hold(1)
    	plot(chisq_time[int(chisq_start):int(chisq_stop)],rsqThreshVect)
    	hold(0)
    xlim(-duration/2., duration/2.)
    xlabel('time (s)',size='x-large')
    ylabel(r'chisq/p',size='x-large')
    grid()
    title(ifoName[0] + ' trigger: ' + gpsTime)
    savefig(ifoName[0] + '_' + str(gpsTime).replace(".","_")  + '_rsq.png')
    	
    
    figure(22)
    plot(chisq_time[int(chisq_start):int(chisq_stop)],rsq_vector[int(chisq_start):int(chisq_stop)])
    if(chisqThreshold < 100.):
    	hold(1)
    	plot(chisq_time[int(chisq_start):int(chisq_stop)],rsqThreshVect)
        hold(0)
    xlim(-duration/20., duration/20.)
    xlabel('time (s)',size='x-large')
    ylabel(r'chisq/p',size='x-large')
    grid()
    title(ifoName[0] + ' trigger: ' + gpsTime + ' Zoom')
    savefig(ifoName[0] + '_' + str(gpsTime).replace(".","_")  + '_rsq_zoom.png')


    # Now plot the normalized chisq time serie !!
    figure(3)
    plot(chisq_time[int(chisq_start):int(chisq_stop)],chisqNorm_vector[int(chisq_start):int(chisq_stop)])
    if(rsqThreshold > 0 & rsqThreshold < 100.):
    	hold(1)
    	plot(chisq_time[int(chisq_start):int(chisq_stop)],chisqThreshVect)
        hold(0)
    xlim(-duration/2., duration/2.)
    xlabel('time (s)',size='x-large')
    ylabel(r'chisq / (p + Delta * snrsq)',size='x-large')
    grid(1)
    title(ifoName[0] + ' trigger: ' + gpsTime)
    savefig(ifoName[0] + '_' + str(gpsTime).replace(".","_")  + '_chisq.png')

    figure(33)
    plot(chisq_time[int(chisq_start):int(chisq_stop)],chisqNorm_vector[int(chisq_start):int(chisq_stop)])
    if(rsqThreshold > 0 & rsqThreshold < 100.):
    	hold(1)
    	plot(chisq_time[int(chisq_start):int(chisq_stop)],chisqThreshVect)
        hold(0)
    xlim(-duration/20., duration/20.)
    xlabel('time (s)',size='x-large')
    ylabel(r'chisq / (p + Delta * snrsq)',size='x-large')
    grid(1)
    title(ifoName[0] + ' trigger: ' + gpsTime + ' Zoom')
    savefig(ifoName[0] + '_' + str(gpsTime).replace(".","_")  + '_chisq_zoom.png')

    
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

parser.add_option("-x","--inspiral-xml-file", action="store",type="string", \
    metavar=" XML",help="use inspiral-file")

parser.add_option("-p","--output-html-file", action="store",type="string", \
    metavar=" XML",help="file to append html tables to")


command_line = sys.argv[1:]
(opts,args) = parser.parse_args()

#################################
# if --version flagged
if opts.version:
  print "$Id: plotsnrchisq_pipe.py,v 1.10 2007/05/14 21:43:51 channa Exp $"
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

doc = utils.load_filename(opts.inspiral_xml_file,None)
proc = table.get_table(doc, lsctables.ProcessParamsTable.tableName)
plotsnrchisq(str(opts.gps),opts.frame_file,opts.output_path,proc,opts.output_html_file)

