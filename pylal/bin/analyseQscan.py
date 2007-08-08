#!/usr/bin/env @PYTHONPROG@
"""
Something

$Id$

This program plots the background distribution of qscan auxialiary channels
"""

__author__ = 'Romain Gouaty <romain@phys.lsu.edu>'
__date__ = '$Date$'
__version__ = '$Revision$'[11:-2]

##############################################################################
# import standard modules and append the lalapps prefix to the python path

#import matplotlib
#matplotlib.use('Agg')

import matplotlib.cm
from matplotlib.patches     import Patch
from matplotlib.axes        import Axes
from matplotlib.collections import PolyCollection
from matplotlib.colors      import normalize, Colormap

import sys, getopt, os, copy, math
import socket, time
import re, string
from optparse import *
from pylab import *
import tempfile
import ConfigParser
import urlparse
from UserDict import UserDict
sys.path.append('@PYTHONLIBDIR@')

# rc('text', usetex=False)

##############################################################################
# import the modules we need to build the pipeline

from glue.ligolw import ligolw
from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw import utils

##############################################################################
# class to read qscan summary files

class readSummaryFiles:

  def __init__(self):
    self.table = {'channel_name':[],'peak_time':[],'peak_frequency':[],'peak_q':[],'peak_significance':[],'peak_amplitude':[],'qscan_id':[]}

  def getAuxChannels(self,inputPath):
    id = 0
    listDir = os.listdir(inputPath)
    for dir in listDir:
      try:
        doc = utils.load_filename(inputPath + "/" + dir + "/summary.xml")
        qscanTable = table.get_table(doc, "qscan:summary:table")
      except: 
        print >> sys.stderr, "failed to read" + inputPath + "/" + dir + "/summary.xml"
        continue
      for channel in qscanTable:
        self.table['channel_name'].append(channel.channelName)
        self.table['peak_time'].append(channel.peakTime)
        self.table['peak_frequency'].append(channel.peakFrequency)
        self.table['peak_q'].append(channel.peakQ)
        self.table['peak_significance'].append(channel.peakSignificance)
        self.table['peak_amplitude'].append(channel.peakAmplitude)
        self.table['qscan_id'].append(id)
      id = id+1
    return self.table      

##############################################################################
# functions to compute and plot significance histograms

def makeHistogramSignificance(table,chan,cp):
  
  z_list = []
  ind_i = 0
  counter = 0
  tableLength = len(table['channel_name'])
  n_chan = table['channel_name'].count(chan)
  
  for i in range(0,n_chan,1):
    ind_i = ind_i + table['channel_name'][ind_i:tableLength].index(chan)
    z = table['peak_significance'][ind_i]
    if z > eval(cp.get('z-distribution','z-threshold')):
      z_list.append(z)
      counter = counter + 1
    ind_i = ind_i + 1

  # set up the bin boundaries for the histogram
  min_val = eval(cp.get('z-distribution','z-min'))
  max_val = eval(cp.get('z-distribution','z-max'))
  nbins = eval(cp.get('z-distribution','z-nbins'))

  step = (float(max_val) - float(min_val))/float(nbins) 
  bins = arange(min_val, max_val, step)
  # bins = range(-1,int(eval(cp.get('z-distribution','z-max')))+1,2)

  if n_chan > 0:
    if len(z_list):
      # compute the histogram values
      [z_dist, bin, info] = hist(z_list,bins,bottom=None,\
      align='edge', orientation='vertical', width=None)
    else:
      print >> sys.stderr, 'no z above threshold for channel ' + chan
  else:
    print >> sys.stderr, 'Channel ' + chan + ' not found'
  
  return z_dist,bin


def saveHistogramSignificance(chan,cp,z_dist,bin):
     
  histoFileName = chan.split(':')[0] + '_' + chan.split(':')[1] + '_z_distribution.txt'
  histoFilePath = string.strip(cp.get('background-output','output-path')) + '/' + histoFileName
  store = open(histoFilePath,'w')
  for i,z in enumerate(z_dist):
    store.write(str(z) + '\t' + str(bin[i]) + '\n')


def readHistogramSignificance(chan,cp):
  
  testOpenFile = 0
  histoFileList = []

  histoFileName = chan.split(':')[0] + '_' + chan.split(':')[1] + '_z_distribution.txt'
  histoFilePath = string.strip(cp.get('background-output','output-path')) + '/' + histoFileName
  try:  histoFile = open(histoFilePath,'r')
  except: 
    print >> sys.stderr, 'could not open ' + histoFilePath
    testOpenFile = 1

  if not testOpenFile:
    histoFileList = histoFile.readlines()
    if not len(histoFileList):
      print >> sys.stderr, 'The file ' + histoFilePath + ' is empty !'
    else:
      histoList = []
      binList = []
      for bin in histoFileList:
        bin_temp = string.strip(bin)
        histoList.append(eval(bin_temp.split('\t')[0]))
        binList.append(eval(bin_temp.split('\t')[1]))

  return histoList,binList


def plotHistogramSignificance(chan,cp,histoList,binList,figNumber):
      
  step = binList[1] - binList[0]
  counter = sum(histoList)
        
  figure(figNumber)
  # semilogy(bins + step/2., z_dist+0.0001, 'r^',markerfacecolor="b",markersize=12)
  # plot(bins + step/2., z_dist)
  bar(binList, histoList, width=step, bottom=0)

  xlabel('Z value',size='large')
  # ylabel(r'#',size='x-large')
  grid()  
  title("Histogram of the Z value for " + chan + ', Statistics = ' + str(counter))

  figName = chan.split(':')[0] + '_' + chan.split(':')[1] + '_z_distribution.png'
  figFileName = string.strip(cp.get('background-output','output-path')) + '/' + figName
  savefig(figFileName) 

##############################################################################
#
#  MAIN PROGRAM
#
##############################################################################

######################## OPTION PARSING  #####################################
usage = """usage: %prog [options]
"""

parser = OptionParser( usage )

parser.add_option("-v", "--version",action="store_true",default=False,\
    help="print version information and exit")

parser.add_option("-f", "--config-file",action="store",type="string",\
    metavar=" FILE",help="ini file")

parser.add_option("-z", "--make-z-distribution",action="store_true",\
    default=False,help="compute the z distributions")

parser.add_option("-Z", "--plot-z-distribution",action="store_true",\
    default=False,help="plot the z distributions")

command_line = sys.argv[1:]
(opts,args) = parser.parse_args()

if opts.version:
  print "$ Id: makeQscanBackground.py,v 1.0 2007/08/02 18:46:00 romain Exp $"
  sys.exit(0)

####################### SANITY CHECKS #####################################

if not opts.config_file:
  print >> sys.stderr, "No configuration file specified."
  print >> sys.stderr, "Use --config-file FILE to specify location"
  sys.exit(1)

if not opts.plot_z_distribution and not opts.make_z_distribution:
  print >> sys.stderr, "No step of the pipeline specified"
  print >> sys.stderr, "Please specify at least one of"
  print >> sys.stderr, "--make-z-distribution, --plot-z-distribution"
  sys.exit(1)

#################### READ IN THE CONFIG (.ini) FILE ########################
cp = ConfigParser.ConfigParser()
cp.read(opts.config_file)

################# NOW START THE REAL WORK ##################################

# Get the list of qscan channels to be analyzed (should become optional later...). The script needs to be improved to follow this behavior: if the text file is not specified in the configuration file, then all the channels found in the summary files should be analyzed...
if len(string.strip(cp.get('qscan-summary','channel-list'))) > 0:
  channelFile = open(string.strip(cp.get('qscan-summary','channel-list')),'r')
  channel_list = []
  channel_list = channelFile.readlines()
  if not len(channel_list):
    print >> sys.stderr, "No channel found in the channel_list file"
    print >> sys.stderr, "Is the first line blank?"
    sys.exit(1)
  channelList = []
  for chan in channel_list:
    channelList.append(string.strip(chan))

# Read the qscan summary files and hold the information in memory
if opts.make_z_distribution:
  summaryFiles = readSummaryFiles()
  qscanTable = summaryFiles.getAuxChannels(string.strip(cp.get('qscan-summary','input-path')))

  # perform a sanity check
  if not (len(qscanTable['channel_name']) == len(qscanTable['qscan_id'])):
    print >> sys.stderr, "the length of channel_name does not match the length of qscan_id"
    print >> sys.stderr, "check for data corruption in the qscan summary files"
    sys.exit(1)

# prepare and plot the distribution of significance
if (opts.make_z_distribution or opts.plot_z_distribution) and len(string.strip(cp.get('qscan-summary','channel-list'))) > 0:
  figNumber = 0
  for channel in channelList:

    if opts.make_z_distribution:
      try: 
        zHisto,zBin = makeHistogramSignificance(qscanTable,channel,cp)
        testMakeHisto = 0
      except: 
        print 'could not make histogram for channel ' + channel
        testMakeHisto = 1
      if not testMakeHisto:
        saveHistogramSignificance(channel,cp,zHisto,zBin)
        if opts.plot_z_distribution:
          figNumber = figNumber + 1
          plotHistogramSignificance(channel,cp,zHisto,zBin,figNumber)
         
    if opts.plot_z_distribution and not opts.make_z_distribution:
      try: 
        zHisto,zBin = readHistogramSignificance(channel,cp)
        testReadHisto = 0
      except: 
        print 'failed to read the txt file containing information for ' + channel + ' background'
        testReadHisto = 1
      if not testReadHisto:
        figNumber = figNumber + 1
        plotHistogramSignificance(channel,cp,zHisto,zBin,figNumber)

