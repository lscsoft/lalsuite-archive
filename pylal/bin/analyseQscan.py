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
# function to check the length of the summary files (for debugging)

def checkSummaryLength(inputPath): # used for debugging only

  outputPath = string.strip(cp.get('output','output-path')) + '/' + 'qscan_length.txt'
  storeLength = open(outputPath,'w')

  listDir = os.listdir(inputPath)
  listDir.sort()
  counter = 0
  for dir in listDir:
    try:
      summary = open(inputPath + "/" + dir + "/summary.txt","r")
      summary_length = len(summary.readlines())
      counter = counter + 1
      storeLength.write(str(counter) + "\t" + dir + "\t" + str(summary_length) + "\n")
    except:
      print >> sys.stderr, "could not check file length for" + inputPath + "/" + dir + "/summary.txt"
      continue
  storeLength.close()

##############################################################################
# class to read qscan summary files

class readSummaryFiles:

  def __init__(self):
    self.table = {'channel_name':[],'peak_time':[],'peak_frequency':[],'peak_q':[],'peak_significance':[],'peak_amplitude':[],'qscan_time':[]}


  def getAuxChannels(self,inputPath):
    #id = 0
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
	self.table['qscan_time'].append(dir)
        #self.table['qscan_id'].append(id)
      #id = id+1
    return self.table      

##############################################################################
# functions to compute and plot histograms

def buildChannelList(table,channel,param1,param2=None):
  ind_j = 0
  dict = {param1:[]}
  if param2: dict[param2]=[]

  tableLength = len(table['channel_name'])
  n_chan = table['channel_name'].count(channel)

  for i in range(0,n_chan,1):
    ind_j = ind_j + table['channel_name'][ind_j:tableLength].index(channel)
    element1 = table[param1][ind_j]
    dict[param1].append(element1)
    if param2:
      element2 = table[param2][ind_j]
      dict[param2].append(element2)
    ind_j = ind_j + 1

  if not n_chan > 0:
    print >> sys.stderr, 'Channel ' + channel + ' not found'

  return dict


def printChannelList(table,channel,param1,param2): # for debugging only

  chan_dict = buildChannelList(table,channel,param1,param2)
  fileName = channel.split(':')[0] + '_' + channel.split(':')[1] + '_' + param1 + '_' + param2 + '.txt'
  filePath = string.strip(cp.get('output','output-path')) + '/' + fileName
  store = open(filePath,'w')
  for i,value in enumerate(chan_dict[param1]):
    store.write(str(value) + '\t' + str(chan_dict[param2][i]) + '\n')
  store.close()


def computeDeltaT(table,chan,cp):
  dt_list = []
  
  if cp.has_option('dt-distribution','t-reference'):
    if len(string.strip(cp.get('dt-distribution','t-reference'))) > 0:
      darm_dict = buildChannelList(table,string.strip(cp.get('dt-distribution','t-reference')),'peak_time','qscan_time')

  chan_dict = buildChannelList(table,chan,'peak_time','qscan_time')  

  for i,time in enumerate(chan_dict['peak_time']):
    if time > 0.0:
      if not cp.has_option('dt-distribution','t-reference') or \
      (cp.has_option('dt-distribution','t-reference') and \
      len(string.strip(cp.get('dt-distribution','t-reference'))) == 0):
        dt = time - float(chan_dict['qscan_time'][i])
        dt_list.append(dt)
      else:
        if len(darm_dict['peak_time']) > 0:
          try: 
            indx_ref = darm_dict['qscan_time'].index(chan_dict['qscan_time'][i])
            if float(darm_dict['peak_time'][indx_ref]) > 0.0:
              dt = time - float(darm_dict['peak_time'][indx_ref])
              dt_list.append(dt)
          except:
            print >> sys.stderr, 'GW channel could not be found in qscan ' + chan_dict['qscan_time'][i]
            continue

  return dt_list


def makeHistogram(list,distribution,chan,cp):

  # set up the bin boundaries for the histogram
  min_val = eval(cp.get(distribution,'min'))
  max_val = eval(cp.get(distribution,'max'))
  nbins = eval(cp.get(distribution,'nbins'))

  step = (float(max_val) - float(min_val))/float(nbins)
  bins = arange(min_val, max_val, step)

  if len(list):
    # compute the histogram values
    [dist, bin, info] = hist(list,bins,bottom=None,\
    align='edge', orientation='vertical', width=None)

  return dist,bin


def selectSignificance(table,chan,cp):

  z_list = []
  chan_dict = buildChannelList(table,chan,'peak_significance')

  for z in chan_dict['peak_significance']:
    if z > eval(cp.get('z-distribution','z-threshold')):
     z_list.append(z)

  return z_list


def saveHistogramSignificance(chan,cp,z_dist,bin):
     
  histoFileName = chan.split(':')[0] + '_' + chan.split(':')[1] + '_z_distribution.txt'
  histoFilePath = string.strip(cp.get('output','output-path')) + '/' + histoFileName
  store = open(histoFilePath,'w')
  for i,z in enumerate(z_dist):
    store.write(str(z) + '\t' + str(bin[i]) + '\n')
  store.close()

def readHistogramSignificance(chan,cp):
  
  testOpenFile = 0
  histoFileList = []

  histoFileName = chan.split(':')[0] + '_' + chan.split(':')[1] + '_z_distribution.txt'
  histoFilePath = string.strip(cp.get('output','output-path')) + '/' + histoFileName
  try:  histoFile = open(histoFilePath,'r')
  except: 
    print >> sys.stderr, 'could not open ' + histoFilePath
    testOpenFile = 1

  if not testOpenFile:
    histoFileList = histoFile.readlines()
    histoFile.close()
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


def plotHistogram(chan,cp,distribution,histoList,binList,figNumber):

  parameter = distribution.split('-')[0]
      
  step = binList[1] - binList[0]
  counter = sum(histoList)
        
  figure(figNumber)
  # semilogy(bins + step/2., z_dist+0.0001, 'r^',markerfacecolor="b",markersize=12)
  # plot(bins + step/2., z_dist)
  bar(binList, histoList, width=step, bottom=0)

  xlabel(parameter + ' value',size='large')
  # ylabel(r'#',size='x-large')
  grid()  
  title("Histogram of the " + parameter + " value for " + chan + ', Statistics = ' + str(counter))

  figName = chan.split(':')[0] + '_' + chan.split(':')[1] + '_' + parameter + '_distribution.png'
  figFileName = string.strip(cp.get('output','output-path')) + '/' + figName
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

parser.add_option("-l", "--check-length",action="store_true",\
    default=False,help="check the length of the summary txt files")

parser.add_option("-c", "--create-param-list",action="store_true",\
    default=False,help="create .txt files containing the list of parameters for each channel (for debugging only)")

parser.add_option("-z", "--make-z-distribution",action="store_true",\
    default=False,help="compute the z distributions")

parser.add_option("-Z", "--plot-z-distribution",action="store_true",\
    default=False,help="plot the z distributions")

parser.add_option("-T", "--plot-dt-distribution",action="store_true",\
    default=False,help="plot the delta t distributions")

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

if not opts.plot_z_distribution and not opts.make_z_distribution and not opts.check_length and not opts.plot_dt_distribution and not opts.make_dt_distribution and not opts.create_param_list:
  print >> sys.stderr, "No step of the pipeline specified"
  print >> sys.stderr, "Please specify at least one of"
  print >> sys.stderr, "--make-z-distribution, --plot-z-distribution"
  print >> sys.stderr, "--plot-dt-distribution"
  print >> sys.stderr, "--check-length, --create-param-list"
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
  channelFile.close()
  if not len(channel_list):
    print >> sys.stderr, "No channel found in the channel_list file"
    print >> sys.stderr, "Is the first line blank?"
    sys.exit(1)
  channelList = []
  for chan in channel_list:
    channelList.append(string.strip(chan))

# Check the summary length
if opts.check_length:
  checkSummaryLength(string.strip(cp.get('qscan-summary','input-path')))


# Read the qscan summary files and hold the information in memory
if opts.make_z_distribution or opts.plot_dt_distribution:
  summaryFiles = readSummaryFiles()
  qscanTable = summaryFiles.getAuxChannels(string.strip(cp.get('qscan-summary','input-path')))

  # perform a sanity check
  if not (len(qscanTable['channel_name']) == len(qscanTable['qscan_time'])):
    print >> sys.stderr, "the length of channel_name does not match the length of qscan_time"
    print >> sys.stderr, "check for data corruption in the qscan summary files"
    sys.exit(1)
 
  # Display the list of parameters for each channel
  if opts.create_param_list:
    for channel in channelList:
      printChannelList(qscanTable,channel,'peak_time','qscan_time')
      printChannelList(qscanTable,channel,'peak_significance','qscan_time')

if len(string.strip(cp.get('qscan-summary','channel-list'))) > 0:
  
  figNumber = 0
  
  # prepare and plot the distribution of delta t
  if (opts.plot_dt_distribution):
    for channel in channelList:
      dtList = computeDeltaT(qscanTable,channel,cp)
      try:
        dtHisto,dtBin = makeHistogram(dtList,'dt-distribution',channel,cp)
      except:
        print >> sys,stderr, 'could not make the dt histogram for channel ' + channel 
        continue
      figNumber = figNumber + 1
      plotHistogram(channel,cp,'dt-distribution',dtHisto,dtBin,figNumber)
 
  # prepare and plot the distribution of significance
  if (opts.make_z_distribution or opts.plot_z_distribution):
    for channel in channelList:
      if opts.make_z_distribution:
        zList = selectSignificance(qscanTable,channel,cp)
        try:
          zHisto,zBin = makeHistogram(zList,'z-distribution',channel,cp)
        except:
          print >> sys,stderr, 'could not make the z histogram for channel ' + channel
          continue
        saveHistogramSignificance(channel,cp,zHisto,zBin)
        if opts.plot_z_distribution:
          figNumber = figNumber + 1
          plotHistogram(channel,cp,'z-distribution',zHisto,zBin,figNumber)
         
      if opts.plot_z_distribution and not opts.make_z_distribution:
        try: 
          zHisto,zBin = readHistogramSignificance(channel,cp)
          testReadHisto = 0
        except: 
          print 'failed to read the txt file containing information for ' + channel + ' z background'
          testReadHisto = 1
        if not testReadHisto:
          figNumber = figNumber + 1
          plotHistogram(channel,cp,'z-distribution',zHisto,zBin,figNumber)

