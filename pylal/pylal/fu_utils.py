import sys
import os
import copy
import re
import exceptions
import glob
import fileinput
import linecache
import string
import random
from optparse import *
from types import *
import matplotlib
matplotlib.use('Agg')

from pylab import *
from glue import segments
from glue import segmentsUtils
from glue.ligolw import ligolw
from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw import utils
from pylal import CoincInspiralUtils

#############################################################################
# function to set up directories
#############################################################################
def setupdirs():

  try:
    os.chdir("clustered")
    os.chdir('..')
  except: os.mkdir("clustered")

  try:
    os.chdir("found")
    os.chdir('..')
  except: os.mkdir("found")

  try:
    os.chdir("missed")
    os.chdir('..')
  except: os.mkdir("missed")

  try:
    os.chdir("followuptrigs")
    os.chdir('..')
  except: os.mkdir("followuptrigs")

  try:
    os.chdir("followupfound")
    os.chdir('..')
  except: os.mkdir("followupfound")

  try:
    os.chdir("followupmissed")
    os.chdir('..')
  except: os.mkdir("followupmissed")


#############################################################################
# function to return the number of slides in a file (as a string)
#############################################################################
def getslidenum(fName):
  command = "grep -m 1 -e '--num-slides' " + fName + " | sed -e 's@\"@\'@g"
  fStr = os.popen(command).readlines()
  if fStr:
    fStrlist = fStr[0].rsplit(",")
    return fStrlist[-2]
  else:
   return "0"

#############################################################################
# function to glob files
#############################################################################
def globxmlandinj(opts):

  #Glob for the files both injections and trigger files
  fList = glob.glob(opts.xml_glob)
  if not fList:
    print >>sys.stderr, "The glob for " + opts.xml_glob + " returned no files"
    sys.exit(1)

  if opts.injection_glob:
    iList = glob.glob(opts.injection_glob)
    if not iList:
      print >> sys.stderr, "The glob for " + opts.injection_glob + " returned no files"
      sys.exit(1)
  else:
    iList = None

  return fList,iList

##############################################################################
# function to coire the input files
##############################################################################
def coire(opts,fList,iList=None):
  """
  do coires
  """
  sims = None
  coincs = None

  # we need to get slide num and data type from the files
  # this is just hard coded now!!!
  for xmls in fList:
    if opts.verbose:
      print "running lalapps_coire on " + xmls
    numslides = getslidenum(xmls)
    command = "lalapps_coire " + \
              "--glob '" + xmls + "'" + \
              " --data-type " + str(opts.data_type) + \
              " --num-slides " + numslides + " --coinc-stat " + \
              opts.statistic + " --cluster-time " + \
              str(1000*opts.cluster_time)

    if iList:
      command += " --injection-file " + iList[0]
      command += " --injection-window " + str(opts.injection_coinc_time)
      command += " --missed-injections missed/MISSED" + xmls
      command += " --output found/FOUND" + xmls
    else:
      command += " --output clustered/CLUSTER"
      command += str(opts.cluster_time) + "s_" + xmls
    if opts.bitten_l_a and opts.bitten_l_b:
      command += " --h1-bittenl-a " + str(opts.bitten_l_a)
      command += " --h2-bittenl-a " + str(opts.bitten_l_a)
      command += " --l1-bittenl-a " + str(opts.bitten_l_a)
      command += " --h1-bittenl-b " + str(opts.bitten_l_b)
      command += " --h2-bittenl-b " + str(opts.bitten_l_b)
      command += " --l1-bittenl-b " + str(opts.bitten_l_b)

    os.system(command)

##############################################################################
# function to extract the statistic information
##############################################################################
def getstatistic(opts):

  if opts.statistic == "effective_snrsq":
    newstat = "effective_snr"
  else:
    newstat = opts.statistic

  statistic=CoincInspiralUtils.coincStatistic( newstat, opts.bitten_l_a,
                                               opts.bitten_l_b)
  return statistic


#############################################################################
# Follow up list class definition
#############################################################################
class followUpList:
  """
  Class to hold gps and ifo pairs to send to subsequent functions
  It also holds an instance of the coinc class which contains
  All of the relevant xml information
  """
  def __init__(self,Coincs = None, Missed = None ):
    self.gpsTime = {"H1" : None, "H2" : None, "L1" : None,
                  "G1" : None, "V1" : None, "T1" : None}
    self.coincs = Coincs
    self.missed = Missed
    self.eventID = None
    self.stat = None
  def add_coincs(self,Coincs):
    setattr(self,"coincs",Coincs)
    self.eventID = Coincs.event_id
    self.statValue = Coincs.stat
  def add_missed(self,Missed):
    setattr(self,"missed",Missed)
  def is_trigs(self):
    if isinstance(self.coincs,CoincInspiralUtils.coincInspiralTable.row):
      return 1
  def is_found(self):
    if isinstance(self.coincs.sim,lsctables.SimInspiral):
      return 1
  def is_missed(self):
    if isinstance(self.missed,lsctables.SimInspiralTable):
      return 1

#############################################################################
# Function to return the follow up list of coinc triggers
#############################################################################
def getfollowuptrigs(opts,coincs=None,missed=None):

  followups = []

  if coincs:
    coincs.sort()
    numTrigs = 0
    for ckey in coincs:
      fuList = followUpList()
      fuList.add_coincs(ckey)
      try:
        getattr(ckey,'H1')
        fuList.gpsTime["H1"] = (float(getattr(ckey,'H1').end_time_ns)/1000000000)+float(getattr(ckey,'H1').end_time)
      except: fuList.gpsTime["H1"] = None
      try:
        getattr(ckey,'H2')
        fuList.gpsTime["H2"] = (float(getattr(ckey,'H2').end_time_ns)/1000000000)+float(getattr(ckey,'H2').end_time)
      except: fuList.gpsTime["H2"] = None
      try:
        getattr(ckey,'L1')
        fuList.gpsTime["L1"] = (float(getattr(ckey,'L1').end_time_ns)/1000000000)+float(getattr(ckey,'L1').end_time)
      except: fuList.gpsTime["L1"] = None
      try:
        getattr(ckey,'G1')
        fuList.gpsTime["G1"] = (float(getattr(ckey,'G1').end_time_ns)/1000000000)+float(getattr(ckey,'G1').end_time)
      except: fuList.gpsTime["G1"] = None
      try:
        getattr(ckey,'V1')
        fuList.gpsTime["V1"] = (float(getattr(ckey,'V1').end_time_ns)/1000000000)+float(getattr(ckey,'V1').end_time)
      except: fuList.gpsTime["V1"] = None
      try:
        getattr(ckey,'T1')
        fuList.gpsTime["T1"] = (float(getattr(ckey,'T1').end_time_ns)/1000000000)+float(getattr(ckey,'T1').end_time)
      except: fuList.gpsTime["T1"] = None
      followups.append(fuList)
      numTrigs += 1
      if numTrigs >= opts.num_trigs:
        break

  # the missed stuff doesnt work yet!!!
  if missed:
    followups

  return followups

#############################################################################
# Class to hold summary HTML information for all of the functions
#############################################################################
class summaryHTMLTable:

  def __init__(self,trig):
    if trig.is_trigs():
      self.eventID = trig.eventID
      self.statValue = trig.statValue
    if trig.is_trigs() and not trig.is_found():
      self.summarypath = "followuptrigs/"
    if trig.is_trigs() and trig.is_found():
      self.summarypath = "followupfound/"
    else:
      self.eventID = None
      self.statValue = None
    self.H1time = trig.gpsTime["H1"]
    self.H2time = trig.gpsTime["H2"]
    self.L1time = trig.gpsTime["L1"]
    self.G1time = trig.gpsTime["G1"]
    self.V1time = trig.gpsTime["V1"]
    self.T1time = trig.gpsTime["T1"]
    self.containers = []

class HTMLcontainer:

  def __init__(self,trig,name):
    # The missed injections dont work yet!!!
    self.name = name.rsplit(".")[-1]

    if trig.is_trigs() and not trig.is_found():    
      os.chdir("followuptrigs")
      try: 
        os.chdir(self.name)
        os.chdir("../../")
      except: 
        os.mkdir(self.name)
        os.chdir("../")
      self.detailpath = "followuptrigs/" + self.name + "/"
    else: 
      self.detailpath = ""

    if trig.is_trigs() and trig.is_found():
      os.chdir("followupfound")
      try:
        os.chdir(self.name)
        os.chdir("../../")
      except:
        os.mkdir(self.name)
        os.chdir("../")
      self.detailpath = "followupfound/" + self.name + "/"
    else:
      self.detailpath = ""

    self.image = self.detailpath + str(trig.statValue) + "_" + str(trig.eventID) + "_" + self.name + ".png" 
    self.text = "click here"
    self.link = self.detailpath + str(trig.statValue) + "_" + str(trig.eventID) + "_" + self.name + ".html"




##############################################################################
# Function to write the HTML tables to pages
##############################################################################
def beginSummaryTable(file, table):
  file.write("<h2>Trigger [" + str(table.eventID) + 
             "] with combined statistic = " + str(table.statValue) + "</h2>")
  file.write('\n<br><table width=1000 border=1>')
  
def endSummaryTable(file, table):
  file.write('\n</table>')

def writeModule(file, container):
  file.write('\n<tr><td width=600><font color="red" size=5> MODULE NAME '+
             container.name + '</font>')
  file.write('\n<br><a href= "' + container.link + '">' + container.text 
             + "</a></td>")
  file.write('\n<td><a href= "' + container.image + '">\n<img src="' + 
              container.image + '" width=400></a></td></tr>')

def writeHTMLTables(summaryHTMLlist):
  for table in summaryHTMLlist:
    tableFile = open(table.summarypath + str(table.statValue) + "_" + 
                      str(table.eventID) + "_summary.html" ,'w')
    beginSummaryTable(tableFile,table)
    for container in table.containers:
      writeModule(tableFile, container)
    endSummaryTable(tableFile,table)
    tableFile.close()

