#!/usr/bin/env @PYTHONPROG@
"""
followup utilities

$Id$

This
"""

__author__ = 'Chad Hanna <channa@phys.lsu.edu>'
__date__ = '$Date$'
__version__ = '$Revision$'[11:-2]



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
import operator
from UserDict import UserDict
import operator

from pylab import *
from glue import segments
from glue import segmentsUtils
from glue.ligolw import ligolw
from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw import utils
from pylal import CoincInspiralUtils
from pylal import itertools # For the 2nd year analysis, this has to be replaced by
# from glue import iterutils
from glue import pipeline
from glue.lal import *
from glue import lal
from lalapps import inspiralutils

########## CLASS TO WRITE LAL CACHE FROM HIPE OUTPUT #########################
class getCache(UserDict):
  """
  An instance of a lal cache
  """
  def __init__(self, options,cp, currentPath):
    UserDict.__init__(self)
    self.dir = os.listdir(os.getcwd())
    self.options = options
    self.types = ['TMPLTBANK', 'TRIGBANK', 'INSPIRAL-', \
                 'INSPIRAL_H', 'THINCA-', 'THINCA_']
    self.iniNames = ['tmpltbank-path', 'trigbank-path', 'first-inspiral-path', \
         'second-inspiral-path', 'first-coinc-path', 'second-coinc-path']
    self.iniNameMaps = map(None, self.iniNames, self.types)
    self.oNames = ['bank.cache', 'trigbank.cache', 'first_inspiral.cache', \
         'second_inspiral.cache', 'first_thinca.cache', 'second_thinca.cache']
    self.nameMaps = map(None, self.oNames, self.types)
    self.ifoTypes = ['H1','H2','L1','H1H2','H1L1','H2L1','H1H2L1']

    if options.generate_fu_cache:
      print >> sys.stderr, "\nOption --generate-fu-cache is specified, it overwrites the hipe cache files  which already exist"
      self.getCacheAll(cp)
      self.writeCacheAll('fu_hipe.cache')
      print >> sys.stderr, "\nHIPE CACHE FILE WRITTEN TO: fu_hipe.cache"
    try:
      os.chdir("cache")
      os.chdir(currentPath)
    except:
      if len(string.strip(cp.get('hipe-cache','hipe-cache-path'))) > 0:
        os.symlink(string.strip(cp.get('hipe-cache','hipe-cache-path')), 'cache')
      else: pass


  def ifoDict(self):
    return {'H1':[],'H2':[],'L1':[],'H1H2':[], \
                    'H1L1':[],'H2L1':[],'H1H2L1':[]}

  def getCacheType(self, iniName, type, cp=None):
    self[type] = []
    p = re.compile(type)
    f = re.compile("FOLLOWUP")
    m = re.compile("-")
    x = re.compile(".xml")
    g = re.compile(".gz")
    try:
      dir = os.listdir(string.strip(cp.get('hipe-cache',iniName)))
      cache_path = os.path.abspath(string.strip(cp.get('hipe-cache',iniName)))
    except:
      dir = self.dir
      cache_path = os.path.abspath(self.options.cache_path)
    for fname in dir:
      if f.search(fname): continue
      if p.search(fname):
        fnamebis =  g.split(fname)[0]
        ifo = m.split(fnamebis)[0]
        tag = m.split(fnamebis)[1]
        start = m.split(fnamebis)[-2]
        dur = x.split(m.split(fnamebis)[-1])
        try:
          # scirun = string.strip(cp.get('hipe-cache','science-run'))
          # tmpentry = ifo+" "+scirun+" "+start+" "+dur[0]+" "+"file://localhost" +cache_path+"/"+fname
          tmpentry = ifo+" "+tag+" "+start+" "+dur[0]+" "+"file://localhost" +cache_path+"/"+fname
          entry = lal.CacheEntry(tmpentry)
          self[type].append(entry)
        except: pass

  def getCacheAll(self,cp=None):
    for iniName, type in self.iniNameMaps:
      self.getCacheType(iniName,type,cp)


  def writeCacheType(self,oName,type,cName):
    #cName = open(oName,'w')
    for fname in self[type]:
      cName.write(str(fname)+"\n")
    #cName.close()

  def writeCacheAll(self,outputfile):
    cName = open(outputfile,'w')
    for oName, type in self.nameMaps:
      self.writeCacheType(str(oName),type,cName)
    cName.close()


  def getListFromCache(self,cache):
    listInstance = []
    for ifo in cache:
      listInstance = listInstance + cache[ifo].pfnlist()
    return listInstance

  def filesMatchingGPSinCache(self, cacheString, time=None, cacheType=None, ifo_tag=None, ifo_in_coinc=None):
    cacheSubSet = self.ifoDict()
    try: 
      cacheList = Cache.fromfile(open(cacheString))
      
      cacheListTest = 0
    except: 
      print >> sys.stderr, "could not open the file " + cacheString
      cacheListTest = 1
    if not cacheListTest:
      for ifo in self.ifoTypes:
        try:
          if time:
            if time[ifo]:
              time_ifo = time[ifo]
            else:
              #if ifo found in the analysed times, but not in coinc...
              if ifo_tag and len(ifo)==2 and re.search(ifo,ifo_tag):
                time_ifo = time[ifo_in_coinc[0]]
              else:
                continue
            seg1 = segments.segment(time_ifo,time_ifo+1)
            seg2 = segments.segment(time_ifo-1,time_ifo)
          else:
            seg1 = None
            seg2 = None          
          listInstance = cacheList.sieve(ifo,cacheType,None,True)
          listInstance = listInstance.sieve(None,None,seg1)
          listInstance = listInstance.sieve(None,None,seg2)
          cacheSubSet[ifo] = listInstance
        except:
          continue
    return(cacheSubSet)


  def getProcessParamsFromCache(self, subCache, tag, time, getInsp = False, ifo_in_coinc=None):

    process = self.ifoDict()
    inspTable = self.ifoDict()
    for ifo in subCache:
      for f in subCache[ifo]:
        # As soon as process[ifo] is filled we can stop iterating in this loop.
        if process[ifo]:
          break

        path = f.path()
        extension = path.split('.')[len(path.split('.'))-1]
        if extension == 'gz': gz = True
        else: gz = False
        doc = utils.load_filename(path,False,gz)
        proc = table.get_table(doc, lsctables.ProcessParamsTable.tableName)
        if getInsp:
          insp = table.get_table(doc, lsctables.SnglInspiralTable.tableName)
        for row in proc:          
          if str(row.param).find("--ifo-tag") >= 0:
             ifoTag = row.value
          # The default is to assume that the file matches if there is no tag
          else: ifoTag = tag
          # this is a hack to handle the "-userTag" string in some xml files...
          if str(row.param).find("-userTag") >= 0:
             row.param = "--user-tag"
          # end of the hack...

        # We check that the ifo names appearing in the tag (input argument) are 
        # found in the ifoTag of the process param table. But we don't require 
        # an exact match. We allow a "tag" being double times, and an "ifoTag"
        #  being triple times. This is to handle the case of candidates which 
        # were identified when one ifo was vetoed while still in science mode.
        # Notice that this method is a hack. However implementing
        # this loose condition on the ifo-tag was required in order to be able
        # to run the followup robustly over the candidates of the 1st calendar
        # of S5 (after the IFO TIME SEPARATION).
        # THIS SHOULD BE IMPROVED IN FUTURE ANALYSES #

        ifoTagTest = True
        for j in range(0,len(tag)-1,2):
          ifostring = tag[j:j+2]
          if ifostring not in ifoTag:
            ifoTagTest = False
            break
        if ifoTagTest:

          if time[ifo]:
            time_ifo = time[ifo]
          else:
            # if time[ifo] is not defined (happens for a double coinc in triple
            #  times) then use the end_time in the first ifo of the coinc
            time_ifo = time[ifo_in_coinc[0]]
          search = table.get_table(doc, lsctables.SearchSummaryTable.tableName)
          for row in search:
            out_start_time = float(row.out_start_time)
            out_start_time_ns = float(row.out_start_time_ns)/1000000000
            out_end_time = float(row.out_end_time)
            out_end_time_ns = float(row.out_end_time_ns)/1000000000
            if ( (time_ifo >= (out_start_time+out_start_time_ns)) and (time_ifo <= (out_end_time+out_end_time_ns)) ):
              process[ifo] = proc
              if getInsp: inspTable[ifo] = insp
              break          

    if getInsp:
      return process, inspTable
    else:
      return process

  def processFollowupCache(self, cp, opts, trig, type="INSPIRAL_"):
    if cp.has_option('hipe-cache','hipe-intermediate-cache'):
      intermediateCache = string.strip(cp.get('hipe-cache','hipe-intermediate-cache'))
    else:
      intermediateCache = ''
    if opts.generate_fu_cache or len(intermediateCache) == 0:
      cacheFile = 'fu_hipe.cache'
    else:
      cacheFile = intermediateCache

    if type == "INSPIRAL_":
      # get the list of interferometers from the ifoTag
      listIfoTime = []
      for j in range(0,len(trig.ifoTag)-1,2):
        ifo = trig.ifoTag[j:j+2]
        listIfoTime.append(ifo)

      # get an ifo tag with wild cards (example: *H1*H2*)
      wildIfoTag = "*" + "*".join(listIfoTime) + "*"
      # now we define a wildIfoTag type, for example '*H1*H2*'
      # This type will be used to sieve (in filesMatchingGPSinCache())
      wildType = type + wildIfoTag

      try:
        inspiral_process_params = self.getProcessParamsFromCache( \
                       self.filesMatchingGPSinCache(cacheFile,\
                       trig.gpsTime, wildType, trig.ifoTag, trig.ifolist_in_coinc), \
                       trig.ifoTag, trig.gpsTime, False, trig.ifolist_in_coinc)
      except:
        print "couldn't get inspiral process params for " + str(trig.eventID)
        inspiral_process_params = []
      return inspiral_process_params
    else:
      try:
        inspiral_process_params, sngl = self.getProcessParamsFromCache( \
                       self.filesMatchingGPSinCache(cacheFile,\
                       trig.gpsTime, type), \
                       trig.ifoTag, trig.gpsTime,True)
      except:
        print "couldn't get "+type+" process params for " + str(trig.eventID)
        inspiral_process_params = []
        sngl = []
      return inspiral_process_params, sngl

  def readTriggerFiles(self,cp):

    # Since we are continuing get useful stuff from the ini file.
    if cp.has_option('triggers','hipe-output-cache'):
      triggerCacheString = string.strip(cp.get('triggers','hipe-output-cache'))
    else:
      triggerCacheString = ''
    if cp.has_option('triggers','triggers-tag'):
      triggerTag = string.strip(cp.get('triggers','triggers-tag'))
    else:
      triggerTag = ''

    if len(triggerCacheString) == 0 or len(triggerTag) == 0:
      xml_glob = string.strip(cp.get('triggers','xml-glob'))
    else:
      triggerCache = self.filesMatchingGPSinCache(triggerCacheString,None,triggerTag)
      triggerList = self.getListFromCache(triggerCache)
      xml_glob = triggerList[0]
      print xml_glob

    numtrigs = string.strip(cp.get('triggers','num-trigs'))
    statistic =  string.strip(cp.get('triggers','statistic'))
    bla =  string.strip(cp.get('triggers','bitten-l-a'))
    blb =  string.strip(cp.get('triggers','bitten-l-b'))
    found, coincs, search = readFiles(xml_glob,getstatistic(statistic,bla,blb))
    return numtrigs, found, coincs, search

##############################################################################
# functions to parse qscan cache files (used in analyseQscan)
##############################################################################

def getPathFromCache(fileName,type,ifo=None,time=None):
  qscanList = []
  cacheList = listFromFile(fileName)
  if len(cacheList) == 0:
    return qscanList
  for line in cacheList:
    test_line = True
    if not re.search(type,line.split(' ')[1]):
      test_line = False
    if ifo:
      if not ifo == line.split(' ')[0]:
        test_line = False
      else: pass
    if time:
      if not time == line.split(' ')[2]:
        test_line = False
      else: pass
    if test_line:
      path_output = line.split(' ')[-1]
      time_output = line.split(' ')[2]
      type_output = line.split(' ')[1]
      ifo_output = line.split(' ')[0]
      qscanList.append([path_output,time_output,type_output,ifo_output])
    else: continue
  return qscanList

##############################################################################
# class for qscan intermediate table
##############################################################################

class QscanIntermediateTable(table.Table):
  tableName = "qscan:intermediate:table"
  validcolumns = {
    "type": "lstring",
    "param": "lstring",
    "ifo": "lstring",
    "qscan_time": "lstring",
    "qscan_dir": "lstring",
  }

##############################################################################
# Function to get qscan background. 
##############################################################################

def getQscanBackgroundTimes(cp, opts, ifo, dq_url_pattern, segFile):
    times = []
  
    if cp.has_option('background-qscan-times',ifo+'range'):
      rangeString = string.strip(cp.get('background-qscan-times',ifo+'range'))
    else:
      rangeString = ''
    if cp.has_option('background-qscan-times',ifo+'segment-list'):
      segmentListFile = string.strip(cp.get('background-qscan-times',ifo+'segment-list'))
    else:
      segmentListFile = ''
    if cp.has_option('background-qscan-times',ifo+'time-list'):
      timeListFile = string.strip(cp.get('background-qscan-times',ifo+'time-list'))
    else:
      timeListFile = ''

    if len(rangeString) == 0 and len(segmentListFile) == 0 and len(timeListFile) == 0:
      print "No qscan background specified for " + ifo
    else:

      # Generate the list of science segments (excluding cat 1 vetoes) if a time range is provided in the ini file
      if not len(rangeString) == 0:
        epochStart = rangeString.split(',')[0]
        epochEnd = rangeString.split(',')[1]
        opts.gps_start_time = int(epochStart)
        opts.gps_end_time = int(epochEnd)
        opts.use_available_data = False
        opts.run_data_quality = False
        # overwrite the ini file if the field "analyze" in section [segments] exist...
        cp.set("segments", "analyze", "Science")

        inspiralutils.findSegmentsToAnalyze(cp,opts,ifo,dq_url_pattern,segFile)
        segmentListFile = segFile[ifo]

      # Use the segment list if provided, and generate a list of random times
      if not len(segmentListFile) == 0:
        segmentList = pipeline.ScienceData()
        segmentMin = cp.getint('background-qscan-times','segment-min-len')
        segmentList.read(segmentListFile,segmentMin)
        segmentListLength = segmentList.__len__()
        segmentListStart = segmentList.__getitem__(1).start()
        segmentListEnd = segmentList.__getitem__(segmentListLength - 1).end()

        seed = cp.getint('background-qscan-times','random-seed')
        statistics = cp.getint('background-qscan-times','background-statistics')

        random.seed(seed)
        counter = 0
        times = []

        while counter < statistics:
          gps = float(segmentListStart) + float(segmentListEnd - segmentListStart)*random.random()
          testList = copy.deepcopy(segmentList)
          secondList = [pipeline.ScienceSegment(tuple([0,int(gps),int(gps)+1,1]))]        
          if testList.intersection(secondList):
            times.append(gps)
            counter = counter + 1
          else:
            continue
        # Save the list of times in a file (for possibly re-using it later)
        timeList = floatToStringList(times)
        fileName = "timeList-" + ifo + "-" + repr(seed) + "-" + repr(statistics) + "-" + repr(segmentListStart) + "-" + repr(segmentListEnd) + ".txt"
        saveRandomTimes(timeList,fileName)

      # Use the time-list file if provided
      if len(segmentListFile) == 0 and not len(timeListFile) == 0:
        timeList = listFromFile(timeListFile)
        if not timeList:
          print >> sys.stderr, "the list of times in file " + timeListFile + " could not be found"
          sys.exit(1)
        times = stringToFloatList(timeList)

    return times

##############################################################################
# function to read/write a list of strings in a file
##############################################################################
def listFromFile(fileName):
  listInstance = []
  try:
    file = open(fileName,"r")
  except:
    print >> sys.stderr, "could not open file " + fileName
    return listInstance
  list_in_file = file.readlines()
  file.close()
  if not len(list_in_file):
    print >> sys.stderr, "No lines found in file " + fileName
    print >> sys.stderr, "Is the first line blank ?"
    return listInstance
  for line in list_in_file:
    if not line[0] == '#':
      listInstance.append(string.strip(line))
    else: pass
  return listInstance

def saveRandomTimes(timeList,fileName):
  file = open(fileName,"w")
  for time in timeList:
    file.write(time + '\n')
  file.close()

def stringToFloatList(listin):
  listout = []
  for line in listin:
    listout.append(float(line))
  return listout

def floatToStringList(listin):
  listout = []
  for line in listin:
    listout.append(repr(line))
  return listout

##############################################################################
# function to read in a list of files and extract the simInspiral tables
# and sngl_inspiral tables
##############################################################################
def readFiles(fileGlob,statistic=None):
  """
  read in the Sngl and SimInspiralTables from a list of files
  if Sngls are found, construct coincs, add injections (if any)
  also return Sims (if any)
  @param fileGlob: glob of input files
  @param statistic: statistic to use in creating coincs
  """
  #if fileGlob is empty return empty structures
  if not fileGlob:
    if opts.verbose:
      print "Warning: No glob specified, returning empty structures..."
    return None, CoincInspiralUtils.coincInspiralTable(), None

  # if there aren't any files globbed exit
  fList = glob.glob(fileGlob)
  if not fList:
    print >>sys.stderr, "The glob for " + fileGlob + " returned no files"
    sys.exit(1)

  sims = None
  coincs = None
  search = None
  for thisFile in fList:
    extension = thisFile.split('.')[len(thisFile.split('.'))-1]
    if extension == 'gz': gz = True
    else: gz = False
    doc = utils.load_filename(thisFile,False,gz)
    # extract the sim inspiral table
    try:
      simInspiralTable = \
          table.get_table(doc, lsctables.SimInspiralTable.tableName)
      if sims: sims.extend(simInspiralTable)
      else: sims = simInspiralTable
    except: simInspiralTable = None

    # extract the sngl inspiral table, construct coincs
    try: snglInspiralTable = \
      table.get_table(doc, lsctables.SnglInspiralTable.tableName)
    except: 
      snglInspiralTable = None
      searchSumTable = None
    if snglInspiralTable:
      coincInspiralTable = \
        CoincInspiralUtils.coincInspiralTable(snglInspiralTable,statistic)
      if simInspiralTable:
        coincInspiralTable.add_sim_inspirals(simInspiralTable)
      # extract the search_summary table only if a sngl inspiral table is found
      searchSumTable = table.get_table(doc,lsctables.SearchSummaryTable.tableName)
      if coincs: 
        coincs.extend(coincInspiralTable)
        search.extend(searchSumTable)         
      else: 
        coincs = coincInspiralTable
        search = searchSumTable    
  return sims,coincs,search


#############################################################################
# function to set up directories
#############################################################################

def createdir(newDir):
  try:
    os.mkdir(newDir)
  except:
    pass


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

  try:
    os.chdir("logs")
    os.chdir('..')
  except: os.mkdir("logs")

  try:
    os.chdir("datafind_cache")
    os.chdir('..')
  except: os.mkdir("datafind_cache")

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
def getstatistic(stat, bla, blb):

  if stat == "effective_snrsq":
    newstat = "effective_snr"
  else:
    newstat = stat

  statistic=CoincInspiralUtils.coincStatistic( newstat, bla, blb )
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
    self.ifos = '' # added to construct a new field in "followups", see function getfollowuptrigs
    self.ifolist_in_coinc = None
    self.ifoTag = ''
    self.coincs = Coincs
    self.missed = Missed
    self.eventID = None
    self.stat = None
    self.page = None
    self.summarydir = None
    self.summarypage = None

  def add_coincs(self,Coincs):
    setattr(self,"coincs",Coincs)
    self.eventID = Coincs.event_id
    self.statValue = Coincs.stat
    self.ifos, self.ifolist_in_coinc = Coincs.get_ifos()
    if self.is_trigs():
      self.summarydir = "followuptrigs"
    if self.is_found():
      self.summarydir = "followupfound"
  
  def add_missed(self,Missed):
    setattr(self,"missed",Missed)
    self.summarydir = "followupmissed"
  
  def is_trigs(self):
    if isinstance(self.coincs,CoincInspiralUtils.coincInspiralTable.row):
      return 1
  
  def is_found(self):
    sim = None
    try: 
      sim = isinstance(self.coincs.sim,lsctables.SimInspiral)
    except: return 0
    if sim:
      return 1
  def is_missed(self):
    if isinstance(self.missed,lsctables.SimInspiralTable):
      return 1
  def add_page(self,page):
    self.page = page

#############################################################################
# Function to generate a trigbank xml file
#############################################################################
def generateXMLfile(ckey,ifo,outputPath=None,table_type='pre-bank-veto'):

  if outputPath:
    try:
      os.mkdir(outputPath) 
    except: pass

  xmldoc = ligolw.Document()
  xmldoc.appendChild(ligolw.LIGO_LW())
  trig = getattr(ckey,ifo)
  # BEFORE WE MAKE A NEW TABLE FIGURE OUT WHAT COLUMNS ARE VALID !!!
  valid_columns = trig.__slots__
  columns = []
  for col in valid_columns:
    try: 
      getattr(trig,col)
      columns.append(col)
    except: pass

  process_params_table = lsctables.New(lsctables.ProcessParamsTable)
  xmldoc.childNodes[-1].appendChild(process_params_table) 

  sngl_inspiral_table = lsctables.New(lsctables.SnglInspiralTable,columns)
  xmldoc.childNodes[-1].appendChild(sngl_inspiral_table)
  sngl_inspiral_table.append(trig)

  fileName = ifo + '-TRIGBANK_FOLLOWUP_' + str(ckey.event_id) + ".xml.gz"
  if outputPath:
    fileName = outputPath + '/' + fileName
  utils.write_filename(xmldoc, fileName, verbose = True, gz = True)   

def generateBankVetoBank(fuTrig, ifo,sngl,subBankSize,outputPath=None):
  
  trig =  getattr(fuTrig.coincs,ifo)
  mass = trig.mass1 + trig.mass2
  if outputPath:
    try:
      os.mkdir(outputPath)  
    except: pass
  xmldoc = ligolw.Document()
  xmldoc.appendChild(ligolw.LIGO_LW())
  process_params_table = lsctables.New(lsctables.ProcessParamsTable)
  xmldoc.childNodes[-1].appendChild(process_params_table)
  valid_columns = lsctables.SnglInspiralTable.validcolumns
  columns = []
  notcolumns = []
  for col in valid_columns:
    try:
      getattr(sngl[0],col)
      columns.append(col)
    except:
      notcolumns.append(col) 
  sngl_inspiral_table = lsctables.New(lsctables.SnglInspiralTable)
  sngl.sort(lambda a, b: cmp(a.mtotal, b.mtotal))
  index = sngl.getColumnByName('mtotal').index(mass)
  fromEnd = len(sngl)-index-int(subBankSize/2)
  if fromEnd < 0:
    sngl_sub = sngl[index+fromEnd-int(subBankSize/2):-1]
  else: 
    sngl_sub = sngl[index-int(subBankSize/2):index+int(subBankSize/2)]
    
  for row in sngl_sub:
    for col in notcolumns:
      setattr(row,col,0)
    sngl_inspiral_table.append(row)
  xmldoc.childNodes[-1].appendChild(sngl_inspiral_table)
 
  fileName = ifo + '-BANKVETO_FOLLOWUP_' + str(fuTrig.eventID) + ".xml.gz"
  if outputPath:
    fileName = outputPath + '/' + fileName
  utils.write_filename(xmldoc, fileName, verbose = True, gz = True)
  return fileName

#############################################################################
# Function to return the follow up list of coinc triggers
#############################################################################
def getfollowuptrigs(numtrigs,page=None,coincs=None,missed=None,search=None,trigbank_test=None):

  followups = []
  if coincs:
    sim = None
    try:
      sim = isinstance(coincs[0].sim,lsctables.SimInspiral)
    except: pass
    if sim: 
      coincs.sort(False) # This does an ascending sort instead for found inj
    else: coincs.sort()
    numTrigs = 0
    for ckey in coincs:
      numTrigs += 1
      if numTrigs > eval(numtrigs):
        break
      fuList = followUpList()
      fuList.add_coincs(ckey)
      if page:
        fuList.add_page(page)
      ifo_list = ['H1','H2','L1','G1','V1','T1']
      for ifo in ifo_list:
        try:
          getattr(ckey,ifo)
          fuList.gpsTime[ifo] = (float(getattr(ckey,ifo).end_time_ns)/1000000000)+float(getattr(ckey,ifo).end_time)
        except: fuList.gpsTime[ifo] = None
        if fuList.gpsTime[ifo] and trigbank_test:
          generateXMLfile(ckey,ifo,'trigTemplateBank')

      # now, find the ifoTag associated with the triggers, 
      # using the search summary tables...
      if fuList.ifolist_in_coinc:
        firstIfo = fuList.ifolist_in_coinc[0]
        triggerTime = fuList.gpsTime[firstIfo]
        if search:  
          for chunk in search:
            out_start_time = float(chunk.out_start_time)
            out_start_time_ns = float(chunk.out_start_time_ns)/1000000000
            out_end_time = float(chunk.out_end_time)
            out_end_time_ns = float(chunk.out_end_time_ns)/1000000000
            if ( (triggerTime >= (out_start_time+out_start_time_ns)) and (triggerTime <= (out_end_time+out_end_time_ns)) ):
              fuList.ifoTag = chunk.ifos
              break 
      followups.append(fuList)
  # the missed stuff doesnt work yet!!!
  if missed:
    followups
  return followups

#############################################################################
# Class to hold summary HTML information for all of the functions
#############################################################################
class summaryHTMLTable:

  def __init__(self,trig):
    if trig.is_trigs() and not trig.is_found():
      self.summarypath =  "followuptrigs/"
    if trig.is_trigs() and trig.is_found():
      self.summarypath = "followupfound/"
    self.eventID = trig.eventID
    self.statValue = trig.statValue
    self.H1time = trig.gpsTime["H1"]
    self.H2time = trig.gpsTime["H2"]
    self.L1time = trig.gpsTime["L1"]
    self.G1time = trig.gpsTime["G1"]
    self.V1time = trig.gpsTime["V1"]
    self.T1time = trig.gpsTime["T1"]
    self.containers = []

class HTMLcontainer:

  def __init__(self,trig,name,alt_web=None):
    # The missed injections dont work yet!!!
    self.name = name.rsplit(".")[-1]
    self.detailpath = ""
    if trig.is_trigs() and not trig.is_found():   
      os.chdir("followuptrigs")
      try: 
        os.chdir(self.name)
        os.chdir("../../")
      except: 
        os.mkdir(self.name)
        os.chdir("../")
      self.detailpath = trig.page + "/followuptrigs/" + self.name + "/"
      self.localdetailpath = "followuptrigs/" + self.name + "/"

    if trig.is_trigs() and trig.is_found():
      os.chdir("followupfound")
      try:
        os.chdir(self.name)
        os.chdir("../../")
      except:
        os.mkdir(self.name)
        os.chdir("../")
      self.detailpath = trig.page + "/followupfound/" + self.name + "/"
      self.localdetailpath = "followupfound/" + self.name + "/"

    self.image = self.detailpath + str(trig.statValue) + "_" + str(trig.eventID) + "_" + self.name + ".png"
    self.localimage = self.localdetailpath + str(trig.statValue) + "_" + str(trig.eventID) + "_" + self.name + ".png" 
    self.text = "click here"
    if alt_web: name = alt_web
    else: name = self.name
    self.link = self.detailpath + str(trig.statValue) + "_" + str(trig.eventID) + "_" + name + ".html"
    self.locallink = self.localdetailpath + str(trig.statValue) + "_" + str(trig.eventID) + "_" + name + ".html"



##############################################################################
# Function to write the HTML tables to pages
##############################################################################
def writeIULHeader(file):
  file.write('<%method title>Follow Up Report</%method><%method headline>Follow Up Report</%method><%method cvsid>$Id$</%method>\n')

def beginSummaryTable(file, table):
  file.write("<h3>Trigger [" + str(table.eventID) + 
             "] with combined statistic = " + str(table.statValue) + "</h3>")
  file.write('\n<br><table width=800 border=1>')
  
def endSummaryTable(file, table):
  file.write('\n</table>')

def writeModule(file, container):
  file.write('\n<tr><td width=400><font color="red" size=5>'+
             container.name + '</font>')
  file.write('\n<br><a href= "' + container.link + '">' + container.text 
             + "</a></td>")
  file.write('\n<td><a href= "' + container.image + '">\n<img src="' + 
              container.image + '" width=400 alt="No Image"></a></td></tr>')

def writeHTMLTables(summaryHTMLlist):
  for table in summaryHTMLlist:
    tableFile = open(table.summarypath + str(table.statValue) + "_" + 
                      str(table.eventID) + "_summary.html" ,'w')
    beginSummaryTable(tableFile,table)
    for container in table.containers:
      writeModule(tableFile, container)
    endSummaryTable(tableFile,table)
    writeIULHeader(tableFile)
    tableFile.close()

class HTMLTable:
  def __init__(self):
    self.columns = []
    self.headers = []

  def add_column(self,rows,header):
    self.columns.append(rows)
    self.headers.append(header)

  def write(self,file):
    file.write('\n<br><table><tr>')
    cnt = 0
    for i in self.columns:
      file.write('<td><b>'+ self.headers[cnt]+'</b><table>')
      cnt +=1
      for j in i:
        file.write('<tr><td>'+str(j)+'</td></tr>\n')
      file.write('</table></td>')
    file.write('\n</tr></table><br>\n')

##############################################################################
# Function to publish the web tree
##############################################################################
def publishOnHydra(page):
  indexFile = open("index.html","w")
  patt = re.compile('.*summary.html')
  
  # First do the found injections
  files = os.chdir('followupfound')
  files = os.listdir('.')
  table = HTMLTable()
  fList = []
  stat = []
  ID = []
  if files:
    for f in files:
      temp = patt.match(f)
      if temp:
        fList.append((float(temp.group().rsplit('_')[0]), \
                  temp.group().rsplit('_')[1]))
    indexFile.write('<h3>Found Injection follow ups</h3>\n')
    sortedF = sorted(fList,key=operator.itemgetter(0),reverse=False)
    for i in sortedF:
      stat.append(str(i[0]))
      ID.append('<a href="' +page+ '/followupfound/'+str(i[0])+'_'+str(i[1])+
                '_summary.html">'+i[1]+'</a>')
    table.add_column(stat,'Stat Value')
    table.add_column(ID,'ID')
    table.write(indexFile)  
  os.chdir('..')

# then do trigs
  files = os.chdir('followuptrigs')
  files = os.listdir('.')
  table = HTMLTable()
  fList = []
  stat = []
  ID = []
  if files:
    for f in files:
      temp = patt.match(f)
      if temp:
        fList.append((float(temp.group().rsplit('_')[0]), \
                  temp.group().rsplit('_')[1]))
    indexFile.write('<h3>Found Trigger follow ups</h3>\n')
    sortedF = sorted(fList,key=operator.itemgetter(0),reverse=True)
    for i in sortedF:
      stat.append(str(i[0]))
      ID.append('<a href="' +page+ '/followuptrigs/'+str(i[0])+'_'+str(i[1])+
                '_summary.html">'+i[1]+'</a>')
    table.add_column(stat,'Stat Value')
    table.add_column(ID,'ID')
    table.write(indexFile)
  os.chdir('..')



#  # Then do the triggers 
#  files = os.chdir('followuptrigs')
#  files = os.listdir('.')
#  table = HTMLTable()
#  fList = []
#  stat = []
#  ID = []
#  if files:
#    for f in files:
#      temp = patt.match(f)
#      if temp:
#        fList.append((float(temp.group().rsplit('_')[0]), \
#                  temp.group().rsplit('_')[1]))
#    indexFile.write('<h3>Trigger follow ups</h3>\n')
#    sortedF = sorted(fList,key=operator.itemgetter(0),reverse=True)
#    for i in sortedF:
#      stat.append(str(i[0]))
#      ID.append('<a href="' +page+ '/followuptrigs/'+str(i[0])+'_'+str(i[1])+
#                '_summary.html">'+i[1]+'</a>')    
#    table.add_column(stat,'Stat Value')
#    table.add_column(ID,'ID')
#    table.write(indexFile)
#  
#  os.chdir('..')

  writeIULHeader(indexFile)
  indexFile.close()
  #This needs to be done so that it doesn't keep writing over this 
  #directory... ;)
  os.system('scp -r index.html followuptrigs followupfound followupmissed ' +
            'hydra.phys.uwm.edu:/home/htdocs/uwmlsc/root/'+
            page + '.')

############################################################
# Class for checking follow up GPS times to Nelson's veto lists.
############################################################
class nVeto:
  """
  Author: Cristina Valeria Torres
  Provides a method for checking S5y1 and S5y2 follow up candidates
  against locally cached static veto lists produced by Nelson Christensen
  """
  def __init__(self):
    """
    Initialize this class.
    """
    self.database=list([["L1","/archive/home/romain/Projects/LowMassCBC/20051104-20061114/triggers/upperlimits_v99_july08/nelsonVetoes/L1/VetoLists",list()],
                        ["H1","/archive/home/romain/Projects/LowMassCBC/20051104-20061114/triggers/upperlimits_v99_july08/nelsonVetoes/H1/VetoLists",list()],
                        ["H2","/archive/home/romain/Projects/LowMassCBC/20051104-20061114/triggers/upperlimits_v99_july08/nelsonVetoes/H2/VetoLists",list()]])
      #Form [[filename,list(gpsIntervals)],[filename2,list(gpsIntervals)]...
    self.tolWin=float(0.01)
    self.filesLoaded=bool(False)
  #End self.__init__()

  def __loadFiles__(self):
    """
    This method is invoked automatically if required.
    It creates the dynamically structures to hold the veto lists.
    """
    #Check that the paths specified exist

    for ifoName,ifoPath,vetoData in self.database:
        if not os.path.exists(ifoPath):
          raise IOError("%s Path Invalid: %s"%(ifoName,ifoPath))
    for index in range(self.database.__len__()):
      self.__processPath__(index)
    self.filesLoaded=True
  #End self.__loadFiles__()

  def __processPath__(self,index=""):
    """
    This method is not to be called explicity.  It is called by
    self.__loadFiles__(). The role of this method is to fill up 
    the instance with the data required to check the veto intervals.
    """
    path=self.database[index][1]
    listOfFiles=os.listdir(path)
    listOfFiles=[os.path.normpath(path+"/"+file) for file in listOfFiles]
    for filename in listOfFiles:
      fileData=list()
      fp=open(str(filename),'r')
      txtData=fp.readlines()
      fp.close()
      for row in txtData:
        a,b,c=map(float,row.split())
        fileData.append([a,b,c])
      self.database[index][2].append([filename,fileData])
  #End self.__processPath__()

  def __checkIntervals__(self,vList,gpsTime):
    """
    This method is not to be called explicity. It is called by
    self.findInterval(IFO,gpsTime).  This method returns a 
    text string giving the veto list names and intervals that
    the gpsTime intersects with.
    """
    vetoMatch=list()
    for vetoName,vetoList in vList:
      for startT,stopT,KWSig in vetoList:
        if ((startT-self.tolWin)<=gpsTime<=(stopT+self.tolWin)):
          vetoMatch.append([vetoName,startT,stopT])
    tmpList=list()
    for a,b,c in vetoMatch:
      tmpList.append("%s %s %s\n"%(str(a),str(b),str(c)))
    outputString=str().join(tmpList)
    if vetoMatch.__len__()==0:
      outputString="%s %s %s\n"%("NONE","0","0")
    return outputString
  #End self.__checkIntervals__()

  def setIfoPath(self,ifoName="",ifoPath=""):
    """
    If you want to use something other than the hardwired path
    defaults to the veto files invoke this method with the IFO
    for which you need to adjust the path to those veto files.
    ie To adjust L1 paths do
    self.setIfopath("L1","/path/to/new/veto/files")
    """
    myIndex=-1
    for i in range(self.database.__len__()):
      if (self.database[i][0].lower()==ifoName.lower()):
        myIndex=i
    if myIndex==-1:
      raise IOError("%s Invalid IFO Name Path Unchanged."%(ifoName))
    self.database[myIndex][1]=ifoPath
  #End self.setIfoPath()

  def findInterval(self,ifoName,gpsTime):
    """
    Given an IFO name and GPS time it searches the veto lists to
    determine which if any of the veto intervals interset this
    gps time.
    """
    if not self.filesLoaded:
      self.__loadFiles__()
    myIndex=-1
    for i in range(self.database.__len__()):
      if (self.database[i][0].lower()==ifoName.lower()):
        myIndex=i
    if myIndex==-1:
      raise IOError("%s Invalid IFO Name"%(ifoName))
    vList=self.database[myIndex][2]
    outputString=self.__checkIntervals__(vList,gpsTime)
    return outputString
  #End self.findInterval()
#
#End nVeto() Class Definition
############################################################
