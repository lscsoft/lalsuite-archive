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

from pylab import *
from glue import segments
from glue import segmentsUtils
from glue.ligolw import ligolw
from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw import utils
from pylal import CoincInspiralUtils
from glue import pipeline
from glue.lal import *
from glue import lal

########## CLASS TO WRITE LAL CACHE FROM HIPE OUTPUT #########################
class getCache(UserDict):
  """
  An instance of a lal cache
  """
  def __init__(self, options):
    UserDict.__init__(self)
    self.dir = os.listdir(os.getcwd())
    self.options = options
    self.types = ['TMPLTBANK', 'TRIGBANK', 'INSPIRAL-', \
                 'INSPIRAL_', 'THINCA-', 'THINCA_']
    self.iniNames = ['tmpltbank-path', 'trigbank-path', 'first-inspiral-path', \
         'second-inspiral-path', 'first-coinc-path', 'second-coinc-path']
    self.iniNameMaps = map(None, self.iniNames, self.types)
    self.oNames = ['bank.cache', 'trigbank.cache', 'first_inspiral.cache', \
         'second_inspiral.cache', 'first_thinca.cache', 'second_thinca.cache']
    self.nameMaps = map(None, self.oNames, self.types)
    self.ifoTypes = ['H1','H2','L1','H1H2','H1L1','H2L1','H1H2L1']

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
    list = []
    for ifo in cache:
      list = list + cache[ifo].pfnlist()
    return list

  def filesMatchingGPSinCache(self, cacheString, time=None, cacheType=None):
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
            seg1 = segments.segment(time[ifo],time[ifo]+1)
            seg2 = segments.segment(time[ifo]-1,time[ifo])
          else:
            seg1 = None
            seg2 = None
          list = cacheList.sieve(ifo,cacheType,seg1)
          list = list.sieve(None,None,seg2)
          cacheSubSet[ifo] = list
        except:
          continue
    return(cacheSubSet)


  def getProcessParamsFromCache(self, subCache, tag, time):
    process = self.ifoDict()
    for ifo in subCache:
      for f in subCache[ifo]:
        path = f.path()
        extension = path.split('.')[len(path.split('.'))-1]
        if extension == 'gz': gz = True
        else: gz = False
        doc = utils.load_filename(path,False,gz)
        proc = table.get_table(doc, lsctables.ProcessParamsTable.tableName)
        for row in proc:          
          if str(row.param).find("--ifo-tag") >= 0:
             ifoTag = row.value
          # this is a hack to handle the "-userTag" string in some xml files...
          if str(row.param).find("-userTag") >= 0:
             row.param = "--user-tag"
          # end of the hack...     
        
        if ifoTag == tag:
          search = table.get_table(doc, lsctables.SearchSummaryTable.tableName)
          for row in search:
            out_start_time = float(row.out_start_time)
            out_start_time_ns = float(row.out_start_time_ns)/1000000000
            out_end_time = float(row.out_end_time)
            out_end_time_ns = float(row.out_end_time_ns)/1000000000
            if ( (time[ifo] >= (out_start_time+out_start_time_ns)) and (time[ifo] <= (out_end_time+out_end_time_ns)) ):
              process[ifo] = proc
              break

    return process
    
##############################################################################
# function to read/write a list of strings in a file
##############################################################################
def listFromFile(fileName):
  list = []
  try:
    file = open(fileName,"r")
  except:
    print >> sys.stderr, "could not open file " + fileName
    return list
  list_in_file = file.readlines()
  if not len(list_in_file):
    print >> sys.stderr, "No lines found in file " + fileName
    print >> sys.stderr, "Is the first line blank ?"
    return list
  for line in list_in_file:
    list.append(string.strip(line))
  return list

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

#  try:
#    os.chdir("hoft_qscan_cache")
#    os.chdir('..')
#  except: os.mkdir("hoft_qscan_cache")

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
def generateXMLfile(ckey,ifo,gpsTime):

  xmldoc = ligolw.Document()
  xmldoc.appendChild(ligolw.LIGO_LW())

  process_params_table = lsctables.New(lsctables.ProcessParamsTable)
  xmldoc.childNodes[-1].appendChild(process_params_table) 

  sngl_inspiral_table = lsctables.New(lsctables.SnglInspiralTable)
  xmldoc.childNodes[-1].appendChild(sngl_inspiral_table)
  trig = getattr(ckey,ifo)
  sngl_inspiral_table.append(trig)

  fileName = ifo + '-TRIGBANK_FOLLOWUP_' + gpsTime + ".xml"
  utils.write_filename(xmldoc, fileName, verbose = True, gz = False)   

#############################################################################
# Function to return the follow up list of coinc triggers
#############################################################################
def getfollowuptrigs(numtrigs,page,coincs=None,missed=None,search=None,trigbank_test=None):

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
      fuList.add_page(page)
      try:
        getattr(ckey,'H1')
        fuList.gpsTime["H1"] = (float(getattr(ckey,'H1').end_time_ns)/1000000000)+float(getattr(ckey,'H1').end_time)
        if trigbank_test:
          try: generateXMLfile(ckey,'H1',str(fuList.gpsTime["H1"])) 
          except: print "the trigBank xml file could not be generated for H1 " + str(fuList.eventID)
      except: fuList.gpsTime["H1"] = None
      try:
        getattr(ckey,'H2')
        fuList.gpsTime["H2"] = (float(getattr(ckey,'H2').end_time_ns)/1000000000)+float(getattr(ckey,'H2').end_time)
        if trigbank_test:
          try: generateXMLfile(ckey,'H2',str(fuList.gpsTime["H2"]))
          except: print "the trigBank xml file could not be generated for H2 " + str(fuList.eventID)
      except: fuList.gpsTime["H2"] = None
      try:
        getattr(ckey,'L1')
        fuList.gpsTime["L1"] = (float(getattr(ckey,'L1').end_time_ns)/1000000000)+float(getattr(ckey,'L1').end_time)
        if trigbank_test:
          try: generateXMLfile(ckey,'L1',str(fuList.gpsTime["L1"]))
          except: print "the trigBank xml file could not be generated for L1 " + str(fuList.eventID)
      except: fuList.gpsTime["L1"] = None
      try:
        getattr(ckey,'G1')
        fuList.gpsTime["G1"] = (float(getattr(ckey,'G1').end_time_ns)/1000000000)+float(getattr(ckey,'G1').end_time)
        if trigbank_test:
          try: generateXMLfile(ckey,'G1',str(fuList.gpsTime["G1"]))
          except: print "the trigBank xml file could not be generated for G1 " + str(fuList.eventID)
      except: fuList.gpsTime["G1"] = None
      try:
        getattr(ckey,'V1')
        fuList.gpsTime["V1"] = (float(getattr(ckey,'V1').end_time_ns)/1000000000)+float(getattr(ckey,'V1').end_time)
        if trigbank_test:
          try: generateXMLfile(ckey,'V1',str(fuList.gpsTime["V1"]))
          except: print "the trigBank xml file could not be generated for V1 " + str(fuList.eventID)
      except: fuList.gpsTime["V1"] = None
      try:
        getattr(ckey,'T1')
        fuList.gpsTime["T1"] = (float(getattr(ckey,'T1').end_time_ns)/1000000000)+float(getattr(ckey,'T1').end_time)
        if trigbank_test:
          try: generateXMLfile(ckey,'T1',str(fuList.gpsTime["T1"]))
          except: print "the trigBank xml file could not be generated for T1 " + str(fuList.eventID)
      except: fuList.gpsTime["T1"] = None

      # now, find the ifoTag associated with the triggers, using the search summary tables...
      if fuList.ifolist_in_coinc:
        firstIfo = fuList.ifolist_in_coinc[0]
        triggerTime = fuList.gpsTime[firstIfo]
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

