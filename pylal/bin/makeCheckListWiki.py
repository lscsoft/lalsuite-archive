#!/usr/bin/python
#
# Copyright (C) 2009 Cristina Valeria Torres
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#
"""
This checklist script is responsible for creating a MoinMoin Wiki
page.  This checklist page is meant to replace the original HTML based
page created by the script makeCheckList.py.  The output of this
script is a text file that can be cut and paste into a MoinMoin Wiki
editor, resulting in a Wikified checklist under MoinMoin version
control.
"""
__author__  = "Cristina Valeria Torres <cristina.torres@ligo.org>"
__date__    = '$Date$'
__version__ = '$Revision$'
__prog__    = 'makeCheckListWiki.py'


####################################################################
# Read in the standard python modules and the LSC python modules
# needed to create the MoinMoin file
####################################################################
import copy
import optparse
import ConfigParser
import os
import random
import socket
import sys
import time
import urllib
import fnmatch
import shutil
import re
from pylal import stfu_pipe
from glue import cbcwebpage
from pylal import git_version
sys.path.append('@PYTHONLIBDIR@')




####################################################################
# Custom methods for building wiki checklists
####################################################################

def scanTreeFnMatch(parentPath='.',levels=int(100),filemask='*'):
  """
  This recurses the subdirectories of parentPath 0,1,...,N levels.  It
  returns a list of files who match FILEMASK.
  """
  matchingFiles=list()
  myTreeLocal=os.walk(os.path.normpath(parentPath))
  parentLevel=parentPath.count(os.sep)
  #Iterate down the levels
  for root,dir,files in myTreeLocal:
    if root.count(os.sep) <= parentLevel+levels:
      for myFile in [os.path.join(root,x) for x in files]:
        if fnmatch.fnmatch(myFile,filemask):
          matchingFiles.append(myFile)
  return matchingFiles

class findFileType(object):
  """
  Initialized with a file structure and coinc data it can return a
  list of files from that structure.
  """
  def __init__(self,fStructure=None,myCoinc=None):
    if fStructure==None or myCoinc==None:
      print "Given None Types FS:%s Coinc:%s"%(type(fStructure),type(myCoinc))
      return None
    else:
      self.fsys=fStructure
      self.coinc=myCoinc

  def __readZranks__(self,myFilename=None):
    """
    Takes a file and returns a structure (list) of information from
    the zvalue files.
    """
    if myFilename == None:
      sys.stderr.write("No filename passed to method __readZranks__()!\n")
      return None
    rawData=[x.rstrip("\n") for x in file(myFilename)]
    #Sort this data by percentiles: Channel Z PercentileSignificance
    try:
      #ignoreOctothorpe
      rawData2=list()
      for index,row in enumerate(rawData):
        if not row.startswith("#") and row != "":
          rawData2.append(row.strip())
      #Split into strings
      listData=[str(a).split() for a in rawData2]
      #Adjust properties to Str,float,float
      listData=[[str(a),float(b),float("%2.3f"%float(c))] for a,b,c in listData]
      tmpListData=[[x[2],x] for x in listData]
      tmpListData.sort(reverse=True)
      listData=[x[1] for x in tmpListData]
    except:
      sys.stderr.write("Error parsing Z value file :%s\n"%(myFilename))
      return []
    return listData
  
  def __readCache__(self,cacheListing=list()):
    """
    Simple mehtod to read in a cache or list of cache files and return
    a list of files to an empty list if nothing found
    """
    #Open the cache entry and search for those entrys
    fileListing=list()
    for entry in cacheListing:
      #Cache files listed themselves comment out following line
      fileListing.append(entry)
      fileListing.extend([x.rstrip("\n") for x in file(entry)])
      #PATCH START to add in the z distribution files
      for fname in fileListing:
        if ".html" in fname:
          zFile=fname.replace(".html",".txt")
          fileListing.append(zFile)
      #PATCH END
    finalList=list()
    for thisFile in fileListing:
      #Search filesystem for file full path
      finalList.extend(fnmatch.filter(self.fsys,"*%s"%thisFile))
      #Look for potential matching thumbnails
      if thisFile.endswith(".png"):
        finalList.extend(fnmatch.filter(self.fsys,"*%s"%thisFile.replace(".png","?thumb?png")))
    if len(finalList) < 1:
      return list()
    else:
      return finalList

    
  def get_hoft_frame(self):
    """
    """
    tmpList=list()
    for sngl in self.coinc.sngls:
      #Determine file type
      frametype,channelName=stfu_pipe.figure_out_type(sngl.time,sngl.ifo,'hoft')
      myMaskIndex="*%s/*/%s/index.html"%(frametype,sngl.time)
      myMaskPNG="*%s/*/%s/*.png"%(frametype,sngl.time)
      tmpList.extend(fnmatch.filter(self.fsys,myMaskIndex))
      tmpList.extend(fnmatch.filter(self.fsys,myMaskPNG))
    return tmpList
    
  def get_RDS_C03_L2(self):
    """
    """
    tmpList=list()
    for sngl in self.coinc.sngls:
      myMaskIndex="*%s_RDS_C03_L2/*/%s/index.html"%(sngl.ifo,sngl.time)
      myMaskPNG="*%s_RDS_C03_L2/*/%s/*.png"%(sngl.ifo,sngl.time)
      tmpList.extend(fnmatch.filter(self.fsys,myMaskIndex))
      tmpList.extend(fnmatch.filter(self.fsys,myMaskPNG))
    return tmpList
    
  def get_RDS_R_L1(self):
    """
    """
    tmpList=list()
    for sngl in self.coinc.sngls:
      myMaskIndex="*/%s_RDS_R_L1/*/%s/index.html"%(sngl.ifo,sngl.time)
      myMaskPNG="*/%s_RDS_R_L1/*/%s/*.png"%(sngl.ifo,sngl.time)
      tmpList.extend(fnmatch.filter(self.fsys,myMaskIndex))
      tmpList.extend(fnmatch.filter(self.fsys,myMaskPNG))
    return tmpList

  def get_RDS_R_L1_SEIS(self):
    """
    """
    tmpList=list()
    for sngl in self.coinc.sngls:
      myMaskIndex="*/%s_RDS_R_L1_SEIS*/%s/*.html"%(sngl.ifo,sngl.time)
      myMaskPNG="*/%s_RDS_R_L1_SEIS*/%s/*.png"%(sngl.ifo,sngl.time)
      tmpList.extend(fnmatch.filter(self.fsys,myMaskIndex))
      tmpList.extend(fnmatch.filter(self.fsys,myMaskPNG))
    return tmpList
      
  def get_findVetos(self):
    tmpList=list()
    #H1,H2,L1-findFlags_H1,H2,L1_831695156.714.wiki
    #instrument,ifos
    ifoString=""
    for i in range(0,len(self.coinc.ifos)/2):ifoString=ifoString+"%s,"%self.coinc.ifos[2*i:2*i+2]
    ifoString=ifoString.rstrip(",")
    insString=""
    for i in range(0,len(self.coinc.instruments)/2):insString=insString+"%s,"%self.coinc.instruments[2*i:2*i+2]
    insString=insString.rstrip(",")
    myMask="*%s*%s-findVetos_%s_%s.wiki"%\
            (self.coinc.type,insString,ifoString,self.coinc.time)
    tmpList.extend(fnmatch.filter(self.fsys,myMask))
    return tmpList
    
  def get_effDRatio(self):
    tmpList=list()
    #H1,H2,L1-findFlags_H1,H2,L1_831695156.714.wiki
    #instrument,ifos
    ifoString=""
    for i in range(0,len(self.coinc.ifos)/2):ifoString=ifoString+"%s,"%self.coinc.ifos[2*i:2*i+2]
    ifoString=ifoString.rstrip(",")
    insString=""
    for i in range(0,len(self.coinc.instruments)/2):insString=insString+"%s,"%self.coinc.instruments[2*i:2*i+2]
    insString=insString.rstrip(",")
    myMask="*%s*%s-effDRatio_%s_%s.wiki"%\
            (self.coinc.type,insString,ifoString,self.coinc.time)
    tmpList.extend(fnmatch.filter(self.fsys,myMask))
    return tmpList
  
  def get_findFlags(self):
    """
    """
    tmpList=list()
    #H1,H2,L1-findFlags_H1,H2,L1_831695156.714.wiki
    #instrument,ifos
    ifoString=""
    for i in range(0,len(self.coinc.ifos)/2):ifoString=ifoString+"%s,"%self.coinc.ifos[2*i:2*i+2]
    ifoString=ifoString.rstrip(",")
    insString=""
    for i in range(0,len(self.coinc.instruments)/2):insString=insString+"%s,"%self.coinc.instruments[2*i:2*i+2]
    insString=insString.rstrip(",")
    myMask="*%s*%s-findFlags_%s_%s.wiki"%\
            (self.coinc.type,insString,ifoString,self.coinc.time)
    tmpList.extend(fnmatch.filter(self.fsys,myMask))
    return tmpList
    
  def get_plotsnrchisq(self):
    """
    """
    tmpList=list()
    insString=""
    for i in range(0,len(self.coinc.instruments)/2):insString=insString+"%s,"%self.coinc.instruments[2*i:2*i+2]
    insString=insString.rstrip(",")
    for sngl in self.coinc.sngls:
      myMask="*%s*/%s-plotsnrchisq_pipe_%s_FOLLOWUP_PLOTSNRCHISQ_%s*.cache"%\
              (self.coinc.type,\
               insString,\
               sngl.ifo,\
               sngl.time)
      tmpList.extend(fnmatch.filter(self.fsys,myMask))
    #Open the cache entry and search for those entrys
    cacheListing=list()
    for entry in tmpList:
      cacheListing.append(entry)
      cacheListing.extend([x.rstrip("\n") for x in file(entry).readlines()])
    finalList=list()
    for thisFile in cacheListing:
      finalList.extend(fnmatch.filter(self.fsys,"*%s"%thisFile))
      #Scan for both thumb types for all PNGs
      if thisFile.endswith(".png"):
        finalList.extend(fnmatch.filter(self.fsys,"*%s"%thisFile.replace(".png","_thumb.png")))
        finalList.extend(fnmatch.filter(self.fsys,"*%s"%thisFile.replace(".png",".thumb.png")))
    if len(finalList) < 1:
      return list()
    else:
      return finalList
    
  def get_plotchiatimeseries(self):
    """
    This is for the coherence based tests.
    """
    tmpList=list()
    myMask="*/%s-plotchiatimeseries_%s_PLOT_CHIA_%s*.cache"%\
            (self.coinc.instruments,\
             self.coinc.ifos,\
             self.coinc.time)
    tmpList.extend(fnmatch.filter(self.fsys,myMask))
    #Open the cache entry and search for those entrys
    cacheListing=list()
    for entry in tmpList:
      cacheListing.append(entry)
      cacheListing.extend([x.rstrip("\n") for x in file(entry).readlines()])
    finalList=list()
    for thisFile in cacheListing:
      finalList.extend(fnmatch.filter(self.fsys,"*%s"%thisFile))
      if thisFile.endswith(".png"):
        finalList.extend(fnmatch.filter(self.fsys,"*%s"%thisFile.replace(".png","_thumb.png")))
        finalList.extend(fnmatch.filter(self.fsys,"*%s"%thisFile.replace(".png",".thumb.png")))
    if len(finalList) < 1:
      return(list())
    else:
      return finalList

  def get_analyzeQscan_SEIS(self):
    """
    This seeks out the html and png files associated with SEIS result
    of an analyzeQscan job.
    """
    cacheList=list()
    cacheFiles=list()
    for sngl in self.coinc.sngls:
      intS=nanS=0
      intS,nanS=str(float(sngl.time)).split(".")
      timeString="%s_%s"%(intS,nanS)
      myCacheMask="*/%s-analyseQscan_%s_%s*_seis_rds*.cache"%(sngl.ifo,sngl.ifo,timeString)
      #Read the cache file or files
      cacheList.extend(fnmatch.filter(self.fsys,myCacheMask))
    cacheFiles=self.__readCache__(cacheList)
    return cacheFiles
      
  def get_analyzeQscan_RDS(self):
    """
    """
    #analyseQscan.py_FG_RDS_full_data/H1-analyseQscan_H1_931176926_116_rds-unspecified-gpstime.cache
    cacheList=list()
    cacheFiles=list()
    for sngl in self.coinc.sngls:
      intS=nanS=0
      intS,nanS=str(float(sngl.time)).split(".")
      timeString="%s_%s"%(intS,nanS)
      myCacheMask="*/%s-analyseQscan_%s_%s*_rds*.cache"%(sngl.ifo,sngl.ifo,timeString)
      #Ignore the files with seis_rds in them
      for x in fnmatch.filter(self.fsys,myCacheMask):
        if not x.__contains__('seis_rds'):
          cacheList.append(x)
    #Read the cache file or files
    cacheFiles=self.__readCache__(cacheList)
    return cacheFiles                   

  def get_analyzeQscan_HT(self):
    """
    """
    #analyseQscan.py_FG_HT_full_data/H1-analyseQscan_H1_931176926_116_ht-unspecified-gpstime.cache
    cacheList=list()
    cacheFiles=list()
    for sngl in self.coinc.sngls:
      intS=nanS=0
      intS,nanS=str(float(sngl.time)).split(".")
      timeString="%s_%s"%(intS,nanS)
      myCacheMask="*/%s-analyseQscan_%s_%s*_ht*.cache"%(sngl.ifo,sngl.ifo,timeString)
      cacheList.extend(fnmatch.filter(self.fsys,myCacheMask))
    #Read the cache file or files
    cacheFiles=self.__readCache__(cacheList)
    return cacheFiles         

  
  def get_all(self):
    """
    """
    globalList=list()
    globalList.extend(self.get_plotsnrchisq())
    globalList.extend(self.get_plotchiatimeseries())
    globalList.extend(self.get_hoft_frame())
    globalList.extend(self.get_analyzeQscan_HT())
    globalList.extend(self.get_RDS_R_L1())
    globalList.extend(self.get_analyzeQscan_RDS())
    globalList.extend(self.get_RDS_R_L1_SEIS())
    globalList.extend(self.get_analyzeQscan_SEIS())
    globalList.extend(self.get_findVetos())
    globalList.extend(self.get_effDRatio())
    globalList.extend(self.get_findFlags())    
    return globalList
                 
def matchFiles(fileList=None,jobString=None,instruments=None,ifos=None,time=None):
  """
  Given a list of file paths is tests this list to select the files
  related to the options given to the method.  A list of matches is
  returned.
  """
  matchList=list()
  for thisFile in fileList:
    tPath,tFile=os.path.split(thisFile.replace(",",""))
    if (tFile.__contains__(jobString) and \
        tFile.__contains__(instruments) and \
        tFile.__contains__(ifos) and \
        tFile.__contains__(time)):
      matchList.append(thisFile)
  return matchList


####################################################################
# Custom wiki class to make writing MoinMoin text simpler
####################################################################
class wiki(object):
  def __init__(self,filename="default.wiki"):
    self.file = open(filename,"w")
    self.content=list()

  def tableOfContents(self,levels="1"):
    """
    Take an integer to determine maxdepth of TOC
    """
    self.content.append("\n <<TableOfContents(%i)>> \n"%(levels))

  def image_link(self,path,webserver):
    thumb = "thumb_" + path
    command = 'convert ' + path + ' -resize 300x300 -antialias ' + thumb
    print command
    popen = subprocess.Popen(command.split())
    popen.communicate()
    status = popen.returncode
    #s = '[[ImageLink('+webserver+'/'+thumb+','+webserver+'/'+path+',width=300][,alt=none])]]'
    s = '[['+webserver+'/'+path+'|{{'+webserver+'/'+thumb+'||width=300}}]]'
    self.file.write(s)

  def linkedRemoteImage(self,image=None,link=None):
    """
    Take a link to an image and display it, then use the image as a
    link to some other web location
    """
    if link==None:
      sys.stdout.write("No link specified.")
      return ""
    if image==None:
      image=str("")
    wikiString="[[%s|{{%s}}|target=\"_blank\"]]"
    return wikiString%(link,image)

  def image_table(self,image_list, webserver):
    if not image_list: return
    for j, i in enumerate(image_list):
      if not (j) % 3: self.file.write("\n\n||")
      self.image_link(i, webserver)
      self.file.write("||")
    self.file.write("\n\n")
  
  def image_glob(self, pat):
    image_list = []
    for image in glob.glob(pat):
      if 'thumb' in image: continue
      else: image_list.append(image)
    image_list.sort()
    return image_list

  def insertHR(self):
    self.content.append("\n----------\n")

  def title(self,text=""):
    s = "= %s =\n"%(text.strip())
    self.content.append(s)

  def section(self,text=""):
    s = "== %s ==\n"%(text.strip())
    self.content.append(s)

  def subsection(self,text=""):
    s = "=== %s ===\n"%(text.strip())
    self.content.append(s)

  def subsubsection(self,text=""):
    s = "==== %s ====\n"%(text.strip())
    self.content.append(s)

  def putText(self,text=""):
    s = "\n%s\n"%(text.strip())
    self.content.append(s)

  def makeExternalLink(self,url="",label=""):
    """
    Returns the MoinMoin specific wikification string
    for an external link, like an HTML page.  This string
    can be inserted into the page by a putText call.
    """
    s = " [[%s|%s]] "%(url.strip(),label.strip())
    return s

  class wikiTable(object):
    """
    Internal class to manipulate a wiki table
    Addressable as
    X[row][col]=Information
    """
    def __init__(self,rows=1,cols=1):
      if rows < 1:
        sys.stderr.write("Error setting number of rows in table. \
Requested %s rows, assuming you meant at least 1.\n"%rows)
        self.rows = 1
      else:
        self.rows=rows
      if cols < 1:
        sys.stderr.write("Error setting number of columns in table. \
Requested %s columns, assuming you meant at least 1.\n"%col)
        self.cols = 1
      else:
        self.cols=cols
      self.tStyle=None
      self.tHeadline=None
      self.data=list()
      #Build headline format string
      self.tHeadlineFormat=""
      for colNum in range(0,self.cols):
        self.tHeadlineFormat+="||"
      self.tHeadlineFormat+="""%s||\n"""
      
      #Create tuple object with number of rows
      for rc in range(0,rows):
        self.data.append(self.__rowbuilder__(self.cols))

    def __rowbuilder__(self,cols):
      return [str(" ") for x in range(0,cols)]

    def setTableHeadline(self,titleString=""):
      """
      Allows you to insert a title row to a table object and place
      text or other table commands. Invoking it with None as the arg
      deletes a previous set table headline.
      """
      #To achieve colspan tabulate extra "||" to place in start of tHeadline
      if titleString==None:
        self.tHeadline=None
      elif self.tStyle != None:
        self.tHeadline=titleString
        
    def setTableStyle(self,fstring=""):
      """
      Allows you to specify table style see MoinMoin help
      Setting arg to NONE removes the current style specified if any.
      """
      if fstring=="NONE":
        self.tStyle=None
      else:
        self.tStyle=fstring

    def __formatAndTitleString__(self,location=None):
      """
      Looks at values for wiki table editing and returns the
      appropriate string for a title row or options for the table to
      cell[0][0]. Valid arguments are the keywords ['TITLE','CELL']
      """
      titleString=""
      cellString=""
      #Title string with specified formatting
      if self.tHeadline and self.tStyle:
        tmpString="""<tablestyle="%s" rowbgcolor="#FFFFE0" > %s"""%(self.tStyle,self.tHeadline)
        titleString = self.tHeadlineFormat%tmpString
      #Title string with no formatting
      if self.tHeadline and not self.tStyle:
        titleString = self.tHeadlineFormat%self.tHeadline
      #No title string with table formatting
      if not self.tHeadline and self.tStyle:
        cellString =  """<tablestyle="%s">"""%self.tStyle
      #Default to return an empty string if neither option set
      if not self.tHeadline and not self.tStyle:
        titleString = ""
        cellString = ""
      if location == None:
        raise Exception, "Argument to self.__formatAndTitleString__() \
expected. Received <None>!\n"
      if location.upper() == "CELL":
        return cellString
      elif location.upper() == "TITLE":
        return titleString
      else:
        return ""
      
  def insertTable(self,obj):
    """
    Pass in a wikiTable object then create the relevant
    wiki markup and place that in to the self.content
    list for writing to the file
    """
    tableContent=""
    #If a title row specified
    tableContent+=obj.__formatAndTitleString__("TITLE")
    oldCell="%s"%obj.data[0][0]
    obj.data[0][0]="%s%s"%(obj.__formatAndTitleString__("CELL"),str(oldCell))
    if type(obj) != type(self.wikiTable()):
      raise Exception,"Expecting nested type instance of WikiTable"
    else:
      for row in range(0,obj.rows):
        for col in range(0,obj.cols):
          try:
            if obj.data[row][col].rstrip().lstrip().__contains__("style"):
              tableContent=tableContent+"||%s "%(obj.data[row][col].rstrip().lstrip())
            else:
              tableContent=tableContent+"||%s "%(obj.data[row][col].rstrip().lstrip())
          except:
            sys.stderr.write("Error creating wiki markup for table. \
R:%i/%i,C:%i/%i,Cells:%i\n"%(row,obj.rows,col,obj.cols,len(obj.data)))
            raise
        tableContent="%s ||\n"%(tableContent)                           
    tableContent="%s\n"%(tableContent)
    self.content.append(tableContent)                      
    obj.data[0][0]=oldCell

  def insertQscanTable(self,images=None,thumbs=None,indexes=None):
    """
    Inserts a table constructured of thumbnails linked to larger
    Qscan plots.  It accounts for the ifo present in the coinc via
    qCoinc.  The method expects a lists of URLs
    Channel naming extraction in the method likely broken see code in
    method insertAnalyzeQscanTable as a fix guide.
    """
    if images.keys() != indexes.keys():
      sys.write.stderr("Error: insertQscanTable ifo keys malformed.\n")
    #Generate Image Labels
    channelNames=list()
    for ifo in images.keys():
      channelNames.extend([os.path.basename(x).split("_",1)[1].rsplit("_",3)[0].split(":",1)[1] \
                       for x in images[ifo]])
    uniqChannelNames=list()
    lastName=None
    channelNames.sort()
    while channelNames:
      myName=channelNames.pop()
      if lastName != myName:
        lastName=myName
        uniqChannelNames.append(myName)
    #Create table object
    rowCount=len(uniqChannelNames)+1
    colCount=len(images.keys())+1
    myTable=self.wikiTable(rowCount,colCount)
    myTable.setTableStyle("text-align:center")
    #Make title row
    myTable.data[0][0]=""
    for i,label in enumerate(indexes.keys()):
      if len(indexes[label]) < 1:
        myTable.data[0][i+1]=" %s "%label
      else:
        myIndexURL="%s"%indexes[label][0]
        myTable.data[0][i+1]="%s"%self.makeExternalLink(myIndexURL,label)

    #Fill in table with thumbnails and links
    for i,channel in enumerate(uniqChannelNames):
      myTable.data[i+1][0]=" %s "%(channel)
      for j,key in enumerate(images.keys()):
        try:
          imageIndex=[x.__contains__(channel) \
                      for x in images[key]].index(True)
          imageURL=images[key][imageIndex]
          thumbIndex=[x.__contains__(channel) \
                      for x in thumbs[key]].index(True)
          thumbURL=thumbs[key][thumbIndex]
          myTable.data[i+1][j+1]=self.linkedRemoteImage(thumbURL,\
                                                        imageURL)
        except:
          myTable.data[i+1][j+1]="Unavailable"
    self.insertTable(myTable)
    
  def insertAnalyzeQscanTable(self,
                              images=None,
                              thumbs=None,
                              indexes=None,
                              imagesAQ=None,
                              thumbsAQ=None,
                              indexesAQ=None,
                              channelRanks=None):
    """
    Insert a multiple IFO table with 5 cols with the AQ underneath
    this depends on the numer of IFO keys in indexes dictionary.
    The option channelRanks is not required to change the plot order!
    Channel ranks is dict similar in shape to other args.
    Cells are shaded light grey if they are top N channels and that
    the trigger is greater in value that 0.5.  Assuming the
    channelRanks dict is not empty. 
    """
    #channelRanks={'ifo':[[chan,Zvalue,rank]...[chan,Zvalue,rank]],'ifo2':[[ ]]}
    #Review the keys for Qscans and analyzeQscans.
    if not images.keys()==thumbs.keys()==indexes.keys():
      sys.stderr.write("Error: Keys for Qscan tables creations inconsistent!\n")
    if not imagesAQ.keys()==thumbsAQ.keys()==indexesAQ.keys():
      sys.stderr.write("Error: Keys for Qscan tables creations inconsistent!\n")

    keyList=indexes.keys()
    if len(keyList) < indexesAQ.keys():
      keyList=indexesAQ.keys()
    for ifo in keyList:
      #Overall loop for each IFO in the dict structures

      #Generate Image Labels
      #H1-analyseQscan_H1_931245125_408_seis_rds_H0_PEM-MY_SEISZ_dt_dist-unspecified-gpstime.png
      #L1-analyseQscan_L1_931182185_246_rds_L0_PEM-LVEA_MAGX_z_dist-unspecified-gpstime.png
      #H1-analyseQscan_H1_931257951_208_ht_H1_LDAS-STRAIN_dt_dist-unspecified-gpstime.png
      #932797512.6862793_H0:PEM-BSC9_ACC1X_16.00_eventgram_autoscaled.png
      #932797512.6862793_H1:LSC-DARM_ERR_512.00_eventgram_raw.thumb.png
      #933259905.03857422_H0:PEM-MX_SEISY_512.00_eventgram_autoscaled.png
      #CHANNEL NAMES in ALL CAPS
      channelNames=list()
      #Extract channel names
      tmpCN=[os.path.basename(x) for x in images[ifo]]
      tmpCN.extend([os.path.basename(x) for x in imagesAQ[ifo]])
      startREG=re.compile('_[H,V,L][0,1,2][:,-,_]')
      stopREG=re.compile('_(?=[0-9,a-z])')
      channelNames=[re.split(stopREG,re.split(startREG,x).pop())[0].strip() \
                    for x in tmpCN]
      uniqChannelNames=list()
      lastName=None
      channelNames.sort()
      while channelNames:
        myName=channelNames.pop()
        if lastName != myName:
          lastName=myName
          uniqChannelNames.append(myName)
      #Check if uniqChannelNames list is empty
      if len(uniqChannelNames) < 1:
        sys.stderr.write("Warning: [%s] No channels available to plot in table!\n"%ifo)
        uniqChannelNames.append("No_Channels_To_Display")
      #Create table object reserve first cell for txt labels
      colCount=3
      fullRows,modRows=divmod(len(uniqChannelNames),colCount)
      if modRows > 0:
        rowCount=fullRows+1
      else:
        rowCount=fullRows
      myTable=self.wikiTable(rowCount,colCount)
      myTable.setTableStyle("text-align:center")
      #Insert HTML links and IFO Label
      contentString=""
      contentString=contentString+" %s  "%(ifo)
      #Add html links
      for newLink in indexes[ifo]:
        contentString=contentString+" %s "%self.makeExternalLink(newLink,"Qscan")
      for newLink in indexesAQ[ifo]:
        contentString=contentString+" %s  "%self.makeExternalLink(newLink,"analyzeQscan")
      #Legend for top N analyzeQscan images
      #ifoColors={'L1':'blue','H1':'orange','H2':'magenta','V1':'grey','DEFAULT':'pink'}
      #Shortlist the channels we will highlight: Sort then cut
      tmpRanks=[[x[2],x] for x in channelRanks[ifo]]
      tmpRanks.sort(reverse=True)
      ifoRanks=[x[1] for x in tmpRanks]
      topN=9
      shortList=ifoRanks[0:min(len(ifoRanks),topN)]
      myTable.setTableHeadline(contentString)
      #Start filling cells with Qscan and analyzeQscan scatter plot
      for cellNum,channel in enumerate(uniqChannelNames):
        #Grab plot info for this channel name
        myName=channel
        try:
          myOmegaIndex=[x.__contains__(myName) for x in images[ifo]].index(True)
        except ValueError:
          myOmegaIndex=None
        try:
          myOmegaIndexT=[x.__contains__(myName) for x in thumbs[ifo]].index(True)
        except ValueError:
          myOmegaIndexT=None
        try:
          myAQIndex=[x.__contains__(myName) for x in imagesAQ[ifo]].index(True)
        except ValueError:
          myAQIndex=None
        try:
          myAQIndexT=[x.__contains__(myName) for x in thumbsAQ[ifo]].index(True)
        except ValueError:
          myAQIndexT=None
        cellString=""
        #Setup shading
        htmlShadeColor='#E0E0FF'
        cutP=0.75
        if [a.__contains__(myName) \
            and c >= cutP for a,b,c in shortList].count(True):
          cellString=cellString+"<%s> "%htmlShadeColor
        #Use indices to get URLs
        if myName:
          cellString=cellString+" %s <<BR>> "%myName
        else:
          cellString=cellString+" Unknown_Channel <<BR>> "
        if myOmegaIndex!=None:
          cellString=cellString+" %s "%self.linkedRemoteImage(thumbs[ifo][myOmegaIndexT],
                                                             images[ifo][myOmegaIndex])
        else:
          cellString=cellString+" Unavailable_Qscan <<BR>> "
        if myAQIndex!=None:
          cellString=cellString+" %s "%self.linkedRemoteImage(thumbsAQ[ifo][myAQIndexT],
                                                              imagesAQ[ifo][myAQIndex])
        else:
          cellString=cellString+" Unavailable_analyzeQScan <<BR>> "
        #Add string to cell
        myRow,myCol=divmod(cellNum,colCount)
        myTable.data[myRow][myCol]=" %s "%cellString
      self.insertTable(myTable)

  def write(self):
    """
    Writes the contents of the wiki object to disk.
    It flushes the current wikification to disk, more
    wikification steps could be done.  If writing is
    not possible an error is raised.
    """
    try:
      self.file.writelines(self.content)
    except:
      raise
    self.content=list()

  def finish(self):
    self.file.close()  
####################################################################
# End Custom wiki class
####################################################################

####################################################################
# Custom method that represents a checklist
# the checklist is written to disk at the end of the method run
####################################################################

def prepareChecklist(wikiFilename=None,wikiCoinc=None,wikiTree=None,file2URL=None):
  """
  Method to prepare a checklist where data products are isolated in
  directory.
  """
  endOfS5=int(875232014)
  wikiFileFinder=findFileType(wikiTree,wikiCoinc)
  #
  # Check to see if wiki file with name already exists
  #
  maxCount=0
  while os.path.exists(wikiFilename) and maxCount < 10:
    sys.stdout.write("File %s already exists.\n"%\
                     os.path.split(wikiFilename)[1])
    wikiFilename=wikiFilename+".wiki"
    maxCount=maxCount+1
  #
  #Create the wikipage object etc
  #
  wikiPage=wiki(wikiFilename)
  #
  # Create top two trigger params tables
  #
  cTable=wikiPage.wikiTable(2,9)
  cTable.data=[
    ["Trigger Type",
     "Rank",
     "FAR",
     "SNR",
     "IFOS(Coinc)",
     "Instruments(Active)",
     "Coincidence Time (s)",
     "Total Mass (mSol)",
     "Chirp Mass (mSol)"
     ],
    ["%s"%(wikiCoinc.type),
     "%s"%(wikiCoinc.rank),
     "%s"%(wikiCoinc.far),
     "%s"%(wikiCoinc.snr),
     "%s"%(wikiCoinc.ifos),
     "%s"%(wikiCoinc.instruments),
     "%s"%(wikiCoinc.time),
     "%s"%(wikiCoinc.mass),
     "%s"%(wikiCoinc.mchirp)
     ]
    ]
  pTable=wikiPage.wikiTable(len(wikiCoinc.sngls_in_coinc())+1,7)
  pTable.data[0]=[
    "IFO",
    "GPS Time(s)",
    "SNR",
    "CHISQR",
    "Mass 1",
    "Mass 2",
    "Chirp Mass"
    ]
  for row,cSngl in enumerate(wikiCoinc.sngls_in_coinc()):
    pTable.data[row+1]=[
      "%s"%(cSngl.ifo),
      "%s"%(cSngl.time),
      "%s"%(cSngl.snr),
      "%s"%(cSngl.chisqr),
      "%s"%(cSngl.mass1),
      "%s"%(cSngl.mass2),
      "%s"%(cSngl.mchirp)
      ]
  #Write the tables into the Wiki object
  wikiPage.putText("Coincident Trigger Event Information: %s\n"\
                   %(stfu_pipe.gpsTimeToReadableDate(wikiCoinc.time)))
  wikiPage.insertTable(cTable)
  wikiPage.putText("Corresponding Coincident Single IFO Trigger Information\n")
  wikiPage.insertTable(pTable)

  #Generate a table of contents to appear after candidate params table
  wikiPage.tableOfContents(3)

  #Begin including each checklist item as section with subsections
  wikiPage.section("Follow-up Checklist")
  #Put each checklist item
  wikiPage.subsection("Checklist Summary")
  wikiPage.subsubsection("Does this candidate pass this checklist?")
  wikiPage.subsubsection("Answer")
  wikiPage.subsubsection("Relevant Information and Comments")
  wikiPage.insertHR()
  #
  #First real checklist item
  wikiPage.subsection("#0 False Alarm Probability")
  wikiPage.subsubsection("Question")
  wikiPage.putText("What is the false alarm rate associated with this candidate?")
  wikiPage.subsubsection("Answer")
  wikiPage.putText("Edit Here")
  wikiPage.subsubsection("Relevant Information")
  farTable=wikiPage.wikiTable(2,1)
  farTable.setTableStyle("background-color: yellow; text-align center;")
  farTable.data[0][0]="False Alarm Rate"
  farTable.data[1][0]="%s"%(wikiCoinc.far)
  wikiPage.insertTable(farTable)
  wikiPage.subsubsection("Investigator Comments")
  wikiPage.putText("Edit Here")
  wikiPage.insertHR()
  #
  #Additional Checklist Item
  #First real checklist item
  wikiPage.subsection("#1 Data Quality Flags")
  wikiPage.subsubsection("Question")
  wikiPage.putText("Can the data quality flags coincident with this candidate be safely disregarded?")
  wikiPage.subsubsection("Answer")
  wikiPage.putText("Edit Here")
  wikiPage.subsubsection("Relevant Information")
  wikiPath=os.path.split(wikiFilename)[0]
  dqFileList=wikiFileFinder.get_findFlags()
  if len(dqFileList) != 1:
    sys.stdout.write("Warning: DQ flags data product import problem.\n")
    print "Found %i files."%len(dqFileList)
    for mf in dqFileList: print mf
  for myFile in dqFileList:
    wikiPage.putText("%s\n"%(file(myFile).read()))
  wikiPage.subsubsection("Investigator Comments")
  wikiPage.putText("Edit Here")
  wikiPage.insertHR()
  #
  #Additional Checklist Item
  #First real checklist item
  wikiPage.subsection("#2 Veto Investigations")
  wikiPage.subsubsection("Question")
  wikiPage.putText("Does the candidate survive the veto investigations performed at its time?")
  wikiPage.subsubsection("Answer")
  wikiPage.putText("Edit Here")
  wikiPage.subsubsection("Relevant Information")
  vetoFileList=wikiFileFinder.get_findVetos()
  if len(vetoFileList) != 1:
    sys.stdout.write("Warning: Veto flags data product import problem.\n")
    for myFile in vetoFileList:print myFile
  for myFile in vetoFileList:
    wikiPage.putText("%s\n"%(file(myFile).read()))
  wikiPage.subsubsection("Investigator Comments")
  wikiPage.putText("Edit Here")
  wikiPage.insertHR()
  #
  #Additional Checklist Item
  #First real checklist item
  wikiPage.subsection("#3 IFO Status")
  wikiPage.subsubsection("Question")
  wikiPage.putText("Are the interferometers operating normally with a reasonable level of sensitivity around the time of the candidate?")
  wikiPage.subsubsection("Answer")
  wikiPage.putText("Edit Here")
  wikiPage.subsubsection("Relevant Information")
  #Add link to Daily Stats
  if wikiCoinc.time <= endOfS5:
    statsLink=wikiPage.makeExternalLink("http://blue.ligo-wa.caltech.edu/scirun/S5/DailyStatistics/",\
                                        "S5 Daily Stats Page")
  else:
    statsLink="This should be a link to S6 Daily Stats!\n"
  wikiPage.putText(statsLink)
  #Link figures of merit
  #Get link for all members of wikiCoinc
  wikiPage.putText("Figures of Merit\n")
  if wikiCoinc.time > endOfS5:
    fomLinks=dict()
    elems=0
    for wikiSngl in wikiCoinc.sngls:
      if not(wikiSngl.ifo.upper().rstrip().lstrip() == 'V1'):
        fomLinks[wikiSngl.ifo]=stfu_pipe.getFOMLinks(wikiCoinc.time,wikiSngl.ifo)
        elems=elems+len(fomLinks[wikiSngl.ifo])
      else:
        for myLabel,myLink,myThumb in stfu_pipe.getFOMLinks(wikiCoinc.time,wikiSngl.ifo):
          wikiPage.putText("%s\n"%(wikiPage.makeExternalLink(myLink,myLabel)))
    if elems%3 != 0:
      sys.stdout.write("Generation of FOM links seems incomplete!\n")
    cols=4
    rows=(elems/3)+1
    fTable=wikiPage.wikiTable(rows,cols)
    fTable.data[0]=["IFO,Shift","FOM1","FOM2","FOM3"]
    currentIndex=0
    for myIFOKey in fomLinks.keys():
      for label,link,thumb in fomLinks[myIFOKey]:
         myRow=currentIndex/int(3)+1
         myCol=currentIndex%int(3)+1
         fTable.data[myRow][0]=label
         thumbURL=thumb
         fTable.data[myRow][myCol]="%s"%(wikiPage.linkedRemoteImage(thumb,link))
         currentIndex=currentIndex+1
    wikiPage.insertTable(fTable)
  else:
    wikiPage.putText("Can not automatically fetch S5 FOM links.")  
  wikiPage.subsubsection("Investigator Comments")
  wikiPage.putText("Edit Here")
  wikiPage.insertHR()
  #
  #Additional Checklist Item
  #First real checklist item
  wikiPage.subsection("#4 Candidate Appearance")
  wikiPage.subsubsection("Question")
  wikiPage.putText("Do the Qscan figures show what we would expect for a gravitational-wave event?")
  wikiPage.subsubsection("Answer")
  wikiPage.putText("Edit Here")
  wikiPage.subsubsection("Relevant Information")
  imageDict=dict()
  indexDict=dict()
  thumbDict=dict()
  for sngl in wikiCoinc.sngls:
    frametype,channelName=stfu_pipe.figure_out_type(sngl.time,sngl.ifo,'hoft')
    indexDict[sngl.ifo]=fnmatch.filter(wikiFileFinder.get_hoft_frame(),\
                                       "*/%s/*/%s/*index.html"%(frametype,sngl.time))
    imageDict[sngl.ifo]=fnmatch.filter(wikiFileFinder.get_hoft_frame(),\
                                       "*%s*_%s_16.00_spectrogram_whitened.png"\
                                       %(sngl.time,channelName))
    thumbDict[sngl.ifo]=fnmatch.filter(wikiFileFinder.get_hoft_frame(),\
                                       "*%s*_%s_16.00_spectrogram_whitened?thumb.png"\
                                       %(sngl.time,channelName))
    #
    #Convert disk locals to URLs
    imageDict[sngl.ifo]=[file2URL.convert(x) for x in imageDict[sngl.ifo]]
    indexDict[sngl.ifo]=[file2URL.convert(x) for x in indexDict[sngl.ifo]]
    thumbDict[sngl.ifo]=[file2URL.convert(x) for x in thumbDict[sngl.ifo]]
    if len(indexDict[sngl.ifo]) < 1:
      wikiPage.putText("GW data channel scans for %s not available.\n"%sngl.ifo)
  enoughImage=[len(imageDict[key])>0 for key in imageDict.keys()].count(True) >= 1
  enoughIndex=[len(indexDict[key])>0 for key in indexDict.keys()].count(True) >= 1
  if enoughImage and enoughIndex:
    wikiPage.insertQscanTable(imageDict,\
                              thumbDict,\
                              indexDict)
  else:
    sys.stdout.write("Warning: Candidate appearance plot import problem.\n")
  wikiPage.subsubsection("Investigator Comments")
  wikiPage.putText("Edit Here")
  wikiPage.insertHR()
  #
  #
  #Additional Checklist Item
  wikiPage.subsection("#5 Seismic Plots")
  wikiPage.subsubsection("Question")
  wikiPage.putText("Is the seismic activity insignificant around the time of the candidate?")
  wikiPage.subsubsection("Answer")
  wikiPage.putText("Edit Here")
  wikiPage.subsubsection("Relevant Information")
  wikiPage.putText("Plots and pipeline data go here!")
  imageDict=dict()
  indexDict=dict()
  thumbDict=dict()
  imageDictAQ=dict()
  indexDictAQ=dict()
  thumbDictAQ=dict()
  zValueDictAQ=dict()
  for sngl in wikiCoinc.sngls_in_coinc():
    indexDict[sngl.ifo]=fnmatch.filter(wikiFileFinder.get_RDS_R_L1_SEIS(),\
                                       "*/%s_RDS_*/%s/*index.html"%(sngl.ifo,sngl.time))
    imageDict[sngl.ifo]=fnmatch.filter(wikiFileFinder.get_RDS_R_L1_SEIS(),\
                                       "*/%s_RDS_*/%s/*SEI*_512.00_spectrogram_whitened.png"%\
                                       (sngl.ifo,sngl.time))
    thumbDict[sngl.ifo]=fnmatch.filter(wikiFileFinder.get_RDS_R_L1_SEIS(),\
                                       "*/%s_RDS_*/%s/*SEI*_512.00_spectrogram_whitened?thumb.png"%\
                                       (sngl.ifo,sngl.time))
    #Search for analyzeQscan files
    intS=nanS=0
    intS,nanS=str(float(sngl.time)).split(".")
    timeString="%s_%s"%(intS,nanS)
    zValueFiles=fnmatch.filter(wikiFileFinder.get_analyzeQscan_SEIS(),\
                                         "*_%s_%s_*.txt"%(sngl.ifo,timeString))
    indexDictAQ[sngl.ifo]=fnmatch.filter(wikiFileFinder.get_analyzeQscan_SEIS(),\
                                         "*_%s_%s_*.html"%(sngl.ifo,timeString))
    thumbDictAQ[sngl.ifo]=fnmatch.filter(wikiFileFinder.get_analyzeQscan_SEIS(),\
                                         "*%s-*_%s_*_z_scat-unspecified-gpstime_thumb.png"\
                                         %(sngl.ifo,timeString))
    imageDictAQ[sngl.ifo]=[x.replace("_thumb.png",".png") for x in thumbDictAQ[sngl.ifo]]
    #Process zValue ranking file if found for IFO
    if len(zValueFiles):
      zValueDictAQ[sngl.ifo]=wikiFileFinder.__readZranks__(zValueFiles[0])
    else:
      zValueDictAQ[sngl.ifo]=list()
      sys.stdout.write("Z ranking file not found for %s. ...skipping...\n"%sngl.ifo)
    #Convert disk locals to URLs
    imageDict[sngl.ifo]=[file2URL.convert(x) for x in imageDict[sngl.ifo]]
    indexDict[sngl.ifo]=[file2URL.convert(x) for x in indexDict[sngl.ifo]]
    thumbDict[sngl.ifo]=[file2URL.convert(x) for x in thumbDict[sngl.ifo]]
    imageDictAQ[sngl.ifo]=[file2URL.convert(x) for x in imageDictAQ[sngl.ifo]]
    indexDictAQ[sngl.ifo]=[file2URL.convert(x) for x in indexDictAQ[sngl.ifo]]
    thumbDictAQ[sngl.ifo]=[file2URL.convert(x) for x in thumbDictAQ[sngl.ifo]]
    if len(indexDict[sngl.ifo]) < 1:
      wikiPage.putText("Seismic scans for %s not available.\n"%sngl.ifo)
  enoughImage=[len(imageDict[key])>0 for key in imageDict.keys()].count(True) >=1
  enoughIndex=[len(indexDict[key])>0 for key in indexDict.keys()].count(True) >=1
  if enoughImage and enoughIndex:
    wikiPage.insertAnalyzeQscanTable(imageDict,
                                     thumbDict,
                                     indexDict,
                                     imageDictAQ,
                                     thumbDictAQ,
                                     indexDictAQ,
                                     zValueDictAQ)
  else:
    sys.stdout.write("Warning: Seismic plots product import problem.\n")
  wikiPage.subsubsection("Investigator Comments")
  wikiPage.putText("Edit Here")
  wikiPage.insertHR()
  #
  #
  #Additional Checklist Item
  wikiPage.subsection("#6 Other environmental causes")
  wikiPage.subsubsection("Question")
  wikiPage.putText("Were the environmental disturbances (other than seismic) insignificant at the time of the candidate?")
  wikiPage.subsubsection("Answer")
  wikiPage.putText("Edit Here")
  wikiPage.subsubsection("Relevant Information")
  wikiPage.putText("Plots and pipeline data go here!")
  imageDict=dict()
  indexDict=dict()
  thumbDict=dict()
  imageDictAQ=dict()
  indexDictAQ=dict()
  thumbDictAQ=dict()
  zValueDictAQ=dict()
  #Select only PEM channels
  for sngl in wikiCoinc.sngls_in_coinc():
    imageDict[sngl.ifo]=list()
    indexDict[sngl.ifo]=list()
    thumbDict[sngl.ifo]=list()
    for myFile in fnmatch.filter(wikiFileFinder.get_RDS_R_L1(),\
                                 "*/%s_RDS_*/%s/*html"%(sngl.ifo,sngl.time)):
      indexDict[sngl.ifo].append(myFile)

    for myFile in fnmatch.filter(wikiFileFinder.get_RDS_R_L1(),\
                                 "*/%s_RDS_*/%s/*_16.00_spectrogram_whitened.png"%\
                                 (sngl.ifo,sngl.time)):
      if myFile.upper().__contains__("PEM"):
        imageDict[sngl.ifo].append(myFile)
        
    for myFile in fnmatch.filter(wikiFileFinder.get_RDS_R_L1(),\
                                 "*/%s_RDS_*/%s/*_16.00_spectrogram_whitened?thumb.png"%\
                                 (sngl.ifo,sngl.time)):
      if myFile.upper().__contains__("PEM"):
        thumbDict[sngl.ifo].append(myFile)
    #Select associated analyzeQscans
    imageDictAQ[sngl.ifo]=list()
    indexDictAQ[sngl.ifo]=list()
    thumbDictAQ[sngl.ifo]=list()
    intS=nanS=0
    intS,nanS=str(float(sngl.time)).split(".")
    timeString="%s_%s"%(intS,nanS)
    for myFile in fnmatch.filter(wikiFileFinder.get_analyzeQscan_RDS(),\
                                 "*%s-*_%s_*html"%(sngl.ifo,timeString)):
      indexDictAQ[sngl.ifo].append(myFile)
    zValuesFiles=fnmatch.filter(wikiFileFinder.get_analyzeQscan_RDS(),\
                                "*%s-*_%s_*.txt"%(sngl.ifo,timeString))
    if len(zValueFiles) > 0:
      zValueDictAQ[sngl.ifo]=wikiFileFinder.__readZranks__(zValueFiles[0])
    else:
      zValueDictAQ[sngl.ifo]=list()
    #H1-analyseQscan_H1_931176926_116_rds_H0_PEM-MY_SEISX_z_scat-unspecified-gpstime_thumb.png
    #H1-analyseQscan_H1_931176926_116_rds_H0_PEM-MY_SEISX_z_scat-unspecified-gpstime.png
    for myFile in fnmatch.filter(wikiFileFinder.get_analyzeQscan_RDS(),\
                                 "*%s-*_%s_*_z_scat-unspecified-gpstime.png"%\
                                 (sngl.ifo,timeString)):
      if myFile.upper().__contains__("PEM"):
        imageDictAQ[sngl.ifo].append(myFile)
        
    for myFile in fnmatch.filter(wikiFileFinder.get_analyzeQscan_RDS(),\
                                 "*%s-*_%s_*_z_scat-unspecified-gpstime?thumb.png"%\
                                 (sngl.ifo,timeString)):
      if myFile.upper().__contains__("PEM"):
        thumbDictAQ[sngl.ifo].append(myFile)
    #Convert disk locals to URLs
    imageDict[sngl.ifo]=[file2URL.convert(x) for x in imageDict[sngl.ifo]]
    indexDict[sngl.ifo]=[file2URL.convert(x) for x in indexDict[sngl.ifo]]
    thumbDict[sngl.ifo]=[file2URL.convert(x) for x in thumbDict[sngl.ifo]]
    imageDictAQ[sngl.ifo]=[file2URL.convert(x) for x in imageDictAQ[sngl.ifo]]
    indexDictAQ[sngl.ifo]=[file2URL.convert(x) for x in indexDictAQ[sngl.ifo]]
    thumbDictAQ[sngl.ifo]=[file2URL.convert(x) for x in thumbDictAQ[sngl.ifo]]
    if len(imageDict[sngl.ifo]) < 1:
      wikiPage.putText("PEM scans for %s not available.\n"%sngl.ifo)
  enoughImage=[len(imageDict[key])>0 for key in imageDict.keys()].count(True) >=1
  enoughIndex=[len(indexDict[key])>0 for key in indexDict.keys()].count(True) >=1
  if enoughImage and enoughIndex:
    wikiPage.insertAnalyzeQscanTable(imageDict,
                                     thumbDict,
                                     indexDict,
                                     imageDictAQ,
                                     thumbDictAQ,
                                     indexDictAQ,
                                     zValueDictAQ)
  else:
    sys.stdout.write("Warning: PEM plots import trouble.\n")
  wikiPage.subsubsection("Investigator Comments")
  wikiPage.putText("Edit Here")
  wikiPage.insertHR()
  #
  #
  #Additional Checklist Item
  wikiPage.subsection("#7 Auxiliary degree of freedom")
  wikiPage.subsubsection("Question")
  wikiPage.putText("Were the auxiliary channel transients coincident with the candidate insignificant?")
  wikiPage.subsubsection("Answer")
  wikiPage.putText("Edit Here")
  wikiPage.subsubsection("Relevant Information")
  imageDict=dict()
  indexDict=dict()
  thumbDict=dict()
  imageDictAQ=dict()
  indexDictAQ=dict()
  thumbDictAQ=dict()
  zValueDictAQ=dict()
  #Select only AUX channels
  for sngl in wikiCoinc.sngls:
    imageDict[sngl.ifo]=list()
    indexDict[sngl.ifo]=list()
    thumbDict[sngl.ifo]=list()
    for myFile in fnmatch.filter(wikiFileFinder.get_RDS_R_L1(),\
                                 "*/%s_RDS_*/%s/*html"%(sngl.ifo,sngl.time)):
      indexDict[sngl.ifo].append(myFile)
    for myFile in fnmatch.filter(wikiFileFinder.get_RDS_R_L1(),\
                                 "*/%s_RDS_*/%s/*_16.00_spectrogram_whitened.png"%\
                                 (sngl.ifo,sngl.time)):
      if not myFile.upper().__contains__("PEM"):
        imageDict[sngl.ifo].append(myFile)
        
    for myFile in fnmatch.filter(wikiFileFinder.get_RDS_R_L1(),\
                                 "*/%s_RDS_*/%s/*_16.00_spectrogram_whitened?thumb.png"%\
                                 (sngl.ifo,sngl.time)):
      if not myFile.upper().__contains__("PEM"):
        thumbDict[sngl.ifo].append(myFile)
    #Select associated analyzeQscans
    imageDictAQ[sngl.ifo]=list()
    indexDictAQ[sngl.ifo]=list()
    thumbDictAQ[sngl.ifo]=list()
    intS=nanS=0
    intS,nanS=str(float(sngl.time)).split(".")
    timeString="%s_%s"%(intS,nanS)
    #H1-analyseQscan_H1_931176926_116_rds-unspecified-gpstime.html
    for myFile in fnmatch.filter(wikiFileFinder.get_analyzeQscan_RDS(),\
                                 "*%s-*_%s_*html"%(sngl.ifo,timeString)):
      indexDictAQ[sngl.ifo].append(myFile)
    zValueFiles=fnmatch.filter(wikiFileFinder.get_analyzeQscan_RDS(),\
                               "*_%s_%s_*.txt"%(sngl.ifo,timeString))
    #Process zValue ranking file if found for IFO
    if len(zValueFiles) > 0:
      zValueDictAQ[sngl.ifo]=wikiFileFinder.__readZranks__(zValueFiles[0])
    else:
      zValueDictAQ[sngl.ifo]=list()
    #H1-analyseQscan_H1_931176926_116_rds_H0_PEM-MY_SEISX_z_scat-unspecified-gpstime_thumb.png
    #H1-analyseQscan_H1_931176926_116_rds_H0_PEM-MY_SEISX_z_scat-unspecified-gpstime.png
    for myFile in fnmatch.filter(wikiFileFinder.get_analyzeQscan_RDS(),\
                                 "*%s-*_%s_*_z_scat-unspecified-gpstime.png"%\
                                 (sngl.ifo,timeString)):
      if not myFile.upper().__contains__("PEM"):
        imageDictAQ[sngl.ifo].append(myFile)
        
    for myFile in fnmatch.filter(wikiFileFinder.get_analyzeQscan_RDS(),\
                                 "*%s-*_%s_*_z_scat-unspecified-gpstime?thumb.png"%\
                                 (sngl.ifo,timeString)):
      if not myFile.upper().__contains__("PEM"):
        thumbDictAQ[sngl.ifo].append(myFile)

    #Convert disk locals to URLs
    imageDict[sngl.ifo]=[file2URL.convert(x) for x in imageDict[sngl.ifo]]
    indexDict[sngl.ifo]=[file2URL.convert(x) for x in indexDict[sngl.ifo]]
    thumbDict[sngl.ifo]=[file2URL.convert(x) for x in thumbDict[sngl.ifo]]
    imageDictAQ[sngl.ifo]=[file2URL.convert(x) for x in imageDictAQ[sngl.ifo]]
    indexDictAQ[sngl.ifo]=[file2URL.convert(x) for x in indexDictAQ[sngl.ifo]]
    thumbDictAQ[sngl.ifo]=[file2URL.convert(x) for x in thumbDictAQ[sngl.ifo]]
    if len(indexDict[sngl.ifo]) < 1:
      wikiPage.putText("Other scans for %s not available.\n"%sngl.ifo)
  enoughImage=[len(imageDict[key])>0 for key in imageDict.keys()].count(True) >=1
  enoughIndex=[len(indexDict[key])>0 for key in indexDict.keys()].count(True) >=1
  if enoughImage and enoughIndex:
    wikiPage.insertAnalyzeQscanTable(imageDict,
                                     thumbDict,
                                     indexDict,
                                     imageDictAQ,
                                     thumbDictAQ,
                                     indexDictAQ,
                                     zValueDictAQ)
  else:
    sys.stdout.write("Warning: AUX plots import trouble.\n")
  wikiPage.subsubsection("Investigator Comments")
  wikiPage.putText("Edit Here")
  wikiPage.insertHR()
  #
  #
  #Additional Checklist Item
  wikiPage.subsection("#8 Electronic Log Book")
  wikiPage.subsubsection("Question")
  wikiPage.putText("Were the instruments behaving normally according to the comments posted by the sci-mons or the operators in the e-log?")
  wikiPage.subsubsection("Answer")
  wikiPage.putText("Edit Here")
  wikiPage.subsubsection("Relevant Information")
  wikiPage.putText("Plots and pipeline data go here!")
  wikiLinkLHOlog=wikiPage.makeExternalLink(stfu_pipe.getiLogURL(myCoinc.time,"H1"),
                                           "Hanford eLog")
  wikiLinkLLOlog=wikiPage.makeExternalLink(stfu_pipe.getiLogURL(myCoinc.time,"L1"),
                                           "Livingston eLog")
  wikiPage.putText("%s\n\n%s\n\n"%(wikiLinkLHOlog,wikiLinkLLOlog))
  wikiPage.subsubsection("Investigator Comments")
  wikiPage.putText("Edit Here")
  wikiPage.insertHR()
  #
  #
  #Additional Checklist Item
  wikiPage.subsection("#9 Glitch Report")
  wikiPage.subsubsection("Question")
  wikiPage.putText("Were the instruments behaving normally according to the weekly glitch report?")
  wikiPage.subsubsection("Answer")
  wikiPage.putText("Edit Here")
  wikiPage.subsubsection("Relevant Information")
  wikiPage.putText("Plots and pipeline data go here!")
  if int(wikiCoinc.time) >= endOfS5:
    wikiLinkGlitch=wikiPage.makeExternalLink(
      "https://www.lsc-group.phys.uwm.edu/twiki/bin/view/DetChar/GlitchStudies",
      "Glitch Reports for S6"
      )
  else:
    wikiLinkGlitch=wikiPage.makeExternalLink(
      "http://www.lsc-group.phys.uwm.edu/glitch/investigations/s5index.html#shift",
      "Glitch Reports for S5"
      )
  wikiPage.putText("%s\n"%(wikiLinkGlitch))
  wikiPage.subsubsection("Investigator Comments")
  wikiPage.putText("Edit Here")
  wikiPage.insertHR()
  #
  #
  #Additional Checklist Item
  wikiPage.subsection("#10 Snr versus time")
  wikiPage.subsubsection("Question")
  wikiPage.putText("Is this trigger significant in a SNR versus time plot of all triggers in its analysis chunk?")
  wikiPage.subsubsection("Answer")
  wikiPage.putText("Edit Here")
  wikiPage.subsubsection("Relevant Information")
  wikiPage.putText("Plots and pipeline data go here!")
  wikiPage.subsubsection("Investigator Comments")
  wikiPage.putText("Edit Here")
  wikiPage.insertHR()
  #
  #
  #Additional Checklist Item
  wikiPage.subsection("#11 Parameters of the candidate")
  wikiPage.subsubsection("Question")
  wikiPage.putText("Does the candidate have a high likelihood of being a gravitational-wave according to its parameters?")
  wikiPage.subsubsection("Answer")
  wikiPage.putText("Edit Here")
  wikiPage.subsubsection("Relevant Information")
  wikiPage.putText("Effective Distance Ratio Test\n")
  effDList=wikiFileFinder.get_effDRatio()
  if len(effDList) != 1:
    sys.stdout.write("Warning: Effective Distance Test import problem.\n")
  for myFile in effDList:
    wikiPage.putText("%s\n"%(file(myFile).read()))
  wikiPage.subsubsection("Investigator Comments")
  wikiPage.putText("Edit Here")
  wikiPage.insertHR()
  #
  #
  #Additional Checklist Item
  wikiPage.subsection("#12 Snr and Chisq")
  wikiPage.subsubsection("Question")
  wikiPage.putText("Are the SNR and CHISQ time series consistent with our expectations for a gravitational wave?")
  wikiPage.subsubsection("Answer")
                                  
  wikiPage.putText("Edit Here")
  wikiPage.subsubsection("Relevant Information")
  #
  #Put plots SNR and Chi sqr
  #
  indexList=fnmatch.filter(wikiFileFinder.get_plotsnrchisq(),"*.html")
  thumbList=fnmatch.filter(wikiFileFinder.get_plotsnrchisq(),"*_snr-*thumb.png")
  thumbList.extend(fnmatch.filter(wikiFileFinder.get_plotsnrchisq(),"*_chisq-*thumb.png"))
  thumbList.sort()
  indexList=[file2URL.convert(x) for x in indexList]
  thumbList=[file2URL.convert(x) for x in thumbList]
  #Two thumb types possible "_thumb.png" or ".thumb.png"
  imageList=[x.replace("_thumb.png",".png").replace(".thumb.png",".png") for x in thumbList]
  ifoCount=len(wikiCoinc.sngls)
  rowLabel={"SNR":1,"CHISQ":2}
  rowCount=len(rowLabel)
  colCount=ifoCount
  if len(indexList) >= 1:
    snrTable=wikiPage.wikiTable(rowCount+1,colCount+1)
    for i,sngl in enumerate(wikiCoinc.sngls):
      myIndex=""
      for indexFile in indexList:
        if indexFile.__contains__("_pipe_%s_FOLLOWUP_"%sngl.ifo):
          myIndex=indexFile
      if myIndex=="":
        snrTable.data[0][i+1]=" %s "%sngl.ifo
      else:
        snrTable.data[0][i+1]=wikiPage.makeExternalLink(myIndex,sngl.ifo)
    for col,sngl in enumerate(wikiCoinc.sngls):
      for row,label in enumerate(rowLabel.keys()):
        snrTable.data[row+1][0]=label
        for k,image in enumerate(imageList):
          if (image.__contains__("_%s-"%label.lower()) \
              and image.__contains__("pipe_%s_FOLLOWUP"%sngl.ifo)):
            snrTable.data[row+1][col+1]=" %s "%(wikiPage.linkedRemoteImage(thumbList[k],thumbList[k]))
    wikiPage.insertTable(snrTable)
  else:
    sys.stdout.write("Warning: SNR and CHISQ plots not found.\n")
    wikiPage.putText("SNR and CHISQ plots not found.\n")
  wikiPage.subsubsection("Investigator Comments")
  wikiPage.putText("Edit Here")
  wikiPage.insertHR()
  #
  #
  #Additional Checklist Item
  wikiPage.subsection("#13 Template bank veto")
  wikiPage.subsubsection("Question")
  wikiPage.putText("Is the bank veto value consistent with our expectations for a gravitational wave?")
  wikiPage.subsubsection("Answer")
  wikiPage.putText("Edit Here")
  wikiPage.subsubsection("Relevant Information")
  wikiPage.putText("Plots and pipeline data go here!")
  wikiPage.subsubsection("Investigator Comments")
  wikiPage.putText("Edit Here")
  wikiPage.insertHR()
  #
  #
  #Additional Checklist Item
  wikiPage.subsection("#14 Coherent studies")
  wikiPage.subsubsection("Question")
  wikiPage.putText("Are the triggers found in multiple interferometers coherent with each other?")
  wikiPage.subsubsection("Answer")
  wikiPage.putText("Edit Here")
  wikiPage.subsubsection("Relevant Information")
  wikiPage.putText("Plots and pipeline data go here!")
  indexList=fnmatch.filter(wikiFileFinder.get_plotchiatimeseries(),"*.html")
  if len(indexList) >= 1:
    myIndex=file2URL.convert(indexList[0])
    wikiPage.putText(wikiPage.makeExternalLink(myIndex,\
                                               "%s Coherence Study Results"%(wikiCoinc.ifos)))
    thumbList=fnmatch.filter(wikiFileFinder.get_plotchiatimeseries(),\
                             "PLOT_CHIA_%s_snr-squared*thumb.png"%(wikiCoinc.time))
    imageList=[x.replace("_thumb.png",".png").replace(".thumb.png",".png") for x in thumbList]
    rowCount=len(imageList)
    colCount=1
    cohSnrTimeTable=wikiPage.wikiTable(rowCount+1,colCount)
    cohSnrTimeTable.data[0][0]="%s Coherent SNR Squared Times Series"%(wikiCoinc.ifos)
    for i,image in enumerate(imageList):
      cohSnrTimeTable.data[i+1][0]=wikiPage.linkedRemoteImage(image,thumbList[i])
    wikiPage.insertTable(cohSnrTimeTable)
  else:
    sys.stdout.write("Warning: Coherent plotting jobs not found.\n")
    wikiPage.putText("Coherent Studies plots not found.\n")
  wikiPage.subsubsection("Investigator Comments")
  wikiPage.putText("Edit Here")
  wikiPage.insertHR()
  #
  #
  #Additional Checklist Item
  wikiPage.subsection("#15 Segmentation Stability")
  wikiPage.subsubsection("Question")
  wikiPage.putText("Is the candidate stable against changes in segmentation?")
  wikiPage.subsubsection("Answer")
  wikiPage.putText("Edit Here")
  wikiPage.subsubsection("Relevant Information")
  wikiPage.putText("Plots and pipeline data go here!")
  wikiPage.subsubsection("Investigator Comments")
  wikiPage.putText("Edit Here")
  wikiPage.insertHR()
  #
  #
  #Additional Checklist Item
  wikiPage.subsection("#16 Calibration Stability")
  wikiPage.subsubsection("Question")
  wikiPage.putText("Is the candidate stable against changes in calibration that are consistent with systematic uncertainties?")
  wikiPage.subsubsection("Answer")
  wikiPage.putText("Edit Here")
  wikiPage.subsubsection("Relevant Information")
  wikiPage.putText("Plots and pipeline data go here!")
  wikiPage.subsubsection("Investigator Comments")
  wikiPage.putText("Edit Here")
  wikiPage.insertHR()
  #
  #
#   #Additional Checklist Item
#   wikiPage.subsection("#17 Frame File Validation")
#   wikiPage.subsubsection("Question")
#   wikiPage.putText("Is the data used in the analysis free from corruption at the time of the candidate?")
#   wikiPage.subsubsection("Answer")
#   wikiPage.putText("Edit Here")
#   wikiPage.subsubsection("Relevant Information")
#   wikiPage.putText("Plots and pipeline data go here!")
#   wikiPage.subsubsection("Investigator Comments")
#   wikiPage.putText("Edit Here")
#   wikiPage.insertHR()

  #
  # A complete separate section in the WIKI page
  #
  wikiPage.section("Parameter Estimation")
  wikiPage.subsection("#1 Parameters of the candidate")
  wikiPage.subsubsection("Question")
  wikiPage.putText("Can we get more accurate information on the parameters of this candidate using MCMC or Bayesian methods?")
  wikiPage.subsubsection("Answer")
  wikiPage.putText("Edit Here")
  wikiPage.subsubsection("Relevant Information")
  wikiPage.putText("Plots and pipeline data go here!")
  wikiPage.subsubsection("Investigator Comments")
  wikiPage.putText("Edit Here")
  wikiPage.insertHR()
  #
  #
  wikiPage.subsection("#2 Coherent follow-up")
  wikiPage.subsubsection("Question")
  wikiPage.putText("Make a followup with coherent multi-detector code.")
  wikiPage.subsubsection("Answer")
  wikiPage.putText("Edit Here")
  wikiPage.subsubsection("Relevant Information")
  wikiPage.putText("Plots and pipeline data go here!")
  wikiPage.subsubsection("Investigator Comments")
  wikiPage.putText("Edit Here")
  wikiPage.insertHR()
  #
  #
  wikiPage.subsection("#3")
  wikiPage.subsubsection("Question")
  wikiPage.putText("Are the results of the Burst analysis astrophysically consistent with a possible detection?")
  wikiPage.subsubsection("Answer")
  wikiPage.putText("Edit Here")
  wikiPage.subsubsection("Relevant Information")
  wikiPage.putText("Plots and pipeline data go here!")
  wikiPage.subsubsection("Investigator Comments")
  wikiPage.putText("Edit Here")
  wikiPage.insertHR()
  #
  #
  wikiPage.subsection("#4")
  wikiPage.subsubsection("Question")
  wikiPage.putText("Are the results of a ringdown search astrophisycally consistent with a possible detection?")
  wikiPage.subsubsection("Answer")
  wikiPage.putText("Edit Here")
  wikiPage.subsubsection("Relevant Information")
  wikiPage.putText("Plots and pipeline data go here!")
  wikiPage.subsubsection("Investigator Comments")
  wikiPage.putText("Edit Here")
  wikiPage.insertHR()
  #
  #
  wikiPage.subsection("#5 EM Triggers")
  wikiPage.subsubsection("Question")
  wikiPage.putText("Are there any EM triggers in coincidence with the candidate? Is the distance estimated from interferometer time-delays or coherent analysis consistent with electro-magnetic observations? Are the distances as measured in several instruments consistent with position information?")
  wikiPage.subsubsection("Answer")
  wikiPage.putText("Edit Here")
  wikiPage.subsubsection("Relevant Information")
  wikiPage.putText("Plots and pipeline data go here!")
  wikiPage.subsubsection("Investigator Comments")
  wikiPage.putText("Edit Here")
  wikiPage.insertHR()
  #
  # Third miscellaneous section
  #
  wikiPage.section("Follow up documentation Miscellaneous Information")
  fuDocLink=wikiPage.makeExternalLink(
      "https://ldas-jobs.ligo.caltech.edu/~ctorres/followUpLivingDoc_LAST.pdf",
      "Living Follow up document"
      )
  wikiPage.putText(fuDocLink)
  wikiPage.insertHR()
  #
  # Wikie page done write out the content
  #
  wikiPage.write()
  #Close the file
  wikiPage.finish()
#
####################################################################
# End custom method prepareChecklist
####################################################################

####################################################################
# Coinc definition
####################################################################
class coinc(object):
  """
  """
  def __init__(self,coincInfoFile=None):
    """
    """
    if os.path.exists(coincInfoFile):
      inputData=open(coincInfoFile).readlines()
    else:
      Raise,"Error coinc info file not found!"
    #First line header
    #Second line is the Coinc Information
    rawCoinc=dict()
    rawCoincKeys=list()
    rawCoincData=list()
    rawCoincKeys=[x.replace("#","") for x in inputData[0].split()]
    rawCoincData=[x.replace(" ","") for x in inputData[1].split()]
    for i in range(0,len(rawCoincKeys)):
      rawCoinc[rawCoincKeys[i]]=rawCoincData[i]
    #Setup Coinc
    self.type=str(rawCoinc["DIR"])
    self.rank=float(rawCoinc["RANK"])
    self.far=float(rawCoinc["FAR"])
    self.snr=float(rawCoinc["SNR"])
    self.ifos=str(rawCoinc["IFOS"])
    self.instruments=str(rawCoinc["INSTRUMENTS"])
    self.time=float(rawCoinc["TIME"])
    self.mass=float(rawCoinc["MASS"])
    self.mchirp=float(rawCoinc["MCHIRP"])
    #Remaining header for sngl information
    rawSngl=list()
    rawSnglKeys=list()
    rawSnglData=list()
    rawSnglKeys=[x.replace("#","") for x in inputData[2].split()]
    rawSnglData=[x.split() for x in inputData[3:]]
    #Setup associated sngl data
    self.sngls=list()
    for rData in rawSnglData:
      tmp=dict()
      for i in range(0,len(rData)):
        tmp[rawSnglKeys[i]]=rData[i]
      self.sngls.append(sngl(tmp["DIR"],tmp["IFO"],tmp["TIME"],tmp["SNR"],tmp["CHISQ"],tmp["MASS1"],tmp["MASS2"],tmp["MCHIRP"]))
      del tmp

  def sngls_in_coinc(self):
    """
    """
    coincSngls=list()
    for sngl in self.sngls:
      if self.ifos.__contains__(sngl.ifo):
        coincSngls.append(sngl)
    return coincSngls

####################################################################
# Sngl definition
####################################################################
class sngl(object):
  """
  """
  def __init__(self,type=None,ifo=None,time=None,snr=None,chisqr=None,mass1=None,mass2=None,mchirp=None):
    """
    """
    self.type=str(type)
    self.ifo=str(ifo)
    self.time=float(time)
    self.snr=float(snr)
    self.chisqr=float(chisqr)
    self.mass1=float(mass1)
    self.mass2=float(mass2)
    self.mchirp=float(mchirp)
    
####################################################################
# Cache file parser
####################################################################

####################################################################
# Main part of script
####################################################################
usage = """usage: %prog [options]"""
userURL="~%s/"%(os.getenv("USER"))
userHOME="%s"%(os.getenv("HOME"))
hostnameURL="http://%s/"%(socket.gethostbyaddr(socket.gethostname())[0])
#
defaultWeblink="%s%sWEBPATH"%(hostnameURL,userURL)
parser = optparse.OptionParser(usage,version=git_version.verbose_msg)
parser.add_option("-f","--followup-directory",\
                  action="store",type="string",\
                  default="./",metavar="FUDIR",\
                  help="Set this to point the the parent directory of \
a follow up pipeline run.  From this location we \
will build checklists for all the events processed \
in this directory structure.")
parser.add_option("-i","--ini-file",\
                  action="store",type="string",\
                  default=None,\
                  help="Set this to point to the INI file \
used to run the follow up pipeline, the default action assumes \
the file is in the top directory of FUDIR.")

(opts,args) = parser.parse_args()
followup_directory=os.path.normpath(opts.followup_directory)
ini_file=opts.ini_file
#Read in the first ini file if none specified
if ini_file == None:
  inifilelist=scanTreeFnMatch(followup_directory,levels=1,filemask="*.ini")
  if len(inifilelist) < 1:
    raise Exception,"No potential ini files seen in %s\n"%(followup_directory)
  ini_file=inifilelist[0]
  sys.stdout.write("Will use file %s as ini config file.\n"%(ini_file))
if not os.path.exists(ini_file):
  raise Exception,"Path to ini file specified not found.\n"
#Extract information about webserver and web publication location
iniOpts=ConfigParser.ConfigParser()
iniOpts.read(ini_file)
publication_directory=None
publication_url=None
if iniOpts.has_option("fu-output","output-dir"):
  publication_directory=iniOpts.get("fu-output","output-dir")
else:
  raise Exception,"Ini file is missing options fu-output,output-dir.\n"
if iniOpts.has_option("fu-output","web-url"):
  publication_url=iniOpts.get("fu-output","web-url")
else:
  raise Exception,"Ini file is missing options fu-output,web-url.\n"
#
#
#Create static listing of pipe directory tree
#Create static listing of publication directory tree
#
pipeTree=scanTreeFnMatch(os.path.abspath(followup_directory),filemask="*")
omega_directory=publication_directory+"/omega/"
omegaTree=scanTreeFnMatch(os.path.abspath(omega_directory),filemask="*")
#
#
coincList=fnmatch.filter(pipeTree,"*coincEvent.info")
listCount=len(coincList)
for listsDone,coincFile in enumerate(coincList):
  #
  #Create directory for checklist in publication location.  We will
  #only move files not already in the html area to the publication
  #location.  This will speed things up since we don't need to redo
  #Qscan stuff.
  #
  myCoinc=coinc(coincFile)
  myChecklistFilename="CHECKLIST_%s_%s_%s_%s.wiki"%(myCoinc.type,
                                                    myCoinc.ifos,
                                                    myCoinc.instruments,
                                                    myCoinc.time)
  sys.stdout.write("Creating list (%i/%i):%s\n"%(listsDone+1,listCount,myChecklistFilename))
                                                    
  mySourcePath=os.path.abspath(followup_directory)
  myDestPath=os.path.abspath(publication_directory+"/"+myChecklistFilename.rstrip(".wiki")+"/")
  sys.stdout.write("Checklist is available at %s\n"%(myDestPath))
  if not os.path.exists(myDestPath):
    os.makedirs(myDestPath)
  #Scan for files required to make checklist.
  myFileFinderPipeTree=findFileType(pipeTree,myCoinc)
  myFileFinderOmegaTree=findFileType(omegaTree,myCoinc)
  allSources={'pipe':list(),
              'omega':list()}
  allSources['pipe'].extend(myFileFinderPipeTree.get_all())
  allSources['omega'].extend(myFileFinderOmegaTree.get_all())
  #Copy the files in allSource to CHECKLIST dir if not in publicationDirectory
  pud=os.path.abspath(publication_directory)
  minFileCount=1
  for key,fileList in allSources.items():
    if len(fileList) > minFileCount:
      commonPath=os.path.commonprefix(fileList)
      for singleFile in fileList:
        if not singleFile.__contains__(pud):
          myDestFile=singleFile.replace(commonPath,myDestPath+"/")
          if not os.path.exists(os.path.split(myDestFile)[0]):
            os.makedirs(os.path.split(myDestFile)[0])
          shutil.copy2(singleFile,myDestFile)
    else:
      sys.stdout.write("Warning: Scanning (%s) found %s files.\n"%\
                       (key,len(fileList)))
  # Create list of files used for checklist generation
  checklistTree=scanTreeFnMatch(myDestPath+"/",filemask="*")
  fileTree=list()
  fileTree.extend(checklistTree)
  fileTree.extend(allSources['omega'])
  mapFileURL=stfu_pipe.filenameToURLMapper(publication_directory,publication_url)
  prepareChecklist(myDestPath+"/"+myChecklistFilename,\
                   myCoinc,\
                   fileTree,\
                   mapFileURL)
  sys.stdout.write("Checklist is prepared.\n\n")
