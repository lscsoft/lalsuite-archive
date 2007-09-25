#!/usr/bin/env @PYTHONPROG@
"""
This module contains condor jobs / node classes for the followup dag

$Id$

This program creates cache files for the output of inspiral hipe
"""

__author__ = 'Chad Hanna <channa@phys.lsu.edu>'
__date__ = '$Date$'
__version__ = '$Revision$'[11:-2]

##############################################################################
# import standard modules and append the lalapps prefix to the python path
import sys, os, copy, math
import socket, time
import re, string
from optparse import *
import tempfile
import ConfigParser 
import urlparse
from UserDict import UserDict
sys.path.append('@PYTHONLIBDIR@')

##############################################################################
# import the modules we need to build the pipeline
from glue import pipeline
from glue import lal
from glue import segments
from glue import segmentsUtils
from pylal.webUtils import *
from lalapps import inspiral
  
###### GENERIC CLASSES TO HELP WEBIFY A CONDOR DAG OUTPUT ####################
##############################################################################

class webTheJob:
  """
  webTheJob is a class intended to be inherited by a class that subclasses
  the condor DAG Job.  It is useful for setting up a standard structure to
  webify the output of a dag.  You'll want to use webTheNode and webTheDAG
  for your condor subclasses too. 
  """

  def __init__(self):
    pass

  def setupJobWeb(self, name, tag_base=None):
    # Give this job a name.  Then make directories for the log files and such
    # This name is important since these directories will be included in
    # the web tree.
    self.name = name
    try:
       os.mkdir(name)
       os.mkdir(name+'/logs')
    except: pass
    # Set up the usual stuff and name the log files appropriately
    self.tag_base = tag_base
    self.add_condor_cmd('environment',"KMP_LIBRARY=serial;MKL_SERIAL=yes")
    self.set_sub_file(name+'.sub')
    self.relPath = name + '/'
    self.outputPath = os.getcwd() + '/' + name + '/'
    self.set_stdout_file(self.outputPath+'/logs/'+name+'-$(macroid).out')
    self.set_stderr_file(self.outputPath+'/logs/'+name+'-$(macroid).err')

class webTheNode:
  """
  webTheNode is a class intended to be inherited by a class that subclasses
  the condor DAG Node .  It is useful for setting up a standard structure to
  webify the output of a dag.  You'll want to use webTheJob and webTheDAG
  for your condor subclasses too. 
  """

  def __init__(self):
    pass

  def setupNodeWeb(self, job, passItAlong=True, content=None, page=None,webOverride=None):
    # setup the node id
    self.add_macro("macroid", self.id)
    # determine the output web file name for the job
    self.webFileName = job.outputPath + self.id + '.html'
    self.jobName = job.name
    if page:
      self.webLink = page+'/'+job.relPath+self.id+'.html'
    if webOverride:
      self.webLink = webOverride
    # standardize the I/O for executables that themselves write web pages.
    # this is great for developers since they can in some sense treat
    # the DAG as a black box for their executable.   They just have to 
    # comply with these few command line arguments to tie everything together
    if passItAlong:
      self.add_var_opt("output-web-file",self.webFileName)
      self.add_var_opt("output-path",job.outputPath)
      self.add_var_opt("page-rel-path",job.relPath)
      self.add_var_opt("page", page)
    if content: self.writeContent(content)

  def writeContent(self,content):
    # The talkBack class is a way for the users job to provide information
    # back to the DAG web.  It is done through the reading and writing of
    # a config file (.ini) that has the same naming convention as the
    # web file
    self.talkBack = talkBack(self.webFileName)
    self.talkBack.read()
    content.appendTable(1,2,0,600)
    content.lastTable.row[0].cell[0].link(self.webLink,self.friendlyName)
    # Each time the dag is generated it checks for the existance of the
    # appropriate config file to include the contents in the web page
    if self.talkBack.summaryText:
      content.lastTable.row[0].cell[0].linebreak()
      content.lastTable.row[0].cell[0].text(self.talkBack.summaryText)
    if self.talkBack.summaryPlot:
      content.lastTable.row[0].cell[1].image(self.talkBack.summaryPlot)
    if self.talkBack.summaryPlotCaption:
      content.lastTable.row[0].cell[1].linebreak()
      content.lastTable.row[0].cell[1].text(self.talkBack.summaryPlotCaption)


class webTheDAG:
  """
  webTheDAG is a class intended to be inherited by a class that subclasses
  the condor DAG Node .  It is useful for setting up a standard structure to
  webify the output of a dag.  You'll want to use webTheJob and webTheDAG
  for your condor subclasses too. 
  """

  def __init__(self):
    pass

  def setupDAGWeb(self,title,filename,root=""):
    self.page = root
    self.webPage = WebPage(title,filename,root)
    self.webDirs = {}
    try:
       os.mkdir('DAGWeb')
    except: pass


  def writeDAGWeb(self,type):
    self.webPage.cleanWrite(type)

  def appendSection(self,name):
    self.webPage.appendSection(name)
    inifile = name.replace(" ","_").replace("@","-").replace("=",'-') + '.ini'
    file = open('DAGWeb/'+inifile,'a')
    file.close()
    talkback = talkBack('DAGWeb/'+inifile)
    talkback.read()
    self.webPage.lastSection.appendTable(1,2,0,600)

    if talkback.summaryText:
      self.webPage.lastSection.lastTable.row[0].cell[0].linebreak()
      self.webPage.lastSection.lastTable.row[0].cell[0].text(talkback.summaryText)
    if talkback.summaryPlot:
      self.webPage.lastSection.lastTable.row[0].cell[1].image(talkback.summaryPlot)
    if talkback.summaryPlotCaption:
      self.webPage.lastSection.lastTable.row[0].cell[1].linebreak()
      self.webPage.lastSection.lastTable.row[0].cell[1].text(talkback.summaryPlotCaption)

  def appendSubSection(self,name):
    self.webPage.lastSection.appendSubSection(name)
    inifile = name.replace(" ","_").replace("@","-").replace("=",'-') + '.ini'
    file = open('DAGWeb/'+inifile,'a')
    file.close()
    talkback = talkBack('DAGWeb/'+inifile)
    talkback.read()
    self.webPage.lastSection.lastSub.appendTable(1,2,0,600)

    if talkback.summaryText:
      self.webPage.lastSection.lastSub.lastTable.row[0].cell[0].linebreak()
      self.webPage.lastSection.lastSub.lastTable.row[0].cell[0].text(talkback.summaryText)
    if talkback.summaryPlot:
      self.webPage.lastSection.lastSub.lastTable.row[0].cell[1].image(talkback.summaryPlot)
    if talkback.summaryPlotCaption:
      self.webPage.lastSection.lastSub.lastTable.row[0].cell[1].linebreak()
      self.webPage.lastSection.lastSub.lastTable.row[0].cell[1].text(talkback.summaryPlotCaption)

  def addNode(self, node,jobType):
    try:
      self.jobsDict[jobType] = self.jobsDict[jobType] + 1
      self.webDirs[node.jobName] = node.jobName
    except:
      self.jobsDict[jobType] = 1
    self.add_node(node)


  def publishToHydra(self):
    dirStr = ''
    for dir in self.webDirs:
      dirStr += dir + ' '
    dirStr = 'rsync -vrz '+dirStr+' DAGWeb index.html '
    print dirStr
    os.system(dirStr+'hydra.phys.uwm.edu:/home/htdocs/uwmlsc/root/.'+self.page)


