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
from pylal.fu_Web import *
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

    

###### WRAPPER FOR CONDOR DAG - TO MAKE THE FOLLOWUP DAG WEBIFIABLE ###########
###############################################################################

class followUpDAG(pipeline.CondorDAG, webTheDAG):

  def __init__(self, config_file, log_path):
    self.basename = re.sub(r'\.ini',r'', config_file) 
    tempfile.tempdir = log_path
    tempfile.template = self.basename + '.dag.log.'
    logfile = tempfile.mktemp()
    fh = open( logfile, "w" )
    fh.close()
    pipeline.CondorDAG.__init__(self,logfile)
    self.set_dag_file(self.basename)
    self.jobsDict = {}
 
  def printNodeCounts(self):
    for jobs in self.jobsDict:
      print "\nFound " + str(self.jobsDict[jobs]) + " " + str(jobs) + " Jobs"
 
  def writeAll(self, type='IUL'):
    self.printNodeCounts()
    print "\n\n.......Writing DAG"
    self.write_sub_files()
    self.write_dag()
    self.writeDAGWeb(type)
    print "\n\n  Created a DAG file which can be submitted by executing"
    print "    condor_submit_dag " + self.get_dag_file()
    print """\n  from a condor submit machine
  Before submitting the dag, you must execute

    export _CONDOR_DAGMAN_LOG_ON_NFS_IS_ERROR=FALSE

  If you are running LSCdataFind jobs, do not forget to initialize your grid
  proxy certificate on the condor submit machine by running the commands

    unset X509_USER_PROXY
    grid-proxy-init -hours 72

  Enter your pass phrase when prompted. The proxy will be valid for 72 hours.
  If you expect the LSCdataFind jobs to take longer to complete, increase the
  time specified in the -hours option to grid-proxy-init. You can check that
  the grid proxy has been sucessfully created by executing the command:

    grid-cert-info -all -file /tmp/x509up_u`id -u`

  This will also give the expiry time of the proxy."""



#### A CLASS TO DO FOLLOWUP INSPIRAL JOBS ####################################
###############################################################################
class followUpInspJob(inspiral.InspiralJob,webTheJob):
  def __init__(self,cp):
    inspiral.InspiralJob.__init__(self,cp)
    self.name = 'followUpInspJob'
    self.setupJobWeb(self.name)
    

class followUpInspNode(inspiral.InspiralNode,webTheNode):
  
  def __init__(self, inspJob, procParams, ifo, trig, cp,opts,dag):    
   
    try:
      inspiral.InspiralNode.__init__(self, inspJob) 
      injFile = self.checkInjections(cp)
      bankFile = ifo + '-TRIGBANK_FOLLOWUP_' + str(trig.gpsTime[ifo]) + '.xml'
      self.add_var_opt("write-snrsq","")
      self.add_var_opt("write-chisq","")
      self.add_var_opt("write-spectrum","")
      self.set_bank(bankFile)
      self.set_user_tag("FOLLOWUP_" + str(trig.gpsTime[ifo]))
      self.__usertag = "FOLLOWUP_" + str(trig.gpsTime[ifo])
      self.set_trig_start( int(trig.gpsTime[ifo]) - 1)
      self.set_trig_end( int(trig.gpsTime[ifo]) + 1 )
      if injFile: self.set_injections( injFile )
      skipParams = ['minimal-match', 'bank-file', 'user-tag', 'injection-file', 'trig-start-time', 'trig-end-time']

      for row in procParams:
        param = row.param.strip("-")
        value = row.value
        if param in skipParams: continue
        self.add_var_opt(param,value)
        if param == 'gps-end-time': self.__end = value
        if param == 'gps-start-time': self.__start = value
        if param == 'ifo-tag': self.__ifotag = value
        if param == 'channel-name': self.inputIfo = value[0:2]

      # the output_file_name is required by the child job (plotSNRCHISQNode)
      self.output_file_name = self.inputIfo + "-INSPIRAL_" + self.__ifotag + "_" + self.__usertag + "-" + self.__start + "-" + str(int(self.__end)-int(self.__start)) + ".xml"
      self.id =  self.inputIfo + "-INSPIRAL_" + self.__ifotag + "_" + self.__usertag + "-" + self.__start + "-" + str(int(self.__end)-int(self.__start))
      self.setupNodeWeb(inspJob)
      if opts.inspiral:
        dag.addNode(self,'inspiral')
        self.validNode = True
      else: self.validNode = False
    except:
      self.validNode = False
      print "couldn't add inspiral job for " + str(ifo) + "@ "+ str(trig.gpsTime[ifo])


  def checkInjections(self,cp):
    if len(string.strip(cp.get('triggers','injection-file'))) > 0:
      injectionFile = string.strip(cp.get('triggers','injection-file'))
    else:
      injectionFile = None
   
########## PLOT SNR  CHISQ TIME SERIES ########################################
###############################################################################

##############################################################################
# jobs class for plot snr chisq 

class plotSNRCHISQJob(pipeline.CondorDAGJob,webTheJob):
  """
  A followup plotting job for snr and chisq time series
  """
  def __init__(self, options, cp, tag_base='PLOT_FOLLOWUP'):
    """
    """
    self.__name__ = 'plotSNRCHISQJob'
    self.__executable = string.strip(cp.get('condor','plotsnrchisq'))
    self.__universe = "vanilla"
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    self.add_condor_cmd('getenv','True')
    self.setupJobWeb(self.__name__,tag_base)

##############################################################################
# node class for plot snr chisq

class plotSNRCHISQNode(pipeline.CondorDAGNode,webTheNode):
  """
  Runs an instance of a plotSNRCHISQ followup job
  """
  def __init__(self,job,ifo,fileName,trig,page,dag, inspiralNode,opts,prev_plotnode):
    """
    job = A CondorDAGJob that can run an instance of plotSNRCHISQ followup.
    """
    time = trig.gpsTime[ifo]
    self.friendlyName = 'Plot SNR/CHISQ/PSD'
    try:
      pipeline.CondorDAGNode.__init__(self,job)
      self.output_file_name = ""
      self.add_var_opt("frame-file",fileName.replace(".xml",".gwf"))
      self.add_var_opt("gps",time)
      self.add_var_opt("inspiral-xml-file",fileName)
      self.id = job.name + '-' + ifo + '-' + str(trig.statValue) + '_' + str(trig.eventID)
      self.setupNodeWeb(job,True, dag.webPage.lastSection.lastSub,page)
      try: 
        if inspiralNode.validNode: self.add_parent(inspiralNode)
      except: pass
      #try: 
      #  if self.idCheck(prev_plotNode): self.add_parent(prev_plotNode)
      #except: pass
      if opts.plots:
        dag.addNode(self,self.friendlyName)
        prev_plotNode = self
        self.validNode = True
      else: self.validNode = False
    except: 
      self.validNode = False
      print "couldn't add plot job for " + str(ifo) + "@ "+ str(trig.gpsTime[ifo])

  
  def idCheck(self, prev_plotNode):
    if (self.id == prev_plotNode.id) and prev_plotNode.validNode: return True
    else: return False



############### QSCAN CLASSES #################################################
###############################################################################

class qscanDataFindJob(pipeline.LSCDataFindJob,webTheJob):
  
  def __init__(self, cache_dir, log_dir, config_file, source):
    pipeline.LSCDataFindJob.__init__(self, cache_dir, log_dir, config_file)
    if source == 'futrig':
      self.setup_cacheconv(config_file)
    self.name = 'qscanDataFindJob'
    self.setupJobWeb(self.name)

  def setup_cacheconv(self,cp):
    # create a shell script to call convertlalcache.pl if the value of $RETURN is 0
    convert_script = open('cacheconv.sh','w')
    convert_script.write("""#!/bin/bash
    if [ ${1} -ne 0 ] ; then
      exit 1
    else
      %s ${2} ${3}
    fi
    """ % string.strip(cp.get('condor','convertcache')))
    convert_script.close()
    os.chmod('cacheconv.sh',0755)

class qscanDataFindNode(pipeline.LSCDataFindNode,webTheNode):
 
  def __init__(self, job, source, type, cp, time, ifo, opts,  prev_dNode, dag):
    try:
      pipeline.LSCDataFindNode.__init__(self,job)
      self.id = str(ifo) + '-' + repr(time) + '-' + str(type) + 'datafind'
      self.setupNodeWeb(job)
      if source == 'futrig':
        self.outputFileName = self.setup_fu_trig(cp, time, ifo, type)
      try:
        if prev_dNode and opts.datafind: dNode.add_parent(prev_dNode)
      except: pass

      if opts.datafind:
        dag.addNode(self,'qscan data find')
        prev_dNode = self
        self.validNode = True
      else: self.validNode = False
    except:
      self.validNode = False
      print >> sys.stderr, "could not set up the datafind jobs for " + type

    
  def setup_fu_trig(self, cp, time, ifo, type):
    # 1s is substracted to the expected startTime to make sure the window
    # will be large enough. This is to be sure to handle the rouding to the
    # next sample done by qscan.
    self.q_time = cp.getint(type,'search-time-range')/2
    self.set_observatory(ifo[0])
    self.set_start(int( time - self.q_time - 1))
    self.set_end(int( time + self.q_time + 1))
    if cp.has_option(type, ifo + '_type'): 
      self.set_type( cp.get(type, ifo + '_type' ))
    else:
      self.set_type( cp.get(type, 'type' ))
    lalCache = self.get_output()
    qCache = lalCache.rstrip("cache") + "qcache"
    self.set_post_script("cacheconv.sh $RETURN %s %s" %(lalCache,qCache) )
    return(qCache)    

##############################################################################
# qscan class for qscan jobs

class qscanJob(pipeline.CondorDAGJob, webTheJob):
  """
  A qscan job
  """
  def __init__(self, opts, cp, tag_base='QSCAN'):
    """
    """
    self.__executable = string.strip(cp.get('condor','qscan'))
    self.__universe = "vanilla"
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    #self.setupJobWeb(self.__name__,tag_base)
    self.setupJobWeb(tag_base)

class qscanLiteJob(pipeline.CondorDAGJob, webTheJob):
  """
  A qscanLite job
  """
  def __init__(self, opts, cp, tag_base='QSCANLITE'):
    """
    """
    self.__executable = string.strip(cp.get('condor','qscanlite'))
    self.__universe = "vanilla"
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    #self.setupJobWeb(self.__name__,tag_base)
    self.setupJobWeb(tag_base)


##############################################################################
# qscan class for qscan Node

class qscanNode(pipeline.CondorDAGNode,webTheNode):
  """
  Runs an instance of a qscan job
  """
  def __init__(self,job,time,cp,qcache,ifo,name, opts, d_node, dag, datafindCommand, qscanCommand, trig=None,qFlag=None):
    """
    job = A CondorDAGJob that can run an instance of qscan.
    """
    page = string.strip(cp.get('output','page'))
    self.friendlyName = name
    try:
      pipeline.CondorDAGNode.__init__(self,job)
      self.add_var_arg(repr(time))
      qscanConfig = string.strip(cp.get(name, ifo + 'config-file'))
      self.add_file_arg(qscanConfig)
      self.add_file_arg(qcache)
      output = string.strip(cp.get(name, ifo + 'output'))
      self.add_var_arg(output)
      self.id = ifo + repr(time)
      #try to extract web output from the ini file, else ignore it
      try: self.setupNodeWeb(job,False,dag.webPage.lastSection.lastSub,page,string.strip(cp.get(name,ifo+'web'))+repr(time))
      except: self.setupNodeWeb(job,False) 
      #get the absolute output path whatever the path might be in the ini file
      currentPath = os.path.abspath('.')
      try:
        os.chdir(output)
        absoutput = os.path.abspath('.')
        os.chdir(currentPath)
      except:
        print >> sys.stderr, 'invalid path for qscan output in the ini file'
        sys.exit(1)
      self.outputName = absoutput + '/' + repr(time) # redirect output name
    
      #prepare the string for the output cache
      self.outputCache = repr(time) + '\t' + name + '\t' + ifo + '\t' + self.outputName
    
      # only add a parent if it exists
      try:
        if d_node.validNode and eval('opts.' + datafindCommand):
          self.add_parent(d_node)
      except: pass

      if eval('opts.' + qscanCommand):
        dag.addNode(self,self.friendlyName)
        self.validNode = True
      else: self.validNode = False
    except: 
      self.validNode = False
      print >> sys.stderr, "could not set up the background qscan jobs"





