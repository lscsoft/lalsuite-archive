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
from pylal.webCondor import *
from lalapps import inspiral

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
      bankFile = 'trigTemplateBank/' + ifo + '-TRIGBANK_FOLLOWUP_' + str(trig.gpsTime[ifo]) + '.xml.gz'
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
        if param == 'write-compress': 
          extension = '.xml.gz'
        else:
          extension = '.xml'

      # the output_file_name is required by the child job (plotSNRCHISQNode)
      self.output_file_name = inspJob.outputPath + self.inputIfo + "-INSPIRAL_" + self.__ifotag + "_" + self.__usertag + "-" + self.__start + "-" + str(int(self.__end)-int(self.__start)) + extension

      self.set_id(self.inputIfo + "-INSPIRAL_" + self.__ifotag + "_" + self.__usertag + "-" + self.__start + "-" + str(int(self.__end)-int(self.__start)))

      self.outputCache = ifo + ' ' + 'INSPIRAL' + ' ' + str(self.__start) + ' ' + str(int(self.__end)-int(self.__start)) + ' ' + self.output_file_name  + '\n' + ifo + ' ' + 'INSPIRAL-FRAME' + ' ' + str(self.__start) + ' ' + str(int(self.__end)-int(self.__start)) + ' ' + self.output_file_name.replace(extension,".gwf") + '\n'

      self.setupNodeWeb(inspJob,False,None,None,None,dag.cache)
      self.add_var_opt("output-path",inspJob.outputPath)

      if opts.inspiral:
        dag.addNode(self,'inspiral')
        self.validate()
      else: self.invalidate()
    except:
      self.invalidate()
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
  def __init__(self,job,ifo,fileName,trig,page,dag,inspiralNode,opts):
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
      self.setupNodeWeb(job,True, dag.webPage.lastSection.lastSub,page,None,dag.cache)
      try: 
        if inspiralNode.validNode: self.add_parent(inspiralNode)
      except: pass
      if opts.plots:
        dag.addNode(self,self.friendlyName)
        self.validate()
      else: self.invalidate()
    except: 
      self.invalidate()
      print "couldn't add plot job for " + str(ifo) + "@ "+ str(trig.gpsTime[ifo])

############### QSCAN CLASSES #################################################
###############################################################################

class qscanDataFindJob(pipeline.LSCDataFindJob,webTheJob):

  def __init__(self, config_file, source):
    self.name = 'qscanDataFindJob'
    # unfortunately the logs directory has to be created before we call LSCDataFindJob
    try:
      os.mkdir(self.name)
      os.mkdir(self.name + '/logs')
    except: pass
    pipeline.LSCDataFindJob.__init__(self, self.name, self.name + '/logs', config_file)
    if source == 'futrig':
      self.setup_cacheconv(config_file)
    self.setupJobWeb(self.name) # this will overwrite the stdout and stderr set up by LSCDataFindJob

  def setup_cacheconv(self,cp):
    # create a shell script to call convertlalcache.pl if the value of $RETURN is 0
    convert_script = open(self.name + '/cacheconv.sh','w')
    convert_script.write("""#!/bin/bash
    if [ ${1} -ne 0 ] ; then
      exit 1
    else
      %s ${2} ${3}
    fi
    """ % string.strip(cp.get('condor','convertcache')))
    convert_script.close()
    os.chmod(self.name + '/cacheconv.sh',0755)

class qscanDataFindNode(pipeline.LSCDataFindNode,webTheNode):
 
  def __init__(self, job, source, type, cp, time, ifo, opts, dag, prev_dNode, datafindCommand):
    try:
      pipeline.LSCDataFindNode.__init__(self,job)
      self.id = str(ifo) + '-' + repr(time) + '-' + str(type)
      self.setupNodeWeb(job,False,None,None,None,dag.cache)
      if source == 'futrig':
        self.outputFileName = self.setup_fu_trig(job, cp, time, ifo, type)
      try:
        if prev_dNode and eval('opts.' + datafindCommand):
          self.add_parent(prev_dNode)
      except: pass

      if eval('opts.' + datafindCommand):
        dag.addNode(self,'qscan data find')
        self.validNode = True
      else: self.validNode = False
    except:
      self.validNode = False
      print >> sys.stderr, "could not set up the datafind jobs for " + type


  def setup_fu_trig(self, job, cp, time, ifo, type):
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
    self.set_post_script(job.name + "/cacheconv.sh $RETURN %s %s" %(lalCache,qCache) )
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
    self.id = ifo + '-' + name + '-' + repr(time)

    pipeline.CondorDAGNode.__init__(self,job)
    self.add_var_arg(repr(time))
    qscanConfig = string.strip(cp.get(name, ifo + 'config-file'))
    self.add_file_arg(qscanConfig)
    self.add_file_arg(qcache)
    output = string.strip(cp.get(name, ifo + 'output'))
    self.add_var_arg(output)

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
    self.outputCache = ifo + ' ' + name + ' ' + repr(time) + ' ' + self.outputName + '\n'

    #try to extract web output from the ini file, else ignore it
    try:
      self.setupNodeWeb(job,False,dag.webPage.lastSection.lastSub,page,string.strip(cp.get(name,ifo+'web'))+repr(time),dag.cache)
    except: 
      self.setupNodeWeb(job,False,None,None,None,dag.cache) 

    # only add a parent if it exists
    try:
      if d_node.validNode and eval('opts.' + datafindCommand):
        self.add_parent(d_node)
    except: pass

    if eval('opts.' + qscanCommand):
      dag.addNode(self,self.friendlyName)
      self.validNode = True
    else: self.validNode = False
 #   except: 
 #     self.validNode = False
 #     print >> sys.stderr, "could not set up the qscan job for " + self.id


##############################################################################
# analyse qscan class: the job

class analyseQscanJob(pipeline.CondorDAGJob, webTheJob):
  """
  A followup analyseQscan job to interprete the qscans
  """
  def __init__(self,options,cp,tag_base='ANALYSE_QSCAN'):
    self.__name__ = 'analyseQscanJob'
    self.__executable = string.strip(cp.get('condor','analyseQscan'))
    self.__universe = "vanilla"
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    self.add_condor_cmd('getenv','True')
    self.setupJobWeb(self.__name__,tag_base)

##############################################################################
# analyse qscan class: the node

class analyseQscanNode(pipeline.CondorDAGNode,webTheNode):
  """
  Runs an instance of a followup analyseQscan job
  """
  def __init__(self,job,time,ifo,name,foregroundCache,backgroundCache,cp,opts,dag,command):
    """
    job = A CondorDAGJob that can run an instance of analyseQscan followup.
    """
    page = string.strip(cp.get('output','page'))
    self.friendlyName = 'analyse ' + name
    self.id = ifo + '-' + name + '-' + repr(time)

    try:
      pipeline.CondorDAGNode.__init__(self,job)
      if cp.has_option('analyse-qscan','generate-qscan-xml'):
        self.add_var_opt('generate-qscan-xml','')
      self.add_var_opt('z-threshold',cp.getfloat('analyse-qscan','z-threshold'))
      if cp.has_option('analyse-qscan','plot-z-distribution'):
        self.add_var_opt('plot-z-distribution','')
        self.add_var_opt('z-min',cp.getfloat('analyse-qscan','z-min'))
        self.add_var_opt('z-max',cp.getfloat('analyse-qscan','z-max'))
        self.add_var_opt('z-bins',cp.getfloat('analyse-qscan','z-bins'))
      if cp.has_option('analyse-qscan','plot-dt-distribution'):
        self.add_var_opt('plot-dt-distribution','')
        self.add_var_opt('dt-min',cp.getfloat('analyse-qscan','dt-min'))
        self.add_var_opt('dt-max',cp.getfloat('analyse-qscan','dt-max'))
        self.add_var_opt('dt-bins',cp.getfloat('analyse-qscan','dt-bins'))
      if cp.has_option('analyse-qscan','plot-z-scattered'):
        self.add_var_opt('plot-z-scattered','')
      if cp.has_option('analyse-qscan','plot-z-scattered') or cp.has_option('analyse-qscan','plot-dt-distribution'):
        self.add_var_opt('ref-channel',cp.get('analyse-qscan','ref-channel'))
      self.add_var_opt('qscan-id',name + '_' + ifo + '_' + repr(time)) 

      self.add_var_opt('qscan-cache-foreground',foregroundCache)
      self.add_var_opt('qscan-cache-background',backgroundCache)

      self.setupNodeWeb(job,True,None,page,None,dag.cache)

      # get the table of the qscan job associated to this trigger
      for node in dag.get_nodes():
        if isinstance(node,qscanNode):
          if node.id == self.id:
            # link the analyseQscan output page to the qscan table
            node.webTable.row[0].cell[0].linebreak()
            node.webTable.row[0].cell[0].link(self.webLink,"qscan background vs qscan foreground")
            break      

      for node in dag.get_nodes():
        if isinstance(node,qscanNode): 
          if node.validNode:
            if node.friendlyName == name or \
            node.friendlyName.replace('background','foreground') == name:
              self.add_parent(node)

      if eval('opts.' + command):
        dag.addNode(self,self.friendlyName)
        self.validNode = True
      else: self.validNode = False

    except:
      self.validNode = False
      print "couldn't add " + name + " analyseQscan job for " + ifo + "@ "+ repr(time)

###############################################################################
# FrCheck Jobs and Nodes

class FrCheckJob(pipeline.CondorDAGJob, webTheJob):
  """
  A followup job for checking frames
  """
  def __init__(self, options, cp, tag_base='FRCHECK'):
    """
    """
    self.__name__ = 'FrCheckJob'
    self.__executable = string.strip(cp.get('condor','frame_check'))
    self.__universe = "vanilla"
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    self.add_condor_cmd('getenv','True')
    self.setupJobWeb(self.__name__,tag_base)


class FrCheckNode(pipeline.CondorDAGNode,webTheNode):
  """
  Runs an instance of a FrCheck followup job
  """
  def __init__(self, FrCheckJob, procParams, ifo, trig, cp,opts,dag):

    for row in procParams:
      param = row.param.strip("-")
      value = row.value
      if param == 'frame-cache': cacheFile = value 

    self.friendlyName = 'Frame Check'

    
    pipeline.CondorDAGNode.__init__(self,FrCheckJob)
    self.add_var_opt("frame-cache", cacheFile)
    self.add_var_opt("frame-check-executable", string.strip(cp.get('frameCheck','executable')))

    self.id = FrCheckJob.name + '-' + ifo + '-' + str(trig.statValue) + '_' + str(trig.eventID)
    self.setupNodeWeb(FrCheckJob,True, dag.webPage.lastSection.lastSub,dag.page,None,dag.cache)
    if opts.frame_check:
      dag.addNode(self,self.friendlyName)
      self.validate()
    else: self.invalidate()
#    except:
#      self.invalidate()
#      print "couldn't add frame check job for " + str(ifo) + "@ "+ str(trig.gpsTime[ifo])


    

