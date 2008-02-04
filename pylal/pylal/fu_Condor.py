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
import math
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
from pylal import fu_utils

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

  def __init__(self,cp,type='plot'):

    inspiral.InspiralJob.__init__(self,cp)
    if type == 'head':
      self.set_executable(string.strip(cp.get('condor','inspiral_head')))
    self.name = 'followUpInspJob' + type
    self.setupJobWeb(self.name)


class followUpInspNode(inspiral.InspiralNode,webTheNode):
  
  def __init__(self, inspJob, procParams, ifo, trig, cp,opts,dag, type='plot',sngl_table = None):    

    try:
      inspiral.InspiralNode.__init__(self, inspJob) 
      injFile = self.checkInjections(cp)      

      if type == 'plot':
        bankFile = 'trigTemplateBank/' + ifo + '-TRIGBANK_FOLLOWUP_' + str(trig.eventID) + '.xml.gz'
        self.add_var_opt("write-snrsq","")
        self.add_var_opt("write-chisq","")
        self.add_var_opt("write-spectrum","")
        self.set_bank(bankFile)
        self.set_trig_start( int(trig.gpsTime[ifo]) - 1)
        self.set_trig_end( int(trig.gpsTime[ifo]) + 1 )

      if injFile: 
        self.set_injections( injFile )
      
      skipParams = ['minimal-match', 'bank-file', 'user-tag', 'injection-file', 'trig-start-time', 'trig-end-time']

      for row in procParams:
        param = row.param.strip("-")
        value = row.value
        if param == 'bank-file':
          bankFile = value
        if param in skipParams: continue
        self.add_var_opt(param,value)
        if param == 'gps-end-time': self.__end = value
        if param == 'gps-start-time': self.__start = value
        if param == 'ifo-tag':
          self.__ifotag = value
        if param == 'channel-name': self.inputIfo = value[0:2]
        if param == 'write-compress': 
          extension = '.xml.gz'
        else:
          extension = '.xml'

      if not ifo == self.inputIfo:
        second_user_tag = "_" + ifo + "tmplt"
      else:
        second_user_tag = ""
      self.set_user_tag("FOLLOWUP_" + str(trig.eventID) + second_user_tag)
      self.__usertag = "FOLLOWUP_" + str(trig.eventID) + second_user_tag


      # THIS IS A HACK FOR NOW, THERE IS PROBABLY A BETTER WAY TO DO THIS
      if (type == 'head'): 
        subBankSize = string.strip(cp.get('inspiral-head','bank-veto-subbank-size'))
        if opts.inspiral_head:
          bankFileName = fu_utils.generateBankVetoBank(trig, ifo, str(trig.gpsTime[ifo]), sngl_table[ifo],int(subBankSize),'BankVetoBank')
        else: bankFileName = 'none'      
        self.add_var_opt("bank-veto-subbank-size", string.strip(cp.get('inspiral-head','bank-veto-subbank-size')))
        self.add_var_opt("order", string.strip(cp.get('inspiral-head','order')))
        self.set_bank(bankFileName)


      # the output_file_name is required by the child job (plotSNRCHISQNode)
      self.output_file_name = inspJob.outputPath + self.inputIfo + "-INSPIRAL_" + self.__ifotag + "_" + self.__usertag + "-" + self.__start + "-" + str(int(self.__end)-int(self.__start)) + extension

      self.set_id(self.inputIfo + "-INSPIRAL_" + self.__ifotag + "_" + self.__usertag + "-" + self.__start + "-" + str(int(self.__end)-int(self.__start)))

      self.outputCache = self.inputIfo + ' ' + 'INSPIRAL' + ' ' + str(self.__start) + ' ' + str(int(self.__end)-int(self.__start)) + ' ' + self.output_file_name  + '\n' + self.inputIfo + ' ' + 'INSPIRAL-FRAME' + ' ' + str(self.__start) + ' ' + str(int(self.__end)-int(self.__start)) + ' ' + self.output_file_name.replace(extension,".gwf") + '\n'

      self.setupNodeWeb(inspJob,False,None,None,None,dag.cache)
      self.add_var_opt("output-path",inspJob.outputPath)

      if type == 'plot':
        if opts.inspiral:
          dag.addNode(self,'inspiral')
          self.validate()
        else: self.invalidate()
      if type == 'head':
        if opts.inspiral_head:
          dag.addNode(self,'inspiral-head')
          self.validate()
        else: self.invalidate()

    except:
      self.invalidate()
      print "couldn't add inspiral job for " + self.inputIfo + "@ "+ str(trig.gpsTime[ifo])


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
  def __init__(self,job,ifo,fileName,trig,page,dag,inspiralNode,opts,ifoString=None):
    """
    job = A CondorDAGJob that can run an instance of plotSNRCHISQ followup.
    """
    if ifoString:
      time = trig.gpsTime[ifoString]
    else:
      time = trig.gpsTime[ifo]
    self.friendlyName = 'Plot SNR/CHISQ/PSD'
    try:
      pipeline.CondorDAGNode.__init__(self,job)
      self.output_file_name = ""
      self.add_var_opt("frame-file",fileName.replace(".xml",".gwf"))
      self.add_var_opt("gps",time)
      self.add_var_opt("inspiral-xml-file",fileName)
      if ifoString:
        self.add_var_opt("user-tag",ifo+"_"+ifoString+'tmplt_'+str(trig.eventID))
        self.id = job.name + '-' + ifo + '-' + ifoString + 'tmplt' + '-' + str(trig.statValue) + '_' + str(trig.eventID)
      else:
        self.add_var_opt("user-tag",ifo+'_'+str(trig.eventID))
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
      print "couldn't add plot job for " + str(ifo) + "@ "+ str(time)

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

    nameList = name.split('-')[1:len(name.split('-'))]
    shortName = ''
    for word in nameList:
      shortName = shortName + word + '-'

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
        self.add_var_opt('dt-min',cp.getfloat('analyse-qscan',shortName + 'dt-min'))
        self.add_var_opt('dt-max',cp.getfloat('analyse-qscan',shortName + 'dt-max'))
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

##############################################################################
# class for h1h2 qevent jobs

class h1h2QeventJob(pipeline.CondorDAGJob, webTheJob):
  """
  A h1h2 qevent job
  """
  def __init__(self, opts, cp):
    """
    """
    self.name = 'h1h2QeventJob'
    self.__executable = string.strip(cp.get('condor','qevent'))
    self.__universe = "vanilla"
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)    
    self.setupJobWeb(self.name)
    self.setup_cachecat()

  def setup_cachecat(self):
    # create a shell script to cat all the required cache files
    cat_script = open(self.name + '/cachecat.sh','w')
    cat_script.write("""#!/bin/bash
    cat ${1} ${2} > ${3}
    """)
    cat_script.close()
    os.chmod(self.name + '/cachecat.sh',0755)

#############################################################################
# class for h1h2 qevent Node

class h1h2QeventNode(pipeline.CondorDAGNode,webTheNode):
  """
  Runs an instance of a qscan job
  """
  def __init__(self,job,dNode,times,ifoList,name,cp,opts,dag,qeventCommand):
    """
    job = A CondorDAGJob that can run an instance of H1H2 qevent.
    """

    ifoString = ''
    for ifo in ifoList:
      ifoString = ifoString + ifo

    page = string.strip(cp.get('output','page'))
    self.friendlyName = name
    self.id = ifoString + '-' + name + '-' + str(times[ifoList[0]])

    pipeline.CondorDAGNode.__init__(self,job)
    self.add_var_arg(repr(times[ifoList[0]]))
    eventDuration = string.strip(cp.get(name, 'duration'))
    self.add_var_arg(eventDuration)
    qeventConfig = string.strip(cp.get(name, ifoString + '-config-file'))
    self.add_file_arg(qeventConfig)

    cache_type_temp = dNode[ifoList[0]].outputFileName.split('-')[1]
    cache_type = cache_type_temp[3:len(cache_type_temp)]
    cache_start = []
    cache_end = []
    for ifo in ifoList:
      cache_temp = dNode[ifo].outputFileName.split('.')[0]
      cache_start.append(cache_temp.split('-')[2])
      cache_end.append(cache_temp.split('-')[-1])
    cache_start_time = max(cache_start)

    qeventcache = job.name + '/' + ifoString + '_' + cache_type + '-' + \
    str(max(cache_start)) + '-' + str(min(cache_end)) + '.qcache'

    self.add_file_arg(qeventcache)

    output = string.strip(cp.get(name, ifoString + '-output'))
    self.add_var_arg(output)
    self.set_pre_script(job.name + "/cachecat.sh %s %s %s" \
    %(dNode[ifoList[0]].outputFileName, dNode[ifoList[1]].outputFileName, \
    qeventcache))

    #get the absolute output path whatever the path might be in the ini file
    currentPath = os.path.abspath('.')
    try:
      os.chdir(output)
      absoutput = os.path.abspath('.')
      os.chdir(currentPath)
    except:
      print >> sys.stderr, 'invalid path for qevent output in the ini file'
      sys.exit(1)
    self.outputName = absoutput + '/' + repr(times[ifoList[0]]) # redirect output name

    #prepare the string for the output cache
    self.outputCache = ifoString + ' ' + name + ' ' + repr(times[ifoList[0]]) + ' ' + self.outputName + '\n'

    #try to extract web output from the ini file, else ignore it
    try:
      self.setupNodeWeb(job,False,dag.webPage.lastSection.lastSub,page,string.strip(cp.get(name,ifoString+'-web'))+repr(times[ifoList[0]]),dag.cache)
    except:
      self.setupNodeWeb(job,False,None,None,None,dag.cache)


    for ifo in ifoList:
      if dNode[ifo].validNode: self.add_parent(dNode[ifo])
      else: pass

    if eval('opts.' + qeventCommand):
      dag.addNode(self,self.friendlyName)
      self.validNode = True
    else: self.validNode = False


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

class IFOstatus_checkJob(pipeline.CondorDAGJob, webTheJob):
  """
  A followup job for downloading summary plots
  """
  def __init__(self, options, cp, tag_base='IFOSTATUS'):
    self.__name__ = 'IFOstatus_checkJob'
    self.__executable = string.strip(cp.get('condor','IFOstatus_check'))
    self.__universe = "local"
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    self.add_condor_cmd('getenv','True')
    self.setupJobWeb(self.__name__,tag_base)
    
class IFOstatus_checkNode(pipeline.CondorDAGNode,webTheNode):
  """
  Runs an instance of a FrCheck followup job
  """
  def __init__(self, IFOstatus_checkJob, ifo, trig, cp,opts,dag):

    self.friendlyName = 'IFO status summary plots'
    pipeline.CondorDAGNode.__init__(self,IFOstatus_checkJob)
    self.add_var_opt("ifo", ifo)
    self.add_var_opt("gps-time", trig.gpsTime[ifo])
    self.id = IFOstatus_checkJob.name + '-' + str(ifo) + '-' + str(trig.statValue) + '_' + str(trig.eventID)
    self.setupNodeWeb(IFOstatus_checkJob,True, dag.webPage.lastSection.lastSub,dag.page,None,dag.cache)
    if opts.ifo_status_check:
      dag.addNode(self,self.friendlyName)
      self.validate()
    else: self.invalidate()



###############################################################################
# MCMC - this is currently only experimental !!!!
# Use for testing only 
###############################################################################

class mcmcJob(pipeline.CondorDAGJob, webTheJob):
  """
  An mcmc job
  """
  def __init__(self, options, cp, tag_base='MCMC'):
    """
    """
    self.__name__ = 'mcmcJob'
    self.__executable = string.strip(cp.get('condor','mcmc'))
    self.__universe = "standard"
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    #self.add_condor_cmd('getenv','True')
    self.setupJobWeb(self.__name__,tag_base)


class mcmcNode(pipeline.CondorDAGNode,webTheNode):
  """
  Runs an instance of an mcmc followup job
  """
  def __init__(self, mcmcJob, procParams, ifo, trig, cp,opts,dag):

    time_margin = string.strip(cp.get('mcmc','prior-coal-time-marg'))

    self.friendlyName = 'MCMC followup'
    pipeline.CondorDAGNode.__init__(self,mcmcJob)

    for row in procParams:
      param = row.param.strip("-")
      value = row.value
      if param == 'frame-cache': cacheFile = value
      if param == 'low-frequency-cutoff': lowCut = value
      if param == 'sample-rate': highCut = str(float(value)/2 - 1.0)
      if param == 'channel-name': channel = value
      if param == 'gps-end-time': chunk_end = value
      if param == 'gps-start-time': chunk_start = value

    # add the arguments that I know now...
    self.add_var_opt("low-cut", lowCut)
    self.add_var_opt("high-cut", str(min(1800,float(highCut))))
    self.add_var_opt("channel-name", channel) # must be STRAIN ?
    # THIS COALESCENCE TIME INFO IS AD HOC AND NEEDS TO BE IMPROVED
    self.add_var_opt("prior-coal-time-mean", str(trig.gpsTime[ifo]))
    self.add_var_opt("prior-coal-time-marg",time_margin)
    self.add_var_opt("random-seed-one", str(trig.gpsTime[ifo]).split('.')[0][5:9])
    self.add_var_opt("random-seed-two", str(trig.gpsTime[ifo]).split('.')[1])
    mass1 = getattr(trig.coincs,ifo).mass1
    mass2 = getattr(trig.coincs,ifo).mass2
    dist = getattr(trig.coincs,ifo).eff_distance
    snr = getattr(trig.coincs,ifo).snr
    duration = getattr(trig.coincs,ifo).template_duration
    # THE DIST LIMITS ARE AD HOC THIS NEEDS TO BE FIXED
    # THESE CAN HAVE LARGE SYSTEMATIC ERRORS THAT WILL NOT
    # BE CAPTURED BY THIS RIGHT??? I am using two sigma instead of 
    # 1.65 which would be right for the 10/90 !!
    #dist10 = 1.0/math.sqrt((snr*snr+2.0*2.0)/snr/snr)*dist
    #dist90 = 1.0/math.sqrt((snr*snr-2.0*2.0)/snr/snr)*dist
    #dist90 = 45.6*math.sqrt(mass1*mass2)/(mass1+mass2)**(1./6.)
    dist90 = 56.5
    dist10 = dist90 + 5. 
    # THE MASS LIMITS ARE AD HOC THIS NEEDS TO BE FIXED
    #if mass1 < mass2:
    #  self.add_var_opt("prior-lower-mass", str(0.9) )
    #  self.add_var_opt("prior-upper-mass", str(10.0) )
    #else:
    #  self.add_var_opt("prior-lower-mass", str(0.9) )
    #  self.add_var_opt("prior-upper-mass", str(10.0) )
    self.add_var_opt("prior-lower-mass", str(0.5) )
    self.add_var_opt("prior-upper-mass", str(34.5) )
    self.add_var_opt("prior-distance-10", str(dist10))
    self.add_var_opt("prior-distance-90", str(dist90))
    self.add_var_opt("before-coal-time", str(duration*1.5))
    #self.add_var_opt("before-coal-time", str(28))
    
    ########################################################################
    # GET THE FRAME FILE INFO - THIS NEEDS TO BE CHANGED !!!
    # THE MCMC CODE SHOULD TAKE THE SAME CACHE FILE AS LALAPPS_INSPIRAL AND
    # COMPUTE THE SAME PSD ETC.  CURRENTLY THE WAY IT HANDLES FRAMES IS A BIT
    # AD HOC.  THIS IS BOUND TO FAIL ON OCCASION SINCE IT ASSUMES ALL FRAMES
    # HAVE THE SPECIFIED DURATION!!!
    ########################################################################
    frameCache = lal.Cache().fromfile(open(cacheFile,'r'))  
    trigSegment = segments.segment(trig.gpsTime[ifo],trig.gpsTime[ifo])
    goodCache = frameCache.sieve(ifo[0],ifo[0], trigSegment)
    frame = goodCache[0].path()
    frameFilePath = str.join('/',frame.split('/')[0:-1])
    frameFile = frame.split('/')[-1]
    frameFilePrefix = str.join('-',frameFile.split('-')[0:2]+[''])
    frameFileSuffix = str.join('-',['',frameFile.split('-')[-1]])
    gpsStartTime = frameFile.split('-')[-2]
    frameFileSize = str(goodCache[0]).split(' ')[3]
    self.add_var_opt("frame-file-path",frameFilePath)
    self.add_var_opt("frame-file-prefix",frameFilePrefix)
    self.add_var_opt("frame-file-suffix",frameFileSuffix)
    self.add_var_opt("frame-file-size",frameFileSize)
    self.add_var_opt("gps-start-time",gpsStartTime)

    if (float(chunk_end) - trig.gpsTime[ifo]) > 512.:
      psdSegment = segments.segment(trig.gpsTime[ifo]+256.,trig.gpsTime[ifo]+256.)
    else:
      psdSegment = segments.segment(trig.gpsTime[ifo]-384.,trig.gpsTime[ifo]-384.)
    psdCache = frameCache.sieve(ifo[0],ifo[0],psdSegment)
    psdFrame = psdCache[0].path()
    psdFilePath = str.join('/',psdFrame.split('/')[0:-1])
    psdFile = psdFrame.split('/')[-1]
    psdFilePrefix = str.join('-',psdFile.split('-')[0:2]+[''])
    psdFileSuffix = str.join('-',['',psdFile.split('-')[-1]])
    psdStartTime = psdFile.split('-')[-2]
    psdFileSize = str(psdCache[0]).split(' ')[3]
    #self.add_var_opt("psd-file-path",psdFilePath)
    #self.add_var_opt("psd-file-prefix",psdFilePrefix)
    #self.add_var_opt("psd-file-suffix",psdFileSuffix)
    #self.add_var_opt("psd-file-size",psdFileSize)
    self.add_var_opt("psd-start-time",psdStartTime)

    self.id = mcmcJob.name + '-' + ifo + '-' + str(trig.statValue) + '_' + str(trig.eventID)
    self.setupNodeWeb(mcmcJob)
    self.add_var_opt("output-file-name",mcmcJob.name+'/'+self.id+'.txt') 
    if opts.mcmc:
      dag.addNode(self,self.friendlyName)
      self.validate()
    else: self.invalidate()
