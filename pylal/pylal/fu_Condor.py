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

###### GENERIC CLASSES TO HELP WEBIFY A CONDOR DAG OUTPUT ####################
##############################################################################

class webTheJob:
  
  def __init__(self):
    pass

  def setupJobWeb(self, name, tag_base=None):
    self.name = name
    try: 
       os.mkdir(name)
       os.mkdir(name+'/logs')
    except: pass
    self.tag_base = tag_base
    self.add_condor_cmd('environment',"KMP_LIBRARY=serial;MKL_SERIAL=yes")
    self.set_sub_file(name+'.sub')
    self.outputPath = os.getcwd() + '/' + name + '/'
    self.set_stdout_file(self.outputPath+'/logs/'+name+'-$(macroid).out')
    self.set_stderr_file(self.outputPath+'/logs/'+name+'-$(macroid).err')

class webTheNode:

  def __init__(self):
    pass

  def setupNodeWeb(self, job, id, passItAlong=True, page=None):
    self.add_macro("macroid", id)
    self.outputFileName = job.outputPath + job.name + '-'+ id
    if passItAlong:
      self.add_var_opt("output-web-file",self.outputFileName)
      self.add_var_opt("output-path",job.outputPath)
      self.add_var_opt("page", page)

class webTheDAG:
 
  def __init__(self):
    pass    

###### WRAPPER FOR CONDOR DAG - TO MAKE THE FOLLOWUP DAG WEBIFIABLE ###########
###############################################################################

class followUpDAG(pipeline.CondorDAG):

  def init(self, config_file, log_path):
    self.basename = re.sub(r'\.ini',r'', config_file) 
    tempfile.tempdir = log_path
    tempfile.template = self.basename + '.dag.log.'
    logfile = tempfile.mktemp()
    fh = open( logfile, "w" )
    fh.close()
    pipeline.CondorDAG.__init__(self,logfile)
    self.set_dag_file(self.basename)

    
    
  
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
  def __init__(self,job,time,fileName,trig,page,plotFlag):
    """
    job = A CondorDAGJob that can run an instance of plotSNRCHISQ followup.
    """
    pipeline.CondorDAGNode.__init__(self,job)
    self.output_file_name = ""
    self.add_var_opt("frame-file",fileName.replace(".xml",".gwf"))
    self.add_var_opt("gps",time)
    self.add_var_opt("inspiral-xml-file",fileName)
    self.id = str(trig.statValue) + '_' + str(trig.eventID)
    self.setupNodeWeb(job,self.id,page)


############### QSCAN CLASSES #################################################
###############################################################################

class qscanDataFindJob(pipeline.LSCDataFindJob):
  
  def __init__(self, cache_dir, log_dir, config_file, source):
    pipeline.LSCDataFindJob.__init__(self, cache_dir, log_dir, config_file)
    if source == 'futrig':
      self.setup_cacheconv()

  def setup_cacheconv():
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

class qscanDataFindNode(pipeline.LSCDataFindNode):
 
  def __init__(self, job, source, type, cp, time, ifo):
    pipeline.LSCDataFindNode.__init__(self,job)
    if source == 'futrig':
      self.outputFileName = self.setup_fu_trig(cp, time, ifo, type)

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
# qscan class for background qscan Node

class qscanBgNode(pipeline.CondorDAGNode,webTheNode):
  """
  Runs an instance of a qscan job
  """
  def __init__(self,job,time,cp,qcache,ifo,name):
    """
    job = A CondorDAGJob that can run an instance of qscan.
    """
    pipeline.CondorDAGNode.__init__(self,job)
    self.add_var_arg(repr(time))
    qscanConfig = string.strip(cp.get(name[1], ifo + 'config-file'))
    self.add_file_arg(qscanConfig)
    self.add_file_arg(qcache)
    self.add_var_arg(string.strip(cp.get(name[1], ifo + 'output')))
    self.id = ifo + repr(time)
    #self.setupNodeWeb(job,self.id,page,False)
    self.setupNodeWeb(job,self.id,False)

##############################################################################
# qscan class for qscan Node

class qscanNode(pipeline.CondorDAGNode,webTheNode):
  """
  Runs an instance of a qscan job
  """
  def __init__(self,job,time,cp,trig,qcache,ifo,name,qFlag):
    """
    job = A CondorDAGJob that can run an instance of qscan.
    """
    pipeline.CondorDAGNode.__init__(self,job)
    self.add_var_arg(repr(time))
    qscanConfig = string.strip(cp.get(name[1], ifo + 'config-file'))
    self.add_file_arg(qscanConfig)
    self.add_file_arg(qcache)
    output = string.strip(cp.get(name[1], ifo + 'output'))
    self.add_var_arg(output)
    #self.setupNodeWeb(job,self.id,page,False)
    self.id = ifo + repr(time)
    self.setupNodeWeb(job,self.id,False)
    self.outputFileName = output + '/' + repr(time) # redirect output name


################ TRIG BANK FROM SIRE FILE CONDOR DAG JOB ######################
def inspiralFollowup(inspJob, procParams, triggerTime, bankFile, injFile):
  """
  Set up an instance of an inspiral followup job
  """
  insp = inspiral.InspiralNode(inspJob)
  insp.add_var_opt("write-snrsq","")
  insp.add_var_opt("write-chisq","")
  insp.add_var_opt("write-spectrum","")
  insp.set_bank(bankFile)

  insp.set_user_tag("FOLLOWUP_" + str(triggerTime))
  insp.__usertag = "FOLLOWUP_" + str(triggerTime)

  insp.set_trig_start( int(triggerTime) - 1)
  insp.set_trig_end( int(triggerTime) + 1 )
  if injFile: insp.set_injections( injFile )

  skipParams = ['minimal-match', 'bank-file', 'user-tag', 'injection-file', 'trig-start-time', 'trig-end-time']

  for row in procParams:
    param = row.param.strip("-")
    value = row.value

    if param in skipParams:
      continue

    insp.add_var_opt(param,value)

    if param == 'gps-end-time':
      #insp.set_end(value)
      insp.__end = value

    if param == 'gps-start-time':
      #insp.set_start(value)
      insp.__start = value

    if param == 'ifo-tag':
      #insp.set_ifo_tag(value)
      insp.__ifotag = value

    if param == 'channel-name':
      #insp.set_ifo(value[0:2])
      insp.inputIfo = value[0:2]

  # the output_file_name is required by the child job (plotSNRCHISQNode)
  insp.output_file_name = insp.inputIfo + "-INSPIRAL_" + insp.__ifotag + "_" + insp.__usertag + "-" + insp.__start + "-" + str(int(insp.__end)-int(insp.__start)) + ".xml"

  return insp




