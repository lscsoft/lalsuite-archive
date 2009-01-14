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
from glue.ligolw import lsctables

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
    # The list remote_nodes will contain the list of nodes run remotely 
    # (such as V1 qscans)
    self.remote_nodes = []

########################################################
#### Methods common to several followup classes ########
########################################################

def checkHipeCachePath(cp):
  try:
    if len(string.strip(cp.get('followup-hipe-cache','hipe-cache-path'))) > 0:
      hipeCachePath = string.strip(cp.get('followup-hipe-cache','hipe-cache-path'))
    else:
      hipeCachePath = None
    return(hipeCachePath)
  except:
    print >> sys.stderr, "ERROR: failure in checkHipeCachePath()"
    return None

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
  
  def __init__(self, inspJob, procParams, ifo, trig, cp,opts,dag, datafindCache, d_node, datafindCommand, type='plot', sngl_table = None):

    try:
      self.output_file_name = ""
      #inspiral.InspiralNode.__init__(self, inspJob) 
      # the use of this class would require some reorganisation in fu_Condor.py
      # and webCondor.py in order to set up the jobs following the same scheme
      # as the way it is done for the Inspiral pipeline...
      pipeline.CondorDAGNode.__init__(self,inspJob)
      injFile = self.checkInjections(cp)      
      hipeCache = checkHipeCachePath(cp)

      if type == "plot" or type == "notrig" or type == "coh":
        bankFile = 'trigTemplateBank/' + ifo + '-TRIGBANK_FOLLOWUP_' + type + str(trig.eventID) + '.xml.gz'
        self.add_var_opt("write-snrsq","")
        self.add_var_opt("write-chisq","")
        self.add_var_opt("write-spectrum","")
        self.add_var_opt("write-cdata","")
        #self.add_var_opt("verbose","")
        self.set_bank(bankFile)
        # Here we define the trig-start-time and the trig-end-time;
        # The difference between these two times should be kept to 2s
        # Otherwise change the clustering window also
        hLengthAnalyzed = 1
	if type == "coh": hLengthAnalyzed = 2
        self.set_trig_start( int(trig.gpsTime[ifo]) - int(hLengthAnalyzed) )
        self.set_trig_end( int(trig.gpsTime[ifo]) + int(hLengthAnalyzed) )

      if injFile: 
        self.set_injections( injFile )

      skipParams = ['minimal-match', 'bank-file', 'user-tag', 'injection-file', 'trig-start-time', 'trig-end-time']
      if not hipeCache:
        skipParams.append('frame-cache')
        self.add_var_opt('frame-cache',datafindCache)        

      # initialize the extension of the output file. If the option 
      # write_compress is found in procParams the extension will be overwritten
      # later as .xml.gz
      extension = ".xml"

      for row in procParams:
        param = row.param.strip("-")
        value = row.value
        if param == 'bank-file':
          bankFile = value
        if type == "notrig" or type == "coh":
        # if forceTrigger is true, we loose the thresholds to
        # make sure to get a trigger
          if param == 'snr-threshold': value = "0.1"
          # rsq veto must be disabled
          if param == 'do-rsq-veto': continue
          if param == 'enable-rsq-veto': continue
          # chisq veto is disabled by loosing its threshold 
          # we still want to generate the chisq time-series
          if param == 'chisq-threshold': value = "1.0e+06"
          # Using a window of 1s for clustering will allow us to always get
          # at least one trigger
          if param == 'cluster-method': value = 'window'
          if param == 'cluster-window': continue
          pass
        if param in skipParams: continue
        self.add_var_opt(param,value)
        # The attributes _AnalysisNode__end, _AnalysisNode__start,
        # _InspiralAnalysisNode__pad_data need to be defined before calling the
        # method "writeAll" of "class webTheDAG". This method calls
        # "write_sub_files()" in pipeline.py, which itself relies on the
        # "finalize()" method of "class InspiralAnalysisNode" in inspiral.py .
        # This is where all these attributes are being used. This hack is 
        # required because "inspiral.InspiralNode.__init__(self, inspJob)"
        # currently does not work within "class followUpInspNoDE"
        if param == 'gps-end-time':
          self.__end = value
          self._AnalysisNode__end = int(value)
        if param == 'gps-start-time':
          self.__start = value
          self._AnalysisNode__start = int(value)
        if param == 'pad-data': 
          self._InspiralAnalysisNode__pad_data = int(value)
        if param == 'ifo-tag':
          self.__ifotag = value
        if param == 'channel-name': self.inputIfo = value[0:2]
        if param == 'write-compress':
          extension = '.xml.gz'

      if type == "notrig" or type == "coh":
        self.add_var_opt('cluster-window',str(hLengthAnalyzed))
        self.add_var_opt('disable-rsq-veto',' ')

      # add the arguments that have been specified in the section 
      # [inspiral-extra] of the ini file (intended for 12-18 month analysis)
      if cp.has_section("followup-inspiral-extra"):
        for (name,value) in cp.items("followup-inspiral-extra"):
          self.add_var_opt(name,value)

      if not ifo == self.inputIfo:
        second_user_tag = "_" + ifo + "tmplt"
      else:
        second_user_tag = ""
      self.set_user_tag("FOLLOWUP_" + str(trig.eventID) + second_user_tag)
      self.__usertag = "FOLLOWUP_" + str(trig.eventID) + second_user_tag


      # THIS IS A HACK FOR NOW, THERE IS PROBABLY A BETTER WAY TO DO THIS
      if (type == 'head'): 
        subBankSize = string.strip(cp.get('followup-inspiral-head','bank-veto-subbank-size'))
        if opts.inspiral_head:
          bankFileName = fu_utils.generateBankVetoBank(trig, ifo, str(trig.gpsTime[ifo]), sngl_table[ifo],int(subBankSize),'BankVetoBank')
        else: bankFileName = 'none'      
        self.add_var_opt("bank-veto-subbank-size", string.strip(cp.get('followup-inspiral-head','bank-veto-subbank-size')))
        self.add_var_opt("order", string.strip(cp.get('followup-inspiral-head','order')))
        self.set_bank(bankFileName)


      # the output_file_name is required by the child job (plotSNRCHISQNode)
      self.output_file_name = inspJob.outputPath + self.inputIfo + "-INSPIRAL_" + self.__ifotag + "_" + self.__usertag + "-" + self.__start + "-" + str(int(self.__end)-int(self.__start)) + extension

      self.set_id(self.inputIfo + "-INSPIRAL_" + self.__ifotag + "_" + self.__usertag + "-" + self.__start + "-" + str(int(self.__end)-int(self.__start)))

      self.outputCache = self.inputIfo + ' ' + 'INSPIRAL' + ' ' + str(self.__start) + ' ' + str(int(self.__end)-int(self.__start)) + ' ' + self.output_file_name  + '\n' + self.inputIfo + ' ' + 'INSPIRAL-FRAME' + ' ' + str(self.__start) + ' ' + str(int(self.__end)-int(self.__start)) + ' ' + self.output_file_name.replace(extension,".gwf") + '\n'

      self.setupNodeWeb(inspJob,False,None,None,None,dag.cache)
      self.add_var_opt("output-path",inspJob.outputPath)

      if not opts.disable_dag_categories:
        self.set_category(inspJob.name.lower())

      try:
        if d_node.validNode and eval('opts.' + datafindCommand):
          self.add_parent(d_node)
      except: 
        print >> sys.stderr, "Didn't find a datafind job, I'll assume I don't need it"

      if type == "plot" or type == "notrig":
        if opts.inspiral:
          dag.addNode(self,'inspiral')
          self.validate()
        else: self.invalidate()

      if type == 'head':
        if opts.inspiral_head:
          dag.addNode(self,'inspiral-head')
          self.validate()
        else: self.invalidate()

      if type == 'coh':
        if opts.coh_inspiral:
          dag.addNode(self,'coh-inspiral')
          self.validate()
        else: self.invalidate()

    except:
      try:
        print "couldn't add inspiral job for " + self.inputIfo + "@ "+ str(trig.gpsTime[ifo])
        # if self.inputIfo does not exist (happens when inspiral cache and xml files not available), then use ifo in the string.
      except:
        print "couldn't add inspiral job for " + ifo + "@ "+ str(trig.gpsTime[ifo])

  def checkInjections(self,cp):
    try:
      if len(string.strip(cp.get('followup-triggers','injection-file'))) > 0:
        injectionFile = string.strip(cp.get('followup-triggers','injection-file'))
      else:
        injectionFile = None
      return(injectionFile)
    except:
      print >> sys.stderr, "ERROR: failure in followUpInspNode.checkInjections()"
      return None

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
      self.add_var_opt("frame-file",fileName.replace(".xml",".gwf").strip(".gz"))
      self.add_var_opt("gps",time)
      self.add_var_opt("inspiral-xml-file",fileName)
      if ifoString:
        self.add_var_opt("user-tag",ifo+"_"+ifoString+'tmplt_'+str(trig.eventID))
        self.id = job.name + '-' + ifo + '-' + ifoString + 'tmplt' + '-' + str(trig.statValue) + '_' + str(trig.eventID)
      else:
        self.add_var_opt("user-tag",ifo+'_'+str(trig.eventID))
        self.id = job.name + '-' + ifo + '-' + str(trig.statValue) + '_' + str(trig.eventID)
      self.setupNodeWeb(job,True, dag.webPage.lastSection.lastSub,page,None,dag.cache)

      if not opts.disable_dag_categories:
        self.set_category(job.name.lower())

      #try: 
      if inspiralNode.validNode: self.add_parent(inspiralNode)
      #except: pass
      if opts.plots:
        dag.addNode(self,self.friendlyName)
        self.validate()
      else: self.invalidate()
    except: 
      self.invalidate()
      print "couldn't add plot job for " + str(ifo) + "@ "+ str(time)

##############################################################################
# job class for producing the skymap

class lalapps_skyMapJob(pipeline.CondorDAGJob,webTheJob):
  """
  Generates sky map data
  """
  def __init__(self, options, cp, tag_base='SKY_MAP'):
    """
    """
    self.__name__ = 'lalapps_skyMapJob'
    self.__executable = string.strip(cp.get('condor','lalapps_skymap'))
    self.__universe = "standard"
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    self.add_condor_cmd('getenv','True')
    self.setupJobWeb(self.__name__,tag_base)

##############################################################################
# job class for producing the skymap

class pylal_skyPlotJob(pipeline.CondorDAGJob,webTheJob):
  """
  Plots the sky map output of lalapps_skymap
  """
  def __init__(self, options, cp, tag_base='SKY_PLOT'):
    """
    """
    self.__name__ = 'pylal_skyPlotJob'
    self.__executable = string.strip(cp.get('condor','pylal_skyPlotJob'))
    self.__universe = "vanilla"
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    self.add_condor_cmd('getenv','True')
    self.setupJobWeb(self.__name__,tag_base)


##############################################################################
# job class for producing the skymap

class lalapps_skyMapNode(pipeline.CondorDAGNode,webTheNode):
  """
  A C code for computing the sky map
  An example command line is:

lalapps_skymap --h1-frame-file H1-INSPIRAL_SECOND_H1H2L1V1_FOLLOWUP_866088314000001908-866088022-2048.gwf --l1-frame-file L1-INSPIRAL_SECOND_H1H2L1V1_FOLLOWUP_866088314000001908-866088022-2048.gwf --v1-frame-file V1-INSPIRAL_SECOND_H1H2L1V1_FOLLOWUP_866088314000001908-866088205-2048.gwf --event-id 866088314000001908 --ra-res 512 --dec-res 256 --h1-xml-file H1-INSPIRAL_SECOND_H1H2L1V1_FOLLOWUP_866088314000001908-866088022-2048.xml.gz --l1-xml-file L1-INSPIRAL_SECOND_H1H2L1V1_FOLLOWUP_866088314000001908-866088022-2048.xml.gz --v1-xml-file V1-INSPIRAL_SECOND_H1H2L1V1_FOLLOWUP_866088314000001908-866088205-2048.xml.gz --output-file chad.txt
  """
  def __init__(self,job,trig,opts):
    self.ifo_list = ["H1","L1","V1"]
    #self.already_added_ifo_list = []
    self.ra_res = 512
    self.dec_res = 256
    pipeline.CondorDAGNode.__init__(self,job)
    self.friendlyName = 'Produce sky map of event'    
    self.id = job.name + '-skymap-' + str(trig.statValue) + '_' + str(trig.eventID)
    self.setupNodeWeb(job)
    # required by pylal_skyPlotNode
    self.output_file_name = job.outputPath + self.id+".txt"
    self.add_var_opt("output-file",self.output_file_name)
    self.add_var_opt("ra-res",self.ra_res)
    self.add_var_opt("dec-res",self.dec_res)
    self.add_var_opt("event-id",trig.eventID)
    self.add_var_opt("h1-frame-file","none");
    self.add_var_opt("h1-xml-file","none");
    self.add_var_opt("l1-frame-file","none");
    self.add_var_opt("l1-xml-file","none");
    self.add_var_opt("v1-frame-file","none");
    self.add_var_opt("v1-xml-file","none");

    if not opts.disable_dag_categories:
      self.set_category(job.name.lower())

  def append_insp_node(self,inspNode,ifo):
    if ifo in self.ifo_list:
      fileName = str(inspNode.output_file_name)
      self.add_var_opt(ifo.lower()+"-frame-file",str(fileName.replace(".xml",".gwf").strip(".gz")))
      self.add_var_opt(ifo.lower()+"-xml-file",str(fileName))
      if inspNode.validNode: self.add_parent(inspNode)
      
    else: print >> sys.stderr, "WARNING: Already added " + ifo


  def add_node_to_dag(self,dag,opts,trig):
    if opts.sky_map:
      dag.addNode(self,self.friendlyName)
      self.validate()
    else: 
      self.invalidate()
      print "couldn't add sky map job for " + str(trig.eventID)



##############################################################################
# job class for producing the skymap

class pylal_skyPlotNode(pipeline.CondorDAGNode,webTheNode):
  """
  A python code for plotting the sky map
  An example command line is

  /pylal_plot_inspiral_skymap --event-id 866088314000001908 --ra-res 512 --dec-res 256 --output-path . --page-rel-path . --output-web-file test.html --page . --injection-right-ascension 0 --injection-declination 0 --map-data-file chad.txt 
  """
  def __init__(self,job,trig,skyMapNode,dag,page,opts):
    # Always initialize the CondorDAGNode
    pipeline.CondorDAGNode.__init__(self,job)
    
    self.friendlyName = 'Produce a plot of the sky map of an event'
    
    self.id = job.name + '-skymap-plot' + str(trig.statValue) + '_' + str(trig.eventID)
    # This node outputs pretty pictures, so we need to tell the setupNodeWeb()
    # method where to put these things.  We'll put it in the last section 
    # not the last subsection since this is an "event" plot not a single ifo
    # trigger plot
    self.setupNodeWeb(job,True, dag.webPage.lastSection,page,None,dag.cache)
    # this is the output of the skyMapNode.  It contains the data to make a
    # sky map plot.  (RA,DEC,Probability)
    # an example sky map plotting command line is:
    #

    self.add_var_opt("map-data-file",skyMapNode.output_file_name)
    self.add_var_opt("event-id",str(trig.eventID))
    self.add_var_opt("ra-res",str(skyMapNode.ra_res))
    self.add_var_opt("dec-res",str(skyMapNode.dec_res))
    # if this is a software injection pass along the information to the
    # plotting code so that it can make a mark where the injection should have
    # been :)
    if trig.is_found():
      inj_ra = trig.coincs.sim.longitude
      inj_dec = trig.coincs.sim.latitude
      self.add_var_opt("injection-right-ascension",str(inj_ra))
      self.add_var_opt("injection-declination",str(inj_dec))

    if not opts.disable_dag_categories:
      self.set_category(job.name.lower())

    try:
      if skyMapNode.validNode: self.add_parent(skyMapNode)
    except: pass
    if opts.sky_map_plot:
      dag.addNode(self,self.friendlyName)
      self.validate()
    else: self.invalidate()
 

############### DATAFIND CLASSES ##############################################
###############################################################################

class followupDataFindJob(pipeline.LSCDataFindJob,webTheJob):

  def __init__(self, config_file, source):

    if source == 'futrig':
      self.name = 'qscanDataFindJob'
    if source == 'inspiral':
      self.name = 'inspiralDataFindJob'

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


class followupDataFindNode(pipeline.LSCDataFindNode,webTheNode):
 
  def __init__(self, job, source, type, cp, time, ifo, opts, dag, datafindCommand, procParams=None):
    try:
      self.outputFileName = ""
      pipeline.LSCDataFindNode.__init__(self,job)
      self.id = str(ifo) + '-' + repr(time) + '-' + str(type)
      self.setupNodeWeb(job,False,None,None,None,dag.cache)
      if source == 'futrig':
        self.outputFileName = self.setup_fu_trig(job, cp, time, ifo, type)
        nodeName = "qscan data find"
      if source == 'inspiral':
        self.outputFileName = self.setup_inspiral(cp,ifo,type,procParams)
        nodeName = "inspiral data find"

      if not opts.disable_dag_categories:
        self.set_category(job.name.lower())

      # if the selected "ifo" needs to be done remotely (this the case for 
      # Virgo qscan datafind) do not add the node to the dag
      if eval('opts.' + datafindCommand) and \
        not( cp.has_option("followup-"+type,"remote-ifo") and \
        cp.get("followup-"+type,"remote-ifo")==ifo ):
          dag.addNode(self,nodeName)
          self.validNode = True
      else: self.validNode = False
    except:
      self.validNode = False
      print >> sys.stderr, "could not set up the datafind jobs for " + type

  def setup_inspiral(self,cp,ifo,type,procParams):
    for row in procParams:
      param = row.param.strip("-")
      value = row.value
      if param == 'gps-start-time': startTime = value
      if param == 'gps-end-time': endTime = value
      if param == 'pad-data': paddataTime = value
    self.set_observatory(ifo[0])
    self.set_start(int(startTime) - int(paddataTime))
    self.set_end(int(endTime) + int(paddataTime))
    self.set_type(cp.get("followup-"+type,ifo + '_type'))
    lalCache = self.get_output()
    return(lalCache)

  def setup_fu_trig(self, job, cp, time, ifo, type):
    # 1s is substracted to the expected startTime to make sure the window
    # will be large enough. This is to be sure to handle the rouding to the
    # next sample done by qscan.
    self.q_time = cp.getint("followup-"+type,'search-time-range')/2
    self.set_observatory(ifo[0])
    self.set_start(int( time - self.q_time - 1))
    self.set_end(int( time + self.q_time + 1))
    if cp.has_option("followup-"+type, ifo + '_type'): 
      self.set_type( cp.get("followup-"+type, ifo + '_type' ))
    else:
      self.set_type( cp.get("followup-"+type, 'type' ))
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
    self.friendlyName = name
    self.id = ifo + '-' + name + '-' + repr(time)

    pipeline.CondorDAGNode.__init__(self,job)
    if name.split('-')[0]=='background':
      self.add_var_arg('scanlite')
    else:
      self.add_var_arg('scan')
    qscanConfig = string.strip(cp.get("followup-"+name, ifo + 'config-file'))
    self.add_var_arg("-c "+qscanConfig)
    self.add_var_arg("-f "+qcache)

    if cp.has_option("followup-"+name, ifo + 'output') and string.strip(cp.get("followup-"+name, ifo + 'output')):
      output = string.strip(cp.get("followup-"+name, ifo + 'output'))
    else:
      output = dag.publish_path + '/' + job.name + '/' + name + '/' + ifo
    if not os.access(output,os.F_OK):
      os.makedirs(output)
    else:
      if not os.access(output,os.W_OK):
        print >> sys.stderr, 'path '+output+' is not writable'
        sys.exit(1)

    self.add_var_arg("-o "+output)
    self.add_var_arg(repr(time))

    #get the absolute output path whatever the path might be in the ini file
    absoutput = os.path.abspath(output)

    self.outputName = absoutput + '/' + repr(time) # redirect output name

    #prepare the string for the output cache
    self.outputCache = ifo + ' ' + name + ' ' + repr(time) + ' ' + self.outputName + '\n'

    #extract web output from the ini file if the job is QSCAN
    if job.name == 'QSCAN':
      if cp.has_option("followup-"+name,ifo+'web') and string.strip(cp.get("followup-"+name,ifo+'web')):
        pageOverride = string.strip(cp.get("followup-"+name,ifo+'web'))+'/'+repr(time)
      else:
        pageOverride = dag.page + '/' + job.name + '/' + name + '/' + ifo + '/' + repr(time)
      self.setupNodeWeb(job,False,dag.webPage.lastSection.lastSub,dag.page,pageOverride,dag.cache)

    else:
      self.setupNodeWeb(job,False,None,None,None,dag.cache) 

    # This command will force Condor to see the qscan jobs successful even
    # they fail. This is useful when the followups are rerun on candidates 
    # already analysed, since when a qscan directory exists, the qscan job
    # will fail. By setting the post_script to true the qscan job will
    # still be reported as successful, so that an analyseQscan job can be run
    # immediately after. 
    # self.set_post_script("/bin/true")

    if not opts.disable_dag_categories:
      self.set_category(job.name.lower())

    # only add a parent if it exists
    try:
      if d_node.validNode and eval('opts.' + datafindCommand):
        self.add_parent(d_node)
    except: pass

    # if the selected "ifo" needs to be done remotely (this the case for 
    # Virgo qscans) do not add the node to the dag
    if eval('opts.' + qscanCommand):
      if not(cp.has_option("followup-"+name,"remote-ifo") and \
      cp.get("followup-"+name,"remote-ifo")==ifo):
        dag.addNode(self,self.friendlyName)
        self.validNode = True
      else:
        dag.remote_nodes.append(self)
    else: self.validNode = False
 #   except: 
 #     self.validNode = False
 #     print >> sys.stderr, "could not set up the qscan job for " + self.id


##############################################################################
# distributeQscanJob class: the job

class distributeQscanJob(pipeline.CondorDAGJob, webTheJob):
  """
  A job to distribute the results of the qscans that have been run remotely (for LV search)
  """
  def __init__(self,cp):
    self.__name__ = 'distributeQscanJob'
    self.__executable = string.strip(cp.get('condor','distribute_q'))
    self.__universe = "vanilla"
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    self.add_condor_cmd('getenv','True')
    self.setupJobWeb(self.__name__)

##############################################################################
# distributeQscanNode class: the node

class distributeQscanNode(pipeline.CondorDAGNode, webTheNode):
  """
  A node to distribute the results of the qscans that have been run remotely (for LV search)
  """
  def __init__(self,job,foregroundCache,backgroundCache,ifo,inputFile,opts,dag):

    self.friendlyName = "distributeQscanResults"

    pipeline.CondorDAGNode.__init__(self,job)
    self.add_var_opt('qscan-input-file',inputFile)
    self.add_var_opt('qscan-cache-background',backgroundCache)
    self.add_var_opt('qscan-cache-foreground',foregroundCache)
    self.add_var_opt('remote-ifo',ifo)

    typeList=""
    for type in ["qscan","seismic-qscan"]:
      typeList += type + ","
    self.add_var_opt('qscan-type-list',typeList.strip(','))

    if opts.distrib_remote_q:
      dag.addNode(self,self.friendlyName)
      self.validNode = True
    else: self.validNode = False

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
    self.friendlyName = 'analyse ' + name
    self.id = ifo + '-' + name + '-' + repr(time)

    nameList = name.split('-')[1:len(name.split('-'))]
    shortName = ''
    for word in nameList:
      shortName = shortName + word + '-'

    try:
      pipeline.CondorDAGNode.__init__(self,job)
      if cp.has_option('followup-analyse-qscan','generate-qscan-xml'):
        self.add_var_opt('generate-qscan-xml','')
      self.add_var_opt('z-threshold',cp.getfloat('followup-analyse-qscan','z-threshold'))
      if cp.has_option('followup-analyse-qscan','plot-z-distribution'):
        self.add_var_opt('plot-z-distribution','')
        self.add_var_opt('z-min',cp.getfloat('followup-analyse-qscan','z-min'))
        self.add_var_opt('z-max',cp.getfloat('followup-analyse-qscan','z-max'))
        self.add_var_opt('z-bins',cp.getfloat('followup-analyse-qscan','z-bins'))
      if cp.has_option('followup-analyse-qscan','plot-dt-distribution'):
        self.add_var_opt('plot-dt-distribution','')
        self.add_var_opt('dt-min',cp.getfloat('followup-analyse-qscan',shortName + 'dt-min'))
        self.add_var_opt('dt-max',cp.getfloat('followup-analyse-qscan',shortName + 'dt-max'))
        self.add_var_opt('dt-bins',cp.getfloat('followup-analyse-qscan','dt-bins'))
      if cp.has_option('followup-analyse-qscan','plot-z-scattered'):
        self.add_var_opt('plot-z-scattered','')
      if cp.has_option('followup-analyse-qscan','plot-z-scattered') or cp.has_option('followup-analyse-qscan','plot-dt-distribution'):
        if not ifo=='V1':
          refChannel = cp.get('followup-analyse-qscan',shortName + 'ref-channel').split(',')[0].strip()
        else:
          refChannel = cp.get('followup-analyse-qscan',shortName + 'ref-channel').split(',')[1].strip()
        self.add_var_opt('ref-channel',refChannel)
      self.add_var_opt('qscan-id',name + '_' + ifo + '_' + repr(time)) 

      self.add_var_opt('qscan-cache-foreground',foregroundCache)
      self.add_var_opt('qscan-cache-background',backgroundCache)

      self.setupNodeWeb(job,True,None,dag.page,None,dag.cache)
      # get the table of the qscan job associated to this trigger
      if not(cp.has_option("followup-"+name,"remote-ifo") and cp.get("followup-"+name,"remote-ifo")==ifo):
        for node in dag.get_nodes():
          if isinstance(node,qscanNode):
            if node.id == self.id:
              # link the analyseQscan output page to the qscan table
              node.webTable.row[0].cell[0].linebreak()
              node.webTable.row[0].cell[0].link(self.webLink,"qscan background vs qscan foreground")
              break
      # if remote-ifo is analysed, find the associated qscan jobs in dag.remote_nodes
      else:
        for node in dag.remote_nodes:
          if isinstance(node,qscanNode):
            if node.id == self.id:
              node.webTable.row[0].cell[0].linebreak()
              node.webTable.row[0].cell[0].link(self.webLink,"qscan background vs qscan foreground")
              break

      if not opts.disable_dag_categories:
        self.set_category(job.name.lower())
      
      # add the parents to this node
      for node in dag.get_nodes():
        # if node distributeQscanNode is valid and remote if is analysed,
        # add distributeQscanNode as parent
        if isinstance(node,distributeQscanNode):
          if cp.has_option("followup-"+name,"remote-ifo") and cp.get("followup-"+name,"remote-ifo")==ifo:
            if node.validNode:
              self.add_parent(node)
        # add all qscan nodes of the same type as parents
        if isinstance(node,qscanNode): 
          if node.validNode:
            if (node.friendlyName == name or \
            node.friendlyName.replace('background','foreground') == name) \
            and node.id.split('-')[0] == ifo:
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
    self.__executable = string.strip(cp.get('condor','qscan'))
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

    self.friendlyName = name
    self.id = ifoString + '-' + name + '-' + str(times[ifoList[0]])

    pipeline.CondorDAGNode.__init__(self,job)

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


    if cp.has_option("followup-"+name, ifoString + '-output') and string.strip(cp.get("followup-"+name, ifoString + '-output')):
      output = string.strip(cp.get("followup-"+name, ifoString + '-output'))
    else:
      output = dag.publish_path + '/' + job.name + '/' + name + '/' + ifoString
    if not os.access(output,os.F_OK):
      os.makedirs(output)
    else:
      if not os.access(output,os.W_OK):
        print >> sys.stderr, 'path '+output+' is not writable'
        sys.exit(1)

    self.add_var_arg('event')
    qeventConfig = string.strip(cp.get("followup-"+name, ifoString + '-config-file'))
    self.add_var_arg('-p '+qeventConfig)
    self.add_file_arg('-f '+qeventcache)
    self.add_var_arg('-o '+output)
    self.add_var_arg(repr(times[ifoList[0]]))
    eventDuration = string.strip(cp.get("followup-"+name, 'duration'))
    self.add_var_arg(eventDuration)

    self.set_pre_script(job.name + "/cachecat.sh %s %s %s" \
    %(dNode[ifoList[0]].outputFileName, dNode[ifoList[1]].outputFileName, \
    qeventcache))

    #get the absolute output path whatever the path might be in the ini file
    absoutput = os.path.abspath(output)
    self.outputName = absoutput + '/' + repr(times[ifoList[0]]) # redirect output name

    #prepare the string for the output cache
    self.outputCache = ifoString + ' ' + name + ' ' + repr(times[ifoList[0]]) + ' ' + self.outputName + '\n'

    if cp.has_option("followup-"+name,ifoString+'-web') and string.strip(cp.get("followup-"+name,ifoString+'-web')):
      pageOverride = string.strip(cp.get("followup-"+name,ifoString+'-web'))+'/'+repr(times[ifoList[0]])
    else:
      pageOverride = dag.page + '/' + job.name + '/' + name + '/' + ifoString + '/' + repr(times[ifoList[0]])
    self.setupNodeWeb(job,False,dag.webPage.lastSection.lastSub,dag.page,pageOverride,dag.cache)

    if not opts.disable_dag_categories:
      self.set_category(job.name.lower())

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
  def __init__(self, FrCheckJob, procParams, ifo, trig, cp, opts, dag, datafindCache, d_node, datafindCommand):

    try:
      hipeCache = checkHipeCachePath(cp)

      if not hipeCache:
        cacheFile = datafindCache
      else:
        for row in procParams:
          param = row.param.strip("-")
          value = row.value
          if param == 'frame-cache': cacheFile = value

      self.friendlyName = 'Frame Check'
    
      pipeline.CondorDAGNode.__init__(self,FrCheckJob)
      self.add_var_opt("frame-cache", cacheFile)
      self.add_var_opt("frame-check-executable", string.strip(cp.get('followup-frameCheck','executable')))

      self.id = FrCheckJob.name + '-' + ifo + '-' + str(trig.statValue) + '_' + str(trig.eventID)
      self.setupNodeWeb(FrCheckJob,True, dag.webPage.lastSection.lastSub,dag.page,None,dag.cache)

      if not opts.disable_dag_categories:
        self.set_category(FrCheckJob.name.lower())

      try:
        if d_node.validNode and eval('opts.' + datafindCommand):
          self.add_parent(d_node)
      except: pass

      if opts.frame_check:
        dag.addNode(self,self.friendlyName)
        self.validate()
      else: self.invalidate()

    except:
      self.invalidate()
      print "couldn't add frame check job for " + str(ifo) + "@ "+ str(trig.gpsTime[ifo])

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

    if not opts.disable_dag_categories:
      self.set_category(IFOstatus_checkJob.name.lower())

    if opts.ifo_status_check:
      dag.addNode(self,self.friendlyName)
      self.validate()
    else: self.invalidate()

##############################################################################

class followupoddsJob(pipeline.CondorDAGJob, webTheJob):
  """
  A model selection job
  """
  def __init__(self,options,cp,tag_base='FOLLOWUPODDS'):
    """
    """
    self.__name__='followupoddsjob'
    self.__executable=string.strip(cp.get('condor','followupodds'))
    self.__universe="standard"
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    self.setupJobWeb(self.__name__,tag_base)

class followupoddsNode(pipeline.CondorDAGNode,webTheNode):
  """
  Runs an instance of the model selection followup job
  """
  def __init__(self,followupoddsJob,procParamsTable,trig,randomseed,cp,opts,dag):
    #try
    if 1:
      IFOs = trig.ifolist_in_coinc
      time_prior = string.strip(cp.get('followup-odds','time_prior'))
      Nlive = string.strip(cp.get('followup-odds','live_points'))
      srate = string.strip(cp.get('followup-odds','sample_rate'))
      Approximant = string.strip(cp.get('followup-odds','approximant'))
      self.friendlyName = 'Odds followup job'
      pipeline.CondorDAGNode.__init__(self,followupoddsJob)
      cacheFiles=[]
      for ifo in IFOs:
        for row in procParamsTable[ifo]:
          param=row.param.strip("-")
          value=row.value
          if param == 'frame-cache': cacheFile=value
          if param == 'gps-start-time': GPSstart=value
          if param == 'gps-end-time': GPSend=value

        self.add_var_arg("--IFO "+str(ifo))
        self.add_var_arg("--cache " +str(cacheFile))

      outputname = followupoddsJob.name + '/'+followupoddsJob.name+'-' \
                   +trig.ifos+'-'+str(trig.statValue)+'_'+str(trig.eventID)+'.dat'
      self.add_var_opt("Nlive",Nlive)
      self.add_var_opt("GPSstart",GPSstart)
      self.add_var_opt("length",str(float(GPSend)-float(GPSstart)))
      self.add_var_opt("approximant",Approximant)
      self.add_var_opt("out",outputname)
      self.add_var_opt("Nsegs",str((int(GPSend)-int(GPSstart))/8))
      self.add_var_opt("dt",time_prior)
      self.add_var_opt("end_time",trig.gpsTime[ifo])
      self.add_var_opt("Mmin",2.8)
      self.add_var_opt("Mmax",30)
      self.add_var_opt("srate",srate)
#      self.add_var_opt("verbose",'')
      self.id = followupoddsJob.name + '-' + trig.ifos + '-' + str(trig.statValue) + '_' + str(trig.eventID)
      self.outputCache = trig.ifos + ' ' + followupoddsJob.name + ' ' +\
                         self.id.split('-')[-1]+' '+outputname+'\n'

      print "Using IFOs " + str(IFOs)

      self.setupNodeWeb(followupoddsJob,False,None,None,None,dag.cache)

      print "Arguments: " + str(self.get_cmd_line())

      if opts.odds:
        dag.addNode(self,self.friendlyName)
        self.validate()
      else: self.invalidate()

#    except:
    else:
      self.invalidate()
      print "Couldn't add followupOdds job for " + str(trig.gpsTime[ifo])
  

            
      

##############################################################################

class followupmcmcJob(pipeline.CondorDAGJob, webTheJob):
  """
  An mcmc job
  """
  def __init__(self, options, cp, tag_base='FOLLOWUPMCMC'):
    """
    """
    self.__name__ = 'followupmcmcJob'
    self.__executable = string.strip(cp.get('condor','followupmcmc'))
    self.__universe = "standard"
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    self.setupJobWeb(self.__name__,tag_base)

###############################################################################

class followupmcmcNode(pipeline.CondorDAGNode,webTheNode):
  """
  Runs an instance of an mcmc followup job
  """
  def __init__(self, followupmcmcJob, procParams, ifo, trig, randomseed, cp,opts,dag):
    try:
      time_margin = string.strip(cp.get('followup-mcmc','prior-coal-time-marg'))
      iterations = string.strip(cp.get('followup-mcmc','iterations'))
      tbefore = string.strip(cp.get('followup-mcmc','tbefore'))
      tafter = string.strip(cp.get('followup-mcmc','tafter'))
      massmin = string.strip(cp.get('followup-mcmc','massmin'))
      massmax = string.strip(cp.get('followup-mcmc','massmax'))
      dist90 = string.strip(cp.get('followup-mcmc','dist90'))
      dist10 = string.strip(cp.get('followup-mcmc','dist10'))

      self.friendlyName = 'MCMC followup'
      pipeline.CondorDAGNode.__init__(self,followupmcmcJob)
      for row in procParams:
        param = row.param.strip("-")
        value = row.value
        if param == 'frame-cache': cacheFile = value
        if param == 'channel-name': channel = value
        if param == 'gps-end-time': chunk_end = value
        if param == 'gps-start-time': chunk_start = value

      self.add_var_opt("template",string.strip(cp.get('followup-mcmc','template')))
      self.add_var_opt("iterations",iterations)
      self.add_var_opt("randomseed","[" + randomseed[0] + "," + randomseed[1] + "]")
      self.add_var_opt("tcenter","%0.3f"%trig.gpsTime[ifo])
      self.add_var_opt("tbefore",tbefore)
      self.add_var_opt("tafter",tafter)

      tmin = trig.gpsTime[ifo] - float(time_margin)
      tmax = trig.gpsTime[ifo] + float(time_margin)
      self.add_var_opt("priorparameters","[" + massmin + "," + massmax + "," + str(tmin) + "," + str(tmax) + "," + dist90 + "," + dist10 + "]")

      param_mchirp = getattr(trig.coincs,ifo).mchirp
      param_eta = getattr(trig.coincs,ifo).eta
      param_distance = getattr(trig.coincs,ifo).eff_distance
      self.add_var_opt("guess","[" + str(param_mchirp) + "," + str(param_eta) + "," + str(trig.gpsTime[ifo]) + "," + str(param_distance) + "]")
      self.add_var_opt("fixed","[altitude=0,azimuth=0,inclination=0,polarisation=0]")

      self.add_var_opt("readdata","")

      ########################################################################
      # GET THE FRAME FILE INFO - THIS NEEDS TO BE CHANGED !!!
      # THE MCMC CODE SHOULD TAKE THE SAME PSD AS LALAPPS_INSPIRAL.
      ########################################################################

      self.add_var_opt("cachefiledata",cacheFile)
      self.add_var_opt("cachefilenoise",cacheFile)

      self.add_var_opt("filechannel",channel)

      datainchunk_before = int(trig.gpsTime[ifo]) - 75 - 64 - int(chunk_start)

      datainchunk_after = int(chunk_end) - 64 - int(trig.gpsTime[ifo]) - 32
      if datainchunk_after > datainchunk_before:
        self.add_var_opt("psdestimatestart",int(trig.gpsTime[ifo]) + 32)
        self.add_var_opt("psdestimateend",int(chunk_end) - 64)
      else:
        self.add_var_opt("psdestimatestart",int(chunk_start) + 64)
        self.add_var_opt("psdestimateend",int(trig.gpsTime[ifo]) - 75)

      self.add_var_opt("importanceresample",10000)

      self.id = followupmcmcJob.name + '-' + ifo + '-' + str(trig.statValue) + '_' + str(trig.eventID) + '_' + randomseed[0] + '_' + randomseed[1]
      outputName = followupmcmcJob.name+'/'+self.id+'.txt'
      self.outputCache = ifo + ' ' + followupmcmcJob.name + ' ' + self.id.split('-')[-1] + ' ' + outputName + '\n'

      self.setupNodeWeb(followupmcmcJob,False,None,None,None,dag.cache)
      self.add_var_opt("logfilename",outputName)

      if not opts.disable_dag_categories:
        self.set_category(followupmcmcJob.name.lower())

      if opts.mcmc:
        dag.addNode(self,self.friendlyName)
        self.validate()
      else: self.invalidate()

    except:
      self.invalidate()
      print "couldn't add followupmcmc job for " + str(ifo) + "@ "+ str(trig.gpsTime[ifo])

###############################################################################

class plotmcmcJob(pipeline.CondorDAGJob, webTheJob):
  """
  A plot mcmc job
  """
  def __init__(self, options, cp, tag_base='PLOTMCMC'):
    """
    """
    self.__name__ = 'plotmcmcJob'
    self.__executable = string.strip(cp.get('condor','plotmcmc'))
    self.__universe = "vanilla"
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    self.setupJobWeb(self.__name__,tag_base)

###############################################################################

class plotmcmcNode(pipeline.CondorDAGNode,webTheNode):
  """
  Runs an instance of  plotmcmc job
  """
  def __init__(self, plotmcmcjob, ifo, trig, mcmcIdList, cp,opts,dag):

    try:
      self.friendlyName = 'plot MCMC'
      pipeline.CondorDAGNode.__init__(self,plotmcmcjob)

      if cp.has_option('followup-plotmcmc','burnin'):
        burnin = string.strip(cp.get('followup-plotmcmc','burnin'))
        if burnin.strip():
          self.add_var_opt("burnin",burnin)

      plot_routine = string.strip(cp.get('followup-plotmcmc','plot_routine'))
      executable = string.strip(cp.get('followup-plotmcmc','executable'))
      sim = None
      try:
        sim = isinstance(trig.coincs.sim,lsctables.SimInspiral)
      except:
        pass
      if sim:
        time = eval("trig.coincs.sim." + ifo[0:1].lower() + "_end_time")
        time_ns = eval("trig.coincs.sim." + ifo[0:1].lower() + "_end_time_ns")
        gps = float(time) + float(time_ns)/1000000000.
        mchirp = trig.coincs.sim.mchirp
        eta = trig.coincs.sim.eta
        distance = trig.coincs.sim.distance
        phi = trig.coincs.sim.phi0
      else:
        gps = trig.gpsTime[ifo]
        mchirp = getattr(trig.coincs,ifo).mchirp
        eta = getattr(trig.coincs,ifo).eta
        distance = getattr(trig.coincs,ifo).eff_distance
        phi = "0.0"

      self.add_var_opt("plot-routine",plot_routine)
      self.add_var_opt("executable",executable)
      self.add_var_opt("reference-time",gps)
      self.add_var_opt("reference-mchirp",mchirp)
      self.add_var_opt("reference-eta",eta)
      self.add_var_opt("reference-distance",distance)
      self.add_var_opt("reference-phi",phi)
      # get the list of MCMC .txt files to be used as input
      mcmcfilelist = ""
      for mcmcId in mcmcIdList:
        mcmcfilelist += mcmcId.split('-')[0]+'/' + mcmcId + '.txt,' # here we assume that the directory name is mcmcId.split('-')[0]
      self.add_var_opt("mcmc-file",mcmcfilelist.strip(','))

      self.id = plotmcmcjob.name + '-' + ifo + '-' + str(trig.statValue) + '_' + str(trig.eventID)
      self.add_var_opt("identity",self.id)

      if cp.has_option('followup-plotmcmc', 'output') and string.strip(cp.get('followup-plotmcmc', 'output')):
        outputpath = string.strip(cp.get('followup-plotmcmc', 'output'))
      else:
        outputpath = dag.publish_path + '/' + plotmcmcjob.name
      if not os.access(outputpath,os.F_OK):
        os.makedirs(outputpath)
      else:
        if not os.access(outputpath,os.W_OK):
          print >> sys.stderr, 'path '+outputpath+' is not writable'
          sys.exit(1)

      if cp.has_option('followup-plotmcmc','web') and string.strip(cp.get('followup-plotmcmc','web')):
        webpath = string.strip(cp.get('followup-plotmcmc','web'))
      else:
        webpath = dag.page + '/' + plotmcmcjob.name

      output_page = webpath + '/' + self.id
      self.outputCache = self.id.replace('-',' ') + " " + os.path.abspath(outputpath) + "/" + self.id + "\n"
      self.setupNodeWeb(plotmcmcjob,False,dag.webPage.lastSection.lastSub,None,output_page,dag.cache)

      self.add_var_opt("output-path",outputpath)

      if not opts.disable_dag_categories:
        self.set_category(plotmcmcjob.name.lower())

      # only add a parent if it exists
      for node in dag.get_nodes():
        if isinstance(node,followupmcmcNode):
          if not node.id.find(ifo + '-' + str(trig.statValue) + '_' + str(trig.eventID)) == -1:
            try:
              if node.validNode: self.add_parent(node)
            except: pass

      if opts.plot_mcmc:
        dag.addNode(self,self.friendlyName)
        self.validate()
      else:
        self.invalidate()

    except:
      self.invalidate()
      print "couldn't add plot mcmc job for " + str(ifo) + "@ "+ str(trig.gpsTime[ifo])


