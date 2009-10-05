"""
This module contains condor jobs / node classes for the SQlite Triggered Follow Up dag

$Id$

"""

__author__ = 'Chad Hanna <channa@phys.lsu.edu>, Cristina Torres <cristina.torres@ligo.org>, Romain Gouaty <gouaty@lapp.in2p3.fr>'
__date__ = '$Date$'
__version__ = '$Revision$'[11:-2]

try:
        import sqlite3
except ImportError:
        # pre 2.5.x
        from pysqlite2 import dbapi2 as sqlite3

import sys, os, copy, math, time, math, subprocess, socket, re, string
from optparse import *
import tempfile
import ConfigParser
import urlparse
from UserDict import UserDict
sys.path.append('@PYTHONLIBDIR@')

from glue import segments
from glue import segmentsUtils
from glue.ligolw import ligolw
from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw import dbtables
from glue.ligolw import utils
from glue import pipeline
from pylal import db_thinca_rings
from lalapps import inspiral

from pylal.xlal.datatypes.ligotimegps import LIGOTimeGPS
dbtables.lsctables.LIGOTimeGPS = LIGOTimeGPS

###############################################################################
##### UTILITY FUNCTIONS #######################################################
###############################################################################


def mkdir(output):
        # MAKE SURE WE CAN WRITE TO THE OUTPUT DIRECTORY
        if not os.access(output,os.F_OK): os.makedirs(output)
        else:
                if not os.access(output,os.W_OK):
                        print >> sys.stderr, 'path '+output+' is not writable'
                        sys.exit(1)

def science_run(time):
	if time >= 815155213 and time <= 875232014: return 's5'
	if time >= 931035296 and time <= 999999999: return 's6'
	print >>sys.stderr, "COULD NOT DETERMINE SCIENCE RUN from %d" % (int(time),)
	sys.exit(1)

def figure_out_type(time, ifo, data_type='hoft'):
        """
        Run boundaries (from CBC analyses):
        VSR1: 863557214 - 875232014
        S5:   815155213 - 875232014
        VSR2/S6: 931035296 - ...
        Frame types for S5/VSR1:
        () For RDS_R_L1 data set:
        type    channel_name
        RDS_R_L1     H1:LSC-DARM_ERR
        RDS_R_L1     H2:LSC-DARM_ERR
        RDS_R_L1     L1:LSC-DARM_ERR
        () For hoft data:
        type    channel_name
        H1_RDS_C03_L2  H1:LSC-STRAIN
        H2_RDS_C03_L2  H2:LSC-STRAIN
        L1_RDS_C03_L2  L1:LSC-STRAIN
        HrecV2_16384Hz      V1:h_16384Hz
        Frame types for S6/VSR2:
        () For RDS_R_L1 data set:
        type    channel_name
        H1_RDS_R_L1   H1:LSC-DARM_ERR
        L1_RDS_R_L1   L1:LSC-DARM_ERR
        () For hoft data:
        H1_LDAS_C00_L2  H1:LDAS-STRAIN
        L1_LDAS_C00_L2  L1:LDAS-STRAIN
        HrecOnline      V1:h_16384Hz
        """
        L1HoftTypes=(
                ("L1_RDS_C03_L2","L1:LSC-STRAIN",815155213,875232014),
                ("L1_LDAS_C00_L2","L1:LDAS-STRAIN",931035296,999999999)
                )
        H1HoftTypes=(
                ("H1_RDS_C03_L2","H1:LSC-STRAIN",815155213,875232014),
                ("H1_LDAS_C00_L2","H1:LDAS-STRAIN",931035296,999999999)
                )
        H2HoftTypes=(
                ("H2_RDS_C03_L2","H2:LSC-STRAIN",815155213,875232014),
                ("H1_LDAS_C00_L2","H1:LDAS-STRAIN",931035296,999999999)
                )
        V1HoftTypes=(
                ("HrecV2_16384Hz","V1:h_16384Hz",863557214,875232014),
                ("HrecOnline","V1:h_16384Hz",931035296,999999999)
                )
        L1RdsTypes=(
                ("RDS_R_L1","L1:LSC-DARM_ERR",815155213,875232014),
                ("L1_RDS_R_L1","L1:LSC-DARM_ERR",931035296,999999999)
                )
        H1RdsTypes=(
                ("RDS_R_L1","H1:LSC-DARM_ERR",815155213,875232014),
                ("H1_RDS_R_L1","H1:LSC-DARM_ERR",931035296,999999999)
                )
        H2RdsTypes=(
                ("RDS_R_L1","H2:LSC-DARM_ERR",815155213,875232014),
                ("H1_RDS_R_L1","H1:LSC-DARM_ERR",931035296,999999999)
                )
        V1RdsTypes=(
                ("raw","V1:Pr_B1_ACp",863557214,875232014),
                ("raw","V1:Pr_B1_ACp",931035296,999999999)
                )
        channelMap={
                "L1":{"hoft":L1HoftTypes,"rds":L1RdsTypes},
                "H1":{"hoft":H1HoftTypes,"rds":H1RdsTypes},
                "H2":{"hoft":H2HoftTypes,"rds":H2RdsTypes},
                "V1":{"hoft":V1HoftTypes,"rds":V1RdsTypes}
                }
        #Use the IFO type to select the channel type
        foundType=""
        foundChannel=""
        for type,channel,start,stop in channelMap[ifo][data_type]:
                if ((start<=time) and (time<=stop)):
                        foundType=type
                        foundChannel=channel
                        break
        if foundType == "":
                print time,ifo + " not found in method stfu_pipe.figure_out_type"
                os.abort()
        return str(foundType), str(foundChannel)

###############################################################################
##### CONDOR JOB CLASSES ######################################################
###############################################################################

# DO SOME STANDARD STUFF FOR FU JOBS
class FUJob(object):
        """
        """

        def __init__(self):
                pass

        def __conditionalLoadDefaults__(self,defaults,cp):
                if not(cp.has_section(defaults["section"])):
                        cp.add_section(defaults["section"])
                for key, val in defaults["options"].iteritems():
                        if not cp.has_option(defaults["section"], key):
                                cp.set(defaults["section"], key, val)

        def setupJob(self, name, dir= '', tag_base=None, cp=None):
                # Give this job a name.  Then make directories for the log files and such
                # This name is important since these directories will be included in
                # the web tree.
		mkdir(dir)
                self.name = name
		self.relPath = dir + '/' + name
		self.outputPath = os.getcwd() + '/' + self.relPath + '/'
		#path = name + '/' + dir
                if not os.path.exists(self.relPath):
                        os.mkdir(self.relPath)
                if not os.path.exists(self.relPath+'/logs'):
                        os.mkdir(self.relPath+'/logs')
                if not os.path.exists(self.relPath+'/Images'):
                        os.mkdir(self.relPath+'/Images')
                if not os.path.exists(self.relPath+'/DataProducts'):
                        os.mkdir(self.relPath+'/DataProducts')
                # Set up the usual stuff and name the log files appropriately
                self.tag_base = tag_base
                self.add_condor_cmd('environment',"KMP_LIBRARY=serial;MKL_SERIAL=yes")
                self.set_sub_file(name+'.sub')
                self.set_stdout_file(self.outputPath+'/logs/'+name+'-$(macroid).out')
                self.set_stderr_file(self.outputPath+'/logs/'+name+'-$(macroid).err')
                if cp:
                        if cp.has_section("condor-memory-requirement") and cp.has_option("condor-memory-requirement",name):
                                requirement = cp.getint("condor-memory-requirement",name)
                                self.add_condor_cmd("Requirements", "(Memory > " + str(requirement) + ")")

# QSCAN JOB CLASS
class qscanJob(pipeline.CondorDAGJob, FUJob):
        """
        A qscan job
        """
        def __init__(self, opts, cp, dir='', tag_base='QSCAN'):
                """
                """
                self.__executable = string.strip(cp.get('fu-condor','qscan'))
                self.__universe = "vanilla"
                pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
                self.setupJob(tag_base,dir=dir,cp=cp,tag_base=tag_base)
                self.setup_checkForDir()

	def is_dax(self):
		return False

        def setup_checkForDir(self):
                # create a shell script to check for the existence of the qscan output directory and rename it if needed
                checkdir_script = open('checkForDir.sh','w')
                checkdir_script.write("""#!/bin/bash
if [ -d $1/$2 ]
then
matchingList=$(echo $(find $1 -name $2.bk*))
COUNTER=1
for file in $matchingList
   do
     let COUNTER=COUNTER+1
   done
mv $1/$2 $1/$2.bk.$COUNTER
fi
                """)
                checkdir_script.close()
                os.chmod('checkForDir.sh',0755)



# A CLASS TO DO FOLLOWUP INSPIRAL JOBS 
class followUpInspJob(inspiral.InspiralJob,FUJob):

        def __init__(self,cp,dir='',type='plot'):

                inspiral.InspiralJob.__init__(self,cp)

                if type == 'head':
                        self.set_executable(string.strip(cp.get('fu-condor','inspiral_head')))
                self.name = 'followUpInspJob' + type
                self.setupJob(self.name, dir=dir, cp=cp)

# JOB CLASS FOR PRODUCING A SKYMAP
class skyMapJob(pipeline.CondorDAGJob,FUJob):
        """
        Generates sky map data
        """
        def __init__(self, options, cp, dir='', tag_base='SKY_MAP'):
                """
                """
                self.__prog__ = 'lalapps_skyMapJob'
                self.__executable = string.strip(cp.get('fu-condor','lalapps_skymap'))
                self.__universe = "standard"
                pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
                self.add_condor_cmd('getenv','True')
                self.setupJob(self.__prog__,dir=dir, tag_base=tag_base)

		self.ra_res = cp.get("fu-skymap", 'ra-res')
		self.dec_res = cp.get("fu-skymap", 'dec-res')
		self.sample_rate = cp.get("fu-skymap", 'sample-rate')
		
# JOB CLASS FOR PRODUCING SKYMAP PLOT
class skyMapPlotJob(pipeline.CondorDAGJob,FUJob):
        """
        Plots the sky map output of lalapps_skymap
        """
        def __init__(self, options, cp, dir='',tag_base='SKY_PLOT'):
                """
                """
                self.__prog__ = 'pylal_skyPlotJob'
                self.__executable = string.strip(cp.get('fu-condor','pylal_skyPlotJob'))
                self.__universe = "vanilla"
                pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
                self.add_condor_cmd('getenv','True')
                self.setupJob(self.__prog__,dir=dir,tag_base=tag_base)

		self.ra_res = cp.get("fu-skymap", 'ra-res')
		self.dec_res = cp.get("fu-skymap", 'dec-res')
		self.sample_rate = cp.get("fu-skymap", 'sample-rate')

# JOB CLASS FOR PLOTTING SNR AND CHISQ
class plotSNRChisqJob(pipeline.CondorDAGJob,FUJob):
        """
        A followup plotting job for snr and chisq time series
        """
        def __init__(self, options, cp, dir='', tag_base='PLOT_FOLLOWUP'):
                """
                """
                self.__prog__ = 'plotSNRCHISQJob'
                self.__executable = string.strip(cp.get('fu-condor','plotsnrchisq'))
                self.__universe = "vanilla"
                pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
                self.add_condor_cmd('getenv','True')
                self.setupJob(self.__prog__,tag_base=tag_base,dir=dir)

class htDataFindJob(pipeline.LSCDataFindJob,FUJob):
        def __init__(self, config_file, dir='', name='datafind'):
    
                self.name = name
    
                # unfortunately the logs directory has to be created before we call LSCDataFindJob
                try:
                        os.mkdir(self.name)
                        os.mkdir(self.name + '/logs')
                except: pass
                pipeline.LSCDataFindJob.__init__(self, self.name, self.name + '/logs', config_file)
                self.setup_cacheconv(config_file)
    
        def setup_cacheconv(self,cp):
                # create a shell script to call convertlalcache.pl if the value of $RETURN is 0
                convert_script = open('cacheconv.sh','w')
		#FIXME changed convert cache script to not fail on previous error?
                convert_script.write("""#!/bin/bash
%s ${1} ${2}
                """ % string.strip(cp.get('fu-condor','convertcache')))
                convert_script.close()
                os.chmod('cacheconv.sh',0755)

#The class responsible for running the data quality flag finding job
class findFlagsJob(pipeline.CondorDAGJob, FUJob):
        """
        A job which queries the ldbd database holding segment
        information about the DQ flags.
        """
        defaults={"section":"fu-condor",
                  "options":{"universe":"local",
                             "dqflags":"followupQueryDQ.py"}
                  }

        def __init__(self, opts, cp, dir='', tag_base="DQFLAGS"):
                """
                """
                self.__conditionalLoadDefaults__(findFlagsJob.defaults,cp)
                self.__prog__ = 'findFlagsJob'
                self.__executable = string.strip(cp.get('fu-condor','dqflags'))
                self.__universe = string.strip(cp.get('fu-condor','universe'))
                pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
                self.add_condor_cmd('getenv','True')
                self.setupJob(self.__prog__,tag_base=tag_base, dir=dir)

#The class responsible for checking for know veto flags
class findVetosJob(pipeline.CondorDAGJob,FUJob):
        """
        A job instance that queries the segment database for known
        types of active veto intervals.
        """
        defaults={"section":"fu-condor",
                  "options":{"universe":"local",
                             "vetoflags":"followupQueryDQ.py"}
                  }
        def __init__(self, opts,cp, dir='', tag_base="VETOS"):
                """
                """
                self.__conditionalLoadDefaults__(findVetosJob.defaults,cp)
                self.__prog__ = 'findVetosJob'
                self.__executable = string.strip(cp.get('fu-condor','vetoflags'))
                self.__universe = string.strip(cp.get('fu-condor','universe'))
                pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
                self.add_condor_cmd('getenv','True')
                self.setupJob(self.__prog__,tag_base=tag_base, dir=dir)

#The class responsible for running one type of parameter consistency check
class effDRatioJob(pipeline.CondorDAGJob,FUJob):
        """
        A job that performs parameter consitency check for a trigger
        being followed up.
        """
        defaults={"section":"fu-condor",
                  "options":{"universe":"local",
                             "effDRatio":"followupRatioTest.py"}
                  }
        def __init__(self, opts, cp, dir='', tag_base="effDRatioTest"):
                """
                """
                self.__conditionalLoadDefaults__(effDRatioJob.defaults,cp)
                self.__prog__ = 'effDRatioTest'
                self.__executable = string.strip(cp.get('fu-condor','effDRatio'))
                self.__universe = string.strip(cp.get('fu-condor','universe'))
                pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
                self.add_condor_cmd('getenv','True')
                self.setupJob(self.__prog__,tag_base=tag_base, dir=dir)


#############################################################################
###### CONDOR NODE CLASSES ##################################################
#############################################################################

class FUNode:
        """
        """

        def __init__(self):
                pass

        def __conditionalLoadDefaults__(self,defaults,cp):
                if not(cp.has_section(defaults["section"])):
                        cp.add_section(defaults["section"])
                for key, val in defaults["options"].iteritems():
                        if not cp.has_option(defaults["section"], key):
                                cp.set(defaults["section"], key, val)

        def setupNodeWeb(self, job, passItAlong=False, content=None, page=None,webOverride=None,cache=None):
                self.jobName = job.name
                if passItAlong:
                        self.add_var_opt("output-path",job.outputPath)
                        self.add_var_opt("enable-output","")
                if cache:
                        cache.appendCache(job.name,job.outputPath)
                try:
                        if self.outputCache:
                                cache.appendSubCache(job.name,self.outputCache)
                except: pass

class fuQscanNode(pipeline.CondorDAGNode):
        """
QScan node.  This node writes its output to the web directory specified in
the inifile + the ifo and gps time.  For example:
 
        /archive/home/channa/public_html/followup/htQscan/H1/999999999.999

The omega scan command line is 

        wpipeline scan -r -c H1_hoft.txt -f H-H1_RDS_C03_L2-870946612-870946742.qcache -o QSCAN/foreground-hoft-qscan/H1/870946677.52929688 870946677.52929688

        """
        def __init__(self, dag, job, cp, opts, time, ifo, frame_cache, p_nodes=[], type="ht", variety="fg"):
                """
                """
                pipeline.CondorDAGNode.__init__(self,job)

                if variety == "bg":
                  self.add_var_arg('scan')
                else:
                  self.add_var_arg('scan -r')
		config = self.fix_config_for_science_run( cp.get('fu-'+variety+'-'+type+'-qscan', ifo+'config').strip(), time )
                self.add_var_arg("-c " + config )

                if cp.has_option('fu-output','output-dir') and cp.get('fu-output','output-dir'):
                  output = cp.get('fu-output','output-dir') + '/' + type + 'Qscan' + '/' + ifo
                else:
                  output = type + 'Qscan' + '/' + ifo

                # CREATE AND MAKE SURE WE CAN WRITE TO THE OUTPUT DIRECTORY
                mkdir(output)

                self.add_var_arg("-o "+output+"/"+repr(time))
		
		# ADD FRAME CACHE FILE
		self.add_var_arg("-f "+frame_cache)
		
                self.add_var_arg(repr(time))

                self.set_pre_script("checkForDir.sh %s %s" %(output, repr(time)))

                if not(cp.has_option('fu-'+variety+'-'+type+'-qscan','remote-ifo') and cp.get('fu-'+variety+'-'+type+'-qscan','remote-ifo').strip()):
                  for node in p_nodes:
                    self.add_parent(node)
                  dag.add_node(self)

	def fix_config_for_science_run(self, config, time):
		run = science_run(time)
		config_path = os.path.split(config)
		out = "/".join([config_path[0], config_path[1].replace('s5',run).replace('s6',run)])
		return out

class fuDataFindNode(pipeline.LSCDataFindNode):
    
        def __init__(self, dag, job, cp, opts, ifo, sngl=None, qscan=False, trigger_time=None, data_type="hoft", p_nodes=[]):

                self.outputFileName = ""
                pipeline.LSCDataFindNode.__init__(self,job)
                if qscan:
                        if sngl: time = sngl.time
                        else: time = trigger_time
                        self.outputFileName = self.setup_qscan(job, cp, time, ifo, data_type)
                else:
                        if not sngl:
                                print >> sys.stderr, "argument \"sngl\" should be provided to class fuDataFindNode"
                                sys.exit(1)
                        self.outputFileName = self.setup_inspiral(job, cp, sngl, ifo)
                if not qscan or not(cp.has_option('fu-q-'+data_type+'-datafind','remote-ifo') and cp.get('fu-q-'+data_type+'-datafind','remote-ifo').strip()):
                    for node in p_nodes:
                        self.add_parent(node)
                    dag.add_node(self)

        def setup_qscan(self, job, cp, time, ifo, data_type):
                # 1s is substracted to the expected startTime to make sure the window
                # will be large enough. This is to be sure to handle the rouding to the
                # next sample done by qscan.
                type, channel = figure_out_type(time,ifo,data_type)
                self.set_type(type)
                self.q_time = float(cp.get("fu-q-"+data_type+"-datafind","search-time-range"))/2.
                self.set_observatory(ifo[0])
                self.set_start(int( time - self.q_time - 1.))
                self.set_end(int( time + self.q_time + 1.))
                lalCache = self.get_output()
                qCache = lalCache.rstrip("cache") + "qcache"
                self.set_post_script(os.getcwd()+"/cacheconv.sh %s %s" %(lalCache,qCache) )
                return(qCache)

        def setup_inspiral(self, job, cp, sngl, ifo):
                # 1s is substracted to the expected startTime to make sure the window
                # will be large enough. This is to be sure to handle the rouding to the
                # next sample done by qscan.
                type, channel = figure_out_type(sngl.get_gps_start_time(),ifo)
                self.set_type(type)
                self.set_observatory(ifo[0])
		#FIXME use proper pad, not hardcode to 64
                self.set_start(sngl.get_gps_start_time()-64)
                self.set_end(sngl.get_gps_end_time()+64)
                lalCache = self.get_output()
                return(lalCache)


class followUpInspNode(inspiral.InspiralNode,FUNode):

  #def __init__(self, inspJob, procParams, ifo, trig, cp,opts,dag, datafindCache, d_node, datafindCommand, type='plot', sngl_table = None):
        def __init__(self, dag, job, cp, opts, sngl, frame_cache, tag, p_nodes=[]):

                tlen = 1.0
                self.output_file_name = ""
                pipeline.CondorDAGNode.__init__(self,job)

                #FIXME HANDLE INJECTION FILES AND datafind cache
                # injFile = self.checkInjections(cp)
                # self.set_injections( injFile )

                self.set_trig_start( int(sngl.time - tlen + 0.5) )
                self.set_trig_end( int(sngl.time + tlen + 0.5) )
                self.add_var_opt("write-snrsq","")
                self.add_var_opt("write-chisq","")
                self.add_var_opt("write-spectrum","")
                self.add_var_opt("write-template","")
                self.add_var_opt("write-cdata","")

                skipParams = ['minimal-match', 'bank-file', 'user-tag', 'injection-file', 'trig-start-time', 'trig-end-time']

                extension = ".xml"
                for row in sngl.process_params:
                        param = row.param.strip("-")
                        value = row.value
                        # override some options
                        if param == 'frame-cache': value = frame_cache
                        if param == 'snr-threshold': value = "0.1"
                        if param == 'do-rsq-veto': continue
                        if param == 'enable-rsq-veto': continue
                        if param == 'chisq-threshold': value = "1.0e+06"
                        if param == 'cluster-method': value = 'window'
                        if param == 'cluster-window': continue
                        if param in skipParams: continue
                        if param == 'injection-file': value = sngl.inj_file_name
                        self.add_var_opt(param,value)
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

                self.add_var_opt('cluster-window',str( tlen / 2.))
                self.add_var_opt('disable-rsq-veto',' ')
                bankFile = self.write_trig_bank(sngl, 'trig_bank/' + sngl.ifo + '-TRIGBANK_FOLLOWUP_' + repr(sngl.time) + '.xml.gz')
                self.set_bank(bankFile)

                self.set_user_tag( tag.upper() + "_FOLLOWUP_" + repr(sngl.time) )
                self.__usertag = tag.upper() + "_FOLLOWUP_" + repr(sngl.time)
      
                self.output_file_name = job.outputPath + sngl.ifo + "-INSPIRAL_" + self.__ifotag + "_" + self.__usertag + "-" + self.__start + "-" + str(int(self.__end)-int(self.__start)) + extension
                self.outputCache = sngl.ifo + ' ' + 'INSPIRAL' + ' ' + str(self.__start) + ' ' + str(int(self.__end)-int(self.__start)) + ' ' + self.output_file_name  + '\n' + sngl.ifo + ' ' + 'INSPIRAL-FRAME' + ' ' + str(self.__start) + ' ' + str(int(self.__end)-int(self.__start)) + ' ' + self.output_file_name.replace(extension,".gwf") + '\n'

                self.add_var_opt("output-path",job.outputPath)

                #add parents and put node in dag
                for node in p_nodes: self.add_parent(node)
                dag.add_node(self)

        def write_trig_bank(self,sngl, name):
                try:
                        os.mkdir('trig_bank')
                except: pass
                xmldoc = ligolw.Document()
                xmldoc.appendChild(ligolw.LIGO_LW())

                process_params_table = lsctables.New(lsctables.ProcessParamsTable)
                xmldoc.childNodes[-1].appendChild(process_params_table)

                sngl_inspiral_table = lsctables.New(lsctables.SnglInspiralTable)
                xmldoc.childNodes[-1].appendChild(sngl_inspiral_table)
                sngl_inspiral_table.append(sngl.row)

                utils.write_filename(xmldoc, name, verbose=False, gz = True)
                return name

class findFlagsNode(pipeline.CondorDAGNode,FUNode):
        """
        This class is resposible for setting up a node to perform a
        query for the DQ flag for the trigger which under review.
        EXAMPLE
        followupQueryDQ.py --window=60,15 --trigger-time=929052945 --output-format=moinmoin --segment-url="ldbd://segdb.ligo.caltech.edu:30015" --output-file=dqResults.wiki
        """
        defaults={"section":"findFlags",
                  "options":{"window":"60,15",
                             "segment-url":"ldbd://segdb.ligo.caltech.edu:30015",
                             "output-format":"moinmoin",
                             "output-file":"dqResults.wiki"}
                  }
        def __init__(self, dag, job, cp, opts, time):
                """
                """
                self.__conditionalLoadDefaults__(findFlagsNode.defaults,cp)
                pipeline.CondorDAGNode.__init__(self,job)
                self.add_var_opt("trigger-time",time)
                #Output filename
                oFilename="%s_%s"%(str(int(float(time))),cp.get('findFlags','output-file'))
                self.add_var_opt("output-file",job.outputPath+'/DataProducts/'+oFilename)
                self.add_var_opt("segment-url",cp.get('findFlags','segment-url'))
                self.add_var_opt("output-format",cp.get('findFlags','output-format'))
                self.add_var_opt("window",cp.get('findFlags','window'))
                dag.add_node(self)

class findVetosNode(pipeline.CondorDAGNode,FUNode):
        """
        This class is responsible for creating a node in the dag which
        queries the segment database for veto segments active around
        the trigger time of the candidate.
        Command line example:
        followupQueryVeto.py --window=60,15 --trigger-time=929052945 --output-format=moinmoin --segment-url="ldbd://segdb.ligo.caltech.edu:30015" --output-file=vetoResults.wiki
        """
        defaults={"section":"findVetoes",
                  "options":{"window":"60,15",
                             "segment-url":"ldbd://segdb.ligo.caltech.edu:30015",
                             "output-format":"moinmoin",
                             "output-file":"vetoResults.wiki"}
                  }
        def __init__(self, dag, job, cp, opts, time):
                """
                """
                self.__conditionalLoadDefaults__(findVetosNode.defaults,cp)
                pipeline.CondorDAGNode.__init__(self,job)
                self.add_var_opt("trigger-time",time)
                #Output filename
                oFilename="%s_%s"%(str(int(float(time))),cp.get('findFlags','output-file'))
                self.add_var_opt("output-file",job.outputPath+'/DataProducts/'+oFilename)
                self.add_var_opt("segment-url",cp.get('findFlags','segment-url'))
                self.add_var_opt("output-format",cp.get('findFlags','output-format'))
                self.add_var_opt("window",cp.get('findFlags','window'))
                dag.add_node(self)

class effDRatioNode(pipeline.CondorDAGNode,FUNode):
        """
        This Node class performs a parameter consistency check using the
        sites claiming to detect the trigger and the observed
        effective distance at each site. A command line example is
        below:
        followupRatioTest.py -R /archive/home/ctorres/public_html/DQstuff/ratioTest.pickle -iL1 -jH1 -kV1 -I10 -J10 -K5 -A 1 -B 1 -C 1.0001 -pmoinmoin -o mmTable.wiki
        """
        defaults={"section":"effDRatio",
                  "options":{"output-file":"effDRatio.wiki",
                             "output-format":"moinmoin",
                             "snr-ratio-test":"/archive/home/ctorres/public_html/DQstuff/ratioTest.pickle"}
                  }
        def __init__(self, dag, job, cp, opts, coincEvent=None):
                """
                """
                self.__conditionalLoadDefaults__(effDRatioNode.defaults,cp)
                pipeline.CondorDAGNode.__init__(self,job)
                oFilename="%s_%s"%(str("%10.3f"%(coincEvent.time)).replace(".","_"),cp.get('effDRatio','output-file'))
                self.add_var_opt("output-file",job.outputPath+'/DataProducts/'+oFilename)
                #Grab Sngl propteries from Coinc object
                index=1
                for ifo,snglEvent in coincEvent.sngl_inspiral.items():
                        myIFO=snglEvent.ifo
                        mySNR=snglEvent.snr
                        myTIME=snglEvent.time
                        self.add_var_opt("ifo%i"%(index),myIFO)
                        self.add_var_opt("snr%i"%(index),mySNR)
                        self.add_var_opt("time%i"%(index),myTIME)
                dag.add_node(self)

##############################################################################
###### CONDOR DAG THINGY #####################################################
##############################################################################

class followUpDAG(pipeline.CondorDAG):
        def __init__(self, config_file, cp):
                log_path = cp.get('fu-output','log-path').strip()
                self.basename = re.sub(r'\.ini',r'', os.path.split(config_file)[1])
                tempfile.tempdir = log_path
                tempfile.template = self.basename + '.dag.log.'
                logfile = tempfile.mktemp()
                fh = open( logfile, "w" )
                fh.close()
                pipeline.CondorDAG.__init__(self,logfile)
                self.set_dag_file(self.basename)
                self.jobsDict = {}
                self.node_id = 0
                # The list remote_nodes will contain the list of nodes run remotely 
                # (such as V1 qscans)
                self.remote_nodes = []
        def add_node(self,node):
                self.node_id += 1
                node.add_macro("macroid", self.node_id)
                pipeline.CondorDAG.add_node(self, node)


###############################################################################
###### CONFIG PARSER WRAPPING #################################################
###############################################################################
class create_default_config(object):
	def __init__(self, config=None):
		cp = ConfigParser.ConfigParser()
		self.cp = cp
		self.time_now = "_".join([str(i) for i in time.gmtime()[0:6]])
		home_base = self.__home_dirs()
		
		# CONDOR SECTION NEEDED BY THINGS IN INSPIRAL.PY
		cp.add_section("condor")
		cp.set("condor","datafind",self.which("ligo_data_find"))
		cp.set("condor","inspiral",self.which("lalapps_inspiral"))
		cp.set("condor","universe","standard")
		
		# DATAFIND SECTION
		cp.add_section("datafind")

		# FU-CONDOR SECTION
		cp.add_section("fu-condor")
		cp.set("fu-condor","plotsnrchisq",self.which("plotsnrchisq_pipe"))
		cp.set("fu-condor","lalapps_skymap",self.which("lalapps_skymap"))
		cp.set("fu-condor","pylal_skyPlotJob",self.which("pylal_plot_inspiral_skymap"))
		cp.set("fu-condor","datafind",self.which("ligo_data_find"))
		cp.set("fu-condor","convertcache",self.which("convertlalcache.pl"))
		cp.set("fu-condor","query_dq",self.which("pylal_query_dq"))
		#FIXME SET THIS TO SOMETHING THAT WORKS
		cp.set("fu-condor","qscan",home_base+"/romain/opt/omega/omega_r2062_glnxa64_binary/bin/wpipeline")


		# fu-q-hoft-datafind SECTION
		cp.add_section("fu-q-hoft-datafind")
		cp.set("fu-q-hoft-datafind","search-time-range","128")

		# fu-q-rds-datafind SECTION
		cp.add_section("fu-q-rds-datafind")
		cp.set("fu-q-rds-datafind","search-time-range","1024")
		
		# fu-fg-ht-qscan SECTION
		cp.add_section("fu-fg-ht-qscan")
		for config in ["H1config","H2config","L1config","V1config"]:
			cp.set("fu-fg-ht-qscan",config,self.__qscan_config("s5_foreground_" + config[:2] + "_hoft_cbc.txt"))

		# FU-SKYMAP SECTION
		cp.add_section("fu-skymap")
		cp.set("fu-skymap","ra-res","1024")
		cp.set("fu-skymap","dec-res","512")
		cp.set("fu-skymap","sample-rate","4096")

		# FU-OUTPUT SECTION
		cp.add_section("fu-output")
		cp.set("fu-output","log-path",self.log_path())
		cp.set("fu-output","output-dir",self.web_dir())
		cp.set("fu-output","web-url", self.web_url())

		# if we have an ini file override the options
		if config: 
			user_cp = ConfigParser.ConfigParser()
			user_cp.read(config)
		else:
			# otherwise see if a file with the standard ini file exists in the directory, the user probably intends to use it
			try: 
				user_cp = ConfigParser.ConfigParser()
				user_cp.read('stfu_pipe.ini')
			except: pass
		# override the default options
		if user_cp: self.overwrite_config(user_cp)

	def write(self):
		self.get_cp().write(open(self.time_now + ".ini","w"))

	def get_cp(self):
		return self.cp

	def set_qscan_executable(self):
		host = self.__get_hostname()
		if 'phy.syr.edu' in host:
			self.cp.set("fu-condor","qscan",self.__home_dirs()+"/rgouaty/opt/omega/omega_r2062_glnxa64_binary/bin/wpipeline")
		else:
			self.cp.set("fu-condor","qscan",self.__home_dirs()+"/romain/opt/omega/omega_r2062_glnxa64_binary/bin/wpipeline")		

	def __qscan_config(self,config):
		#FIXME why isn't there an environment variable for things in lalapps share?
		path = self.which('lalapps_inspiral')
		if path: path = os.path.split(path)[0]
		else: 
			print >>sys.stderr, "COULD NOT FIND QSCAN CONFIG FILE %s IN %s, ABORTING" % (config, path)
			raise ValueError
			sys.exit(1)
		out = path.replace('bin','share/lalapps') + '/' + config
		if not os.path.isfile(out):
			print >>sys.stderr, "COULD NOT FIND QSCAN CONFIG FILE %s IN %s, ABORTING" % (config, out)
			raise ValueError
			sys.exit(1)
		return out

	def web_dir(self):
		host = self.__get_hostname()
		#FIXME add more hosts as you need them
		if 'caltech.edu' in host: return os.environ['HOME'] + '/public_html/followups/' + self.time_now
		if 'phys.uwm.edu' in host: return os.environ['HOME'] + '/public_html/followups/' + self.time_now
		if 'phy.syr.edu' in host: return os.environ['HOME'] + '/public_html/followups/' + self.time_now
		if 'aei.uni-hannover.de' in host: return os.environ['HOME'] + '/WWW/LSC/' + self.time_now
		print sys.stderr, "WARNING: could not find web directory, returning empty string"
		return ''

	def web_url(self):
		host = self.__get_hostname()
		#FIXME add more hosts as you need them
		if 'ligo.caltech.edu' in host: return "https://ldas-jobs.ligo.caltech.edu/~" +os.environ['USER'] + '/followups/' + self.time_now
		if 'ligo-la.caltech.edu' in host: return "https://ldas-jobs.ligo-la.caltech.edu/~" +os.environ['USER'] + '/followups/' + self.time_now
		if 'ligo-wa.caltech.edu' in host: return "https://ldas-jobs.ligo-wa.caltech.edu/~" +os.environ['USER'] + '/followups/' + self.time_now
		if 'phys.uwm.edu' in host: return "https://ldas-jobs.phys.uwm.edu/~" + os.environ['USER'] + '/followups/' + self.time_now
		if 'phy.syr.edu' in host: return "https://sugar-jobs.phy.syr.edu/~" + os.environ['USER'] + '/followups/' + self.time_now
		if 'aei.uni-hannover.de' in host: return "https://atlas.atlas.aei.uni-hannover.de/~" + os.environ['USER'] + '/LSC/' + self.time_now
		print sys.stderr, "WARNING: could not find web server, returning empty string"
		return ''

	def log_path(self):
		host = self.__get_hostname()
		#FIXME add more hosts as you need them
		if 'caltech.edu' in host: return '/usr1/' + os.environ['USER']
		if 'phys.uwm.edu' in host: return '/people/' + os.environ['USER']
		if 'aei.uni-hannover.de' in host: return '/local/user/' + os.environ['USER']

	def __get_hostname(self):
		host = socket.getfqdn()
		return host
		
	def __home_dirs(self):
		return os.path.split(os.environ['HOME'])[0]
		
	def which(self,prog):
		which = subprocess.Popen(['which',prog], stdout=subprocess.PIPE)
		out = which.stdout.read().strip()
		if not out: print >>sys.stderr, "WARNING: could not find %s in your path, unless you have an ini file to overide the path to %s the DAG will fail" % (prog,prog)
		return out

	def overwrite_config(self,config):
		for section in config.sections():
			if not cp.has_section(section): cp.add_section(section)
			for option in config.options(section):
				cp.set(section,option,config.get(section,option))

