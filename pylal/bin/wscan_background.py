#!/usr/bin/env python

__version__ = "$Revision$"
__date__ = "$Date$"
__prog__ = "wscan_background.py"
__Id__ = "$Id$"
__title__ = "Generate bakground of omega scans"

##############################################################################

import os, sys
from optparse import *
import ConfigParser
import time

from glue import pipeline
from glue import gpstime
from glue.ligolw import dbtables
from pylal import fu_utils
from pylal import date
from pylal import stfu_pipe

from pylal import xlal

##############################################################################
# Useful methods

def create_default_config(home_base):
    cp = ConfigParser.ConfigParser()

    cp.add_section("condor")
    cp.set("condor","datafind","ligo_data_find")

    cp.add_section("fu-condor")
    cp.set("fu-condor","convertcache","convertlalcache.pl")
    cp.set("fu-condor","qscan",home_base+"/cbc/opt/omega/omega_r2062_glnxa64_binary/bin/wpipeline")
    cp.set("fu-condor","query_dq","pylal_query_dq")

    cp.add_section("datafind")

    cp.add_section("fu-q-rds-datafind")
    cp.set("fu-q-rds-datafind","search-time-range","1024")
    cp.set("fu-q-rds-datafind","remote-ifo","V1")

    cp.add_section("fu-q-hoft-datafind")
    cp.set("fu-q-hoft-datafind","search-time-range","128")

    cp.add_section("followup-background-qscan-times")
    cp.set("followup-background-qscan-times","H1range","")
    cp.set("followup-background-qscan-times","L1range","")
    cp.set("followup-background-qscan-times","V1range","")
    cp.set("followup-background-qscan-times","segment-min-len","2048")
    cp.set("followup-background-qscan-times","segment-pading","64")
    cp.set("followup-background-qscan-times","random-seed","1")
    cp.set("followup-background-qscan-times","background-statistics","20")

    cp.add_section("fu-bg-rds-qscan")
    cp.set("fu-bg-rds-qscan","L1config-file",home_base+"/cbc/FOLLOWUP_QSCAN_STUFF_S6/wscan/configurations/background/standard_configuration/L0L1-RDS_R_L1-cbc.txt")
    cp.set("fu-bg-rds-qscan","L1config-file",home_base+"/cbc/FOLLOWUP_QSCAN_STUFF_S6/wscan/configurations/background/standard_configuration/H0H1-RDS_R_L1-cbc.txt")
    cp.set("fu-bg-rds-qscan","remote-ifo","V1")

    cp.add_section("fu-bg-ht-qscan")
    cp.set("fu-bg-ht-qscan","L1config-file",home_base+"/cbc/FOLLOWUP_QSCAN_STUFF_S6/wscan/configurations/background/hoft_configuration/L1_hoft_cbc.txt")
    cp.set("fu-bg-ht-qscan","H1config-file",home_base+"/cbc/FOLLOWUP_QSCAN_STUFF_S6/wscan/configurations/background/hoft_configuration/H1_hoft_cbc.txt")
    cp.set("fu-bg-ht-qscan","V1config-file",home_base+"/cbc/FOLLOWUP_QSCAN_STUFF_S6/wscan/configurations/background/hoft_configuration/L1_hoft_cbc.txt")

    cp.add_section("fu-bg-seismic-qscan")
    cp.set("fu-bg-seismic-qscan","L1config-file",home_base+"/cbc/FOLLOWUP_QSCAN_STUFF_S6/wscan/configurations/background/seismic_configuration/L0-RDS_R_L1-seismic-cbc.txt")
    cp.set("fu-bg-seismic-qscan","H1config-file",home_base+"/cbc/FOLLOWUP_QSCAN_STUFF_S6/wscan/configurations/background/seismic_configuration/H0-RDS_R_L1-seismic-cbc.txt")
    cp.set("fu-bg-seismic-qscan","remote-ifo","V1")

    cp.add_section("fu-output")
    cp.set("fu-output","log-path","/usr1/" + os.getenv("USER"))

    return cp

def overwrite_config(cp,config):
  for section in config.sections():
    if not cp.has_section(section): cp.add_section(section)
    for option in config.options(section):
      cp.set(section,option,config.get(section,option))
  return cp

def get_times():

  # determine the start time : 00:00:00 UTC from the day before
  # and the end time, 00:00:00 UTC the current day

  gps = xlal.date.LIGOTimeGPS(gpstime.GpsSecondsFromPyUTC(time.time()))
  end_gps = int(date.utc_midnight(gps))
  start_gps = end_gps - 86400

  print "Start time : "+str(start_gps)+"   End Time : "+str(end_gps)
  return str(start_gps)+","+str(end_gps)


##############################################################################
#MAIN PROGRAM
##############################################################################

home_dir = os.getenv("HOME")
home_base = "/".join(home_dir.split("/")[0:-1])

######################## OPTION PARSING  #####################################
usage = """usage: %prog [options]
"""

parser = OptionParser( usage )

parser.add_option("-v", "--version",action="store_true",default=False,\
    help="print version information and exit")

parser.add_option("-f","--config-file",action="store",type="string",\
    default="",help="configuration file is optional")

parser.add_option("-m", "--datafind",action="store_true",\
    default=False, help="use datafind to get qscan/trends data")

parser.add_option("-M", "--hoft-datafind",action="store_true",\
    default=False, help="use datafind to get hoft data (for qscan)")

parser.add_option("-Q", "--background-qscan",action="store_true",\
    default=False, help="do qscans over a list of times")

parser.add_option("-N", "--background-hoft-qscan",action="store_true",\
    default=False, help="do hoft qscans over a list of times")

parser.add_option("-S", "--background-seis-qscan",action="store_true",\
    default=False, help="do seismic qscans over a list of times")

command_line = sys.argv[1:]
(opts,args) = parser.parse_args()

if opts.version:
  print "$Id$"
  sys.exit(0)

#############################################################################

cp = create_default_config(home_base)
if opts.config_file: 
 config = ConfigParser.ConfigParser()
 config.read(opts.config_file)
 cp = overwrite_config(cp,config)

ifos_list = ['H1','H2','L1','G1','V1','T1']

#Initialize dag
if opts.config_file:
  dag = stfu_pipe.followUpDAG(opts.config_file,cp)
else:
  dag = stfu_pipe.followUpDAG("wscan_background.ini",cp)

# CONDOR JOB CLASSES
dataJob         = stfu_pipe.htDataFindJob(cp,'qdatafind')
qscanBgJob      = stfu_pipe.qscanJob(opts,cp,'QSCANLITE')

ifo_range = get_times()

for ifo in ifos_list:

    if cp.has_option("followup-background-qscan-times",ifo+"range") and not cp.get("followup-background-qscan-times",ifo+"range"):
      cp.set("followup-background-qscan-times",ifo+"range",ifo_range)

    # FIX ME: input argument segFile is not needed any more
    segFile = {}
    times, timeListFile = fu_utils.getQscanBackgroundTimes(cp,opts,ifo,segFile)

    for qtime in times:
      # SETUP DATAFIND JOBS FOR BACKGROUND QSCANS (REGULAR DATA SET)
      dNode = stfu_pipe.fuDataFindNode(dag,dataJob,cp,opts,ifo,sngl=None,qscan=True,trigger_time=qtime,data_type='rds')

      # SETUP DATAFIND JOBS FOR BACKGROUND QSCANS (HOFT)
      dHoftNode = stfu_pipe.fuDataFindNode(dag,dataJob,cp,opts,ifo,sngl=None,qscan=True,trigger_time=qtime)

      # SETUP BACKGROUND QSCAN JOBS
      qBgNode = stfu_pipe.fuQscanNode(dag,qscanBgJob,cp,opts,qtime,ifo,p_nodes=[dHoftNode],type="ht",variety="bg")

      qBgNode = stfu_pipe.fuQscanNode(dag,qscanBgJob,cp,opts,qtime,ifo,p_nodes=[dNode],type="rds",variety="bg")

      qBgNode = stfu_pipe.fuQscanNode(dag,qscanBgJob,cp,opts,qtime,ifo,p_nodes=[dNode],type="seismic",variety="bg")

#### ALL FINNISH ####
dag.write_sub_files()
dag.write_dag()
dag.write_script()


