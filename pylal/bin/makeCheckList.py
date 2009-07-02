#!/usr/bin/env python
"""
Something

$Id$

This program generates a detection checklist for a candidate.
"""

__author__ = 'Romain Gouaty <gouaty@lapp.in2p3.fr>'
__date__ = '$Date$'
__version__ = '$Revision$'[11:-2]
__prog__ = "makeChecklist"

##############################################################################
# import standard modules and append the lalapps prefix to the python path
import sys, os, copy, math, random
import socket, time
import re, string
import commands
from optparse import *
import tempfile
import ConfigParser
import urlparse
import urllib
import numpy
import fnmatch
from UserDict import UserDict
sys.path.append('@PYTHONLIBDIR@')

##############################################################################
# import the modules we need to build the pipeline
from glue import lal
from glue import segments
from glue import segmentsUtils
from glue.ligolw import ligolw
from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw import utils
from glue import markup
from pylal import InspiralUtils
from pylal import CoincInspiralUtils
from pylal import fu_utils
from pylal.fu_writeXMLparams import *
from pylal import Fr
from pylal.scrapeHtmlUtils import scrapePage
from lalapps import inspiral


##############################################################################
# A few new methods 

def getFileMatchingTrigger(jobname,string_id,file_list=False):
  if os.access(jobname,os.F_OK):
    filesInDir = os.listdir(jobname)
    fileList = []
    for paramFile in filesInDir:
      if fnmatch.fnmatch(paramFile, "*"+string_id+"*.html"):
        if file_list:
          fileList.append("../"+jobname+"/"+paramFile)
        else:
          return "../"+jobname+"/"+paramFile
    if fileList:
      return fileList
    else:
      return False
  else: return False

######################## OPTION PARSING  #####################################
usage = """usage: %prog [options]
"""

parser = OptionParser( usage )

parser.add_option("-v", "--version",action="store_true",default=False,\
    help="print version information and exit")

parser.add_option("","--trigger-id",action="store",type="string",\
    metavar=" STRING",help="eventid of the analysed trigger")

parser.add_option("","--trigger-gps",action="store",type="string",\
    metavar=" STRING",help="list of gps times for each ifo where the \
coincidence is found. The gps times must be separated by a coma, for example \
trigger-gps=\"860308882.71533203,860308882.74438477\"")

parser.add_option("","--hoft-channel-ref",action="store",type="string",\
    metavar=" STRING",default="DMT-STRAIN,h_16384Hz",help="this indicates"\
    " the hoft channel name used by LIGO and Virgo ifos")

parser.add_option("","--ifolist-in-coinc",action="store",type="string",\
    metavar=" STRING",help="string cointaing the ifo names found in coincidence, for example: \"H1H2L1\"")

parser.add_option("-o","--output-path",action="store",type="string",\
    default="", metavar=" PATH",\
    help="path where the figures would be stored")

parser.add_option("-O","--enable-output",action="store_true",\
    default="false",  metavar="OUTPUT",\
    help="enable the generation of the html and cache documents")

parser.add_option("-T","--user-tag", action="store",type="string", \
    default=None, metavar=" USERTAG",help="user tag for the output file name")

parser.add_option("","--ifo-times",action="store",\
    type="string", default=None, metavar=" IFOTIMES",\
    help="provide ifo times for naming figure")

parser.add_option("","--ifo-tag",action="store",\
    type="string",  metavar=" IFOTAG",\
    help="ifo tag gives the information about ifo times and stage")

parser.add_option("-C","--cumulhisto-page",action="store",type="string",\
    metavar=" STRING",help="url to the cumulative histogram of combined statistics")

parser.add_option("-H","--histo-page",action="store",type="string",\
    metavar=" STRING",help="url to the histogram of combined statistics")

parser.add_option("-I","--ifar-page",action="store",type="string",\
    metavar=" STRING",help="url to the ifar plot")

parser.add_option("","--ifar-combined-page",action="store",type="string",\
    metavar=" STRING",help="url to the combined ifar plot")

parser.add_option("-X","--data-quality-url",action="store",type="string",\
                  metavar="URL",default=None, dest="defaultldbd",\
                  help="Using this argument specify a URL the LDBD \
server that you want to query DQ Veto segment information from for\
example ldbd://metaserver.phy.syr.edu:30015")

parser.add_option("-Q","--data-quality-database",action="store",type="string",\
    metavar=" PATH2FILE",default=None, dest="defaultSQL",\
    help="This is the disk location of\
the data quality sqlite database to use for DQ information queries.\
Omission of this option will cause a default search for \
~/followupDQ.sqlite rebuilding it if needed.")

parser.add_option("-R","--SNR-ratio-test",action="store",type="string",\
metavar=" PATH2FILE", default=None, dest="defaultRatioTestPickle", \
help="Set the location of the data (pickle) file used to perform the\
 ratio check on the candidate file.")

parser.add_option("","--remote-qscan-web",action="store",type="string",\
    metavar=" STRING",default=None,help="this option should be provided " \
    "only if an ifo is processed remotely (for instance V1 qscans). The " \
    "format should be an ifo string followed by a coma and a string " \
    "that contains the url of the remote qscan result page. For instance: " \
    "--remote-qscan-web V1,http://virgo.in2p3.fr/followups/romain/tests/lv" \
    "/local_followup_043009/871147814-873567014/CAT3/qscan/V1/")

parser.add_option("","--remote-seismic-qscan-web",action="store",type="string",\
    metavar=" STRING",default=None,help="this option should be provided " \
    "only if an ifo is processed remotely (for instance V1 qscans). The " \
    "format should be an ifo string followed by a coma and a string " \
    "that contains the url of the remote qscan result page. For instance: " \
    "--remote-seismic-qscan-web V1,http://virgo.in2p3.fr/followups/romain/" \
    "tests/lv/local_followup_043009/871147814-873567014/CAT3/seismic-qscan/V1/")

command_line = sys.argv[1:]
(opts,args) = parser.parse_args()

#################################
# if --version flagged
if opts.version:
  print "$Id: generate_checklist.py, v 1.0 2008/05/20 07:00:00 romain Exp"
  sys.exit(0)

##############################################################################
# main program

# List of ifos (used to check for Nelson KW vetoes)
#ifoList = ['H1','H2','L1']

opts = InspiralUtils.initialise(opts,__prog__,__version__)

page = markup.page(mode="strict_html")
page._escape = False
doctype="""<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">"""
doctype+="""\n<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">"""

title = "Detection Checklist for candidate " + str(opts.trigger_id)
page.init(title=title, doctype=doctype)
#page.init(title=title)
page.h1()
page.add("Detection Checklist for Candidate " + str(opts.trigger_id))
page.h1.close()

page.h2()
page.add("Inspiral triggers found by CBC search:")
page.h2.close()

# Check if PARAM_TABLES directory exists 
# and look for a file name matching opts.trigger_id
if os.access("PARAM_TABLES",os.F_OK):
  filesInDir = os.listdir("PARAM_TABLES")
  for paramFile in filesInDir:
    if fnmatch.fnmatch(paramFile, "table*" + opts.trigger_id + "*"):
      # Copy the table of parameters inside the checklist
      paramFile = os.path.normpath("PARAM_TABLES/" + paramFile)
      paramTable = scrapePage()
      paramTable.setContextKeys("<table bgcolor=cyan border=1px>","</table>")
      paramTable.readfile(paramFile)
      page.add(paramTable.buildTableHTML("border=1 bgcolor=yellow"))
      page.hr()
      break
    else: pass


#n_veto = nVeto()
nelsonVeto = []
dailyStat = []
hoft_qscan = []
rds_qscan = []
seis_qscan = []
analyse_hoft_qscan = []
analyse_rds_qscan = []
analyse_seismic_qscan = []
snrchisq = []
coherent_qscan = []
framecheck = []
chia = []
skymap = []
singlemcmc = []
coherentmcmc = []
fu_triggers = []

# get channel name for hoft data
LIGO_channel = opts.hoft_channel_ref.split(",")[0].strip()
Virgo_channel = opts.hoft_channel_ref.split(",")[1].strip()

# prepare strings containing information on Nelson's DQ investigations
#for ifo in ifoList:
#  if ifo in trig.ifolist_in_coinc:
#    nelsonVeto.append(n_veto.findInterval(ifo,trig.gpsTime[ifo]))
#  else:
#    nelsonVeto.append(n_veto.findInterval(ifo,trig.gpsTime[trig.ifolist_in_coinc[0]]))

ifolist = []
for j in range(0,len(opts.ifolist_in_coinc)-1,2):
  ifolist.append(opts.ifolist_in_coinc[j:j+2])

for ifo_index,ifo in enumerate(ifolist):
  gpstime = opts.trigger_gps.split(",")[ifo_index].strip()

  # links to daily stats
  dailyStatFile = getFileMatchingTrigger("IFOstatus_checkJob",ifo+"_"+opts.trigger_id)
  if dailyStatFile:
    dailyStat.append(dailyStatFile)
  else:
    dailyStat.append("")

  # links to qscans
  hoft_qscan.append("../QSCAN/foreground-hoft-qscan/" + ifo + "/" + gpstime)
  if opts.remote_qscan_web and ifo == opts.remote_qscan_web.split(",")[0]:
    rds_qscan.append(opts.remote_qscan_web.split(",")[1] + "/qscan/" + gpstime)
  else:
    rds_qscan.append("../QSCAN/foreground-qscan/" + ifo + "/" + gpstime)
  if opts.remote_seismic_qscan_web and ifo == opts.remote_seismic_qscan_web.split(",")[0]:
    seis_qscan.append(opts.remote_seismic_qscan_web.split(",")[1] + "/seismic_qscan/" + gpstime)
  else:
    seis_qscan.append("../QSCAN/foreground-seismic-qscan/" + ifo + "/" + gpstime)

  # links to analyse qscans
  analyseSeismicQscanFile = getFileMatchingTrigger("analyseQscanJob",ifo+"_"+gpstime.replace(".","_")+"_seismic_qscan")
  if analyseSeismicQscanFile:
    analyse_seismic_qscan.append(analyseSeismicQscanFile) 
  else:
    analyse_seismic_qscan.append("")

  analyseQscanFile = getFileMatchingTrigger("analyseQscanJob",ifo+"_"+gpstime.replace(".","_")+"_qscan")
  if analyseQscanFile:
    analyse_rds_qscan.append(analyseQscanFile)
  else:
    analyse_rds_qscan.append("")

  analyseHoftQscanFile = getFileMatchingTrigger("analyseQscanJob",ifo+"_"+gpstime.replace(".","_")+"_hoft_qscan")
  if analyseHoftQscanFile:
    analyse_hoft_qscan.append(analyseHoftQscanFile)
  else:
    analyse_hoft_qscan.append("")

  # links to snrchisq plots
  snrchisqFile = getFileMatchingTrigger("plotSNRCHISQJob",ifo+"_"+opts.trigger_id)
  if snrchisqFile:
    snrchisq.append(snrchisqFile)
  else:
    snrchisq.append("")

  # links to frame checks
  framecheckFile = getFileMatchingTrigger("FrCheckJob",ifo+"_"+opts.trigger_id)
  if framecheckFile:
    framecheck.append(framecheckFile)
  else:
    framecheck.append("")

  # links to single MCMC chains
  if os.access("plotmcmcJob",os.F_OK):
    filesInDir = os.listdir("plotmcmcJob")
    for element in filesInDir:
      if fnmatch.fnmatch(element, "*" + ifo + "-*" + opts.trigger_id):
        singlemcmc.append("../plotmcmcJob/" + element)
        break
      else: pass
  if len(singlemcmc) < ifo_index+1:
    singlemcmc.append("")

# loop over ifos not found in coincidence (though in the analysed times)
gpstime0 = opts.trigger_gps.split(",")[0].strip()
for j in range(0,len(opts.ifo_times)-1,2):
  ifo = opts.ifo_times[j:j+2]
  if not ifolist.count(ifo):
     # links to qscans
     hoft_qscan.append("../QSCAN/foreground-hoft-qscan/" + ifo + "/" + gpstime0)
     if opts.remote_qscan_web and ifo == opts.remote_qscan_web.split(",")[0]:
       rds_qscan.append(opts.remote_qscan_web.split(",")[1] + "/qscan/" + gpstime0)
     else:
       rds_qscan.append("../QSCAN/foreground-qscan/" + ifo + "/" + gpstime0)
     # links to snrchisq plots
     for ifo_ref in ifolist:
       snrchisqFile = getFileMatchingTrigger("plotSNRCHISQJob",ifo+"_"+ifo_ref+"tmplt_"+opts.trigger_id)
       if snrchisqFile:
         snrchisq.append(snrchisqFile)
       else:
         snrchisq.append("")

# link to coherent qscan
try:
  ifolist.index("H1")
  ifolist.index("H2")
  coherent_qscan.append("../h1h2QeventJob/qevent/H1H2/" + gpstime0)
except: pass

# links to coherent MCMC
if os.access("plotmcmcJob",os.F_OK):
  filesInDir = os.listdir("plotmcmcJob")
  for element in filesInDir:
    if fnmatch.fnmatch(element, "*" + opts.ifolist_in_coinc + "-*" + opts.trigger_id):
      coherentmcmc.append("../plotmcmcJob/" + element)
      break
    else: pass


# link to inspiral coherent followups
coherentInspiralFile = getFileMatchingTrigger("plotChiaJob",opts.ifo_times+"_"+opts.trigger_id)
if coherentInspiralFile:
  coherentParamFile = coherentInspiralFile.strip("../")
  coherentParamTable = scrapePage()
  coherentParamTable.setContextKeys("<table bgcolor=pink border=1px>","</table>")
  coherentParamTable.readfile(coherentParamFile)
  chia.append(coherentInspiralFile)
else:
  pass

# links to skymaps
skymapFile = getFileMatchingTrigger("pylal_skyPlotJob",opts.trigger_id)
if skymapFile:
  skymap.append(skymapFile)

# links to followup of triggers in the chunk 
followupTriggerFile = getFileMatchingTrigger("followUpTriggers",opts.trigger_id,True)
if followupTriggerFile:
  fu_triggers.append(followupTriggerFile)

# build the checklist table
page.h2()
page.add("Follow-up tests")
page.h2.close()

page.add("<table bgcolor=wheat border=1px>")
page.tr()
page.td("<b>ID</b>")
page.td("<b>Questions</b>")
page.td("<b>Answers</b>")
page.td("<b>Relevant information (flags, plots and links)</b>")
page.td("<b>Comments</b>")
page.tr.close()


#file.write("<tr bgcolor=red>\n")
page.add("<tr bgcolor=red>")
page.td()
page.td("<b>Is this candidate a possible gravitational-wave ?</b>")
page.td("<b>YES/NO</b>")
page.td("Main arguments:")
page.tr.close()


# Row #0
page.tr()
page.td("#0 False alarm probability")
page.td("What is the false alarm rate associated with this candidate ?")
page.td()

far_links = ""
if opts.cumulhisto_page:
  far_links += "<a href=\"" + opts.cumulhisto_page + "\">Cumulative histogram</a><br>\n"
if opts.histo_page:
  far_links += "<a href=\"" + opts.histo_page + "\">Non-cumulative histogram</a><br>\n"
if opts.ifar_page:
  far_links += "<a href=\"" + opts.ifar_page + "\">IFAR plot</a><br>\n"
if opts.ifar_combined_page:
  far_links += "<a href=\"" + opts.ifar_combined_page + "\">Combined IFAR plot</a><br>\n"
page.td(far_links)

page.td()
page.tr.close()

# Row #1
###############################Cristina Tue-Feb-10-2009:200902101541 
page.tr()
page.td("#1 DQ flags")
page.td("Can the data quality flags coincident with this candidate be safely disregarded ?")
page.td()
#Only use the SQL file way of selecting DQ IFF ini
#code called with flag --data-quality-database
if opts.defaultSQL != None:
  #Create database object without initialization
  #Assume there is a parser option setting the variable below?
  preBuiltDB=opts.defaultSQL
  if opts.defaultSQL == None:
    preBuiltDB=""
  if os.path.isfile(preBuiltDB):
    checklistDB=fu_utils.followupdqdb(True)
    checklistDB.setDB(preBuiltDB)
  else:
    checklistDB=fu_utils.followupdqdb()
  results=dict()
  results=checklistDB.queryDB(int(float(gpstime0)),600)
  checklistDB.close()
  htmlStart="<table bgcolor=grey border=1px><tr><th>IFO</th><th>Flag</th><th>Start(s)</th><th>Offset(s)</th><th>Stop(s)</th><th>Offset(s)</th></tr>"
  htmlRows=list()
  htmlRowsVeto=list()
  htmlStop="</table>"
  trigGPS=float(gpstime0)
  for ifo in checklistDB.ifoList:
    for segment in results[ifo]:
      token=segment
      if token.__len__() >= 5:
        rowText="<tr bgcolor=%s><td>%s</td><td>%s</td><td>%i</td><td>%i</td><td>%i</td><td>%i</td></tr>"
        if ((token[2]-trigGPS)<0 and (token[3]-trigGPS)<0):
          mycolor="yellow"
        elif ((token[2]-trigGPS)>0 and (token[3]-trigGPS)>0):
          mycolor="green"
        else:
          mycolor="red"
        if token[0].split("_")[0] != "VETO":
          htmlRows.append(rowText%(mycolor,ifo,\
                                   token[0],\
                                   token[2],\
                                   token[2]-trigGPS,\
                                   token[3],\
                                   token[3]-trigGPS))
        else:
          htmlRowsVeto.append(rowText%(mycolor,ifo,\
                                   token[0],\
                                   token[2],\
                                   token[2]-trigGPS,\
                                   token[3],\
                                   token[3]-trigGPS))
  htmlMiddle=""
  for row in htmlRows:
    htmlMiddle="%s %s"%(htmlMiddle,row)
  dqTable=htmlStart+htmlMiddle+htmlStop
  htmlMiddleVeto=""
  for row in htmlRowsVeto:
    htmlMiddleVeto="%s %s"%(htmlMiddleVeto,row)
  vetoTable=htmlStart+htmlMiddleVeto+htmlStop
else:
  dqTable=""
  vetoTable=""
  
if opts.defaultldbd != None:
  defaultServer="ldbd://metaserver.phy.syr.edu:30015"
  windowSize=int(600)
  versionNumber=int(1)
  x=followupDQV(defaultServer)
  x.fetchInformation(float(gpstime0),windowSize,versionNumber)
  dqTable=x.generateHTMLTable()
#
# Insert the new text string of a table using markup.py functions
page.td(dqTable)
page.td()
page.tr.close()
###############################Cristina Tue-Feb-10-2009:200902101541 


# Row #2
page.tr()
page.td("#2 Veto investigations")
page.td("Does the candidate survive the veto investigations performed at its time ?")
page.td()
page.td(vetoTable)
#file.write("  <td>")
#for j,ifo in enumerate(ifoList):
#  file.write("<table>")
#  file.write("\n<b>" + ifo + ":</b>\n")
#  file.write('<tr><td>' + nelsonVeto[j].strip("\n").replace("\n","</td></tr><tr><td>").replace(" ", "</td><td>") + '</td></tr>')
#  file.write("</table>")
#file.write("</td>\n")
page.td()
page.tr.close()

# Row #3
page.tr()
page.td("#3 Ifo status")
page.td("Are the interferometers operating normally with a reasonable level of sensitivity around the time of the candidate ?")
page.td()
ifoStatusLinks = "<a href=\"http://blue.ligo-wa.caltech.edu/scirun/S5/DailyStatistics/\">Daily Stats pages</a>:"
for j,ifo in enumerate(ifolist):
  ifoStatusLinks += " <a href=\"" + dailyStat[j] + "\">" + ifo + "</a>"
#file.write("\n" + ScSegTable.buildTableHTML("border=1 bgcolor=green").replace("\n","") + "<br>" + dateScSeg)
page.td(ifoStatusLinks)
page.td()
page.tr.close()


# Row #4
page.tr()
page.td("#4 Candidate appearance")
page.td("Do the Qscan figures show what we would expect for a gravitational-wave event ?")
page.td()
hoftQscanLinks = "h(t) Qscans:<br>"
for j,ifo in enumerate(ifolist):
  if ifo == "V1": strain = Virgo_channel
  else: strain = LIGO_channel
  gpstime = opts.trigger_gps.split(",")[j].strip()
  hoftQscanLinks += " <a href=\"" + hoft_qscan[j] + "\">" + ifo + "</a><br>"
  hoftQscanLinks += " <a href=\"" + analyse_hoft_qscan[j] + "\"> Background information for " + ifo + "</a>"
  hoftQscanLinks += " <img src=\"" + hoft_qscan[j] + "/" + gpstime + "_" + ifo + ":" + strain + "_1.00_spectrogram_whitened_thumbnail.png\" width=\"50%\">"
  hoftQscanLinks += " <img src=\"" + hoft_qscan[j] + "/" + gpstime + "_" + ifo + ":" + strain + "_16.00_spectrogram_whitened_thumbnail.png\" width=\"50%\">"
i=0
for k in range(0,len(opts.ifo_times)-1,2):
  ifo = opts.ifo_times[k:k+2]
  if ifo == "V1": strain = Virgo_channel
  else: strain = LIGO_channel
  if not ifolist.count(ifo):
    i=i+1
    hoftQscanLinks += " <a href=\"" + hoft_qscan[i + len(ifolist) - 1] + "\">" + ifo + "</a><br>"
    hoftQscanLinks += " <img src=\"" + hoft_qscan[i + len(ifolist) - 1] + "/" + gpstime0 + "_" + ifo + ":" + strain + "_1.00_spectrogram_whitened_thumbnail.png\" width=\"50%\"><br>"
    hoftQscanLinks += " <img src=\"" + hoft_qscan[i + len(ifolist) - 1] + "/" + gpstime0 + "_" + ifo + ":" + strain + "_16.00_spectrogram_whitened_thumbnail.png\" width=\"50%\"><br>"
page.td(hoftQscanLinks)
page.td()
page.tr.close()


# Row #5
page.tr()
page.td("#5 Seismic plots")
page.td("Is the seismic activity insignificant around the time of the candidate ?")
page.td()
seismicQscanLinks = "Seismic Qscans:"
for j,ifo in enumerate(ifolist):
  seismicQscanLinks += " <a href=\"" + seis_qscan[j] + "\">" + ifo + "</a>"
seismicQscanLinks += "<br>Background information on qscans:"
for j,ifo in enumerate(ifolist):
  seismicQscanLinks += " <a href=\"" + analyse_seismic_qscan[j] + "\">" + ifo + "</a>"
page.td(seismicQscanLinks)
page.td()
page.tr.close()


# Row #6
page.tr()
page.td("#6 Other environmental causes")
page.td("Were the environmental disturbances (other than seismic) insignificant at the time of the candidate ?")
page.td()
qscanLinks = "RDS Qscans:"
for j,ifo in enumerate(ifolist):
  qscanLinks += " <a href=\"" + rds_qscan[j] + "\">" + ifo + "</a>"
i=0
for k in range(0,len(opts.ifo_times)-1,2):
  ifo = opts.ifo_times[k:k+2]
  if not ifolist.count(ifo):
    i=i+1
    qscanLinks += " <a href=\"" + rds_qscan[i + len(ifolist) - 1] + "\">" + ifo + "</a>"
qscanLinks += "<br>Background information on qscans:"
for j,ifo in enumerate(ifolist):
  qscanLinks += " <a href=\"" + analyse_rds_qscan[j] + "\">" + ifo + "</a>"
page.td(qscanLinks)
page.td()
page.tr.close()

# Row #7
page.tr()
page.td("#7 Auxiliary degree of freedom")
page.td("Were the auxiliary channel transients coincident with the candidate insignificant ?")
page.td()
qscanLinks = "RDS Qscans:"
for j,ifo in enumerate(ifolist):
  qscanLinks += " <a href=\"" + rds_qscan[j] + "\">" + ifo + "</a>"
i=0
for k in range(0,len(opts.ifo_times)-1,2):
  ifo = opts.ifo_times[k:k+2]
  if not ifolist.count(ifo):
    i=i+1
    qscanLinks += " <a href=\"" + rds_qscan[i + len(ifolist) - 1] + "\">" + ifo + "</a>"
qscanLinks += "<br>Background information on qscans:"
for j,ifo in enumerate(ifolist):
  qscanLinks += " <a href=\"" + analyse_rds_qscan[j] + "\">" + ifo + "</a>"
page.td()
page.tr.close()


# Row #8
page.tr()
page.td("#8 Elog")
page.td("Were the instruments behaving normally according to the comments posted by the sci-mons or the operators in the e-log ?")
page.td()
elogLinks = "<a href=\"http://ilog.ligo-wa.caltech.edu/ilog/pub/ilog.cgi?group=detector\">Hanford elog</a><br>\n"
elogLinks += "<a href=\"http://ilog.ligo-la.caltech.edu/ilog/pub/ilog.cgi?group=detector\">Livingston elog</a>"
page.td(elogLinks)
page.td()
page.tr.close()


# Row #9
page.tr()
page.td("#9 Glitch report")
page.td("Were the instruments behaving normally according to the weekly glitch report ?")
page.td()
page.td("<a href=\"http://www.lsc-group.phys.uwm.edu/glitch/investigations/s5index.html#shift\">Glitch reports</a><br>")
page.td()
page.tr.close()

# Row #10
page.tr()
page.td("#10 Snr versus time")
page.td("Is this trigger significant in a SNR versus time plot of all triggers in its analysis chunk ?")
page.td()
fuTriggerLinks = ""
if fu_triggers:
  for fu_link in fu_triggers[0]:
    window = fu_link.split("-")[-1].strip(".html")
    fuTriggerLinks += " <a href=\"" + fu_link + "\">Triggers versus time (" + window + ")</a><br>"
page.td(fuTriggerLinks)
page.td()
page.tr.close()

# Row #11
######################
#Code to perform test
resultString=(" <table border=1px>\
 <tr><th>IFO:IFO</th><th>ToF</th><th>Deff Ratio</th><th>Prob</th><th>Figure</th></tr>")
#Text insert into page giving the SNR ratio probabilities
preBuiltPickle=opts.defaultRatioTestPickle
if opts.defaultRatioTestPickle == None:
  preBuiltPickle=""
ratioTest=fu_utils.ratioTest()
if os.path.isfile(preBuiltPickle):
  ratioTest.setPickleLocation(preBuiltPickle)
for index1,ifo1 in enumerate(ifolist):
  for index2,ifo2 in enumerate(ifolist):
    ifoA=ratioTest.mapToObservatory(ifo1)
    ifoB=ratioTest.mapToObservatory(ifo2)
    if ifoA != ifoB:
      gpsA=numpy.float64(opts.trigger_gps.split(",")[index1].strip())
      gpsB=numpy.float64(opts.trigger_gps.split(",")[index2].strip())
      snrA=float(str(paramTable.getColumnByText(ifo1,9)).strip().strip("<td>").strip("</td>"))
      snrB=float(paramTable.getColumnByText(ifo2,9).strip().strip("<td>").strip("</td>"))
      try:
        snrRatio=snrA/snrB
      except:
        snrRatio=0
      gpsDiff=gpsA-gpsB
      result=ratioTest.testRatio(ifoA,ifoB,gpsDiff,snrRatio)
      pairURL=ratioTest.findURL(ifoA,ifoB)
      myURL=str('<a href="%s"><img height=150px src="%s"></a>'%(pairURL,pairURL))
      myString="<tr><td>%s:%s</td><td>%2.4f</td><td>%5.2f</td><td>%1.3f</td><td>%s</td></tr>"%\
          (ifoA,ifoB,gpsDiff,snrRatio,result,myURL)
      resultString="%s %s"%(resultString,myString)
imageURL='<a href="https://ldas-jobs.ligo.caltech.edu/~ctorres/DQstuff/delayRatio_090504.png"><img height=200px src="https://ldas-jobs.ligo.caltech.edu/~ctorres/DQstuff/delayRatio_090504.png"></a>'
resultString=" %s </table>"%(resultString)
##############
page.tr()
page.td("#11 Parameters of the candidate")
page.td("Does the candidate have a high likelihood of being a gravitational-wave according to its parameters ?")
page.td()
page.td(resultString)
page.td()
page.tr.close()

# Row #12
page.tr()
page.td("#12 Snr and Chisq")
page.td("Are the SNR and CHISQ time series consistent with our expectations for a gravitational wave ?")
page.td()
snrchisqLinks = ""
for j,ifo in enumerate(ifolist):
  snrchisqLinks += " <a href=\"" + snrchisq[j] + "\">" + ifo + "</a>"
snrchisqLinks += "<br>\n"
i=0
for k in range(0,len(opts.ifo_times)-1,2):
  ifo = opts.ifo_times[k:k+2]
  if not ifolist.count(ifo):
    for ifo_ref in ifolist:
      i=i+1
      snrchisqLinks += " <a href=\"" + snrchisq[i + len(ifolist) - 1] + "\">" + ifo + " with " + ifo_ref + " template" + "</a>"
page.td(snrchisqLinks)
page.td()
page.tr.close()

# Row #13
page.tr()
page.td("#13 Template bank veto")
page.td("Is the bank veto value consistent with our expectations for a gravitational wave ?")
page.td()
page.td()
page.td()
page.tr.close()

# Row #14
page.tr()
page.td("#14 Coherent studies")
page.td("Are the triggers found in multiple interferometers coherent with each other ?")
page.td()
coherentLinks = ""
if coherent_qscan:
  coherentLinks += " <a href=\"" + coherent_qscan[0] + "\">H1H2 coherent qevent</a><hr />"
if chia:
  coherentLinks += coherentParamTable.buildTableHTML("border=1 bgcolor=yellow")
  coherentLinks += " <a href=\"" + chia[0] + "\">Coherent inspiral followup</a><hr />"
if skymap:
  coherentLinks += " <a href=\"" + skymap[0] + "\">Sky Map</a><hr />"
page.td(coherentLinks)
page.td()
page.tr.close()

# Row #15
page.tr()
page.td("#15")
page.td("Is the candidate stable against changes in segmentation ?")
page.td()
page.td()
page.td()
page.tr.close()

# Row #16
page.tr()
page.td("#16")
page.td("Is the candidate stable against changes in calibration that are consistent with systematic uncertainties ?")
page.td()
page.td()
page.td()
page.tr.close()

# Row #17
page.tr()
page.td("#17")
page.td("Is the data used in the analysis free from corruption at the time of the candidate ?")
page.td()
frameCheckLinks = "Frame checks: "
for j,ifo in enumerate(ifolist):
  frameCheckLinks += " <a href=\"" + framecheck[j] + "\">" + ifo + "</a>"
page.td(frameCheckLinks)
page.td()
page.tr.close()

page.table.close()
page.hr()


# Write parameter estimation table
page.h2("Parameter estimation")
page.add("<table bgcolor=chartreuse border=1px>")
page.tr()
page.td("ID")
page.td("Questions")
page.td("Answers")
page.td("Relevant information (flags, plots and links)")
page.td("Comments")
page.tr.close()

# Row #1
page.tr()
page.td("#1 Parameters of the candidate")
page.td("Can we get more accurate information on the parameters of this candidate using MCMC or Bayesian methods ?")
page.td()
singlemcmcLinks = "Single MCMC:"
for j,ifo in enumerate(ifolist):
  singlemcmcLinks += " <a href=\"" + singlemcmc[j] + "\">" + ifo + "</a>"
if coherentmcmc:
  coherentmcmcLink = "<br><br>Coherent MCMC:"
  coherentmcmcLink += " <a href=\"" + coherentmcmc[0] + "\">" + opts.ifolist_in_coinc + "</a>"
else:
  coherentmcmcLink = ""
page.td(singlemcmcLinks+coherentmcmcLink)
page.td()
page.tr.close()

# Row #2
page.tr()
page.td("#2 Coherent follow-up")
page.td("Make a followup with coherent multi-detector code.")
page.td()
coherentLinks = ""
if skymap:
  coherentLinks += " <a href=\"" + skymap[0] + "\">Sky Map</a><hr />"
page.td(coherentLinks)
page.td()
page.tr.close()

# Row #3
page.tr()
page.td("#3")
page.td("Are the results of the Burst analysis astrophysically consistent with a possible detection ?")
page.td()
page.td()
page.td()
page.tr.close()

# Row #4
page.tr()
page.td("#4")
page.td("Are the results of a ringdown search astrophisycally consistent with a possible detection ?")
page.td()
page.td()
page.td()
page.tr.close()

# Row #5
page.tr()
page.td("#5 EM triggers")
page.td("Are there any EM triggers in coincidence with the candidate ?<br>Is the distance estimated from interferometer time-delays or coherent analysis consistent with electro-magnetic observations?<br>? Are the distances as measured in several instruments consistent with position information?<br>")
page.td()
page.td("<a href=\"http://www.uoregon.edu/~ileonor/ligo/s5/grb/online/S5grbs_list.html\">List of GRBs during S5</a><br><a href=\"http://ldas-jobs.ligo.caltech.edu/~dietz/pages/s5/GRB/CVS/overviewS5.html\">CBC compiled list of GRBs</a>")
page.td()
page.tr.close()

page.table.close()
page.h2()
page.add("Follow up documentation")
page.add("Miscellaneous Information")
page.h2.close()
timeString=str(time.gmtime()[0:6]).replace(" ","").replace(",","-")
page.add("Checklist compiled:%s <br>\n"%(timeString))
page.add("<a\
 href=\"https://ldas-jobs.ligo.caltech.edu/~ctorres/followUpLivingDoc_LAST.pdf\">Living\
 follow up document</a>")

if opts.enable_output:
  if not os.access(opts.output_path,os.F_OK):
    os.mkdir(opts.output_path)
  html_filename = opts.output_path + opts.prefix + opts.suffix + ".html"
  html_file = file(html_filename, "w")
  html_file.write(page(False))
  html_file.close()

