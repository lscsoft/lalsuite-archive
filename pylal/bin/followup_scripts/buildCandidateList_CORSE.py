#!/usr/bin/env python
"""
Something

$Id$

This program generates an html page containing the list of candidates to be followed-up. There is one html page for each mass bin, each ifo time, and each coincidence type.
"""

__author__ = 'Romain Gouaty <romain@phys.lsu.edu>'
__date__ = '$Date$'
__version__ = '$Revision$'[11:-2]

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
from pylal import CoincInspiralUtils
from pylal.fu_utils import *
from pylal.fu_writeXMLparams import *
from pylal.webUtils import *
from pylal import Fr
from lalapps import inspiral

##############################################################################
#
#  MAIN PROGRAM
#
##############################################################################

######################## OPTION PARSING  #####################################
usage = """usage: %prog [options]
"""

parser = OptionParser( usage )

parser.add_option("-v", "--version",action="store_true",default=False,\
    help="print version information and exit")

parser.add_option("-l", "--log-path",action="store",type="string",\
    metavar=" PATH",help="directory to write condor log file")

parser.add_option("-p","--page",action="store",type="string",\
    metavar=" PATH",help="use web link PATH")

parser.add_option("","--web-output-file",action="store",type="string",\
    metavar=" PATH",help="use web output PATH")

parser.add_option("-o","--output-path",action="store",type="string",\
    metavar=" PATH",help="use output path PATH")

parser.add_option("","--header",action="store",type="string",\
    metavar=" STRING",help="This string will be the title of the web page")

parser.add_option("-s","--statistic",action="store",type="string",\
    metavar=" STRING",help="use statistic STRING to sort triggers (ex: effective_snrsq)")

parser.add_option("-a","--bla",action="store",type="float",\
    metavar=" VALUE",help="bitten l a")

parser.add_option("-b","--blb",action="store",type="float",\
    metavar=" VALUE",help="bitten l b")

parser.add_option("-n","--num-triggers",action="store",type="int",\
    metavar=" VALUE",help="number of triggers to followup")

parser.add_option("","--far-normalization",action="store",type="int",\
    metavar=" VALUE",help="constant used to normalize the FAR, to get the combined FAR. If there are three mass bins and 3 coincidence types (for candidates in triple times), this number would then be equal to 9. If there are three mass bins and 1 coincidence type (for candidates in double times), this number would be equal to 3.")


parser.add_option("","--zerolag-cat1",action="store",type="string",\
    default=False,metavar=" PATH",help="input path PATH for cat1 zerolag file")

parser.add_option("","--zerolag-cat12",action="store",type="string",\
    default=False,metavar=" PATH",help="input path PATH for cat12 zerolag file")

parser.add_option("","--zerolag-cat123",action="store",type="string",\
    default=False,metavar=" PATH",help="input path PATH for cat123 zerolag file")

parser.add_option("","--zerolag-cat1234",action="store",type="string",\
    default=False,metavar=" PATH",help="input path PATH for cat1234 zerolag file")


command_line = sys.argv[1:]
(opts,args) = parser.parse_args()


#################################
# if --version flagged
if opts.version:
  print "$Id$"
  sys.exit(0)

#################################
# Sanity check of input arguments

# To be done

#################################
# Here starts the real stuff...


stat = opts.statistic
if stat == "effective_snrsq":
  stat = "effective_snr"

# The purpose of the following lines of code is to read the xml trigger files that contain the lists of candidates. This program should be run once for each of the mass bins, each of the ifo times, and each of the coincidence types. There are four list of times because there are four DQ vetoes categories.
found1, coincs1, search1 = readFiles(opts.zerolag_cat1,getstatistic(opts.statistic,opts.bla,opts.blb))
found12, coincs12, search12 = readFiles(opts.zerolag_cat12,getstatistic(opts.statistic,opts.bla,opts.blb))
found123, coincs123, search123 = readFiles(opts.zerolag_cat123,getstatistic(opts.statistic,opts.bla,opts.blb))
found1234, coincs1234, search1234 = readFiles(opts.zerolag_cat1234,getstatistic(opts.statistic,opts.bla,opts.blb))

# Now we need to get lists of candidates sorted according to their combined effective SNR. We will need them later.
followuptrigs1 = getfollowuptrigs(str(opts.num_triggers),None,coincs1,None,None,None)
followuptrigs12 = getfollowuptrigs(str(opts.num_triggers),None,coincs12,None,None,None)
followuptrigs123 = getfollowuptrigs(str(opts.num_triggers),None,coincs123,None,None,None)
followuptrigs1234 = getfollowuptrigs(str(opts.num_triggers),None,coincs1234,None,None,None)

followuptrigsList = []
followuptrigsList.append(followuptrigs1)
followuptrigsList.append(followuptrigs12)
followuptrigsList.append(followuptrigs123)
followuptrigsList.append(followuptrigs1234)

# Now that we have the four list of triggers for the four DQ vetoes categories, we will merge all the lists into a single one. And we will also sort it according to the combined effective SNR of the candidates.
coincs = coincs1 
coincs.extend(coincs12)
coincs.extend(coincs123)
coincs.extend(coincs1234)

# Notice that when calling the method "getfollowuptrigs", the number of triggers toanalyse is set to "opts.num_triggers*4" instead of "opts.num_triggers". Because we have merged all the veto categories together, we might have replicated trigger in this list. Therefore we need take more than "opts.num_triggers" from this list in order to make sure we will get "opts.num_triggers" different triggers... The fact 4 comes from the fact that there are 4 veto categories, so that each trigger should not be replicated more than 4 times.
followuptrigs = getfollowuptrigs(str(opts.num_triggers*4),None,coincs,None,None,None)

webpage = WebPage(opts.header,opts.output_path + opts.web_output_file, opts.page)

trig_prev = None

# Now we will loop on all the triggers contained in the list "followuptrigs"
iteration = 0
for trig in followuptrigs:
  
  if iteration == opts.num_triggers:
    break

  #First let's make sure that this trigger is not duplicated several times in the list by checking that it is a different trigger than the one in the previous iteration (we compare the end_time and the end_time_ns to check for the identity between triggers). "trig_prev" is defined at the end of the loop. If it happens to be the same trigger as in the previous iteration we will simply ignore it.
  if trig_prev:
    if getattr(trig.coincs,trig.ifolist_in_coinc[0]).end_time == getattr(trig_prev.coincs,trig_prev.ifolist_in_coinc[0]).end_time and getattr(trig.coincs,trig.ifolist_in_coinc[0]).end_time_ns == getattr(trig_prev.coincs,trig_prev.ifolist_in_coinc[0]).end_time_ns:
      continue

  iteration+=1

  # get the GPS time of the trigger (we will use the time of the trigger in the first ifo appearing in "ifolist_in_coinc" as reference)
  gpsTime = trig.gpsTime[trig.ifolist_in_coinc[0]]
  # get the number of ifos involved in this coincidence
  n_ifos = len(trig.ifolist_in_coinc)
  
  # add a new candidate to the web page
  webpage.appendSection("candidate " + str(gpsTime) + " / " + "statistic = " + str(trig.statValue**2))
  #define the tables for this candidate
  webpage.lastSection.appendTable(n_ifos + 1,9,1)
  webpage.lastSection.appendTable(2,4,1)

  #Fill the first row of the first table
  webpage.lastSection.table[0].row[0].cell[0].text("IFO")
  webpage.lastSection.table[0].row[0].cell[1].text("End Time")
  webpage.lastSection.table[0].row[0].cell[2].text("SNR")
  webpage.lastSection.table[0].row[0].cell[3].text("CHISQ")
  webpage.lastSection.table[0].row[0].cell[4].text("Chirp Mass")
  webpage.lastSection.table[0].row[0].cell[5].text("Eta")
  webpage.lastSection.table[0].row[0].cell[6].text("mass1")
  webpage.lastSection.table[0].row[0].cell[7].text("mass2")
  webpage.lastSection.table[0].row[0].cell[8].text("Eff Dist (Mpc)")

  #Fill the first row of the second table
  webpage.lastSection.table[1].row[0].cell[0].text("FAR after cat1")
  webpage.lastSection.table[1].row[0].cell[1].text("FAR after cat12")
  webpage.lastSection.table[1].row[0].cell[2].text("FAR after cat123")
  webpage.lastSection.table[1].row[0].cell[3].text("FAR after cat1234")

  #Fill the first table
  i = 0
  for ifo in trig.ifolist_in_coinc:
    i = i+1
    webpage.lastSection.table[0].row[i].cell[0].text(ifo)
    webpage.lastSection.table[0].row[i].cell[1].text(repr(trig.gpsTime[ifo]))
    webpage.lastSection.table[0].row[i].cell[2].text(str(getattr(trig.coincs,ifo).snr))
    webpage.lastSection.table[0].row[i].cell[3].text(str(getattr(trig.coincs,ifo).chisq))
    webpage.lastSection.table[0].row[i].cell[4].text(str(getattr(trig.coincs,ifo).mchirp))
    webpage.lastSection.table[0].row[i].cell[5].text(str(getattr(trig.coincs,ifo).eta))
    webpage.lastSection.table[0].row[i].cell[6].text(str(getattr(trig.coincs,ifo).mass1))
    webpage.lastSection.table[0].row[i].cell[7].text(str(getattr(trig.coincs,ifo).mass2))
    webpage.lastSection.table[0].row[i].cell[8].text(str(getattr(trig.coincs,ifo).eff_distance))

  # Fill the second table
  # We loop over the four separated lists of triggers (one for each of the DQ vetoes categories)
  j = -1
  for trigger_cat_list in followuptrigsList:
    j+=1
    for trigger_cat in trigger_cat_list:

      # find out whether the trigger "trig" is in the list "trigger_cat", and if so write the FAR value in the table
      if getattr(trig.coincs,trig.ifolist_in_coinc[0]).end_time == getattr(trigger_cat.coincs,trigger_cat.ifolist_in_coinc[0]).end_time and getattr(trig.coincs,trig.ifolist_in_coinc[0]).end_time_ns == getattr(trigger_cat.coincs,trigger_cat.ifolist_in_coinc[0]).end_time_ns:
        FAR = getattr(trigger_cat.coincs,trigger_cat.ifolist_in_coinc[0]).alpha
        webpage.lastSection.table[1].row[1].cell[j].text("FAR = " + "%0.3f"%FAR + "\n")
        # compute and print the combined FAR (called "FARc"), unless it is a H1H2 coincidence
        try: 
          trig.ifolist_in_coinc.index('H1')
          trig.ifolist_in_coinc.index('H2')
          if len(trig.ifolist_in_coinc)==2:
            pass
          else:
            FARc = FAR * opts.far_normalization
            webpage.lastSection.table[1].row[1].cell[j].text("FARc = " + "%0.3f"%FARc + "\n")
        except:
          FARc = FAR * opts.far_normalization
          webpage.lastSection.table[1].row[1].cell[j].text("FARc = " + "%0.3f"%FARc + "\n")
        break
 
      else:
        continue

  trig_prev = trig

webpage.cleanWrite('IUL')
