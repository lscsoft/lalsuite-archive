#!/usr/bin/python
#
# Copyright (C) 2009  Steve Fairhurst, based on glitch-page.sh by Duncan
# Brown, ligolw_glitch_page.py by Larne Pekowsky
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#

from optparse import OptionParser

try:
    import sqlite3
except ImportError:
    # pre 2.5.x
    from pysqlite2 import dbapi2 as sqlite3

import sys
import copy
import os
import glob

import glue.segments

from glue.ligolw import ligolw
from glue.ligolw import lsctables
from glue.ligolw import utils
from glue.ligolw.utils import ligolw_sqlite
from glue.ligolw import dbtables

from glue.segmentdb import query_engine

from glue import gpstime

import exceptions
import commands
import matplotlib
matplotlib.use('Agg')
from types import *
from pylab import *
  
##############################################################################
# redefine the SnglInspiral columns of interest
##############################################################################
sngl_cols = [
    "ifo",
    "end_time",
    "end_time_ns",
    "eff_distance",
    "mass1",
    "mass2",
    "mchirp",
    "snr"]

sim_cols = [
    "geocent_end_time",
    "geocent_end_time_ns",
    "distance",
    "mass1",
    "mass2",
    "mchirp",
    "simulation_id"]

lsctables.SnglInspiralTable.loadcolumns = sngl_cols

##############################################################################
def readInjLog(fname):
  """
  read the injection log, and keep those rows corresponding to
  successful injections
  """
  f = open(fname,'r')
  injections = f.readlines()
  f.close()
  injections.pop(0)
  inj_details = []
  for line in injections:
    details = line.strip().split('\t')
    if details[3].strip() == 'Successful':
      # injection performed successfully, so we want to record it
      if (opts.gps_start_time < int(details[0])) and \
          (opts.gps_end_time > int(details[0])):
        inj_details.append(details)

  inj_details.reverse()
  return inj_details

##############################################################################
def readSchedule(fname):
  """
  read the injection schedule, keep only inspiral injections
  """
  f = open(fname,'r')
  injections = f.readlines()
  f.close()
  inj_details = []
  for line in injections:
    details = line.strip().split(' ')
    if "inspiral" in details[3]:
      if (opts.gps_end_time < int(details[0])) and \
          ((opts.gps_end_time + 86400*opts.future_schedule) > int(details[0])):
        inj_details.append(details)

  return inj_details

##############################################################################
# generate sim inspiral info
def generateInjDetails(source,inj_details,max_mass=None,file_type="injected"):
  """
  parse the injection details from the text file, making use of the 
  injection xml
  @param source: the source xml describing the hardware injections
  @param inj_details: the parsed text file
  @param max_mass: the maximum mass to consider
  @param file_type: whether this is an "injected" or "scheduled" file
  """
  doc = utils.load_filename(source)
  simTable = table.get_table(doc, lsctables.SimInspiralTable.tableName)


  # move rows to a holding location
  simTemplates = simTable[:]
  simTable = lsctables.New(lsctables.SimInspiralTable)

  # parse injections performed
  for details in inj_details:
    if (file_type == "scheduled"):
      inj_number = int(details[3].split('_')[3]) - 1
    else:  
      inj_number = int(details[1].split('_')[3]) - 1
    inj = copy.deepcopy(simTemplates[inj_number])
    # set the correct times
    gpsTime = int(details[0])
    inj.geocent_end_time += gpsTime
    inj.h_end_time += gpsTime
    inj.l_end_time += gpsTime
    # add injection
    if (not max_mass) or (inj.mass1 <= max_mass and inj.mass2 <= max_mass):
      simTable.append(inj)

  # construct a valid document with only that table
  doc = ligolw.Document()
  doc.appendChild(ligolw.LIGO_LW())
  doc.childNodes[0].appendChild(simTable)

  return doc, simTable


##############################################################################
def findTriggers(inj):
    """
    query the database for triggers around the time of an injection
  
    @param inj: the injection of interest
    @param xmlfile: the xml file name to be returned by the query
    """

    # we query for triggers in 1 seconds around the time of the injection
    gps_end_time = inj.geocent_end_time + 1
    gps_start_time = inj.geocent_end_time 
    # XXX need to do better and query within ~100ms, by my SQL isn't up to it


    connection = setup_files(options.trigger_dir, gps_start_time, gps_end_time)

    rows = connection.cursor().execute("""SELECT sngl_inspiral.ifo, 
             sngl_inspiral.end_time,
             sngl_inspiral.end_time_ns,
             sngl_inspiral.snr, 
             sngl_inspiral.eff_distance,
             sngl_inspiral.mass1, 
             sngl_inspiral.mass2, 
             sngl_inspiral.mchirp, 
             sngl_inspiral.ttotal
      FROM sngl_inspiral 
      WHERE (cast(sngl_inspiral.end_time as char) || " -- " || cast(sngl_inspiral.snr as char)) IN
             (SELECT cast(end_time as char) || " -- " || cast(MAX(snr) as char)
              FROM sngl_inspiral
              WHERE end_time >= ? AND
                    end_time <  ? AND
                    ifo = ? AND
              GROUP BY end_time
              ORDER BY snr desc)
      AND sngl_inspiral.ifo = ? 
      ORDER BY snr DESC""",  (gps_start_time, gps_end_time, options.ifo, options.ifo) )

    return rows 


##############################################################################
def findRange(inj):
  """
  Query the database to check whether the data was analyzed, 
  and also to find the range of the detector at the time.
  Not sure where the info will be stored -- currently just returns 10.0

  @param inj: the injection of interest
  """
  # XXX need to work out how to query for this info, for now hardcode to 10.
  range = 10.0

  return range

##############################################################################
class InjTable:
  def __init__(self):
    self.num_trigs = 0
    self.num_missed = 0
    self.text = None
    self.missed = None
    self.start_table()
    self.found_doc = ligolw.Document()
    self.found_doc.appendChild(ligolw.LIGO_LW())
    self.inj_table = lsctables.New(lsctables.SimInspiralTable,\
        columns=sim_cols)
    self.found_table = lsctables.New(lsctables.SnglInspiralTable,\
        columns=sngl_cols)
    self.found_doc.childNodes[0].appendChild(self.inj_table)
    self.found_doc.childNodes[0].appendChild(self.found_table)

  def add_found(self,inj,trig):
    """
    check if a trigger was found.  If so, add inj and trig to the found list
    """
    if not trig: return
    self.inj_table.append(copy.deepcopy(inj))
    self.found_table.append(copy.deepcopy(trig))

  def get_chirp_dist(self):
    """
    calculate the chirp distance from the injection table
    """
    return self.inj_table.get_column("distance") * \
        (1.2/self.inj_table.get_column("mchirp"))**(5./6)

  def plot_mchirp(self):
    """
    make a plot of the chirp mass accuracy
    """
    figure()
    chirp_diff = (self.found_table.get_column("mchirp") - \
        self.inj_table.get_column("mchirp")) / \
        self.inj_table.get_column("mchirp")
    mark_size = (self.inj_table.get_column("geocent_end_time") - \
        opts.gps_start_time)/ 6e5
    mark_size = mark_size.astype('l')
    scatter(asarray(self.get_chirp_dist()), asarray(chirp_diff), \
        5*mark_size, mark_size)
    colorbar()
    xlim(xmin=0)
    xlabel("Chirp Distance (Mpc)",size="x-large")
    ylabel("Fractional error in chirp mass",size="x-large")
    grid(True)
    savefig( opts.ifo + "-mchirp_accuracy.png")

  def plot_timing(self):
    """
    make a plot of the timing accuracy
    """
    figure()
    time_diff = 1e-6 * (self.found_table.get_column("end_time_ns") - \
        self.inj_table.get_column("geocent_end_time_ns") )
    distance = self.inj_table.get_column("distance")
    mark_size = (self.inj_table.get_column("geocent_end_time") - \
        opts.gps_start_time)/ 6e5
    mark_size = mark_size.astype('l')
    scatter(asarray(self.get_chirp_dist()), asarray(time_diff), \
        5*mark_size,mark_size)
    colorbar()
    xlim(xmin=0)
    xlabel("Chirp Distance (Mpc)",size="x-large")
    ylabel("Error in end time (ms)",size="x-large")
    grid(True)
    savefig( opts.ifo + "-timing_accuracy.png")

  def plot_distance(self):
    """
    make a plot of the distance accuracy
    """
    figure()
    dist_diff = self.found_table.get_column("eff_distance")/ \
        self.inj_table.get_column("distance") - 1
    snr = self.found_table.get_column("snr")
    mark_size = (self.inj_table.get_column("geocent_end_time") - \
        opts.gps_start_time)/ 6e5
    mark_size = mark_size.astype('l')
    scatter(asarray(self.get_chirp_dist()), asarray(dist_diff), \
        5*mark_size,mark_size)
    colorbar()
    xlim(xmin=0)
    xlabel("Chirp Distance (Mpc)",size="x-large")
    ylabel("Fractional error in distance",size="x-large")
    grid(True)
    savefig( opts.ifo + "-distance_accuracy.png")

  def start_table(self):
    """
    start the table with all the injection details
    """
    self.text ="<table id=\"table\" border=1>\n" + \
               "<tr >" + \
               "<th bgcolor=\"#9999ff\" colspan=\"6\">" + \
               "Hardware Injections</th>" + \
               "<th bgcolor=\"#ffff99\" colspan=\"2\">" + \
               "Instrumental Sensitivity</th>" + \
               "<th bgcolor=\"#99ff9\" colspan=\"6\">" + \
               "Inspiral Triggers</th>" + \
               "<th bgcolor=\"#ffff99\" colspan=\"3\">" + \
               "Parameter Accuracy <br>" + \
               "(detected - injected)</th></tr>\n"\

    # the injection detail header
    self.text+= "<tr><th>GPS End Time</th>" + \
                "<th>End Time</th>" + \
                "<th>Dist (Mpc)</th>" + \
                "<th>Mass 1</th>" + \
                "<th>Mass 2</th>" + \
                "<th>Chirp Mass</th>\n"

    # the instrumental details header
    self.text+= "<th>BNS Horizon</th>" + \
                "<th>Expected SNR</th>\n"

    # the inspiral trigger header
    self.text+= "<th>GPS End Time</th>" + \
                "<th>SNR</th>" + \
                "<th>Dist (Mpc)</th>" + \
                "<th>Mass 1</th>" + \
                "<th>Mass 2</th>" + \
                "<th>Chirp Mass</th>\n" 

    # the parameter accuracy header
    self.text+= "<th>Time (ms)</th>" + \
                "<th>Chirp Mass</th>" + \
                "<th>Frac Dist</th></tr>\n" 


  def add_injection(self,inj,range,trig):
    """
    add an injection to the table
    @param inj:   the injection detail (sim_inspiral table)
    @param range: the instrumental sensitivity (summ_value table)
    @param trig:  the trigger details (sngl_inspiral table)
    """
    # determine type:
    if not range:
      type = "not analyzed"
      range_col = "bgcolor=\"#cccccc\""
      trig_col  = "bgcolor=\"#cccccc\""
      acc_col = range_col
    else:
      horizon = range.value
      expected_snr = 8. * horizon * (inj.mchirp / 1.21)**(5./6) / inj.distance

      if expected_snr < 7:
        type = "near thresh"
        range_col = "bgcolor=\"#ffcc99\""
        trig_col  = "bgcolor=\"#ffcc99\""
        acc_col = range_col
      elif not trig:
        type = "missed"
        range_col = "bgcolor=\"#ffffdd\""
        trig_col = "bgcolor=\"#ff0000\""
        acc_col = trig_col
      else:
        type = "found"
        range_col = "bgcolor=\"#ffffdd\""
        trig_col  = "bgcolor=\"#ddffdd\"" 
        acc_col   = range_col

    inj_col   = "bgcolor=\"#ccccff\""
    
    # start the row
    self.text += "<tr class=\"" + type + "\">"

    # add injection details
    end_time = inj.geocent_end_time + 1e-9 * inj.geocent_end_time_ns
    utc_time = commands.getoutput('tconvert ' + str(inj.geocent_end_time))
    self.text += "<td " + inj_col + " >%.3f</td>" % end_time + \
                 "<td " + inj_col + " >" + utc_time  + "</td>" + \
                 "<td " + inj_col + " >%.2f</td>" % inj.distance + \
                 "<td " + inj_col + " >%.3f</td>" % inj.mass1 + \
                 "<td " + inj_col + " >%.3f</td>" % inj.mass2 + \
                 "<td " + inj_col + " >%.3f</td>" % inj.mchirp

    # add instrument details
    if type == "not analyzed":
      self.text += "<td " + range_col + " colspan=\"2\">" + \
                   "Data Not Analyzed</td>"
    else:
      self.text += "<td " + range_col + " >%.3f</td>" % horizon + \
                   "<td " + range_col + " >%.2f</td>" % expected_snr

    # add trigger details
    if not trig:
      self.text += "<td " + trig_col + " colspan=\"6\">" + \
                   "No Triggers within 10 ms of Injection</td>"
      self.text += "<td " + acc_col + " colspan=\"3\">No Triggers</td>"

    else:
      end_time = trig.end_time + 1e-9 * trig.end_time_ns
      self.text  += "<td " + trig_col + " >%.3f</td>" % end_time + \
                    "<td " + trig_col + " >%.2f</td>" % trig.snr + \
                    "<td " + trig_col + " >%.2f</td>" % trig.eff_distance + \
                    "<td " + trig_col + " >%.3f</td>" % trig.mass1 + \
                    "<td " + trig_col + " >%.3f</td>" % trig.mass2 + \
                    "<td " + trig_col + " >%.3f</td>" % trig.mchirp
      t_diff = 1e-6 * (trig.end_time_ns - inj.geocent_end_time_ns)
      m_diff = trig.mchirp - inj.mchirp
      d_diff = trig.eff_distance/inj.distance - 1
      self.text  += "<td " + range_col + " >%.2f</td>" % t_diff + \
                    "<td " + range_col + " >%.3f</td>" % m_diff + \
                    "<td " + range_col + " >%.3f</td>" % d_diff

    # end the row
    self.text += "</tr>\n"

  def end_table(self):
    """
    end the table with all the injection details
    """
    self.text += "</table>\n"

  def get_table(self):
    """
    return the table with all the injection details
    """
    return self.text


##############################################################################
class SchedTable:
  def __init__(self,schedInj,senseMon):
    self.num_sched = 0
    self.text = None
    self.sensemon = senseMon
    self.sched_doc = ligolw.Document()
    self.sched_doc.appendChild(ligolw.LIGO_LW())
    self.inj_table = lsctables.New(lsctables.SimInspiralTable,\
        columns=sim_cols)
    self.inj_table.extend(schedInj)
    self.sched_doc.childNodes[0].appendChild(self.inj_table)
    self.start_table()
    for inj in schedInj:
      self.add_injection(inj)
    self.end_table()

  def get_chirp_dist(self):
    """
    calculate the chirp distance from the injection table
    """
    return self.inj_table.get_column("distance") * \
        (1.2/self.inj_table.get_column("mchirp"))**(5./6)

  def start_table(self):
    """
    start the table with all the injection details
    """
    # the injection detail header
    self.text = "<table border=1>\n" + \
                "<tr><th>GPS End Time</th>" + \
                "<th>End Time</th>" + \
                "<th>Dist (Mpc)</th>" + \
                "<th>Mass 1</th>" + \
                "<th>Mass 2</th>" + \
                "<th>Chirp Mass</th>\n"

    # the instrumental details header
    for range in self.sensemon:
      self.text+= "<th>Expected SNR for range of " + str(range) + \
          " Mpc</th>\n"


  def add_injection(self,inj):
    """
    add an injection to the table
    @param inj:   the injection detail (sim_inspiral table)
    @param range: the instrumental sensitivity (summ_value table)
    @param trig:  the trigger details (sngl_inspiral table)
    """
    self.num_sched +=1
    if (inj.mass1 > opts.max_mass) or (inj.mass2 > opts.max_mass):
      self.text += "<tr bgcolor=\"#999999\">"
    elif (self.num_sched % 2):
      self.text += "<tr bgcolor=\"#ccccff\">"
    else:
      self.text += "<tr bgcolor=\"#ddddff\">"

    # add injection details
    end_time = inj.geocent_end_time + 1e-9 * inj.geocent_end_time_ns
    utc_time = commands.getoutput('tconvert ' + str(inj.geocent_end_time))
    self.text += "<td>%.3f</td>" % end_time + \
                 "<td>" + utc_time  + "</td>" + \
                 "<td>%.2f</td>" % inj.distance + \
                 "<td>%.3f</td>" % inj.mass1 + \
                 "<td>%.3f</td>" % inj.mass2 + \
                 "<td>%.3f</td>" % inj.mchirp

    # add instrument details
    if (inj.mass1 > opts.max_mass) or (inj.mass2 > opts.max_mass):
      self.text += "<td colspan=\"" + str(len(self.sensemon)) + \
          "\">Masses outside online bank range</td>"
    else:
      for range in self.sensemon:
        expected_snr = 8. * 2.3 * range * \
            (inj.mchirp / 1.21)**(5./6) / inj.distance
        self.text += "<td>%.2f</td>" % expected_snr

    # end the row
    self.text += "</tr>\n"

  def end_table(self):
    """
    end the table with all the injection details
    """
    self.text += "</table></p>\n"

  def get_table(self):
    """
    return the table with all the injection details
    """
    return self.text

##############################################################################
usage = """ %prog [options]
Program to parse the inspiral injection log
"""

parser = OptionParser( usage )

parser.add_option("-s","--gps-start-time",action="store",type="int",\
    default=815119213, metavar="START",
    help="start of GPS time range (default = 815119213)" )

parser.add_option("-e","--gps-end-time",action="store",type="int",\
    default=None, metavar="END",
    help="end of GPS time range (default: read from $HOME/local/etc/IFO-glitch-page-gps)" )

parser.add_option("-d","--days-delay",action="store",type="int",\
    default=2, metavar="DELAY",
    help="Number of days to look back for new triggers (default = 2)")

parser.add_option("-i","--ifo",action="store",type="string",\
    default=None, metavar="IFO",help="name of ifo (default = L1)" )

parser.add_option("-E","--server",action="store",type="string",\
    default="ldas.ligo-wa.caltech.edu:30021", metavar="SERVER",
    help="name of server (default = ldas.ligo-la/wa.caltech.edu:30021)" )

parser.add_option("-M","--max-mass",action="store",type="float",\
    default=3.0, metavar=" MAX_MASS",
    help="maximum mass to consider (default = 3.0)")

parser.add_option("-p","--path",action="store",type="string",\
    default=".",metavar="PATH",help="path to location of output")

parser.add_option("-x","--source-xml",action="store",type="string",\
    default="s5_sources.xml", metavar="IN_XML",\
    help="input xml file of injections (default = s5_sources.xml)" )

parser.add_option("-f","--future-schedule",action="store",type="int",\
    default=0, metavar="DAYS",\
    help="produce injection schedule for coming DAYS")

parser.add_option("-S","--sensemon-range",action="append",type="float",\
    default=None, metavar="RANGE",\
    help="the approximate sensemon range of the instrument")

(opts,args) = parser.parse_args()


##############################################################################
# Sanity check on the arguments                                              #
##############################################################################
if not opts.ifo:
  print >>sys.stderr, "Must specify an ifo"
  sys.exit(1)


if opts.future_schedule and not opts.sensemon_range:
  print >>sys.stderr, \
      "Must specify at least one sensemon range for future schedule"
  sys.exit(1)

##############################################################################
# cd to right directory
os.chdir(opts.path)

##############################################################################
# Update the hardware injection details:                                     #
##############################################################################
if opts.ifo == "L1":
  opts.server = "ldas.ligo-la.caltech.edu:30021"
  os.system("curl -O http://london.ligo-la.caltech.edu/scirun/S5/HardwareInjection/Details/biinj/L1inspirallist.txt")
  fname = 'L1inspirallist.txt'
elif opts.ifo == "H1":
  os.system("curl -O http://blue.ligo-wa.caltech.edu/scirun/S5/HardwareInjection/Details/biinj/H1inspirallist.txt")
  fname = 'H1inspirallist.txt'
elif opts.ifo == "H2":
  os.system("curl -O http://blue.ligo-wa.caltech.edu/scirun/S5/HardwareInjection/Details/biinj/H2inspirallist.txt")
  fname = 'H2inspirallist.txt'
else:
  print >>sys.stderr, "Invalid IFO, must be one of H1, H2, L1"
  sys.exit(1)

if not opts.gps_end_time:
  f = open('/archive/home/inspiralbns/local/etc/' + opts.ifo + \
      '-glitch-page-gps')
  opts.gps_end_time = int(f.read())
  f.close()

# then parse the log:
inj_details = readInjLog(fname)

# generate injection details
doc, simTable = generateInjDetails(opts.source_xml,inj_details,opts.max_mass)
f_xml = file(opts.ifo + "-HARDWARE-INJECTIONS.xml","w")
doc.write(f_xml)
##############################################################################
# Parse the injections looking for instrumental range and triggers           #
##############################################################################
files = os.listdir(".")

# initialize the injection table:
injections = InjTable()

for inj in simTable:

  # first check the range of the instrument.  
  # Also tells us if data has been analyzed.
  xmlfile = opts.ifo + "-SUMM_VALUE-" + str(inj.geocent_end_time) + "-1.xml"
  
  # if no file, or empty recent file, check range:
  if (not xmlfile in files) or  \
      ( (inj.geocent_end_time > opts.gps_end_time - 86400 * opts.days_delay)
      and not (commands.getoutput("lwtprint " + xmlfile)) ):
    # query the database:
    findRange(inj,xmlfile)

  try:
    doc = utils.load_filename(xmlfile)
    rangeTable = table.get_table(doc, lsctables.SummValueTable.tableName)
    if len(rangeTable):
      range = rangeTable[0]
    else:
      range = None
  except: 
    # not a valid LIGOLwxml file, delete
    os.remove(xmlfile)
    range = None
    print 'bad query'

  # now, look for triggers:  
  xmlfile = opts.ifo + "-INSPIRAL-" + str(inj.geocent_end_time) + "-1.xml"

  # get triggers from database if necessary:
  if (not xmlfile in files) or \
      ( (inj.geocent_end_time > opts.gps_end_time - 86400 * opts.days_delay)
      and not (commands.getoutput("lwtprint " + xmlfile)) ):  
    findTrigger(inj,xmlfile)

  # read in trigger
  try:
    doc = utils.load_filename(xmlfile)
    trigTable = table.get_table(doc, lsctables.SnglInspiralTable.tableName)
    if len(trigTable):
      trig = trigTable[0]
    else:
      trig = None
  except: 
    # not a valid LIGOLwxml file, delete
    os.remove(xmlfile)
    trig = None
    print 'bad query'

  # add the injection to the list
  injections.add_injection(inj,range,trig)
  injections.add_found(inj,trig)

##############################################################################
# generate the web page:
injections.end_table()

injections.plot_mchirp()
injections.plot_timing()
injections.plot_distance()

# write an xml file with found injections
f_xml = file(opts.ifo + "-FOUND-HARDWARE-INJECTIONS.xml","w")
injections.found_doc.write(f_xml)

# write out the html file with all injections
fname = opts.ifo + "-HARDWARE_INJECTIONS-" + str(opts.gps_start_time) + \
         "-" + str(opts.gps_end_time - opts.gps_start_time) + ".html"
f = open(fname,"w")

start_time = commands.getoutput('tconvert ' + str(opts.gps_start_time))
end_time = commands.getoutput('tconvert ' + str(opts.gps_end_time))

f.write("<html>\n" + \
    "<head>\n" + \
    "<title>" + opts.ifo + " Hardware Injection Details</title>" + \
    "<script src=\"hwinj_script.js\"></script>" + \
    "</head>\n" + \
    "<body>\n" + \
    "\n" + \
    "<center><h2>" + opts.ifo + " Inspiral Hardware Injections between "+ \
    start_time + " and " + end_time + "</h2></center>\n \n" + \
    "<table border=1><tr><td>Display:</td>\n" + \
    "<td bgcolor=\"#ff0000\"><input id=\"missed\" type=\"checkbox\" name==\"default\" checked onClick=\"changeVisible('table')\">Missed</td>" + \
    "<td bgcolor=\"#ffcc99\"><input id=\"near thresh\" type=\"checkbox\" name=\"default\" checked onClick=\"changeVisible('table')\">Weak</td>"
    "<td bgcolor=\"#ddffdd\"><input id=\"found\" type=\"checkbox\" name=\"default\" checked onClick=\"changeVisible('table')\">Found</td>" + \
    "<td bgcolor=\"#cccccc\"><input id=\"not analyzed\" type=\"checkbox\" name=\"default\" checked onClick=\"changeVisible('table')\">Not Analyzed</td>" + \
    "</tr></table>\n")

f.write(injections.get_table())
f.close()

os.system('rm ' + opts.ifo.lower() + '-hwinj.html')
os.system('ln -s ' + fname + ' ' + opts.ifo.lower() + '-hwinj.html')

##############################################################################
# Write the table with the future hardware injections
##############################################################################
if opts.future_schedule:
  os.system("curl -O http://blue.ligo-wa.caltech.edu/scirun/S5/HardwareInjection/Details/biinj/schedule")

  future_inj = readSchedule("schedule")
  doc, simTable = generateInjDetails(opts.source_xml, future_inj,\
      file_type="scheduled")
  schedDetails = SchedTable(simTable,opts.sensemon_range)

  fsched = opts.ifo + "-SCHEDULED-HARDWARE-INJECTIONS-" + \
      str(opts.gps_end_time) + \
      "-" + str(opts.gps_end_time + 86400 * opts.future_schedule) + ".html"

  sched_start = commands.getoutput('tconvert ' + str(opts.gps_end_time))
  sched_end   = commands.getoutput('tconvert ' + \
      str(opts.gps_end_time + 86400 * opts.future_schedule))

  f = open(fsched,"w")
  f.write("<center><h2> Scheduled Inspiral Hardware injections between " \
      + sched_start + " and " + sched_end + "</h2></center>\n \n")
  f.write(schedDetails.get_table())
  f.close()

  os.system('rm ' + opts.ifo.lower() + '-hwinj-schedule.html')
  os.system('ln -s ' + fsched + ' ' + opts.ifo.lower() + '-hwinj-schedule.html')


