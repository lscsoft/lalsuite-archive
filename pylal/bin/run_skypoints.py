#!/usr/bin/python

from glue import git_version

__author__ = "Larry Price <larry.price@ligo.org> and Patrick Brady <patrick.brady@ligo.org>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date

import os
import sys
import glob
from optparse import *
import cPickle
from pylal import skylocutils
from glue.ligolw import ligolw, lsctables

##############################################################################
#
#          options
#
##############################################################################

usage = """
usage: %prog [options] 

Estimate the sky position from a coincident trigger.

"""


def parse_command_line():
  """
  Parser function dedicated
  """
  parser = OptionParser( usage=usage, version="%prog "+__version__ )

  # options related to input and output
  parser.add_option("-g","--glob",action="store",type="string",\
      default=None, metavar=" GLOB",help="GLOB of thinca files to read" )
  parser.add_option("-G","--grids",action="store",type="string",\
      default=None, metavar=" GRID",help="pickled sky grids")
  parser.add_option("-R","--ranks",action="store",type="string",\
      default=None, metavar=" RANKS",help="pickled ranking object")
#plotting functionality will be moved elsewhere and is disabled for now
#  parser.add_option("-p","--plotpoints",action="store_true",\
#      default=False, help="make a color coded plot of the sky" )
  parser.add_option("-u","--usecatalog",action="store",type="string",\
      default=None, metavar=" CATALOG_NAME", help="galaxy catalog to use; must be specified if --listgalaxies option is used")
  parser.add_option("-f","--reference-frequency",action="store",type="float", default=0.0, metavar=" REFERENCE_FREQUENCY", \
    help="reference frequency for signal timing" )
  parser.add_option("-s","--snr-threshold",action="store_true",\
    default=False, help="use an snr-dependent quantity for the timing ranking" )
  parser.add_option("-o","--output-prefix",action="store",type="string",default='',\
                    help="appends ouput-prefix to output file names")
  parser.add_option("-z","--input-type",action="store",default="coinctable",\
                    help="specify the type of input in the glob.  valid options are coinctable (DEFAULT) and coire")
  parser.add_option("-y","--timing-only",action="store_true",default=False,\
                    help="only use timing information for sky localization")
  parser.add_option("-d","--debug",action="store_true",default=False,\
                    help="write out debugging info in addition to normal output")

  (options,args) = parser.parse_args()

  return options, sys.argv[1:]



##############################################################################
#
#          i/o setup
#
##############################################################################

opts, args = parse_command_line()

#dump the args into a single string to write with the grids
argstring = " ".join(args)


#deal with the glob 
files = []
if opts.glob is not None:
  for gl in opts.glob.split(" "):
    files.extend(glob.glob(gl))
  if len(files) < 1:
    print >>sys.stderr, "The glob for " + opts.glob + " returned no files" 
    sys.exit(1)
else:
  print >>sys.stderr, "Need to specify a glob"
  sys.exit(1)
#put the files into the coinc data structure
coincs = skylocutils.Coincidences(files,opts.input_type)

#set the reference frequency if requested
if opts.reference_frequency:
  ref_freq = opts.reference_frequency
else:
  ref_freq = None

#setup the output xml tables
xmldoc = ligolw.Document()
xmldoc.appendChild(ligolw.LIGO_LW())
skyloctable = lsctables.New(skylocutils.SkyLocTable)
skylocinjtable = lsctables.New(skylocutils.SkyLocInjTable)
xmldoc.childNodes[0].appendChild(skyloctable)
xmldoc.childNodes[0].appendChild(skylocinjtable)
#FIXME: insert a process params table!


#make it work with the XSL stylesheet
ligolw.Header += u"""\n\n"""\
                 +u"""<?xml-stylesheet type="text/xsl" href="ligolw.xsl"?>"""\
                 +u"""\n\n"""

#setup the output filenames
base_name = 'SKYPOINTS' + opts.output_prefix
grid_fname = base_name + '_grid_GPSTIME.txt.gz'
outfile = base_name + '_GPSTIME.xml'

##############################################################################
#
#          convenience functions
#
##############################################################################

def get_unique_filename(name):
  """
  use this to avoid name collisions
  """
  counter = 1
  base_name, ext = os.path.splitext(name)
  while os.path.isfile(base_name):
    base_name = base_name + '_' + str(counter) + ext
    counter += 1

  return base_name + ext

##############################################################################
#
#          main program
#
##############################################################################

#open up the pickled grids
skygrids = open(opts.grids,'r')
grid = cPickle.load(skygrids)
skygrids.close()

#open up the pickled rankings
rankfile = open(opts.ranks,'r')
rankings = cPickle.load(rankfile)
rankfile.close()
dtr = rankings['dt']
dDr = rankings['dD']
dtdDr = rankings['dtdD']


#the area of each pixel on the fine grid in square degrees
#this gets recorded for each point and makes computing areas simple
#FIXME: record this information in the grid pickles!
fine_area = 0.25
coarse_area = 16

for coinc in coincs:
  sp = skylocutils.SkyPoints()
  
  #for gathering information about area on the sky
  dt90_area = 0.0
  dt60_area = 0.0
  r90_area = 0.0
  r60_area = 0.0

  #compute combined snr if snr dependent thresholds are specified
  if opts.snr_threshold:
    rhosquared = 0.0
    for ifo in coinc.ifo_list:
      rhosquared += coinc.snr[ifo]*coinc.snr[ifo]
    snrfac = sqrt(rhosquared)/10.0
  else:
  #otherwise just multiply by unity
    snrfac = 1.0

  #main loop over the coarse grid
  for coarse_pt in grid.keys():

    #use timing alone to determine if we should move to the fine grid
    dtrss_coarse = skylocutils.get_delta_t_rss(coarse_pt,coinc,ref_freq)
    coarse_rank = dtr.get_rank(dtrss_coarse*snrfac)
    if coarse_rank >= 0.05:
      
      #loop over points on the fine grid
      for fine_pt in grid[coarse_pt]:
        dtrss_fine = skylocutils.get_delta_t_rss(fine_pt,coinc,ref_freq)
        dtrank = dtr.get_rank(dtrss_fine*snrfac)
        dDrss_fine = skylocutils.get_delta_D_rss(fine_pt,coinc)
        dDrank = dDr.get_rank(dDrss_fine)
        ranking = dtrank*dDrank
        sp.append((fine_pt[0],fine_pt[1],ranking,dtrss_fine*snrfac,dDrss_fine,fine_area))
        #compute relevant areas
        #note that 1-rank is the percantage of injections with the same or lower rss
        if 0.1 <= dtrank < 0.4:
          dt90_area += fine_area
        elif dtrank >= 0.4:
          dt90_area += fine_area
          dt60_area += fine_area
        if 0.1 <= ranking < 0.4:
          r90_area += fine_area
        elif ranking >=0.4:
          r90_area += fine_area
          r60_area += fine_area

        #FIXME: put galaxy catalog stuff here!!!
     
    else:
      sp.append((coarse_pt[0],coarse_pt[1],None,coarse_rank,dtrss_coarse*snrfac,None,coarse_area))
  
  #check for name collisions and then write the grid
  #use seconds of the smallest gpstime to label the event
  grid_file = get_unique_filename(grid_fname.replace('GPSTIME',str(coinc.time.seconds)))
  sp.write(grid_file,argstring)

  #populate the output tables
  #list of points has been sorted so the best one is at the top
  #FIXME: replace None with a link to the skymap file name!!!
  skylocutils.populate_SkyLocTable(skyloctable,coinc,dt60_area,dt90_area,r60_area,\
                       r90_area,sp[0],grid_file,None)
  if coinc.is_injection:
    #NB: using the *recovered* snr for the snr dependent threshold
    inj_pt = (coinc.latitude_inj,coinc.longitude_inj)
    #compute combined snr if snr dependent thresholds are specified
    if opts.snr_threshold:
      rhosquared = 0.0
      for ifo in coinc.ifo_list:
        rhosquared += coinc.snr[ifo]*coinc.snr[ifo]
      snrfac = sqrt(rhosquared)/10.0
    else:
      #otherwise just multiply by unity
      snrfac = 1.0
    
    dtrss_inj = skylocutils.get_delta_t_rss(inj_pt,coinc,ref_freq)
    dtrank_inj = dtr.get_rank(dtrss_inj*snrfac)
    dDrss_inj = skylocutils.get_delta_D_rss(inj_pt,coinc)
    dDrank_inj = dDr.get_rank(dDrss_inj)
    rank_inj = dtrank_inj*dDrank_inj
    dt_area = 0.0
    rank_area = 0.0
    for pt in sp:
      if pt[2] == None:
        pass
      elif pt[3] <= dtrss_inj:
        dt_area += pt[5]
      #FIXME: utlimately the next line is the place where we compute area according to rank
      if pt[2] == None:
        pass
      elif pt[2] >= rank_inj:
        rank_area += pt[5]
    skylocutils.populate_SkyLocInjTable(skylocinjtable,coinc,rank_inj,dt_area,rank_area,\
                                        dtrss_inj,dDrss_inj)

#name the xml file
if len(coincs) > 1:
  tmin = min([min(c.gps.values()) for c in coincs]).seconds
  tmax = max([max(c.gps.values()) for c in coincs]).seconds
  ofname=outfile.replace('GPSTIME',str(tmin)+'-'+str(tmax))
else:
  tmin = min([c for c in coincs[0].gps.values()])
  ofname=outfile.replace('GPSTIME',str(tmin))
output = get_unique_filename(ofname)
#write the xml file and we're done
f = open(output,'w')
xmldoc.write(f)
f.close()

