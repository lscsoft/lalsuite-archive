#!/usr/bin/python

import git_version

__author__ = "Larry Price <larry.price@ligo.org> and Patrick Brady <patrick.brady@ligo.org>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date

import glob
from optparse import *
from pylal import skylocutils


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
  parser = OptionParser( usage=usage, version="%prog CVS $Id$ " )

  # options related to input and output
  parser.add_option("-g","--glob",action="store",type="string",\
      default=None, metavar=" GLOB",help="GLOB of thinca files to read" )
  parser.add_option("-d","--plotpoints",action="store_true",\
      default=False, help="make a color coded plot of the sky" )
  parser.add_option("-u","--usecatalog",action="store",type="string",\
      default=None, metavar=" CATALOG_NAME", help="galaxy catalog to use; must be specified if --listgalaxies option is used")
  parser.add_option("-f","--reference-frequency",action="store",type="float", default=0.0, metavar=" REFERENCE_FREQUENCY", \
    help="reference frequency for signal timing" )
  parser.add_option("-n","--dt90",action="store",type="float",\
      default=4, metavar=" DT90",help="90% timing threshold" )
  parser.add_option("-w","--dt60",action="store",type="float",\
      default=15, metavar=" DT60",help="60% timing threshold" )
  parser.add_option("-N","--snr-dt90",action="store",type="float",\
      default=None, metavar=" SNRDT90",help="90% snr-dependent timing threshold" )
  parser.add_option("-W","--snr-dt60",action="store",type="float",\
      default=None, metavar=" SNRDT60",help="60% snr-dependent timing threshold" )
  parser.add_option("-D","--dD90",action="store",type="float",\
      default=None, metavar=" DD90",help="90% effective distance threshold" )
  parser.add_option("-E","--dD60",action="store",type="float",\
      default=None, metavar=" DD60",help="60% effective distance threshold" )
  parser.add_option("-o","--output",action="store",type="string",default=None,\
                    help="name of the xml file to store output in")
  parser.add_option("-z","--input-type",action="store",default="coinctable",\
                    help="specify the type of input in the glob.  valid options are coinctable (DEFAULT) and coire")
  parser.add_option("-y","--timing-only",action="store_true",default=False,\
                    help="only use timing information for sky localization")
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

#put the thresholds here for easy access later
thresholds_60 = (opts.dt60,opts.snr_dt60,opts.dD60)
thresholds_90 = (opts.dt90,opts.snr_dt90,opts.dD90)

#deal with the glob 
files = []
for gl in opts.glob.split(" "):
  files.extend(glob.glob(gl))
if len(files) < 1:
  print >>sys.stderr, "The glob for " + opts.glob + " returned no files" 
  sys.exit(1)

#put the files into the coinc data structure
coincs = Coincidences(files,opts.input_type)

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

#make it work with the XSL stylesheet
ligolw.Header += u"""\n\n"""\
                 +u"""<?xml-stylesheet type="text/xsl" href="ligolw.xsl"?>"""\
                 +u"""\n\n"""


##############################################################################
#
#          convenience functions
#
##############################################################################


##############################################################################
#
#          main program
#
##############################################################################

#open up the pickled grids

#the area of each pixel on the fine grid in square degrees
pixel_area = 0.25

for coinc in coincs:
  sp = SkyPoints()

  #for gathering information about area on the sky
  dt90_area = 0.0
  dt60_area = 0.0
  dt90dD90_area = 0.0
  dt60dD90_area = 0.0

  #compute combined snr if snr dependent thresholds are specified
  if thresholds_90[1] and thresholds_60[1]:
    rhosquared = 0.0
    for ifo in coinc.ifo_list:
      rhosquared += coinc.snr[ifo]*coinc.snr[ifo]
    snr = sqrt(rhosquared)
  else:
    snr = None

  #main loop over the coarse grid
  for coarse_pt in grid.keys():

    #use timing alone to determine if we should move to the fine grid
    dtrss_coarse = skylocutils.get_delta_t_rss(coarse_pt,coinc,ref_freq)
    if dtrss_coarse <= 1.1*thresholds_90[0]:

      #loop over points on the fine grid
      for fine_pt in grid[coarse_pt]:
        dtrss_fine = skylocutils.get_delta_t_rss(fine_pt,coinc,ref_freq)
        dDrss_fine = skylocutils.get_delta_D_rss(fine_pt,coinc)

        #now compute relevant sky areas
        if snr:
          #compute the quanitity for the snr-dependent threshold
          snr_dtrss = dtrss_fine*snr/10.0
          sp.append(fine_pt[0],fine_point[1],snr_dtrss,dDrss_fine)
          if snr_dtrss <= thresholds_60[1]:
            dt60_area += pixel_area
            dt90_area += pixel_area
            if dDrss_fine <= thresholds_60[2]:
              dt60dD60_area += pixel_area
              dt90dD90_area += pixel_area
          elif snr_dt_rss <= thresholds_90[1]:
            dt90_area += pixel_area
            if dDrss_fine <= thresholds_90[2]:
              dt90dD90_area += pixel_area
        else:
          #just timing
          sp.append(fine_pt[0],fine_point[1],dtrss_fine,dDrss_fine)
          if dtrss_fine <= thresholds_60[0]:
            dt60_area += pixel_area
            dt90_area += pixel_area
            if dDrss_fine <= thresholds_60[2]:
              dt60dD60_area += pixel_area
              dt90dD90_area += pixel_area
          elif dtrss_fine <= thresholds_90[0]:
            dt90_area += pixel_area
            if dDrss_fine <= thresholds_90[2]:
              dt90dD90_area += pixel_area


        #FIXME: put galaxy catalog stuff here!!!
     
    else:
      sp.append(coarse_pt[0],coarse_pt[1],dtrss_coarse,None)
  
  #write the grid
  sp.write(grid_file,argstring)

  #populate the output tables
  #list of points has been sorted so the best one is at the top
  populate_SkyLocTable(skyloctable,coinc,dt60_area,dt90_area,dt60dD60_area,\
                       dt90dD90_area,sp[0],grid_file,skymap_file)
  if coinc.is_injection:
    #NB: using the *recovered* snr for the snr dependent threshold
    inj_pt = (coinc.latitude_inj,coinc.longitude_inj)
    if snr:
      dtrss_inj = snr*skylocutils.get_delta_t_rss(inj_pt,coinc,ref_freq)/10.0
    else:
      dtrss_inj = skylocutils.get_delta_t_rss(inj_pt,coinc,ref_freq)
    dDrss_inj = skylocutils.get_delta_D_rss(inj_pt,coinc)
    area_inj = 0.0
    for pt in sp:
      if pt[2] <= dtrss_inj and pt[3] <= dDrss_inj:
        area_inj += pixel_area
    populate_SkyLocInjTable(skylocinjtable,coinc,area_inj,L_inj,dtrss_inj,\
                            dDrss_inj)


