#!/usr/bin/python

from pylal import git_version

__author__ = "Larry Price <larry.price@ligo.org>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date


import sys
import cPickle
import glob
from math import sqrt
from optparse import *
from pylal import skylocutils

usage = """
usage: %prog [options]

Create a pickle of rankings for use by run_skypoints.py.
See ligovirgocvs/cbc/protected/projects/s6/sky_localization for an example.

"""

def parse_command_line():
  """
  Parser function dedicated
  """
  parser = OptionParser(usage=usage)
  #FIXME: make this work with the xml files
  parser.add_option("-g","--glob",action="store",type="string",\
    default=None, metavar=" GLOB",help="GLOB of coire or coinc (xml) files to read" )
  parser.add_option("-z","--input-type",action="store",default="coire",\
    help="specify the type of input in the glob.  valid options are coinctable (DEFAULT) and coire")
  parser.add_option("-s","--snr-threshold",action="store_true",\
    default=False, help="use an snr-dependent quantity for the timing ranking" )
  parser.add_option("-f","--reference-frequency",action="store",type="float", default=0.0, metavar=" REFERENCE_FREQUENCY", \
    help="reference frequency for signal timing" )


  (options,args) = parser.parse_args()

  return options, sys.argv[1:]

opts, args = parse_command_line()

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

dt = []
dD = []

for coinc in coincs:
  inj_pt = (coinc.latitude_inj,coinc.longitude_inj)
  
  if len(coinc.ifo_list) < 3:
    continue
  else:
    if opts.snr_threshold:
      rhosquared = 0.0
      for ifo in coinc.ifo_list:
        rhosquared += coinc.snr[ifo]*coinc.snr[ifo]
      dtrss_inj = sqrt(rhosquared)*skylocutils.get_delta_t_rss(inj_pt,coinc,opts.reference_frequency)/10.0
    else:
      dtrss_inj = skylocutils.get_delta_t_rss(inj_pt,coinc,opts.reference_frequency)
    dDrss_inj = skylocutils.get_delta_D_rss(inj_pt,coinc)
    dt.append(dtrss_inj)
    dD.append(dDrss_inj)

dtrankings = skylocutils.Ranking(dt)
dDrankings = skylocutils.Ranking(dD)

#by making a ranking of the combined (multiplied) dt and dD rankings
#we can easily read off sky areas since a dtdD ranking of 0.x will 
#reflect an x% probability, i.e., we're just approximating a cdf from the 
#joint probabilities
dtdD = [dtrankings.get_rank(i)*dDrankings.get_rank(j) for i,j in zip(dt,dD)]
dtdDrankings = skylocutils.Ranking(dtdD,smaller_is_better=False)

rankings = {}
rankings['dt'] = dtrankings
rankings['dD'] = dDrankings
rankings['dtdD'] = dtdDrankings

f = open('rankings.pkl','w')
cPickle.dump(rankings,f,protocol=2)
f.close


