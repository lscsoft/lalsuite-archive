#!/usr/bin/python

from glue import git_version

__author__ = "Larry Price <larry.price@ligo.org>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date

import sys
import gzip
from optparse import *
import glob
from matplotlib import pyplot
import matplotlib.cm as cm
from mpl_toolkits.basemap import Basemap
import numpy as np

usage = """
usage: %prog [options] 

Make a skymap from the gzipped output of run_skypoints.py
NOTE: Requires the python basemap module!

"""

def parse_command_line():
  """
  Parser function dedicated
  """
  parser = OptionParser( usage=usage, version="%prog "+__version__ )

  # options related to input and output
  parser.add_option("-g","--glob",action="store",type="string",\
      default=None, metavar=" GLOB",help="GLOB of files to read" )
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

for file in files:
  
  latitude = []
  longitude = []
  color = []
  f = gzip.open(file,'r')
  #set limits for the colors
  #these are chosen by trial and error and likely needs more tweaking
  if 'posterior' in file:
    cutoff = .00010
    minc = cutoff
    maxc = 0.003
  elif 'probability' in file:
    cutoff = 0.01
    minc = cutoff
    maxc = 0.7
  else:
    cutoff = 0.0
    minc = None
    maxc = None
  for line in f:
    if '#' in line:
      continue
    elif line.strip():
      lon, lat, c = line.strip().split()
      #remove everything with a probability below the cutoff
      #this makes the skymaps cleaner
      if float(c) >= cutoff:
        latitude.append(float(lat)*180/np.pi)
        longitude.append(float(lon)*180/np.pi)
        color.append(float(c))

  f.close()

  #make sure we've got something to plot!
  if not latitude or not longitude or not color:
    print >>sys.stderr, "Unable to plot " + file
    continue

  #sort them from least to most likely to avoid 
  #strange plotting artifacts
  sortme = zip(color,latitude,longitude)
  sortme.sort()
  color, latitude, longitude = zip(*sortme)

  pyplot.clf()
  m = Basemap(projection='moll', lon_0=180.0, lat_0 = 0.0)
  plx,ply = m(np.asarray(longitude),np.asarray(latitude))
  pyplot.scatter(plx,ply,s=5,c=np.asarray(color),faceted=False,cmap=cm.jet,vmin=minc,vmax=maxc)
  m.drawmapboundary()
  m.drawparallels(np.arange(-90.,120.,45.),labels=[1,0,0,0],labelstyle='+/-') # draw parallels
  m.drawmeridians(np.arange(0.,420.,90.),labels=[0,0,0,1],labelstyle='+/-') # draw meridians
  pyplot.title("Skymap for "+str(file)) # add a title
  pyplot.colorbar()
  pyplot.savefig(file.replace('gz','png'))
