#!/usr/bin/env python
import sys, getopt
from lgen import metaio

class readSnglInspiralTable:
  def __init__(self, triggerfile):
    self.filename = triggerfile
    M = metaio.read_sngl_inspiral('%s' % triggerfile)
    self.end_time     = M[:,0]+1.0e-9*M[:,1]
    self.eff_distance = M[:,2]
    self.coa_phase    = M[:,3]
    self.mass1        = M[:,4]
    self.mass2        = M[:,5]
    self.mchirp       = M[:,6]
    self.eta          = M[:,7]
    self.snr          = M[:,8]
    self.chisq        = M[:,9]

class readSnglBurstTable: 
  def __init__(self, triggerfile):
    self.filename = triggerfile
    M = metaio.read_sngl_burst('%s' % triggerfile)
    self.start_time   = M[:,0]+1.0e-9*M[:,1]
    self.peak_time    = M[:,2]+1.0e-9*M[:,3]
    self.duration     = M[:,4]
    self.central_freq = M[:,5]
    self.bandwidth    = M[:,6]
    self.snr          = M[:,7]
    self.confidence   = M[:,8]

class readSimInspiralTable: 
  def __init__(self, triggerfile):
    self.filename = triggerfile
    M = metaio.read_sim_inspiral('%s' % triggerfile)
    self.geo_end_time = M[:,0]+1.0e-9*M[:,1]
    self.distance     = M[:,2]
    self.f_final      = M[:,3]
    self.mass1        = M[:,4]
    self.mass2        = M[:,5]
    self.mchirp       = M[:,6]
    self.eta          = M[:,7]
    self.eff_dist_h   = M[:,8]
    self.eff_dist_l   = M[:,9]

def usage():
        print "readMeta.py -x xml file"
        sys.exit(0)

def main():
  # parse options
  try:
    opts,args=getopt.getopt(sys.argv[1:],"x:h",["xmlfile=","help"])
  except getopt.GetoptError:
    usage()
    sys.exit(2)
  for o,a in opts:
    if o in ("-x","--xmlfile"):
             xmlfile = a
    if o in ("-h", "--help"):
             usage()
             sys.exit()
  print 'main() does nothing at the moment' 

# execute main if this module is explicitly invoked by the user
if __name__=="__main__":
        main()
