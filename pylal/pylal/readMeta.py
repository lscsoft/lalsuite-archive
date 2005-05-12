#!/usr/bin/env python
import sys, getopt
from lgen import metaio
from pylab    import *

class searchSummaryTable: 
  def readfiles(self, triggerfile):
    self.table = metaio.read_search_summary(triggerfile)

class snglInspiralTable:
  def readfiles(self, triggerfile):
    self.table = metaio.read_sngl_inspiral(triggerfile)

  def nevents(self):
    return len(self.table)

  def mkarray(self, colname):
    myarray = asarray( [ self.table[i][colname] for i in range(self.nevents())] )
    return myarray


class simInspiralTable: 
  def __init__(self, triggerfile):
    self.filename = triggerfile
    self.table = metaio.read_sim_inspiral(triggerfile)

class snglBurstTable: 
  def __init__(self, triggerfile):
    self.filename = triggerfile
    self.table = metaio.read_sngl_burst(triggerfile)
    
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
