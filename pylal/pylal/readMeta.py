#!/usr/bin/env python
import sys, getopt
import re
from lgen import metaio
from pylab    import *

def uniq(list):
  l = []
  for m in list:
    if m not in l:
      l.append(m)
  return l

class metaDataTable:
  def __init__(self, triggerfile, tabletype):
    self.tabletype = tabletype
    if ( triggerfile ):
      self.readfiles(triggerfile, tabletype)
    else:
      self.table = []

  def readfiles(self, triggerfile, tabletype):
    if tabletype == "search_summary":
      self.table = metaio.read_search_summary(triggerfile)
    if tabletype == "summ_value":
      self.table = metaio.read_summ_value(triggerfile)
    if tabletype == "sngl_inspiral":
      self.table = metaio.read_sngl_inspiral(triggerfile)
    if tabletype == "sngl_burst":
      self.table = metaio.read_sngl_burst(triggerfile)

  def nevents(self):
    return len(self.table)

  def mkarray(self, colname):
    myarray = asarray( [ self.table[i][colname] for i in range(self.nevents())] )
    return myarray


class coincInspiralTable:
  def __init__(self, inspTriggers): 
    # an instance of the snglInspiralTable
    h1triggers = metaDataTable("", "sngl_inspiral")
    h1triggers.table = [e for e in inspTriggers.table if (re.match("H1",e["ifo"]))]
    h2triggers = metaDataTable("", "sngl_inspiral")
    h2triggers.table = [e for e in inspTriggers.table if (re.match("H2",e["ifo"]))]
    l1triggers = metaDataTable("", "sngl_inspiral")
    l1triggers.table = [e for e in inspTriggers.table if (re.match("L1",e["ifo"]))]

    # use the supplied method to convert these columns into numarrays
    eventidlist = asarray(uniq(inspTriggers.mkarray("event_id")))
    self.table = []
    for m in eventidlist: 
      self.table.append({ "event_id":m, "H1":{}, "H2":{}, "L1":{} })
  
    for m in self.table:
      for k in h1triggers.table:
        if m["event_id"] == k["event_id"]:
          m["H1"] = k
          break
  
    for m in self.table:
      for k in h2triggers.table:
        if m["event_id"] == k["event_id"]:
          m["H2"] = k
          break
  
    for m in self.table:
      for k in l1triggers.table:
        if m["event_id"] == k["event_id"]:
          m["L1"] = k
          break

  def mkarray(self, colname, ifoname):
    mylist = []
    for i in range(len(self.table)):
      tmpvalue = 0;
      if self.table[i][ifoname].has_key(colname):
        tmpvalue += (self.table[i][ifoname][colname])
      mylist.append( tmpvalue )

    return asarray( mylist )



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
