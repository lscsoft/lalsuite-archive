from pylab import *

import argparse
from optparse import OptionParser
import os
import sys
import random as rand
py_version = sys.version_info[:2]
np_version = np.__version__

'''
parser = OptionParser()


parser.add_option('-i', '--ifos', nargs=3, type=str, metavar='FILE', dest='ifos', help='List of IFOs to be considered') # FIXME: Variable length???
parser.add_option('-s', '--segfiles', nargs=3, type=str, metavar='FILE', dest='segf', help='path to the segment files for IFOS in the order entered above')
parser.add_option('-v', '--vetofiles', nargs=3, type=str, metavar='FILE', dest='vetof', help='path to the veto files for H1, L1, V1 in this order')
parser.add_option('-o', '--outfile', nargs=1, type=str, metavar='FILE', dest='outf', help='path to the output file')
parser.add_option('-l', '--length', nargs=1, type=float, metavar='FLOAT', dest='length', help='length of signal segments in seconds', default=45.0)
#parser.add_option('-d', '--doubles', type=

(options, args) = parser.parse_args()

outfile = str(options.outf)
segfiles = options.segf
vetofiles = options.vetof
ifos = options.ifos
nifos = len(ifos)
seglen = options.length


#parser = argparse.ArgumentParser(description='Find data segments on single/double/triple HLV time free of vetoes.')
#parser.add_argument('integers', metavar='N', type=int, nargs='+', help='an integer for the accumulator')
#parser.add_argument('--sum', dest='accumulate', action='store_const', const=sum, default=max, help='sum the integers (default: find the max)')

#args = parser.parse_args()

# Read segments data in rawlist and only keep long enough segments in fitlist for each of the IFOs

rawlist = []
fitlist = []
for i in arange(nifos):
  ifoname = ifos[i]
  tmpdata = genfromtxt(segfiles[i], dtype=int)
  rawlist.append(tmpdata) # use the IFO class!
  fitlist.append(tmpdata[where(tmpdata[:,2] - tmpdata[:,1] > seglen)[0]])

singlelist = fitlist

'''

class IFO:
  '''An interferometer. Can also refer to multiple interferometers to host doubles, triples etc.'''
  def __init__(self, name, segments, vetoes, minlen=1):
    self._name = name
    self._minlen = minlen

    if type(segments) is str:
      print "Reading segments for " + self._name + ""
      self._segments = self.readSegments(segments)
    elif type(segments) is ndarray:
      self._segments = segmentData(list(segments))
    elif  isinstance(segments, segmentData):
      self._segments = segments
    else:
      print "Cannot recognize segments!"
      return -1

    if type(vetoes) is str:
      print "Reading veto segments for " + self._name + ""
      self._vetoes = self.readSegments(vetoes)
    elif type(vetoes) is ndarray:
      self._vetoes = segmentData(list(segments))
    elif isinstance(vetoes, segmentData):
      self._vetoes = vetoes
    else:
      print "Cannot recognize veto segments!"
      return -1

    self.setUnvetoed()
    print "Number of unvetoed segments that fit: " + str(len(self._unvetoed._seglist))

  def setUnvetoed(self):
    '''This is getting a list of unvetoed segments that fit the minimum length'''
    self._unvetoed = self._segments.fits(self._minlen).getUnvetoed(self._vetoes).fits(self._minlen)

  def readSegments(self, segfname):
    print "Reading " + segfname
    segdata = genfromtxt(segfname, dtype=int)
    return segmentData(list(segdata))

  def printUnvetoedToFile(self, outname):
    self._unvetoed.printToFile(outname)

  def getTrigTimes(self, n, outfile=None):
    '''Returns a list of gps times on which injections can be made'''
    trigtimes = []
    l = self._minlen
    for seg in self._unvetoed._seglist:
      t = seg[1]
      while t + l <= seg[2]:
        trigtimes.append(t + l - 2)
        t += l
    
    if outfile is None:
      return trigtimes
    else:
      savetxt(outfile, array(trigtimes), fmt='%i')

class segment:
  '''A time segment in an IFO'''

  def __init__(self, data, gpsstart=None, gpsend=None):
    '''Creates a segment'''
    if gpsstart is not None and gpsend is not None:
      self._id = data
      self._start = gpsstart
      self._end = gpsend
      self._len = max(gpsend-gpsstart, 0)
    elif gpsstart is None and gpsend is None:
      '''Overloading to generate segment from an array (id, start, end, length)'''
      if len(data) is not 4:
        print "Segment data doesn't have the correct length!"
        return -1
      self._id = data[0]
      self._start = data[1]
      self._end = data[2]
      self._len = data[3]
    else:
      print "Wrong arguments to create segment!"
      return -1
#    if self._len != (self._end - self._start):
#      print "Error in segment data: inconsistent length! " + str(self._len) + " " + str(self._end - self._start)

  def intersectWithList(self, id0, other):
    '''Intersect segment with list of (non-overlapping) segments'''
    newlist = []
    id = id0
    for oseg in other:
      newseg = self.intersectWith(segment(oseg))
      if newseg._start < newseg._end:
        newseg._id = id
        id += 1
        newlist.append(newseg.toArray())
    return newlist


  def intersectWith(self, other):
    '''Intersects with another segment'''
    newstart = max(self._start, other._start)
    newend = min(self._end, other._end)
    newseg = segment(self._id, newstart, newend)
    return newseg

  def toArray(self):
    a = array([self._id, self._start, self._end, self._len])
    return a
    
class segmentData:   
  '''Data that holds segments for one or more IFOs.'''

  def __init__(self, seglist, gpsstart=None, gpsend=None):
    if gpsstart is None:
      gpsstart = seglist[0][1]
    if gpsend is None:
      gpsend = seglist[-1][2]
    sortedlist = list(array(seglist)[argsort(array(seglist)[:,1])])
    self._seglist = sortedlist
    self._start = gpsstart
    self._end = gpsend

  def fits(self, minlen):
    '''Returns the list of segments that fit a minimum required length.'''
    data = array(self._seglist)
    fitted = data[where(data[:,3] >= minlen)[0]]
    fitted = vstack( (arange(len(fitted)),fitted[:, 1:].T) ).T
    return segmentData(list(fitted), self._start, self._end)

  def intersectSegments(self, other):
    '''Returns the intersection with another list of segments.'''
    newdata = []
    id0 = self._seglist[0][0]
    for s in self._seglist:
      seg = segment(s)
      id = id0 + len(newdata)
      newdata += seg.intersectWithList(id, other._seglist)
    return segmentData(newdata)

  def notSegments(self, start, end):
    '''Get the complement of the (non-overlapping) segments in this data'''
    notlist = []
    c = 0
    times = array(self._seglist)[:,[1,2]]
    print shape(times)
    times = times[argsort(times[:,1])]
    t1 = start
    for i in (arange(len(times))):
      t2 = times[i,0]
      if t2 > t1:
        newseg = array([c, t1, t2, t2-t1]) #FIXME: check c initial value
        notlist.append(newseg)
        c += 1
      t1 = times[i,1]
    if t1 < end:
      newseg = array([c, t1, end, end-t1])
      notlist.append(newseg)
    return segmentData(notlist, start, end)

  def getUnvetoed(self, vetodata):
    '''Returns a segmentData object with unvetoed segments'''
    newdata = []
    noveto = vetodata.notSegments(self._start, self._end)
    unvetoed = self.intersectSegments(noveto)
    return unvetoed

  def unionSegments(self, other):
    '''Take the union of segment lists. This is used for combining veto segments of different detectors.'''
    '''Here I use AUB = not(notA^notB)'''
    start = min(self._start, other._start)
    end = max(self._end, other._end)
    united = self.notSegments(start, end).intersectSegments(other.notSegments(start, end)).notSegments(start, end)
    return united

  def printToFile(self, outfile):
    print "Printing segment list to file " + outfile
    savetxt(outfile, array(self._seglist), fmt='%i')
    


def getDoubles(ifo1, ifo2):
  '''Combines 2 interferometer objects into a new one, taking intersection of segments and union of vetoes'''
  name1 = ifo1._name
  name2 = ifo2._name
  seg1 = ifo1._segments
  seg2 = ifo2._segments
  vet1 = ifo1._vetoes
  vet2 = ifo2._vetoes
  unv1 = ifo1._unvetoed
  unv2 = ifo2._unvetoed
  minl1 = ifo1._minlen
  minl2 = ifo2._minlen
  
  name12 = ifo1._name + ifo2._name
  seg12 = seg1.intersectSegments(seg2)
  vet12 = vet1.unionSegments(vet2)
#  unv12 = unv1.intersectSegments(unv2)
  minl12 = max(minl1, minl2)

  ifo12 = IFO(name12, seg12, vet12, minl12)
#  if ifo12
  return ifo12
  

def getTriples(ifo1, ifo2, ifo3):
  '''Combines 3 IFO objects into one'''
  ifo12 = getDoubles(ifo1, ifo2)
  ifo123 = getDoubles(ifo12, ifo3)
  return ifo123
