import matplotlib as mpl
import argparse
from optparse import OptionParser
import os
import sys
import random as rd
mpl.use('Agg')
from pylab import *
py_version = sys.version_info[:2]
np_version = np.__version__


class IFO:
  '''An interferometer. Can also refer to multiple interferometers to host doubles, triples etc.'''
  def __init__(self, name, segments, vetoes, minlen=1):
    self._name = name
    self._minlen = minlen

    if type(segments) is str:
      print "Reading segments for " + self._name + ""
      self._allsegments = self.readSegments(segments)
    elif type(segments) is ndarray:
      self._allsegments = segmentData(list(segments))
    elif  isinstance(segments, segmentData):
      self._allsegments = segments
    else:
      print "Cannot recognize segments!"
      return -1
    self._segments = self._allsegments.fits(self._minlen)

    if type(vetoes) is str:
      print "Reading veto segments for " + self._name + ""
      self._vetoes = self.readSegments(vetoes)
    elif type(vetoes) is ndarray:
      self._vetoes = segmentData(list(vetoes))
    elif isinstance(vetoes, segmentData):
      self._vetoes = vetoes
    else:
      print "Cannot recognize veto segments!"
      return -1

    self.setUnvetoed()
    print "Number of unvetoed segments that fit " + str(self._minlen) + ": " + str(self._unvetoed.length())

  def setUnvetoed(self):
    '''This is getting a list of unvetoed segments that fit the minimum length'''
    self._unvetoed = self._segments.getUnvetoed(self._vetoes).fits(self._minlen)
    #self._unvetoed = self._segments.fits(self._minlen).getUnvetoed(self._vetoes).fits(self._minlen)

  def readSegments(self, segfname):
    print "Reading " + segfname
    segdata = genfromtxt(segfname, dtype=int)
    return segmentData(list(segdata))

  def printUnvetoedToFile(self, outname):
    self._unvetoed.printToFile(outname)

  def getTrigTimes(self, interval=None, rmargin=2, n=None, outfile=None):
    '''Returns a list of gps times on which injections can be made'''
    trigtimes = []
    if interval is None:
      l = self._minlen
    else:
      l = interval
    for seg in self._unvetoed._seglist:
      #print seg
      t = int(seg[1])
      while t + l <= seg[2]:
#        trigtimes.append(t + seg[3]/2)
        trigtimes.append(t + l - rmargin)
        t += l
    
    if n is not None:
      trigtimes = trigtimes[:n]
    if outfile is None:
      return array(trigtimes)
    else:
      savetxt(outfile, array(trigtimes), fmt='%i')

  def plotCumulativeDurations(self, outfile, maxdur=None):
    print "Plotting segment lengths distribution to file " + outfile
    fig = figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel("segment length")
    ax.set_ylabel("# segments")
    sdur = array(self._segments._seglist)[:,-1]
    udur = array(self._unvetoed._seglist)[:,-1]
    if maxdur is not None:
      sdur = sdur[where(sdur <= maxdur)[0]]
      udur = udur[where(udur <= maxdur)[0]]
    ax.set_ylim((0,max([len(udur),len(sdur)])))
    ax.set_title("Segment length distribution for " + self._name)
    ax.hist(sdur, bins=sort(sdur), cumulative=True, histtype='stepfilled', alpha=0.3, label="All segments (" + str(self._segments.length()) + ")")
    ax.hist(udur, bins=sort(udur), cumulative=True, histtype='stepfilled', alpha=0.3, label="Unvetoed segments (" + str(self._unvetoed.length()) + ")")
    ax.legend()
    fig.savefig(outfile)
    


class segment:
  '''A time segment in an IFO'''

  def __init__(self, data, gpsstart=None, gpsend=None):
    '''Creates a segment'''
    if gpsstart is not None and gpsend is not None:
      self._id = data[0]
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

  def hasTime(self, time):
    '''Checks if time is in (open) segment'''
    return (time > self._start and time < self._end)

  def intersectWith(self, other):
    '''Intersects with another segment'''
    newstart = max(self._start, other._start)
    newend = min(self._end, other._end)
    newseg = segment([self._id], newstart, newend)
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

  def length(self):
    return len(self._seglist)
  
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
      if t2 > t1 and t1 >= start and t2 <= end:
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

  def hasTime(self, time):
    '''Checks whether time is in the segments'''
    has = False
    for s in self._seglist:
      has = has or s.hasTime(time)
      if has:
        return True
    return False

  def printToFile(self, outfile):
    print "Printing segment list to file " + outfile
    savetxt(outfile, array(self._seglist), fmt='%i')
    
  def plotCumulativeDurations(self, outfile, maxdur=None):
    print "Plotting segment lengths distribution to file " + outfile
    fig = figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel("segment length")
    ax.set_ylabel("# segments")
    dur = array(self._seglist)[:,-1]
    if maxdur is not None:
      dur = dur[where(dur <= maxdur)[0]]
    ax.hist(dur, bins=sort(dur), cumulative=True, histtype='stepfilled', alpha=0.3)
    fig.savefig(outfile)
    
    


def getDoubles(ifo1, ifo2, unvetoed=False):
  '''
  Combines 2 interferometer objects into a new one, taking intersection of segments and union of vetoes
  Setting unvetoed=True will speed up the calculations by ignoring small/vetoed segments. Use with caution!
  '''
  name1 = ifo1._name
  name2 = ifo2._name
  seg1 = ifo1._allsegments
  seg2 = ifo2._allsegments
  vet1 = ifo1._vetoes
  vet2 = ifo2._vetoes
  unv1 = ifo1._unvetoed
  unv2 = ifo2._unvetoed
  minl1 = ifo1._minlen
  minl2 = ifo2._minlen
  
  name12 = name1 + name2
  if unvetoed:
    seg12 = unv1.intersectSegments(unv2)
  else:
    seg12 = seg1.intersectSegments(seg2)
  vet12 = vet1.unionSegments(vet2)
  minl12 = max(minl1, minl2)

  ifo12 = IFO(name12, seg12, vet12, minl12)
#  if ifo12
  return ifo12
  

def getTriples(ifo1, ifo2, ifo3, unvetoed=False):
  '''Combines 3 IFO objects into one'''
  ifo12 = getDoubles(ifo1, ifo2, unvetoed)
  ifo123 = getDoubles(ifo12, ifo3, unvetoed)
  return ifo123
  
  
def generateTimeslides(tdict, n, ref=None, outfolder=None):
  '''
     Generate a list of injection times and timeslides, 
     given a dictionary of (single) injection times
  '''
  ifos = tdict.keys()
  nifo = len(ifos)
  if ref is None:
    ref = ifos[0]
  injtimes = []
  injdict = {k: [] for k in tdict.keys()}
  slidedict = {k: [] for k in tdict.keys()}
  for i in arange(n):
    for ifo in ifos:
      injdict[ifo].append(rd.sample(tdict[ifo], 1)[0]) 
  
  injtimes = injdict[ref]
  for ifo in ifos:
    slidedict[ifo] = array(injdict[ifo]) - array(injdict[ref])
  
  if outfolder is None:
    return injtimes, slidedict
  else:
    label=''.join(ifos) + '_' + str(n)
    injtimesfile = os.path.join(outfolder,'injtimes_'+label+'.dat')
    slidefile = open(os.path.join(outfolder,'timeslides_'+label+'.dat'), 'w')
    header = " ".join(ifos)
    slidefile.write('# '+header+'\n')
    slidedata = array(slidedict.values()).T
    print injtimes
    print slidedict
    savetxt(injtimesfile, array(injtimes), fmt='%i')
    savetxt(slidefile, slidedata, delimiter=' ', fmt='%i')
    slidefile.close()

if __name__ == "__main__":

###########################################
#
#  Parse arguments (argparser)
#
###########################################

  #parser = OptionParser()
  parser = argparse.ArgumentParser(description="Reads segment files and outputs unvetoed injection times.")
  
  
  parser.add_argument('-i', '--ifos', nargs=3, type=str, metavar='IFO', dest='ifos', help='List of IFOs to be considered') # FIXME: Variable length???
  parser.add_argument('-l', '--length', type=int, metavar='INT', dest='length', help='length of signal segments in seconds', default=45.0)
  parser.add_argument('-n', '--notriple', action="store_true", dest="notriple", help="restrict to doubles", default=False)
  parser.add_argument('-N', '--Ninj', type=int, dest="Ninj", help="# injections", default=1000)
  parser.add_argument('--nodouble', action="store_true", dest="nodouble", help="restrict to sinles", default=False)
  parser.add_argument('-o', '--outfile', type=str, metavar='FILE', dest='outf', help='path to the output file')
  parser.add_argument('-P', '--psdlength', type=int, metavar='INT', dest='psdlength', help='minimum length for calculating PSD in seconds', default=1024)
  parser.add_argument('-p', '--plot', action="store_true", dest="plotsegdist", help="plot cumulative segment length distribution", default=False)
  parser.add_argument('-s', '--segfiles', nargs=3, type=str, metavar='FILE', dest='segf', help='path to the segment files for IFOS in the order entered above')
  parser.add_argument('-t', '--timeslides', action="store_true", dest="timeslides", help="enable timeslides", default=False)
  parser.add_argument('-v', '--vetofiles', nargs=3, type=str, metavar='FILE', dest='vetof', help='path to the veto files for IFOS in the order entered above')
  
  args = parser.parse_args()
  
  outfile = str(args.outf)
  segfiles = args.segf
  vetofiles = args.vetof
  ifos = args.ifos
  nifos = len(ifos)
  Ninj = args.Ninj
  seglen = args.length
  psdlen = args.psdlength
  plotsegdist = args.plotsegdist
  nodouble = args.nodouble
  notriple = args.notriple
  timeslides = args.timeslides
  maxplot = 5000
  
  #parser = argparse.ArgumentParser(description='Find data segments on single/double/triple HLV time free of vetoes.')
  #parser.add_argument('integers', metavar='N', type=int, nargs='+', help='an integer for the accumulator')
  #parser.add_argument('--sum', dest='accumulate', action='store_const', const=sum, default=max, help='sum the integers (default: find the max)')
  
  #args = parser.parse_args()
  
  # Read segments data in rawlist and only keep long enough segments in fitlist for each of the IFOs
  timesdict = {}
  minlen = max(psdlen, seglen)
  print minlen, psdlen, seglen
  
  IFOlist = []
  for i in arange(nifos):
    IFOlist.append(IFO(ifos[i], segfiles[i], vetofiles[i], minlen))
    if timeslides:
      timesdict[ifos[i]] = IFOlist[-1].getTrigTimes(interval=seglen)
      print i, ifos[i]
      print shape(timesdict[ifos[i]])
  
  if not nodouble:
    doubleIFOs = []
    doubleIFOs.append( getDoubles(IFOlist[0], IFOlist[1], unvetoed=True) )
    doubleIFOs.append( getDoubles(IFOlist[1], IFOlist[2], unvetoed=True) )
    doubleIFOs.append( getDoubles(IFOlist[0], IFOlist[2], unvetoed=True) )
    if not notriple:     #FIXME: put this 2 lines above if necessary
      tripleIFO = getDoubles(doubleIFOs[0], IFOlist[2])
  
    for double in doubleIFOs:
      double.getTrigTimes(outfile='injtimes_' + double._name + '_' + str(double._minlen) +'.dat')
      if plotsegdist:
          double.plotCumulativeDurations('doubleseg_' + double._name +'.png', maxplot)
    if not notriple: 
      tripleIFO.getTrigTimes(outfile='injtimes_' + tripleIFO._name + '_' + str(tripleIFO._minlen) +'.dat')
      if plotsegdist:
        tripleIFO.plotCumulativeDurations('tripleseg_' + tripleIFO._name +'.png', maxplot)

  
  if timeslides:
    generateTimeslides(timesdict, Ninj, ref='H1', outfolder='./timeslides')
    

