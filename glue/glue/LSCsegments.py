"""
This modules contains Classes that make it simple to handle GPS segments
"""

__author__ = "Greg Mendell: segment Classes taken from Duncan Brown's pipline.py under lalapps/src/lalapps"
__date__ = '$Date$'
__version__ = '$Revision$'[0:0]

import os
import sys
import string, re
import exceptions
import time
import random
import md5
import math

class ScienceSegment:
  """
  A ScienceSegment is a period of time where the experimenters determine
  that the inteferometer is in a state where the data is suitable for 
  scientific analysis. A science segment can have a list of AnalysisChunks
  asscociated with it that break the segment up into (possibly overlapping)
  smaller time intervals for analysis.
  """
  def __init__(self,segment):
    """
    @param segment: a tuple containing the (segment id, gps start time, gps end
    time, duration) of the segment.
    """
    self.__id = segment[0]
    self.__start = segment[1]
    self.__end = segment[2]
    self.__dur = segment[3]
    self.__chunks = []
    self.__unused = self.dur()
    self.__ifo = None

  def __getitem__(self,i):
    """
    Allows iteration over and direct access to the AnalysisChunks contained
    in this ScienceSegment.
    """
    if i < 0: raise IndexError, "list index out of range"
    return self.__chunks[i]
    
  def __len__(self):
    """
    Returns the number of AnalysisChunks contained in this ScienceSegment.
    """
    return len(self.__chunks)

  def __repr__(self):
    return '<ScienceSegment: id %d, start %d, end %d, dur %d, unused %d>' % (
    self.id(),self.start(),self.end(),self.dur(),self.__unused)

  def __cmp__(self,other):
    """
    ScienceSegments are compared by the GPS start time of the segment.
    """
    return cmp(self.start(),other.start())

  def make_chunks(self,length=0,overlap=0,play=0,sl=0,excl_play=0):
    """
    Divides the science segment into chunks of length seconds overlapped by
    overlap seconds. If the play option is set, only chunks that contain S2
    playground data are generated. If the user has a more complicated way
    of generating chunks, this method should be overriden in a sub-class.
    Any data at the end of the ScienceSegment that is too short to contain a 
    chunk is ignored. The length of this unused data is stored and can be
    retrieved with the unused() method.
    @param length: length of chunk in seconds.
    @param overlap: overlap between chunks in seconds.
    @param play: 1 : only generate chunks that overlap with S2 playground data.
                 2 : as play = 1 plus compute trig start and end times to coincide
                 with the start/end of the playground
    @param sl: slide by sl seconds before determining playground data.
    @param excl_play: exclude the first excl_play second from the start and end
    of the chunk when computing if the chunk overlaps with playground.
    """
    time_left = self.dur()
    start = self.start()
    increment = length - overlap
    while time_left >= length:
      end = start + length
      if (not play) or (play and (((end-sl-excl_play-729273613) % 6370) < 
        (600+length-2*excl_play))):
	if (play == 2):
	  # calculate the start of the playground preceeding the chunk end
	  play_start = 729273613 + 6370 * \
		math.floor((end-sl-excl_play-729273613) / 6370)
 	  play_end = play_start + 600
	  trig_start = 0
	  trig_end = 0
	  if ( (play_end - 6370) > start ):
	    print "Two playground segments in this chunk:",
	    print "  Code to handle this case has not been implemented"
	    sys.exit(1)
          else:
	    if play_start > start:
	      trig_start = int(play_start)
	    if play_end < end:
	      trig_end = int(play_end)
	  self.__chunks.append(AnalysisChunk(start,end,trig_start,trig_end))
        else:
          self.__chunks.append(AnalysisChunk(start,end))
      start += increment
      time_left -= increment
    self.__unused = time_left - overlap

  def add_chunk(self,start,end,trig_start=0,trig_end=0):
    """
    Add an AnalysisChunk to the list associated with this ScienceSegment.
    @param start: GPS start time of chunk.
    @param end: GPS end time of chunk.
    @param trig_start: GPS start time for triggers from chunk
    """
    self.__chunks.append(AnalysisChunk(start,end,trig_start,trig_end))

  def unused(self):
    """
    Returns the length of data in the science segment not used to make chunks.
    """
    return self.__unused

  def set_unused(self,unused):
    """
    Set the length of data in the science segment not used to make chunks.
    """
    self.__unused = unused

  def id(self):
    """
    Returns the ID of this ScienceSegment.
    """
    return self.__id
    
  def start(self):
    """
    Returns the GPS start time of this ScienceSegment.
    """
    return self.__start

  def end(self):
    """
    Returns the GPS end time of this ScienceSegment.
    """
    return self.__end

  def set_start(self,t):
    """
    Override the GPS start time (and set the duration) of this ScienceSegment.
    @param t: new GPS start time.
    """
    self.__dur += self.__start - t
    self.__start = t

  def set_end(self,t):
    """
    Override the GPS end time (and set the duration) of this ScienceSegment.
    @param t: new GPS end time.
    """
    self.__dur -= self.__end - t
    self.__end = t

  def dur(self):
    """
    Returns the length (duration) in seconds of this ScienceSegment.
    """
    return self.__dur




class ScienceData:
  """
  An object that can contain all the science data used in an analysis. Can
  contain multiple ScienceSegments and has a method to generate these from
  a text file produces by the LIGOtools segwizard program.
  """
  def __init__(self):
    self.__sci_segs = []
    self.__file = None

  def __getitem__(self,i):
    """
    Allows direct access to or iteration over the ScienceSegments associated
    with the ScienceData.
    """
    return self.__sci_segs[i]

  def __repr__(self):
    return '<ScienceData: file %s>' % self.__file

  def __len__(self):
    """
    Returns the number of ScienceSegments associated with the ScienceData.
    """
    return len(self.__sci_segs)

  def read(self,file,min_length,slide_sec=0,buffer=0):
    """
    Parse the science segments from the segwizard output contained in file.
    @param file: input text file containing a list of science segments generated by
    segwizard.
    @param min_length: only append science segments that are longer than min_length.
    @param slide_sec: Slide each ScienceSegment by::

      delta > 0:
        [s,e] -> [s+delta,e].
      delta < 0:
        [s,e] -> [s,e-delta].

    @param buffer: shrink the ScienceSegment::

      [s,e] -> [s+buffer,e-buffer]
    """
    self.__file = file
    octothorpe = re.compile(r'\A#')
    for line in open(file):
      if not octothorpe.match(line) and int(line.split()[3]) >= min_length:
        (id,st,en,du) = map(int,line.split())

        # slide the data if doing a background estimation
        if slide_sec > 0:
          st += slide_sec
        elif slide_sec < 0:
          en += slide_sec
        du -= abs(slide_sec)

        # add a buffer
        if buffer > 0:
          st += buffer
          en -= buffer
          du -= 2*abs(buffer)

        x = ScienceSegment(tuple([id,st,en,du]))
        self.__sci_segs.append(x)

  def make_chunks(self,length,overlap=0,play=0,sl=0,excl_play=0):
    """
    Divide each ScienceSegment contained in this object into AnalysisChunks.
    @param length: length of chunk in seconds.
    @param overlap: overlap between segments.
    @param play: if true, only generate chunks that overlap with S2 playground data.
    @param sl: slide by sl seconds before determining playground data.
    @param excl_play: exclude the first excl_play second from the start and end
    of the chunk when computing if the chunk overlaps with playground.
    """
    for seg in self.__sci_segs:
      seg.make_chunks(length,overlap,play,sl,excl_play)

  def make_chunks_from_unused(self,length,trig_overlap,play=0,min_length=0,sl=0,excl_play=0):
    """
    Create an extra chunk that uses up the unused data in the science segment.
    @param length: length of chunk in seconds.
    @param trig_overlap: length of time start generating triggers before the
    start of the unused data.
    @param play: 
                - 1 : only generate chunks that overlap with S2 playground data.
                - 2 : as 1 plus compute trig start and end times to coincide
                        with the start/end of the playground
    @param min_length: the unused data must be greater than min_length to make a
    chunk.
    @param sl: slide by sl seconds before determining playground data.
    @param excl_play: exclude the first excl_play second from the start and end
    of the chunk when computing if the chunk overlaps with playground.
    """
    for seg in self.__sci_segs:
      # if there is unused data longer than the minimum chunk length
      if seg.unused() > min_length:
        start = seg.end() - length
        end = seg.end()
        middle = start + length / 2
        if (not play) or (play and (((end-sl-excl_play-729273613)%6370) < 
          (600+length-2*excl_play))):
          trig_start = end - seg.unused() - trig_overlap
	  if (play == 2):
            # calculate the start of the playground preceeding the chunk end
	    play_start = 729273613 + 6370 * \
              math.floor((end-sl-excl_play-729273613) / 6370)
 	    play_end = play_start + 600
	    trig_end = 0
	    if ( (play_end - 6370) > start ):
	      print "Two playground segments in this chunk"
	      print "  Code to handle this case has not been implemented"
	      sys.exit(1)
            else:
	      if play_start > trig_start:
	        trig_start = int(play_start)
	      if (play_end < end):
	        trig_end = int(play_end)
	      if (trig_end == 0) or (trig_end > trig_start):
                seg.add_chunk(start, end, trig_start, trig_end)
          else:
            seg.add_chunk(start, end, trig_start)
        seg.set_unused(0)

  def make_short_chunks_from_unused(
    self,min_length,overlap=0,play=0,sl=0,excl_play=0):
    """
    Create a chunk that uses up the unused data in the science segment
    @param min_length: the unused data must be greater than min_length to make a
    chunk.
    @param overlap: overlap between chunks in seconds.
    @param play: if true, only generate chunks that overlap with S2 playground data.
    @param sl: slide by sl seconds before determining playground data.
    @param excl_play: exclude the first excl_play second from the start and end
    of the chunk when computing if the chunk overlaps with playground.
    """
    for seg in self.__sci_segs:
      if seg.unused() > min_length:
        start = seg.end() - seg.unused() - overlap
        end = seg.end()
        length = start - end
        if (not play) or (play and (((end-sl-excl_play-729273613)%6370) < 
        (600+length-2*excl_play))):
          seg.add_chunk(start, end, start)
        seg.set_unused(0)

  def intersection(self, other):
    """
    Replaces the ScienceSegments contained in this instance of ScienceData
    with the intersection of those in the instance other. Returns the number
    of segments in the intersection.
    @param other: ScienceData to use to generate the intersection
    """

    # we only deal with the case of two lists here
    length1 = len(self)
    length2 = len(other)

    # initialize list of output segments
    ostart = -1
    outlist = []
    iseg2 = -1
    start2 = -1
    stop2 = -1

    for seg1 in self:
      start1 = seg1.start()
      stop1 = seg1.end()
      id = seg1.id()

      # loop over segments from the second list which overlap this segment
      while start2 < stop1:
        if stop2 > start1:
          # these overlap

          # find the overlapping range
          if start1 < start2:
            ostart = start2
          else:
            ostart = start1
          if stop1 > stop2:
            ostop = stop2
          else:
            ostop = stop1

          x = ScienceSegment(tuple([id, ostart, ostop, ostop-ostart]))
          outlist.append(x)

          if stop2 > stop1:
            break

        # step forward
        iseg2 += 1
        if iseg2 < len(other):
          seg2 = other[iseg2]
          start2 = seg2.start()
          stop2 = seg2.end()
        else:
          # pseudo-segment in the far future
          start2 = 2000000000
          stop2 = 2000000000

    # save the intersection and return the length
    self.__sci_segs = outlist
    return len(self)

  

  def union(self, other):
    """
    Replaces the ScienceSegments contained in this instance of ScienceData
    with the union of those in the instance other. Returns the number of
    ScienceSegments in the union.
    @param other: ScienceData to use to generate the intersection
    """

    # we only deal with the case of two lists here
    length1 = len(self)
    length2 = len(other)

    # initialize list of output segments
    ostart = -1
    seglist = []

    i1 = -1
    i2 = -1
    start1 = -1
    start2 = -1
    id = -1
    
    while 1:
      # if necessary, get a segment from list 1
      if start1 == -1:
        i1 += 1
        if i1 < length1:
          start1 = self[i1].start()
          stop1 = self[i1].end()
          id = self[i1].id()
        elif i2 == length2:
          break

      # if necessary, get a segment from list 2
      if start2 == -1:
        i2 += 1
        if i2 < length2:
          start2 = other[i2].start()
          stop2 = other[i2].end()
        elif i1 == length1:
          break

      # pick the earlier segment from the two lists
      if start1 > -1 and ( start2 == -1 or start1 <= start2):
        ustart = start1
        ustop = stop1
        # mark this segment has having been consumed
        start1 = -1
      elif start2 > -1:
        ustart = start2
        ustop = stop2
        # mark this segment has having been consumed
        start2 = -1
      else:
        break

      # if the output segment is blank, initialize it; otherwise, see
      # whether the new segment extends it or is disjoint
      if ostart == -1:
        ostart = ustart
        ostop = ustop
      elif ustart <= ostop:
        if ustop > ostop:
          # this extends the output segment
          ostop = ustop
        else:
          # This lies entirely within the current output segment
          pass
      else:
         # flush the current output segment, and replace it with the
         # new segment
         x = ScienceSegment(tuple([id,ostart,ostop,ostop-ostart]))
         seglist.append(x)
         ostart = ustart
         ostop = ustop

    # flush out the final output segment (if any)
    if ostart != -1:
      x = ScienceSegment(tuple([id,ostart,ostop,ostop-ostart]))
      seglist.append(x)

    self.__sci_segs = seglist
    return len(self)


  def coalesce(self):
    """
    Coalesces any adjacent ScienceSegments. Returns the number of 
    ScienceSegments in the coalesced list.
    """

    # check for an empty list
    if len(self) == 0:
      return 0

    # sort the list of science segments
    self.__sci_segs.sort()

    # coalesce the list, checking each segment for validity as we go
    outlist = []
    ostop = -1

    for seg in self:
      start = seg.start()
      stop = seg.end()
      id = seg.id()
      if start > ostop:
        # disconnected, so flush out the existing segment (if any)
        if ostop >= 0:
          x = ScienceSegment(tuple([id,ostart,ostop,ostop-ostart]))
          outlist.append(x)
        ostart = start
        ostop = stop
      elif stop > ostop:
        # extend the current segment
        ostop = stop

    # flush out the final segment (if any)
    if ostop >= 0:
      x = ScienceSegment(tuple([id,ostart,ostop,ostop-ostart]))
      outlist.append(x)

    self.__sci_segs = outlist
    return len(self)


  def invert(self):
    """
    Inverts the ScienceSegments in the class (i.e. set NOT).  Returns the
    number of ScienceSegments after inversion.
    """

    # check for an empty list
    if len(self) == 0:
      # return a segment representing all time
      self.__sci_segs = ScienceSegment(tuple(0,0,1999999999,1999999999))

    # go through the list checking for validity as we go
    outlist = []
    ostart = 0
    for seg in self:
      start = seg.start()
      stop = seg.end()
      if start < 0 or stop < start or start < ostart:
        raise SegmentError, "Invalid list"
      if start > 0:
        x = ScienceSegment(tuple([0,ostart,start,start-ostart]))
        outlist.append(x)
      ostart = stop

    if ostart < 1999999999:
      x = ScienceSegment(tuple([0,ostart,1999999999,1999999999-ostart]))
      outlist.append(x)

    self.__sci_segs = outlist
    return len(self)

  def play(self):
    """
    Keep only times in ScienceSegments which are in the playground
    """

    length = len(self)

    # initialize list of output segments
    ostart = -1
    outlist = []
    begin_s2 = 729273613
    play_space = 6370
    play_len = 600

    for seg in self:
      start = seg.start()
      stop = seg.end()
      id = seg.id()
     
      # select first playground segment which ends after start of seg
      play_start = begin_s2+play_space*( 1 + 
	int((start - begin_s2 - play_len)/play_space) )

      while play_start < stop:
	if play_start > start:
	  ostart = play_start
	else:
	  ostart = start
        
        
	play_stop = play_start + play_len

	if play_stop < stop:
	  ostop = play_stop
	else:
	  ostop = stop

        x = ScienceSegment(tuple([id, ostart, ostop, ostop-ostart]))
        outlist.append(x)

        # step forward
	play_start = play_start + play_space 

    # save the playground segs and return the length
    self.__sci_segs = outlist
    return len(self)
