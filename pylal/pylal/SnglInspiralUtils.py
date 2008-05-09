# $Id$
#
# Copyright (C) 2006  Duncan A. Brown
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#

from pylal.date import LIGOTimeGPS
from glue.ligolw import ligolw
from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw import utils
from glue.ligolw.utils import ligolw_add
from glue import segments

#
# =============================================================================
#
#                                   Input
#
# =============================================================================
#


def ReadSnglInspiralFromFiles(fileList, mangle_event_id=False, verbose=False):
  """
  Read the SnglInspiralTables from a list of files

  @param fileList: list of input files
  @param mangle_event_id: ID remapping is necessary in cases where multiple
    files might have event_id collisions (ex: exttrig injections)
  @param verbose: print ligolw_add progress
  """
  # turn on ID remapping if necessary
  if mangle_event_id:
    lsctables.SnglInspiralTable.next_id = lsctables.SnglInspiralID_old(0)

  # ligolw_add will merge all tables, which is overkill, but merge time is
  # much less than I/O time.
  xmldoc = ligolw_add.ligolw_add(ligolw.Document(), fileList, verbose=verbose)

  # extract the SnglInspiral table
  try:
    snglInspiralTriggers = table.get_table(xmldoc, \
      lsctables.SnglInspiralTable.tableName)
  except:
    snglInspiralTriggers = None

  # return ID remapping to its normal state (off)
  lsctables.SnglInspiralTable.next_id = None

  return snglInspiralTriggers


def ReadSnglInspiralSlidesFromFiles(fileList, shiftVector, vetoFile=None,
  mangleEventId=False, verbose=False):
  """
  Function for reading time-slided single inspiral triggers
  with automatic resliding the times, given a list of input files.

  @param fileList: List of files containing single inspiral time-slided
                   triggers
  @param shiftVector: Dictionary of time shifts to apply to triggers
                      keyed by IFO  
  @param vetoFile: segwizard formatted file used to veto all triggers
  @param mangleEventId: ID remapping is necessary in cases where multiple
    files might have event_id collisions (ex: exttrig injections)
  @param verbose: print ligolw_add progress
  """

  # read raw triggers
  inspTriggers = ReadSnglInspiralFromFiles(\
    fileList, verbose=verbose, mangle_event_id=mangleEventId )

  # get the rings
  segDict = SearchSummaryUtils.GetSegListFromSearchSummaries(fileList)
  rings = segments.segmentlist(iterutils.flatten(segDict.values()))
  rings.sort()

  # perform the veto
  if vetoFile is not None:
    segList = segmentsUtils.fromsegwizard(open(vetoFile))
    inspTriggers = inspTriggers.veto(segList)

  # now slide all the triggers within their appropriate ring
  slideTriggersOnRingWithVector(inspTriggers, rings, shiftVector)

  # return the re-slided triggers
  return inspTriggers


#
# =============================================================================
#
#                                  Clustering
#
# =============================================================================
#

def CompareSnglInspiralByEndTime(a, b):
  """
  Orders a and b by peak time.
  """
  return cmp(a.get_end(), b.get_end())


def CompareSnglInspiralBySnr(a, b):
  """
  Orders a and b by peak time.
  """
  return cmp(a.snr, b.snr)


def CompareSnglInspiral(a, b, twindow = LIGOTimeGPS(0)):
  """
  Returns 0 if a and b are less than twindow appart.
  """
  tdiff = abs(a.get_end() - b.get_end())
  if tdiff < twindow:
    return 0
  else:
    return cmp(a.get_end(), b.get_end())

#
# =============================================================================
#
#                              Time Sliding on Ring
#
# =============================================================================
#

def slideTimeOnRing(time, shift, ring):
  """
  Return time after adding shift, but constrained to lie along the ring
  """
  newTime = time + shift
  if shift > 0:
    newTime = ring[0] + (newTime - ring[0]) % abs(ring)
  if shift < 0:
    newTime = ring[1] - (ring[1] - newTime) % abs(ring)

  return newTime


def slideTriggersOnRings(triggerList, rings, shifts):
  """
   In-place modify trigger_list so that triggers are slid by appropriate value
   of shifts along their enclosing ring segment by the algorithm given in XXX.
   This function calls the function slideTimeOnRing

   @param triggerList: a SnglInspiralTable
   @param rings:       sorted segment list of possible rings
   @param shifts:      a dictionary of the time-shifts keyed by IFO
  """
  for trigger in triggerList:
    oldTime = LIGOTimeGPS(trigger.end_time,trigger.end_time_ns)
    ringIdx = rings.find(oldTime)
    ring = rings[ringIdx]

    # check that trigger was represented in rings, actually
    if ringIdx == len(rings):
      print >> sys.stderr, "DireErrorAboutRounding"
      sys.exit(1)
    # algorithm sanity check: trigger in chosen ring
    assert trigger.end_time in ring

    newTime = slideTimeOnRing(oldTime, shifts[trigger.ifo], ring)
    trigger.end_time = newTime.seconds
    trigger.end_time_ns = newTime.nanoseconds


def unslideTriggersOnRings(triggerList, rings, shifts):
  """
   In-place modify trigger_list so that triggers are unslid by appropriate
   value of shifts along their enclosing ring segment by the algorithm given in
   XXX.
   This function calls the function slideTriggersOnRing

   @param triggerList: a SnglInspiralTable
   @param rings:       sorted segment list of possible rings
   @param shifts:      a dictionary of the time-shifts keyed by IFO
  """
  negativeShifts = dict([(ifo, -shift) for ifo,shift in shifts.itervalues()])
  slideTriggersOnRings(triggerList, rings, negativeShifts)


def slideTriggersOnRingWithVector(triggerList, shiftVector, rings):
   """
   In-place modify trigger_list so that triggers are slid by
   along their enclosing ring segment by the algorithm given in XXX.
   Slide numbers are extracted from the event_id of each trigger,
   and multiplied by the corresponding (ifo-keyed) entry in shift_vector
   to get the total slide amount.
   This function is called by ReadSnglInspiralSlidesFromFiles and
   calls the function slideTimeOnRing

   @param triggerList: a SnglInspiralTable
   @param shiftVector: a dictionary of the unit time-shift vector,
                       keyed by IFO
   @param rings:       sorted segment list of possible rings
   """
   for trigger in triggerList:
       # get shift
       shift = trigger.get_slide_number() * shiftVector[trigger.ifo]

       # get ring
       ringIndex = rings.find(trigger.end_time)
       ring = rings[ringIndex]

       # check that trigger was represented in rings, actually
       if ringIndex == len(rings):
         print >> sys.stderr, "DireErrorAboutRounding"
         sys.exit(1)
       # algorithm sanity check: trigger in chosen ring
       assert trigger.end_time in ring

       # perform shift
       oldTime = date.LIGOTimeGPS(trigger.end_time, trigger.end_time_ns)
       newTime = slide_time_on_ring(oldTime, shift, ring)
       trigger.end_time = newTime.seconds
       trigger.end_time_ns = newTime.nanoseconds


def slideSegListDictOnRing(ring, seglistdict, shifts):
  """
   Return seglistdict with segments that are slid by appropriate values of
   shifts along the ring segment by the algorithm given in XXX.
   This function calls the function slideTimeOnRing

   @param ring:        segment on which to cyclicly slide segments in
                       seglistdict
   @param seglistdict: segments to be slid on ring
   @param shifts:      a dictionary of the time-shifts keyed by IFO
  """
  # calculate start and end of ring as well as duration
  start,end = ring
  dur = end - start

  # create a new seglistdict so we don't modify the old one
  slidseglistdict = seglistdict.copy()

  # loop over ifos
  for key in slidseglistdict.keys():
    # shift seglistdict to have segment start at 0
    slidseglistdict[key] = slidseglistdict[key].shift(shifts[key] - start)

    # keep track of whether a segment needs to be split at the end
    splitseg = 0

    tmpslidseglist = segments.segmentlist([])
    for seg in slidseglistdict[key]:
      if shifts[key] >= 0:
        while seg[0] >= dur:
          # segment start has slid past end of ring, subtract duration of
          # ring until start of segment is in ring
          seg = seg.shift(-dur)

        if seg[1] > dur:
          # segment end is after end of ring, we need to split this segment
          # so keep track of how much overflow there is and contract/shift
          # segment so that it ends at the end of the ring
          splitseg = seg[1] - dur
          seg = seg.contract(splitseg/2.)
          seg = seg.shift(-splitseg/2.)

        # append seg to new seglist
        tmpslidseglist.append(seg)

      if shifts[key] < 0:
        while seg[1] <= 0:
          # segment end has slid past start of ring, add duration of
          # ring until end of segment is in ring
          seg = seg.shift(dur)

        if seg[0] < 0:
          # segment start is before start of ring, we need to split this
          # segment so keep track of how much underflow there is and
          # contract/shift segment so that it starts at the start of the ring
          splitseg = seg[0]
          seg = seg.contract(-splitseg/2.)
          seg = seg.shift(-splitseg/2.)

        # append seg to new seglist
        tmpslidseglist.append(seg)

    # set this ifos seglist to the new seglist
    slidseglistdict[key] = tmpslidseglist

    # if there was a segment needing to be split, append the amount
    # that was overflow/underflow as well
    if splitseg > 0:
      slidseglistdict[key].append(segments.segment(0,splitseg))
    if splitseg < 0:
      slidseglistdict[key].append(segments.segment(dur+splitseg,dur))

  # coalesce the new seglistdict
  slidseglistdict = slidseglistdict.coalesce()

  # shift new seglistdict to start of ring and return
  for key in slidseglistdict.keys():
    slidseglistdict[key] = slidseglistdict[key].shift(start)

  return slidseglistdict

