# $Id$
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

import copy

from pylal import date
from pylal import SearchSummaryUtils
from glue.ligolw import ligolw
from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw import utils
from glue.ligolw.utils import ligolw_add
from glue import iterutils
from glue import segments

try:
  all
except NameError:
  from iterutils import all

#
# =============================================================================
#
#                                   Input
#
# =============================================================================
#


def ReadSnglInspiralFromFiles(fileList, mangle_event_id=False, verbose=False, non_lsc_tables_ok=False, old_document=False):
  """
  Read the SnglInspiralTables from a list of files

  @param fileList: list of input files
  @param mangle_event_id: ID remapping is necessary in cases where multiple
    files might have event_id collisions (ex: exttrig injections)
  @param verbose: print ligolw_add progress
  """
  # turn on ID remapping if necessary
  if mangle_event_id or old_document:
    next_id_orig = lsctables.SnglInspiralTable.next_id
    lsctables.SnglInspiralTable.next_id = lsctables.SnglInspiralID_old(0)
  if old_document:
    event_id_orig = lsctables.SnglInspiralTable.validcolumns["event_id"]
    lsctables.SnglInspiralTable.validcolumns["event_id"] = "int_8s"

  # ligolw_add will merge all tables, which is overkill, but merge time is
  # much less than I/O time.
  xmldoc = ligolw_add.ligolw_add(ligolw.Document(), fileList, non_lsc_tables_ok=non_lsc_tables_ok,  verbose=verbose)

  # extract the SnglInspiral table
  try:
    snglInspiralTriggers = table.get_table(xmldoc, \
      lsctables.SnglInspiralTable.tableName)
  except:
    snglInspiralTriggers = None

  # return ID remapping to its normal state (off)
  if mangle_event_id or old_document:
    lsctables.SnglInspiralTable.next_id = next_id_orig
  if old_document:
    lsctables.SnglInspiralTable.validcolumns["event_id"] = event_id_orig

  return snglInspiralTriggers


def ReadSnglInspiralSlidesFromFiles(fileList, shiftVector, vetoFile=None,
  mangleEventId=False, verbose=False, non_lsc_tables_ok=False):
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
    fileList, verbose=verbose, mangle_event_id=mangleEventId, non_lsc_tables_ok=non_lsc_tables_ok )
  if inspTriggers:
    # get the rings
    segDict = SearchSummaryUtils.GetSegListFromSearchSummaries(fileList)
    rings = segments.segmentlist(iterutils.flatten(segDict.values()))
    rings.sort()

    # perform the veto
    if vetoFile is not None:
      segList = segmentsUtils.fromsegwizard(open(vetoFile))
      inspTriggers = inspTriggers.veto(segList)

    # now slide all the triggers within their appropriate ring
    slideTriggersOnRingWithVector(inspTriggers, shiftVector, rings)

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


def CompareSnglInspiral(a, b, twindow = date.LIGOTimeGPS(0)):
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
  assert time in ring
  # use ring[0] as an epoch, do arithmetic using floats relative to epoch
  return ring[0] + (float(time - ring[0]) + shift) % float(abs(ring))


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
    end_time = trigger.get_end()
    trigger.set_end(slideTimeOnRing(end_time, shifts[trigger.ifo], rings[rings.find(end_time)]))

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
  slideTriggersOnRings(triggerList, rings, dict((ifo, -shift) for ifo, shift in shifts.items()))

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
     end_time = trigger.get_end()
     trigger.set_end(slideTimeOnRing(end_time, trigger.get_slide_number() * shiftVector[trigger.ifo], rings[rings.find(end_time)]))

def slideSegListDictOnRing(ring, seglistdict, shifts):
  """
   Return seglistdict with segments that are slid by appropriate values of
   shifts along the ring segment by the algorithm given in XXX.

   @param ring:        segment on which to cyclicly slide segments in
                       seglistdict
   @param seglistdict: segments to be slid on ring
   @param shifts:      a dictionary of the time-shifts keyed by IFO
  """
  # don't do this in loops
  ring_duration = float(abs(ring))
  
  # automate multi-list arithmetic
  ring = segments.segmentlistdict.fromkeys(seglistdict.keys(), segments.segmentlist([ring]))

  # make a copy, and extract the segments that are in the ring.
  seglistdict = seglistdict & ring

  # apply the shift vector.  the shift vector is first normalized so that
  # all the shifts are non-negative and less than the duration of the ring.
  seglistdict.offsets.update(dict((instrument, shift % ring_duration) for instrument, shift in shifts.items()))

  # split the result into the pieces that are still in the ring, and the
  # pieces that have fallen off the edge.  both retain the shift vector in
  # the offsets attribute.
  extra = seglistdict - ring
  seglistdict &= ring

  # wrap the fallen-off pieces around the ring.
  for instrument in extra.keys():
    extra.offsets[instrument] -= ring_duration

  # return the union of the results.  the shifts vector is retained in the
  # offsets attribute of the result
  return seglistdict | extra


def compute_thinca_livetime(on_instruments, off_instruments, rings, vetoseglistdict, offsetvectors):
  """
  @on_instruments is an iterable of the instruments that must be on.

  @off_instruments is an iterable of the instruments that must be off.

  on_instruments and off_instruments must be disjoint.

  @rings is a list of segments defining the analysis ring boundaries.  They
  can overlap, and do not need to be ordered.

  @vetoseglistdict is a coalesced glue.segments.segmentlistdict object
  providing the veto segments for whatever instruments have vetoes defined
  for them.  This can include veto lists for instruments other than those
  listed in on_ and off_instruments (extra veto lists will be ignored), and
  it need not provide lists for all instruments (instruments for which
  there are no veto segment lists are assumed to be on at all times).

  @offsetvectors is an iterable of dictionaries of instrument-offset pairs.
  Each dictionary must contain entries for all instruments in the union of
  on_instruments and off_instruments (it is allowed to name others as well,
  but they will be ignored).  An example of one dictionary of
  instrument-offset pairs:  {"H1": 0.0, "H2": 5.0, "L1": 10.0}.

  The return value is a float giving the livetime in seconds.
  """
  # local copies so they can be modified and iterated over more than once
  # (in case generator expressions have been passed in)
  on_instruments = set(on_instruments)
  off_instruments = set(off_instruments)

  # check that the on and off instruments are disjoint
  if on_instruments & off_instruments:
    raise ValueError, "on_instruments and off_instruments not disjoint"

  # instruments that are not vetoed are assumed to be on
  on_instruments &= set(vetoseglistdict.keys())

  # performance aid:  only need offsets for instruments whose state is
  # important
  all_instruments = on_instruments | off_instruments
  offsetvectors = tuple(dict((key, value) for key, value in offsetvector.items() if key in all_instruments) for offsetvector in offsetvectors)

  # check that each offset vector provides values for all instruments of
  # interest
  for offsetvector in offsetvectors:
    if not set(offsetvector.keys()).issuperset(all_instruments):
      raise ValueError, "incomplete offset vector %s;  missing instrument(s) %s" % (repr(offsetvector), ", ".join(all_instruments - set(offsetvector.keys())))

  # the livetime is trivial if an instrument that must be off is never
  # vetoed
  if not set(vetoseglistdict.keys()).issuperset(off_instruments):
    return 0.0

  # performance aid:  don't need veto segment lists for instruments whose
  # state is unimportant, nor veto segments that don't intersect the rings
  coalesced_rings = segments.segmentlist(rings).coalesce()
  vetoseglistdict = segments.segmentlistdict((key, segments.segmentlist(seg for seg in seglist if coalesced_rings.intersects_segment(seg))) for key, seglist in vetoseglistdict.items() if key in all_instruments)

  # tot up the time when exactly the instruments that must be on are on
  live_time = 0.0
  for ring in rings:
    # don't do this in loops
    ring = segments.segmentlist([ring])

    # performance aid:  this is done in the loop, inside
    # slideSegListDictOnRing(), but we can make that go faster by doing it
    # here first
    clipped_vetoseglistdict = segments.segmentlistdict((key, seglist & ring) for key, seglist in vetoseglistdict.items())

    # performance aid:  if an instrument that must be vetoed is never
    # vetoed in this ring, the livetime is zero
    if not all(clipped_vetoseglistdict[key] for key in off_instruments):
      continue

    # iterate over offset vectors
    for offsetvector in offsetvectors:
      # apply the offset vector to the vetoes, wrapping around the ring
      slidvetoes = slideSegListDictOnRing(ring[0], clipped_vetoseglistdict, offsetvector)

      # slidvetoes = times when instruments are vetoed,
      # slidvetoes.union(on_instruments) = times when an instrument that
      # must be on is vetoed
      #
      # ~slidvetoes = times when instruments are not vetoed,
      # (~slidvetoes).union(off_instruments) = times when an instrument
      # that must be off is not vetoed
      live_time += float(abs(ring - slidvetoes.union(on_instruments) - (~slidvetoes).union(off_instruments)))

  # done
  return live_time
