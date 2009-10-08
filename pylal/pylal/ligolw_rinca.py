# $Id$
#
# Copyright (C) 2008  Kipp C. Cannon, Drew G. Keppel
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


import bisect
import math
import sys


from glue import iterutils
from glue.ligolw import lsctables
from glue.ligolw.utils import process as ligolw_process
from pylal import llwapp
from pylal import snglcoinc
from pylal.xlal import tools as xlaltools
from pylal.xlal.datatypes.ligotimegps import LIGOTimeGPS
try:
	all
except NameError:
	# Python < 2.5.x
	from glue.iterutils import all as all


__author__ = "Kipp Cannon <kipp@gravity.phys.uwm.edu>"
__version__ = "$Revision$"[11:-2]
__date__ = "$Date$"[7:-2]


#
# =============================================================================
#
#                                 Speed Hacks
#
# =============================================================================
#


#
# Use C row classes for memory and speed
#


lsctables.CoincMapTable.RowType = lsctables.CoincMap = xlaltools.CoincMap


#
# Construct a subclass of the C sngl_ringdown row class with the methods
# that are needed
#


class SnglRingdown(xlaltools.SnglRingdownTable):
	__slots__ = ()

	def get_start(self):
		return LIGOTimeGPS(self.start_time, self.start_time_ns)

	def set_start(self, gps):
		self.start_time, self.start_time_ns = gps.seconds, gps.nanoseconds

	def __cmp__(self, other):
		# compare self's start time to the LIGOTimeGPS instance
		# other.  allows bisection searches by GPS time to find
		# ranges of triggers quickly
		return cmp(self.start_time, other.seconds) or cmp(self.start_time_ns, other.nanoseconds)


#
# Use C LIGOTimeGPS type
#


lsctables.LIGOTimeGPS = LIGOTimeGPS


#
# Use C segments module
#


def use___segments(modulename):
	from glue import __segments
	modulename.segments.infinity = __segments.infinity
	modulename.segments.NegInfinity = __segments.NegInfinity
	modulename.segments.PosInfinity = __segments.PosInfinity
	modulename.segments.segment = __segments.segment
	modulename.segments.segmentlist = __segments.segmentlist
use___segments(llwapp)
use___segments(lsctables)


#
# =============================================================================
#
#                           Add Process Information
#
# =============================================================================
#


process_program_name = "ligolw_rinca"


def append_process(xmldoc, comment = None, force = None, ds_sq_threshold = None, save_small_coincs = None, vetoes_name = None, verbose = None):
	process = llwapp.append_process(xmldoc, program = process_program_name, version = __version__, cvs_repository = u"lscsoft", cvs_entry_time = __date__, comment = comment)

	params = [
		(u"--ds-sq-threshold", u"real_8", ds_sq_threshold)
	]
	if comment is not None:
		params += [(u"--comment", u"lstring", comment)]
	if force is not None:
		params += [(u"--force", None, None)]
	if save_small_coincs is not None:
		params += [(u"--save-small-coincs", None, None)]
	if vetoes_name is not None:
		params += [(u"--vetoes-name", u"lstring", vetoes_name)]
	if verbose is not None:
		params += [(u"--verbose", None, None)]

	ligolw_process.append_process_params(xmldoc, process, params)

	return process


#
# =============================================================================
#
#                          CoincTables Customizations
#
# =============================================================================
#


#
# The sngl_ringdown <--> sngl_ringdown coinc type.
#


RingdownCoincDef = lsctables.CoincDef(search = u"ring", search_coinc_type = 0, description = u"sngl_ringdown<-->sngl_ringdown coincidences")


#
# Custom snglcoinc.CoincTables subclass.
#


class RingdownCoincTables(snglcoinc.CoincTables):
	def __init__(self, xmldoc, vetoes = None, program = u"ring"):
		snglcoinc.CoincTables.__init__(self, xmldoc)

		#
		# create a string uniquifier
		#

		self.uniquifier = {}

		#
		# find the coinc_ringdown table or create one if not found
		#

		try:
			self.coinc_ringdown_table = lsctables.table.get_table(xmldoc, lsctables.CoincRingdownTable.tableName)
		except ValueError:
			self.coinc_ringdown_table = lsctables.New(lsctables.CoincRingdownTable)
			xmldoc.childNodes[0].appendChild(self.coinc_ringdown_table)

		#
		# extract the coalesced out segment lists from lalapps_ring
		#

		self.seglists = llwapp.segmentlistdict_fromsearchsummary(xmldoc, program = program).coalesce()
		if vetoes is not None:
			self.seglists -= vetoes

	def append_coinc(self, process_id, time_slide_id, coinc_def_id, events):
		#
		# populate the coinc_event and coinc_event_map tables
		#

		coinc = snglcoinc.CoincTables.append_coinc(self, process_id, time_slide_id, coinc_def_id, events)

		#
		# populate the coinc_ringdown table:
		#
		# - start_time is the SNR-weighted mean of the start times
		# - frequency and quality are the SNR-weighted means of the
		#   frequencies and qualities
		# - snr is root-sum-square of SNRs
		# - false-alarm rate is blank
		#

		coinc_ringdown = self.coinc_ringdown_table.RowType()
		coinc_ringdown.coinc_event_id = coinc.coinc_event_id
		coinc_ringdown.snr = sum(event.snr**2. for event in events)**.5
		coinc_ringdown.false_alarm_rate = None
		# do time arithmetic using floats relative to epoch
		coinc_ringdown.set_start(events[0].get_start() + sum(event.snr * float(event.get_start() - events[0].get_start()) for event in events) / sum(event.snr for event in events))
		coinc_ringdown.set_ifos(event.ifo for event in events)
		coinc_ringdown.frequency = sum(event.snr * event.frequency for event in events) / sum(event.snr for event in events)
		coinc_ringdown.quality = sum(event.snr * event.quality for event in events) / sum(event.snr for event in events)
		self.coinc_ringdown_table.append(coinc_ringdown)

		#
		# record the instruments that were on at the time of the
		# coinc.  note that the start time of the coinc must be
		# unslid to compare with the instrument segment lists
		#

		tstart = coinc_ringdown.get_start()
		coinc.set_instruments(instrument for instrument, segs in self.seglists.items() if tstart - self.time_slide_index[time_slide_id][instrument] in segs)

		#
		# save memory by re-using strings
		#

		coinc.instruments = self.uniquifier.setdefault(coinc.instruments, coinc.instruments)
		coinc_ringdown.ifos = self.uniquifier.setdefault(coinc_ringdown.ifos, coinc_ringdown.ifos)

		#
		# done
		#

		return coinc


#
# =============================================================================
#
#                            Event List Management
#
# =============================================================================
#


class RingdownEventList(snglcoinc.EventList):
	"""
	A customization of the EventList class for use with the ringdown
	search.
	"""
	def make_index(self):
		"""
		Sort events by start time so that a bisection search can
		retrieve them.  Note that the bisection search relies on
		the __cmp__() method of the SnglRingdown row class having
		previously been set to compare the event's start time to a
		LIGOTimeGPS.
		"""
		self.sort(lambda a, b: cmp(a.start_time, b.start_time) or cmp(a.start_time_ns, b.start_time_ns))

	def set_dt(self, dt):
		"""
		If an event's start time differs by more than this many
		seconds from the start time of another event then it is
		*impossible* for them to be coincident.
		"""
		# add 1% for safety, and pre-convert to LIGOTimeGPS to
		# avoid doing type conversion in loops
		self.dt = LIGOTimeGPS(dt * 1.01)

	def _add_offset(self, delta):
		"""
		Add an amount to the start time of each event.
		"""
		for event in self:
			event.set_start(event.get_start() + delta)

	def get_coincs(self, event_a, ds_sq_threshold, comparefunc):
		#
		# event_a's start time
		#

		start = event_a.get_start()

		#
		# extract the subset of events from this list that pass
		# coincidence with event_a (use bisection searches for the
		# minimum and maximum allowed start times to quickly identify
		# a subset of the full list)
		#

		return [event_b for event_b in self[bisect.bisect_left(self, start - self.dt) : bisect.bisect_right(self, start + self.dt)] if not comparefunc(event_a, event_b, ds_sq_threshold)]


#
# =============================================================================
#
#                              Coincidence Tests
#
# =============================================================================
#


def ringdown_max_dt(events, ds_sq_threshold):
	"""
	Given a ds_sq threshold and a list of sngl_ringdown events,
	return the greatest \Delta t that can separate two events and they
	still be considered coincident.
	"""
	# for each instrument present in the event list, compute the
	# largest \Delta t interval for the events from that instrument,
	# and return the sum of the largest two such \Delta t's.

	# FIXME: get these from somewhere else
	LAL_REARTH_SI = 6.378140e6 # m
	LAL_C_SI = 299792458 # m s^-1

	return sum(sorted(max(xlaltools.XLALRingdownTimeError(event, ds_sq_threshold) for event in events if event.ifo == instrument) for instrument in set(event.ifo for event in events))[-2:]) + 2. * LAL_REARTH_SI / LAL_C_SI


def ringdown_coinc_compare(a, b, ds_sq_threshold):
	"""
	Returns False (a & b are coincident) if they pass the metric
	rinca test.
	"""
	try:
		# FIXME:  should it be ">" or ">="?
		return xlaltools.XLAL3DRinca(a, b) > ds_sq_threshold
	except ValueError:
		# ds_sq test failed to converge == events are not
		# coincident
		return True


#
# =============================================================================
#
#                                 Library API
#
# =============================================================================
#


def replicate_threshold(ds_sq_threshold, instruments):
	"""
	From a single threshold and a list of instruments, return a
	dictionary whose keys are every instrument pair (both orders), and
	whose values are all the same single threshold.

	Example:

	>>> replicate_threshold(6, ["H1", "H2"])
	{("H1", "H2"): 6, ("H2", "H1"): 6}
	"""
	instruments = sorted(instruments)
	thresholds = dict((pair, ds_sq_threshold) for pair in iterutils.choices(instruments, 2))
	instruments.reverse()
	thresholds.update(dict((pair, ds_sq_threshold) for pair in iterutils.choices(instruments, 2)))
	return thresholds


def ligolw_rinca(
	xmldoc,
	process_id,
	EventListType,
	CoincTables,
	coinc_definer_row,
	event_comparefunc,
	thresholds,
	ntuple_comparefunc = lambda events: False,
	small_coincs = False,
	veto_segments = None,
	verbose = False
):
	#
	# prepare the coincidence table interface.
	#

	if verbose:
		print >>sys.stderr, "indexing ..."
	coinc_tables = CoincTables(xmldoc, vetoes = veto_segments)
	coinc_def_id = llwapp.get_coinc_def_id(xmldoc, coinc_definer_row.search, coinc_definer_row.search_coinc_type, create_new = True, description = coinc_definer_row.description)
	sngl_index = dict((row.event_id, row) for row in lsctables.table.get_table(xmldoc, lsctables.SnglRingdownTable.tableName))

	#
	# build the event list accessors, populated with events from those
	# processes that can participate in a coincidence.  apply vetoes by
	# removing events from the lists that fall in vetoed segments
	#

	eventlists = snglcoinc.make_eventlists(xmldoc, EventListType, lsctables.SnglRingdownTable.tableName)
	if veto_segments is not None:
		for eventlist in eventlists.values():
			iterutils.inplace_filter((lambda event: event.ifo not in veto_segments or event.get_start() not in veto_segments[event.ifo]), eventlist)

	#
	# set the \Delta t parameter on all the event lists
	#

	max_dt = ringdown_max_dt(lsctables.table.get_table(xmldoc, lsctables.SnglRingdownTable.tableName), thresholds)
	if verbose:
		print >>sys.stderr, "event bisection search window will be %.16g s" % max_dt
	for eventlist in eventlists.values():
		eventlist.set_dt(max_dt)

	#
	# replicate the ds_sq threshold for every possible instrument
	# pair
	#

	thresholds = replicate_threshold(thresholds, set(eventlists))

	#
	# construct offset vector assembly graph
	#

	time_slide_graph = snglcoinc.TimeSlideGraph(coinc_tables.time_slide_index, verbose = verbose)

	#
	# loop over the items in time_slide_graph.head, producing all of
	# those n-tuple coincidences
	#

	if verbose:
		print >>sys.stderr, "constructing coincs for target offset vectors ..."
	for n, node in enumerate(time_slide_graph.head):
		if verbose:
			print >>sys.stderr, "%d/%d: %s" % (n + 1, len(time_slide_graph.head), ", ".join(("%s = %+.16g s" % x) for x in sorted(node.offset_vector.items())))
		for coinc in node.get_coincs(eventlists, event_comparefunc, thresholds, verbose):
			ntuple = [sngl_index[id] for id in coinc]
			if not ntuple_comparefunc(ntuple):
				coinc_tables.append_coinc(process_id, node.time_slide_id, coinc_def_id, ntuple)
		if small_coincs:
			for coinc in node.unused_coincs:
				ntuple = [sngl_index[id] for id in coinc]
				if not ntuple_comparefunc(ntuple):
					coinc_tables.append_coinc(process_id, node.time_slide_id, coinc_def_id, ntuple)

	#
	# remove time offsets from events
	#

	eventlists.remove_offsetdict()

	#
	# done
	#

	return xmldoc
