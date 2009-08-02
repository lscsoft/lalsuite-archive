# $Id$
#
# Copyright (C) 2006  Kipp C. Cannon
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


"""
Generic coincidence engine for use with time-based event lists in LIGO
Light Weight XML documents.
"""


import copy
import sys
# Python 2.3 compatibility
try:
	set
except NameError:
	from sets import Set as set


from glue import iterutils
from glue.ligolw import table
from glue.ligolw import lsctables
from pylal import llwapp


__author__ = "Kipp Cannon <kipp@gravity.phys.uwm.edu>"
__version__ = "$Revision$"[11:-2]
__date__ = "$Date$"[7:-2]


#
# =============================================================================
#
#                             Time Slide Utilities
#
# =============================================================================
#


def time_slide_consideration_order(time_slide_table):
	"""
	Given a time_slide table, return a list of the unique
	time_slide_id's in the table sorted by number of participating
	instruments, then by inter-instrument deltas (then by ID string if
	there are duplicate time slide vectors in the table).

	If coincs are constructed for time slides in this order, each
	offset vector will always be considered after any offset vectors it
	contains.  For example the vector {"H1": 0, "H2": 10, "L1": 20}
	will be considered after {"H1": 0, "H2": 10} because the latter is
	contained in the former.

	If the coincidence test is such that an H1,H2,L1 triple as above
	would also be found as the H1,H2 double above, then if the offset
	vectors are processed in the order returned by this function a
	performance improvement can be achieved in the coincidence
	algorithm by constructing triples by scanning
	previously-constructed doubles for those which match the third
	instrument instead of constructing the triples from scratch from
	the raw trigger list (which involves paying the price of
	constructing the doubles again).
	"""
	#
	# sorted list of the unique instruments appearing in the offset
	# vectors
	#

	instruments = sorted(set(time_slide_table.getColumnByName("instrument")))

	#
	# all possible pairs
	#

	instruments = iterutils.choices(instruments, 2)

	#
	# the offset vectors indexed by time_slide_id
	#

	offsets = time_slide_table.as_dict()

	#
	# order IDs by number of participating instruments then by
	# inter-instrument deltas, then by ID
	#

	def deltas(offset_vector):
		for a, b in instruments:
			try:
				yield offset_vector[b] - offset_vector[a]
			except KeyError:
				yield None

	return sorted(offsets.keys(), lambda a, b: cmp(len(offsets[a]), len(offsets[b])) or cmp(tuple(deltas(offsets[a])), tuple(deltas(offsets[b]))) or cmp(a, b))


#
# =============================================================================
#
#                              Document Interface
#
# =============================================================================
#


class CoincTables(object):
	"""
	A convenience interface to the XML document's coincidence tables,
	allowing for easy addition of coincidence events.
	"""
	def __init__(self, xmldoc):
		# find the coinc table or create one if not found
		try:
			self.coinctable = table.get_table(xmldoc, lsctables.CoincTable.tableName)
		except ValueError:
			self.coinctable = lsctables.New(lsctables.CoincTable)
			xmldoc.childNodes[0].appendChild(self.coinctable)
		self.coinctable.sync_next_id()

		# find the coinc_map table or create one if not found
		try:
			self.coincmaptable = table.get_table(xmldoc, lsctables.CoincMapTable.tableName)
		except ValueError:
			self.coincmaptable = lsctables.New(lsctables.CoincMapTable)
			xmldoc.childNodes[0].appendChild(self.coincmaptable)

		# find the time_slide table
		self.time_slide_table = table.get_table(xmldoc, lsctables.TimeSlideTable.tableName)

	def get_ordered_time_slides(self):
		"""
		Returns a list of (time_slide_id, offset vector) tuples for
		the document's time slide vectors in the order defined by
		the time_slide_consideration_order() function.
		"""
		time_slides = self.time_slide_table.as_dict()
		return [(id, time_slides[id]) for id in time_slide_consideration_order(self.time_slide_table)]

	def append_coinc(self, process_id, time_slide_id, coinc_def_id, events):
		"""
		Takes a process ID, a time slide ID, and a list of events,
		and adds the events as a new coincidence to the coinc_event
		and coinc_map tables.
		"""
		coinc = self.coinctable.RowType()
		coinc.process_id = process_id
		coinc.coinc_def_id = coinc_def_id
		coinc.coinc_event_id = self.coinctable.get_next_id()
		coinc.time_slide_id = time_slide_id
		coinc.set_instruments(None)
		coinc.nevents = len(events)
		coinc.likelihood = None
		self.coinctable.append(coinc)
		for event in events:
			coincmap = self.coincmaptable.RowType()
			coincmap.coinc_event_id = coinc.coinc_event_id
			coincmap.table_name = event.event_id.table_name
			coincmap.event_id = event.event_id
			self.coincmaptable.append(coincmap)
		return coinc


#
# =============================================================================
#
#                                Process Filter
#
# =============================================================================
#


def coincident_process_ids(xmldoc, max_segment_gap, program):
	"""
	Take an XML document tree and determine the set of process IDs
	that will participate in coincidences identified by the time slide
	table therein.  It is OK for xmldoc to contain time slides
	involving instruments not represented in the list of processes:
	these time slides are ignored.  max_segment_gap is the largest gap
	(in seconds) that can exist between two segments and it still be
	possible for the two to provide events coincident with one another.
	"""
	# get the list of all process IDs for the given program
	proc_ids = table.get_table(xmldoc, lsctables.ProcessTable.tableName).get_ids_by_program(program)
	if not proc_ids:
		# hmm... that program's output is not in this file.
		raise KeyError, program

	# extract a segmentlistdict
	search_summ_table = table.get_table(xmldoc, lsctables.SearchSummaryTable.tableName)
	seglistdict = search_summ_table.get_out_segmentlistdict(proc_ids).coalesce()

	# fast path:  if the largest gap anywhere in the lists is smaller
	# than max_segment_gap then all process_ids participate.  NOTE:
	# this also handles the case of max_segment_gap being passed in as
	# float("inf") (which would otherwise break the LIGOTimeGPS
	# arithmetic).
	#
	# this is checking the gaps between segments in the *same
	# instrument*.  this assumption is that if max_segment_gap can
	# close all the gaps within each instrument's segment list then,
	# finally, all processes will be found to be required for the
	# coincidence analysis.  this is not really true, but it's safe.
	if max(b[0] - a[1] for segs in seglistdict.values() for a, b in zip(segs[:-1], segs[1:])) <= max_segment_gap:
		return proc_ids

	# protract by half the largest coincidence window so as to not miss
	# edge effects
	seglistdict.protract(max_segment_gap / 2)

	# determine what time slides are possible given the instruments in
	# the search summary table
	avail_instruments = set(seglistdict.keys())
	timeslides = [offset_vector for offset_vector in table.get_table(xmldoc, lsctables.TimeSlideTable.tableName).as_dict().values() if set(offset_vector.keys()).issubset(avail_instruments)]

	# determine the coincident segments for each instrument
	seglistdict = llwapp.get_coincident_segmentlistdict(seglistdict, timeslides)

	# find the IDs of the processes that contributed to the coincident
	# segments
	coinc_proc_ids = set()
	for row in search_summ_table:
		if row.process_id in proc_ids and row.process_id not in coinc_proc_ids and seglistdict.intersection(row.get_ifos()).intersects_segment(row.get_out()):
			coinc_proc_ids.add(row.process_id)
	return coinc_proc_ids


#
# =============================================================================
#
#                             Event List Interface
#
# =============================================================================
#


class EventList(list):
	"""
	A parent class for managing a list of events:  applying time
	offsets, and retrieving subsets of the list selected by time
	interval.  To be useful, this class must be subclassed with
	overrides provided for certain methods.  The only methods that
	*must* be overridden in a subclass are the _add_offset() and
	get_coincs() methods.  The make_index() method can be overridden if
	needed.  None of the other methods inherited from the list parent
	class need to be overridden, indeed they probably should not be
	unless you know what you're doing.
	"""
	def __init__(self, instrument):
		self.offset = lsctables.LIGOTimeGPS(0)
		self.instrument = instrument

	def make_index(self):
		"""
		Provided to allow for search-specific look-up tables or
		other indexes to be constructed for use in increasing the
		speed of the get_coincs() method.  This will be called
		after all events have been added to the list, and again if
		the list is ever modified, and before get_coincs() is ever
		called.
		"""
		pass

	def _add_offset(self, delta):
		"""
		Add an amount to the time of each event.  This method is
		used internally by the set_offset() method, and is provided
		to simplify the construction of search-specific subclasses.
		Typically, the _add_offset() method will be overridden with
		a search-specific implementation, while the set_offset()
		method is not modified (which does some additional book
		keeping, such as updating the offset attribute).
		"""
		raise NotImplementedError

	def set_offset(self, offset):
		"""
		Set an offset on the times of all events in the list.
		"""
		# cast offset to LIGOTimeGPS to avoid repeated conversion
		# when applying the offset to each event.  also prevents
		# round-off errors since addition and subtraction of
		# LIGOTimeGPS objects are exact, so by recording the offset
		# as a LIGOTimeGPS it should be possible to return the
		# events to exactly their original times before exiting.
		offset = lsctables.LIGOTimeGPS(offset)

		# check for no-op
		if offset != self.offset:
			# apply the delta
			self._add_offset(offset - self.offset)

			# record the new offset
			self.offset = offset

	def get_coincs(self, event_a, threshold, comparefunc):
		"""
		Return a list of the events from this list that are
		coincident with event_a.

		The threshold argument will be the thresholds appropriate
		for "instrument_a, instrument_b", in that order, where
		instrument_a is the instrument for event_a, and
		instrument_b is the instrument for the events in this
		EventList.
		"""
		raise NotImplementedError


class EventListDict(dict):
	"""
	A dictionary of EventList objects, indexed by instrument,
	initialized from an XML trigger table and a list of process IDs
	whose events should be included.
	"""
	def __new__(self, *args, **kwargs):
		# wrapper to shield dict.__new__() from our arguments.
		return dict.__new__(self)

	def __init__(self, EventListType, event_table, process_ids = None):
		"""
		Initialize a newly-created instance.  EventListType is a
		subclass of EventList (the subclass itself, not an instance
		of the subclass).  event_table is a list of events (e.g.,
		an instance of a glue.ligolw.table.Table subclass).  If the
		optional process_ids arguments is not None, then it is
		assumed to be a list or set or other thing implementing the
		"in" operator which is used to define the set of
		process_ids whose events should be considered in the
		coincidence analysis, otherwise all events are considered.
		"""
		for event in event_table:
			if (process_ids is None) or (event.process_id in process_ids):
				if event.ifo not in self:
					self[event.ifo] = EventListType(event.ifo)
				self[event.ifo].append(event)
		for l in self.values():
			l.make_index()

	def set_offsetdict(self, offsetdict):
		"""
		Set the event list offsets to those in the dictionary of
		instrument/offset pairs.  Instruments not in offsetdict are
		not modified.  KeyError is raised if the dictionary of
		instrument/offset pairs contains a key (instrument) that
		this dictionary does not.
		"""
		for instrument, offset in offsetdict.items():
			self[instrument].set_offset(offset)

	def remove_offsetdict(self):
		"""
		Remove the offsets from all event lists (reset them to 0).
		"""
		for l in self.values():
			l.set_offset(0)


def make_eventlists(xmldoc, EventListType, event_table_name, max_segment_gap, program_name):
	"""
	Convenience wrapper for constructing a dictionary of event lists
	from an XML document tree, the name of a table from which to get
	the events, a maximum allowed time window, and the name of the
	program that generated the events.
	"""
	return EventListDict(EventListType, table.get_table(xmldoc, event_table_name), process_ids = coincident_process_ids(xmldoc, max_segment_gap, program_name))


#
# =============================================================================
#
#                            Coincidence Iterators
#
# =============================================================================
#


def CoincidentNTuples(eventlists, comparefunc, instruments, thresholds, verbose = False):
	"""
	Given an EventListDict object, a list (or iterator) of instruments,
	and a dictionary of instrument pair thresholds, generate a sequence
	of tuples of mutually coincident events.  Each tuple returned by
	this generator will contain exactly one event from each of the
	instruments in the instrument list.
	"""
	# retrieve the event lists for the requested instrument combination

	eventlists = [(eventlists[instrument], instrument) for instrument in instruments]
	if not eventlists:
		# no instruments == nothing to do
		return

	# sort the list of event lists from longest to shortest, extract
	# the shortest list, the name of its instrument, and record its length

	eventlists.sort(lambda a, b: cmp(len(a[0]), len(b[0])), reverse = True)
	shortestlist, shortestinst = eventlists.pop()
	length = len(shortestlist)

	try:
		# pre-construct the sequence of thresholds for pairs of the
		# remaining instruments in the order in which they will be
		# used by the inner loop

		threshold_sequence = tuple([thresholds[pair] for pair in iterutils.choices([instrument for eventlist, instrument in eventlists], 2)])

		# retrieve the thresholds to be used in comparing events
		# from the shortest list to those in each of the remaining
		# event lists

		eventlists = tuple([(eventlist, thresholds[(shortestinst, instrument)]) for eventlist, instrument in eventlists])
	except KeyError, e:
		raise KeyError, "no coincidence thresholds provided for instrument pair %s" % str(e)

	# for each event in the shortest list

	for n, event in enumerate(shortestlist):
		if verbose and not (n % 500):
			print >>sys.stderr, "\t%.1f%%\r" % (100.0 * n / length),

		head = (event,)

		# iterate over n-tuples of events from the other lists that
		# are coincident with the outer event.  each n-tuple
		# contains one event from each list.

		for tail in iterutils.MultiIter(*[eventlist.get_coincs(event, threshold, comparefunc) for eventlist, threshold in eventlists]):

			# if the events in the inner n-tuple are
			# mutually-coincident, combine with the outer event
			# and report as a coincident n-tuple.

			for (a, b), threshold in zip(iterutils.choices(tail, 2), threshold_sequence):
				if comparefunc(a, b, threshold):
					# not coincident
					break
			else:
				# made it through the loop, all pairs are
				# coincident
				yield head + tail
	if verbose:
		print >>sys.stderr, "\t100.0%"
