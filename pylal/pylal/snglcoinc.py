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

import sys
# Python 2.3 compatibility
try:
	set
except NameError:
	from sets import Set as set


from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw import ilwd
from pylal import itertools
from pylal import llwapp
from pylal.date import LIGOTimeGPS

__author__ = "Kipp Cannon <kipp@gravity.phys.uwm.edu>"
__version__ = "$Revision$"[11:-2]
__date__ = "$Date$"[7:-2]


#
# =============================================================================
#
#                             Command Line Helpers
#
# =============================================================================
#


def parse_thresholds(thresholdstrings):
	"""
	Turn a list of strings of the form
	inst1,inst2=threshold1[,threshold2,...] into a dictionary with
	(inst1, inst2) 2-tuples as keys and the thresholds parsed into
	lists of individual strings as the values, with their order
	preserved.

	For each pair of instruments represented among the input strings,
	the two possible orders are considered independent:  the input
	strings are allowed to contain one set of thresholds for (inst1,
	inst2), and a different set of thresholds for (inst2, inst1).  Be
	aware that no input checking is done to ensure the user has not
	provided duplicate, incompatible, thresholds;  this is considered
	the responsibility of the application program to verify.

	The output dictionary contains threshold sets for both instrument
	orders regardless of whether or not they were supplied,
	independently, by the input strings.  If the input strings
	specified thresholds for only one of the two orders, the thresholds
	for the other order are copied from the one that was provided.

	Whitespace is removed from the start and end of all strings.

	A typical use for this function is in parsing command line
	arguments or configuration file entries.

	Example:

	>>> from pylal.snglcoinc import parse_thresholds
	>>> parse_thresholds(["H1,H2=0.1,100", "H1,L1=.2,100"])
	{('H1', 'H2'): ['0.1', '100'], ('H1', 'L1'): ['.2', '100'], ('H2', 'H1'): ['0.1', '100'], ('L1', 'H1'): ['.2', '100']}
	"""
	thresholds = {}
	for pair, delta in map(lambda w: str.split(w, "="), thresholdstrings):
		try:
			A, B = map(str.strip, pair.split(","))
		except Exception:
			raise ValueError, "cannot parse instruments %s" % pair
		thresholds[(A, B)] = map(str.strip, delta.split(","))
	for (A, B), value in thresholds.items():
		if (B, A) not in thresholds:
			thresholds[(B, A)] = value
	return thresholds


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
	def __init__(self, xmldoc, contributor_table_names):
		# get the coinc_def_id for coincidences involving the given
		# list of contributing tables
		self.coinc_def_id = llwapp.get_coinc_def_id(xmldoc, contributor_table_names)

		# find the coinc table or create one if not found
		try:
			self.coinctable = table.get_table(xmldoc, lsctables.CoincTable.tableName)
		except ValueError:
			self.coinctable = lsctables.New(lsctables.CoincTable)
			xmldoc.childNodes[0].appendChild(self.coinctable)
		self.coinctable.sync_ids()

		# find the coinc_map table or create one if not found
		try:
			self.coincmaptable = table.get_table(xmldoc, lsctables.CoincMapTable.tableName)
		except ValueError:
			self.coincmaptable = lsctables.New(lsctables.CoincMapTable)
			xmldoc.childNodes[0].appendChild(self.coincmaptable)

		# find the time_slide table
		self.time_slide_table = table.get_table(xmldoc, lsctables.TimeSlideTable.tableName)


	def time_slide_ids(self):
		"""
		Return a list of the time slide IDs.  The list is sorted in
		increasing order by ID number.
		"""
		ids = list(set([row.time_slide_id for row in self.time_slide_table]))
		ids.sort(lambda a, b: cmp(ilwd.ILWDID(a), ilwd.ILWDID(b)))
		return ids


	def get_time_slide(self, id):
		"""
		Return the time slide with the given ID as a dictionary of
		instrument/offset pairs.
		"""
		return self.time_slide_table.get_offset_dict(id)


	def append_coinc(self, process_id, time_slide_id, events):
		"""
		Takes a process ID, a time slide ID, and a list of events,
		and adds the events as a new coincidence to the coinc_event
		and coinc_map tables.
		"""
		coinc = lsctables.Coinc()
		coinc.process_id = process_id
		coinc.coinc_def_id = self.coinc_def_id
		coinc.coinc_event_id = self.coinctable.ids.next()
		coinc.time_slide_id = time_slide_id
		coinc.nevents = len(events)
		coinc.likelihood = 1.0
		self.coinctable.append(coinc)
		for event in events:
			coincmap = lsctables.CoincMap()
			coincmap.coinc_event_id = coinc.coinc_event_id
			coincmap.table_name = ilwd.ILWDTableName(event.event_id)
			coincmap.event_id = event.event_id
			self.coincmaptable.append(coincmap)


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
	interval.  Each search must provide a subclass of this class, and
	set the EventListType attribute of the EventListDict class to their
	own private subclass.  The only methods that *must* be overridden
	in a subclass are the _add_offset() and get_coincs() methods.  The
	make_index() method can be overridden if needed.  None of the other
	methods inherited from the list parent class need to be overridden,
	indeed they probably should not be unless you know what you're
	doing.
	"""
	def __init__(self):
		self.offset = LIGOTimeGPS(0)

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
		if offset != self.offset:
			self._add_offset(offset - self.offset)
			self.offset = offset

	def get_coincs(self, event_a, threshold, comparefunc):
		"""
		Return a list of the events from this list that are
		coincident with event_a.  The events must be coincident
		with event_a, not merely be likely to be coincident with
		event_a given more careful scrutiny, because the events
		returned by this method will never again been compared to
		event_a.  However, it is not necessary for the events in
		the list returned to be themselves mutually coincident in
		any way (that might not even make sense, since each list
		contains only events from a single instrument).

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
	whose events should be included.  The EventListType attribute must
	be set to a subclass of the EventList class.
	"""
	EventListType = None

	def __new__(self, *args):
		return dict.__new__(self)

	def __init__(self, event_table, process_ids):
		"""
		Initialize a newly-created instance.  event_table is a list
		of events (e.g., an instance of a glue.ligolw.table.Table
		subclass), and process_ids is a list or set of the
		process_ids whose events should be considered in the
		coincidence analysis.
		"""
		for event in event_table:
			if event.process_id in process_ids:
				if event.ifo not in self:
					self[event.ifo] = self.EventListType()
				self[event.ifo].append(event)
		for l in self.itervalues():
			l.make_index()

	def set_offsetdict(self, offsetdict):
		"""
		Set the event list offsets to those in the dictionary of
		instrument/offset pairs.  Instruments not in offsetdict are
		not modified.  KeyError is raised if the dictionary of
		instrument/offset pairs contains a key (instrument) that
		this dictionary does not.
		"""
		for instrument, offset in offsetdict.iteritems():
			self[instrument].set_offset(offset)

	def remove_offsetdict(self):
		"""
		Remove the offsets from all event lists (reset them to 0).
		"""
		for l in self.itervalues():
			l.set_offset(LIGOTimeGPS(0))


def make_eventlists(xmldoc, event_table_name, max_segment_gap, program_name):
	"""
	Convenience wrapper for constructing a dictionary of event lists
	from an XML document tree, the name of a table from which to get
	the events, a maximum allowed time window, and the name of the
	program that generated the events.
	"""
	return EventListDict(table.get_table(xmldoc, event_table_name), coincident_process_ids(xmldoc, max_segment_gap, program_name))


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

	# extract a segmentlistdict;  protract by half the largest
	# coincidence window so as to not miss edge effects
	search_summ_table = table.get_table(xmldoc, lsctables.SearchSummaryTable.tableName)
	seglistdict = search_summ_table.get_out_segmentlistdict(proc_ids).protract(max_segment_gap / 2)

	# determine which time slides are possible given the instruments in
	# the search summary table
	time_slide_table = table.get_table(xmldoc, lsctables.TimeSlideTable.tableName)
	timeslides = time_slide_table.get_offsets().values()
	for i in xrange(len(timeslides) - 1, -1, -1):
		for instrument in timeslides[i].iterkeys():
			if instrument not in seglistdict:
				del timeslides[i]
				break

	# determine the coincident segments for each instrument
	seglistdict = llwapp.get_coincident_segmentlistdict(seglistdict, timeslides)

	# find the IDs of the processes which contributed to the coincident
	# segments
	coinc_proc_ids = set()
	for row in search_summ_table:
		if row.process_id not in proc_ids or row.process_id in coinc_proc_ids:
			continue
		if seglistdict[row.ifos].intersects_segment(row.get_out()):
			coinc_proc_ids.add(row.process_id)
	return coinc_proc_ids


#
# =============================================================================
#
#                            Coincidence Iterators
#
# =============================================================================
#


def Level1Iterator(eventlists, comparefunc, instruments, thresholds):
	"""
	First-pass coincidence generator.  Generates a sequence of tuples
	whose elements are, in order:

	tick:  a progress indicator whose value is in the range [0, ticks)
	ticks:  the upper bound for tick
	event:  a burst event
	ntuples:  see below

	ntuples is (yet) another generator that produces a sequence of
	lists of burst events.  Each list of burst events, when event is
	added to it, constitutes a potential coincidence with exactly one
	event from each instrument.  Each event in the list is guaranteed
	to be coincident with event, but the mutual coincidence of the
	events in the list has not yet been established.
	"""
	instruments = list(instruments)	# so we can safely modify it
	eventlists = map(eventlists.__getitem__, instruments)
	lengths = map(len, eventlists)
	length = min(lengths)
	shortest = lengths.index(length)
	shortestlist = eventlists.pop(shortest)
	shortestinst = instruments.pop(shortest)
	thresholds = map(lambda inst: thresholds[(shortestinst, inst)], instruments)
	for n, event in enumerate(shortestlist):
		yield n, length, event, itertools.MultiIter(map(lambda eventlist, threshold: eventlist.get_coincs(event, threshold, comparefunc), eventlists, thresholds))


def mutually_coincident(events, comparefunc, thresholds):
	"""
	Return True if the all the events in the list are mutually
	coincident.
	"""
	try:
		for a, b in itertools.choices(events, 2):
			if comparefunc(a, b, thresholds[(a.ifo, b.ifo)]):
				return False
	except KeyError, e:
		raise KeyError, "no coincidence thresholds provided for instrument pair %s" % str(e)
	return True


def CoincidentNTuples(eventlists, comparefunc, instruments, thresholds, verbose = False):
	"""
	Given an EventListDict object, a list (or iterator) of instruments,
	and a dictionary of instrument pair thresholds, generate a sequence
	of lists of mutually coincident events.  Each list of mutually
	coincident events yielded by this generator will contain exactly
	one event from each of the instruments in the instrument list.
	"""
	for tick, ticks, event, ntuples in Level1Iterator(eventlists, comparefunc, instruments, thresholds):
		if verbose and not (tick % 500):
			print >>sys.stderr, "\t%.1f%%\r" % (100.0 * tick / ticks),
		for ntuple in ntuples:
			if mutually_coincident(ntuple, comparefunc, thresholds):
				ntuple.append(event)
				yield ntuple
	if verbose:
		print >>sys.stderr, "\t100.0%"
