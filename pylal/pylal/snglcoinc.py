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

import bisect
import sys

from glue.ligolw import table
from glue.ligolw import lsctables
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
	Turn a list of strings of the form "inst1,inst2=delta" into a
	dictionary with (inst1, inst2) 2-tuples as keys and the deltas as
	the values as strings.  In the input, each pair of instruments is
	allowed two entries, once for each order.  The output will contain
	exactly two entries for each pair of instruments;  if the input
	strings do not contain a set for a particular instrument order, the
	value for the missing pair will be copied from the one provided.
	"""
	thresholds = {}
	for [pair, delta] in map(lambda w: str.split(w, "="), thresholdstrings):
		try:
			[A, B] = pair.split(",")
		except ValueError:
			raise ValueError, "incorrect number of instruments"
		thresholds[(A, B)] = delta
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

		# initialize the coinc_event_id iterator
		self.coincids = lsctables.NewILWDs(self.coinctable, "coinc_event_id")

		# find the coinc_map table or create one if not found
		try:
			self.coincmaptable = table.get_table(xmldoc, lsctables.CoincMapTable.tableName)
		except ValueError:
			self.coincmaptable = lsctables.New(lsctables.CoincMapTable)
			xmldoc.childNodes[0].appendChild(self.coincmaptable)

		# find the time_slide table, and cast all offsets to
		# LIGOTimeGPS.
		self.tisitable = table.get_table(xmldoc, lsctables.TimeSlideTable.tableName)
		for row in self.tisitable:
			row.offset = LIGOTimeGPS(row.offset)


	def time_slide_ids(self):
		"""
		Return a list of the time slide IDs.
		"""
		return self.tisitable.dict.keys()


	def get_time_slide(self, id):
		"""
		Return the time slide with the given ID as a dictionary of
		instrument/offset pairs.
		"""
		return self.tisitable.get_offset_dict(id)


	def append_coinc(self, process_id, time_slide_id, events):
		"""
		Takes a process ID, a time slide ID, and a list of events,
		and adds the events as a new coincidence to the coinc_event
		and coinc_map tables.
		"""
		coinc = lsctables.Coinc()
		coinc.process_id = process_id
		coinc.coinc_def_id = self.coinc_def_id
		coinc.coinc_event_id = self.coincids.next()
		coinc.time_slide_id = time_slide_id
		coinc.nevents = len(events)
		self.coinctable.append(coinc)
		for event in events:
			coincmap = lsctables.CoincMap()
			coincmap.coinc_event_id = coinc.coinc_event_id
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
		for event in event_table:
			if llwapp.bisect_contains(process_ids, event.process_id):
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


def make_eventlists(xmldoc, event_table_name, max_delta_t, program):
	"""
	Convenience wrapper for constructing a dictionary of event lists
	from an XML document tree, the name of a table from which to get
	the events, a maximum allowed time window, and the name of the
	program that generated the events.
	"""
	return EventListDict(table.get_table(xmldoc, event_table_name), coincident_process_ids(xmldoc, max_delta_t, program))


#
# =============================================================================
#
#                                Process Filter
#
# =============================================================================
#

def coincident_process_ids(xmldoc, max_delta_t, program):
	"""
	Take an XML document tree and determine the list of process IDs
	that will participate in coincidences identified by the time slide
	table therein.  It is OK for xmldoc to contain time slides
	involving instruments not represented in the list of processes,
	these time slides are ignored.  max_delta_t is the largest time
	window that will be considered in the coincidence tests;  that is
	after applying a time slide, two segments can have a gap this large
	between them and still be considered coincident.
	"""
	# extract a segmentlistdict;  protract by half the largest
	# coincidence window so as to not miss edge effects
	seglistdict = llwapp.segmentlistdict_fromsearchsummary(xmldoc, program).protract(max_delta_t / 2)

	# determine which time slides are possible given the instruments in
	# the search summary table
	tisitable = table.get_table(xmldoc, lsctables.TimeSlideTable.tableName)
	timeslides = map(tisitable.get_offset_dict, tisitable.dict.keys())
	i = 0
	while i < len(timeslides):
		for instrument in timeslides[i].keys():
			if instrument not in seglistdict:
				del timeslides[i]
				break
		else:
			i += 1

	# determine the coincident segments for each instrument
	seglistdict = llwapp.get_coincident_segmentlistdict(seglistdict, timeslides)

	# get the list of all process IDs for the given program
	proc_ids = llwapp.get_process_ids_by_program(xmldoc, program)

	# find the IDs of the processes which contributed to the coincident
	# segments
	coinc_proc_ids = []
	for row in table.get_table(xmldoc, lsctables.SearchSummaryTable.tableName):
		if (not llwapp.bisect_contains(proc_ids, row.process_id)) or llwapp.bisect_contains(coinc_proc_ids, row.process_id):
			continue
		if seglistdict[row.ifos].intersects_segment(row.get_out()):
			bisect.insort_left(coinc_proc_ids, row.process_id)
	return coinc_proc_ids


#
# =============================================================================
#
#                            Coincidence Iterators
#
# =============================================================================
#

CompareFunc = None


def Level1Iterator(eventlists, instruments, thresholds):
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
	thresholds = map(lambda inst: thresholds[(inst, shortestinst)], instruments)
	for n, event in enumerate(shortestlist):
		yield n, length, event, itertools.MultiIter(map(lambda eventlist, threshold: eventlist.get_coincs(event, threshold, CompareFunc), eventlists, thresholds))


def mutually_coincident(events, thresholds):
	"""
	Return True if the all the events in the list are mutually
	coincident.
	"""
	try:
		for [a, b] in itertools.choices(events, 2):
			if CompareFunc(a, b, thresholds[(a.ifo, b.ifo)]):
				return False
	except KeyError, e:
		raise KeyError, "no coincidence window provided for instrument pair %s" % str(e)
	return True


def CoincidentNTuples(eventlists, instruments, thresholds, verbose = False):
	"""
	Given an EventListDict object, a list (or iterator) of instruments,
	and a dictionary of instrument pair thresholds, generate a sequence
	of lists of mutually coincident events.  Each list of mutually
	coincident events yielded by this generator will contain exactly
	one event from each of the instruments in the instrument list.
	"""
	for tick, ticks, event, ntuples in Level1Iterator(eventlists, instruments, thresholds):
		if verbose and not (tick % (ticks / 200 or 1)):
			print >>sys.stderr, "\t%.1f%%\r" % (100.0 * tick / ticks),
		for ntuple in ntuples:
			if mutually_coincident(ntuple, thresholds):
				ntuple.append(event)
				yield ntuple
	if verbose:
		print >>sys.stderr, "\t100.0%"
