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

import bisect
import sys

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

def parse_windows(windowstrings):
	"""
	Turn a list of strings of the form "inst1,inst2=delta" into a
	dictionary with (inst1, inst2) 2-tuples as keys and the deltas as
	the values.  Each pair of instruments gets two entries, once for
	each order.
	"""
	windows = {}
	for [pair, delay] in map(lambda w: str.split(w, "="), windowstrings):
		AB = tuple(pair.split(","))
		if len(AB) != 2:
			raise ValueError, "incorrect number of instruments"
		BA = (AB[1], AB[0])
		if (AB in windows) or (BA in windows):
			raise ValueError, "duplicate instrument pair"
		windows[AB] = windows[BA] = LIGOTimeGPS(delay)
	return windows


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
			self.coinctable = llwapp.get_table(xmldoc, lsctables.CoincTable.tableName)
		except ValueError:
			self.coinctable = lsctables.New(lsctables.CoincTable)
			xmldoc.childNodes[0].appendChild(self.coinctable)

		# initialize the coinc_event_id iterator
		self.coincids = lsctables.NewILWDs(self.coinctable, "coinc_event_id")

		# find the coinc_map table or create one if not found
		try:
			self.coincmaptable = llwapp.get_table(xmldoc, lsctables.CoincMapTable.tableName)
		except ValueError:
			self.coincmaptable = lsctables.New(lsctables.CoincMapTable)
			xmldoc.childNodes[0].appendChild(self.coincmaptable)

		# find the time_slide table, and cast all offsets to
		# LIGOTimeGPS.
		self.tisitable = llwapp.get_table(xmldoc, lsctables.TimeSlideTable.tableName)
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

def sngl_burst_cmp(self, other):
	"""
	For sorting sngl_burst events by start time.
	"""
	return cmp(self.start_time, other.start_time) or cmp(self.start_time_ns, other.start_time_ns)


def sngl_burst_get_start(self):
	# get_start() override to use pylal.date.LIGOTimeGPS instead of
	# glue.lal.LIGOTimeGPS
	return LIGOTimeGPS(self.start_time, self.start_time_ns)

lsctables.SnglBurst.get_start = sngl_burst_get_start


class EventList(list):
	"""
	A class for managing a list of events:  applying time offsets, and
	retrieving subsets of the list selected by time interval.
	"""
	def __init__(self):
		self.offset = LIGOTimeGPS(0)

	def make_index(self):
		"""
		Build look-up tables for the events.  Must be called after
		events have been added to the list, and before the
		get_coincs() method is called.  Must be called again if the
		list is modified.
		"""
		self.sort(sngl_burst_cmp)
		self.start_times = [event.get_start() for event in self]
		self.max_duration = LIGOTimeGPS(max([event.duration for event in self]))

	def set_offset(self, offset):
		"""
		Add an offset to the times of all events in the list.
		"""
		if offset != self.offset:
			delta = offset - self.offset
			for event in self:
				event.set_start(event.get_start() + delta)
			self.offset = offset

	def get_coincs(self, event_a, window, comparefunc):
		"""
		Return a list of the events coincident with event_a.
		"""
		min_start = event_a.get_start() - window - self.max_duration - self.offset
		max_start = event_a.get_start() + event_a.duration + window - self.offset
		return [event_b for event_b in self[bisect.bisect_left(self.start_times, min_start) : bisect.bisect_right(self.start_times, max_start)] if not comparefunc(event_a, event_b, window)]


class EventListDict(dict):
	"""
	A dictionary of EventList objects, indexed by instrument,
	initialized from an XML sngl_burst table and a list of process IDs
	whose events should be included.
	"""
	EventListType = EventList

	def __new__(self, *args):
		return dict.__new__(self)

	def __init__(self, event_table, process_ids):
		for event in event_table:
			if llwapp.bisect_contains(process_ids, event.process_id):
				try:
					self[event.ifo].append(event)
				except KeyError:
					self[event.ifo] = self.EventListType()
					self[event.ifo].append(event)
		for l in self.itervalues():
			l.make_index()

	def set_offsetdict(self, offsetdict):
		"""
		Set the events list offsets to those in the dictionary of
		instrument/offset pairs.  Instruments not in offsetdict are
		not modified.
		"""
		for instrument, offset in offsetdict.iteritems():
			self[instrument].set_offset(offset)

	def remove_offsetdict(self):
		"""
		Remove the offsets from all event lists.
		"""
		for l in self.itervalues():
			l.set_offset(LIGOTimeGPS(0))


def make_eventlists(xmldoc, event_table_name, windows, program):
	"""
	Convenience wrapper for constructing a dictionary of event lists
	from an XML document tree, the name of a table from which to get
	the events, a dictionary of inter-instrument time windows, and the
	name of the program that generated the events.
	"""
	return EventListDict(llwapp.get_table(xmldoc, event_table_name), coincident_process_ids(xmldoc, windows, program))


#
# =============================================================================
#
#                                Process Filter
#
# =============================================================================
#

def coincident_process_ids(xmldoc, windows, program):
	"""
	Take an XML document tree and determine the list of process IDs
	that will participate in coincidences identified by the time slide
	table therein.
	"""
	# find the largest coincidence window
	try:
		halfmaxwindow = max(windows.itervalues()) / 2
	except:
		halfmaxwindow = LIGOTimeGPS(0)

	# extract a segmentlistdict;  protract by half the largest
	# coincidence window so as to not miss edge effects
	seglistdict = llwapp.segmentlistdict_fromsearchsummary(xmldoc, program).protract(halfmaxwindow)

	# determine the coincident segments for each instrument
	tisitable = llwapp.get_table(xmldoc, lsctables.TimeSlideTable.tableName)
	seglistdict = llwapp.get_coincident_segmentlistdict(seglistdict, map(tisitable.get_offset_dict, tisitable.dict.keys()))

	# get the list of all process IDs for the given program
	proc_ids = llwapp.get_process_ids_by_program(xmldoc, program)

	# find the IDs of the processes which contributed to the coincident
	# segments
	coinc_proc_ids = []
	for row in llwapp.get_table(xmldoc, lsctables.SearchSummaryTable.tableName):
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


def Level1Iterator(eventlists, instruments, windows):
	"""
	First-pass coincidence generator.  Generates a sequence of tuples
	whose elements are, in order:

	n:  a progress indicator whose value is in the range [0, length)
	length:  the upper bound for n
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
	windows = map(lambda inst: windows[(inst, shortestinst)], instruments)
	for n, event in enumerate(shortestlist):
		yield n, length, event, itertools.MultiIter(map(lambda eventlist, window: eventlist.get_coincs(event, window, CompareFunc), eventlists, windows))


def mutually_coincident(events, windows):
	"""
	Return True if the all the events in the list are mutually
	coincident.
	"""
	try:
		for [a, b] in itertools.choices(events, 2):
			if CompareFunc(a, b, windows[(a.ifo, b.ifo)]):
				return False
	except KeyError, e:
		raise KeyError, "no coincidence window provided for instrument pair %s" % str(e)
	return True


def CoincidentNTuples(eventlists, instruments, windows, verbose = False):
	"""
	Given a EventListDict object, a list of instruments, and a
	dictionary of instrument pair time windows, generate a sequence of
	lists of mutually coincident events.  Each list has exactly one
	event from each instrument.
	"""
	for n, length, event, ntuples in Level1Iterator(eventlists, instruments, windows):
		if verbose and not (n % (length / 200 or 1)):
			print >>sys.stderr, "\t%.1f%%\r" % (100.0 * n / length),
		for ntuple in ntuples:
			if mutually_coincident(ntuple, windows):
				ntuple.append(event)
				yield ntuple
	if verbose:
		print >>sys.stderr, "\t100.0%"
