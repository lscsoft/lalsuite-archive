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


from glue import iterutils
from glue.ligolw import table
from glue.ligolw import lsctables
from pylal import llwapp
from pylal import ligolw_tisi


__author__ = "Kipp Cannon <kipp@gravity.phys.uwm.edu>"
__version__ = "$Revision$"[11:-2]
__date__ = "$Date$"[7:-2]


#
# =============================================================================
#
#                               Time Slide Graph
#
# =============================================================================
#


class TimeSlideGraphNode(object):
	def __init__(self, offset_vector, time_slide_id = None):
		self.time_slide_id = time_slide_id
		self.offset_vector = offset_vector
		self.deltas = frozenset(ligolw_tisi.offset_vector_to_deltas(offset_vector).items())
		self.components = None
		self.coincs = None
		self.unused_coincs = set()


	def get_coincs(self, verbose = False):
		#
		# has this node already been visited?  if so, return the
		# answer we already know
		#

		if self.coincs is not None:
			if verbose:
				print >>sys.stderr, "\treusing %s" % (", ".join(("%s = %+.16g s" % x) for x in sorted(self.offset_vector.items())))
			return self.coincs

		#
		# is this a leaf node?  construct the coincs explicitly
		#

		if self.components is None:
			# FIXME:  move the code from ligolw_thinca to here
			raise NotImplementedError

		#
		# is this a head node, or some other node that magically
		# has only one component?  copy coincs from component
		#

		if len(self.components) == 1:
			if verbose:
				print >>sys.stderr, "\tcopying from %s" % (", ".join(("%s = %+.16g s" % x) for x in sorted(self.components[0].offset_vector.items())))
			self.coincs = self.components[0].get_coincs(verbose = verbose)
			self.unused_coincs = self.components[0].unused_coincs

			#
			# done.  unlink the graph as we go to release
			# memory
			#

			self.components = None
			return self.coincs

		#
		# len(self.components) == 2 is impossible
		#

		assert len(self.components) > 2

		#
		# this is a regular node in the graph.  use coincidence
		# synthesis algorithm to populate its coincs
		#

		if verbose:
			print >>sys.stderr, "\tassembling %s ..." % (", ".join(("%s = %+.16g s" % x) for x in sorted(self.offset_vector.items())))

		self.unused_coincs = reduce(lambda a, b: a & b, (component.unused_coincs for component in self.components)) | reduce(lambda a, b: a | b, (set(component.coincs) for component in self.components))

		self.coincs = []
		allcoincs0 = self.components[0].get_coincs(verbose = verbose)
		allcoincs1 = self.components[1].get_coincs(verbose = verbose)
		allcoincs2 = self.components[-1].get_coincs(verbose = verbose)
		length = len(allcoincs0)
		for n, coinc0 in enumerate(allcoincs0):
			if verbose and not (n % 200):
				print >>sys.stderr, "\t%.1f%%\r" % (100.0 * n / length),
			coincs1 = allcoincs1[bisect.bisect_left(allcoincs1, coinc0[:-1]):bisect.bisect_left(allcoincs1, coinc0[:-2] + (coinc0[-2] + 1,))]
			coincs2 = allcoincs2[bisect.bisect_left(allcoincs2, coinc0[1:]):bisect.bisect_left(allcoincs2, coinc0[1:-1] + (coinc0[-1] + 1,))]
			for coinc1 in coincs1:
				i = bisect.bisect_left(coincs2, coinc0[1:] + coinc1[-1:])
				if i < len(coincs2) and coincs2[i] == coinc0[1:] + coinc1[-1:]:
					new_coinc = coinc0[:1] + coincs2[i]
					self.unused_coincs -= set(iterutils.choices(new_coinc, len(new_coinc) - 1))
					self.coincs.append(new_coinc)
		if verbose:
			print >>sys.stderr, "\t100.0%"
		self.coincs.sort()
		self.coincs = tuple(self.coincs)

		#
		# done.  unlink the graph as we go to release memory
		#

		self.components = None
		return self.coincs


class TimeSlideGraph(object):
	def __init__(self, offset_vector_dict):
		#
		# populate the graph head nodes.  these represent the
		# target offset vectors requested by the calling code.
		#

		self.head = tuple(TimeSlideGraphNode(offset_vector, id) for id, offset_vector in sorted(offset_vector_dict.items()))

		#
		# populate the graph generations.  generations[n] is a
		# tuple of the nodes in the graph representing all unique
		# n-instrument offset vectors to be constructed as part of
		# the analysis (including normalized forms of any
		# n-instrument target offset vectors).
		#

		self.generations = {}
		n = max(len(offset_vector) for offset_vector in offset_vector_dict.values())
		self.generations[n] = tuple(TimeSlideGraphNode(offset_vector) for offset_vector in ligolw_tisi.time_slide_component_vectors((node.offset_vector for node in self.head), n))
		for n in range(n, 2, -1):	# [n, ..., 3]
			offset_vectors = [node.offset_vector for node in self.head if len(node.offset_vector) == n]
			if n in self.generations:
				offset_vectors += [node.offset_vector for node in self.generations[n]]
			self.generations[n - 1] = tuple(TimeSlideGraphNode(offset_vector) for offset_vector in ligolw_tisi.time_slide_component_vectors(offset_vectors, n - 1))

		#
		# link each n-instrument node to the n-1 instrument nodes
		# from which it will be constructed.  NOTE:  the components
		# are sorted according to the alphabetically-sorted tuples
		# of instrument names involved in each component;  this is
		# a critical part of the coincidence synthesis algorithm
		#

		for node in self.head:
			node.components = tuple(sorted((component for component in self.generations[len(node.offset_vector)] if node.deltas == component.deltas), lambda a, b: cmp(sorted(a.offset_vector), sorted(b.offset_vector))))
		for n, nodes in self.generations.items():
			if n <= 2:
				continue
			for node in nodes:
				component_deltas = set(frozenset(ligolw_tisi.offset_vector_to_deltas(offset_vector).items()) for offset_vector in ligolw_tisi.time_slide_component_vectors([node.offset_vector], n - 1))
				node.components = tuple(sorted((component for component in self.generations[n - 1] if component.deltas in component_deltas), lambda a, b: cmp(sorted(a.offset_vector), sorted(b.offset_vector))))

		#
		# done
		#


	def write(self, fileobj):
		"""
		Write a DOT graph representation of the time slide graph
		object to fileobj.
		"""
		vectorstring = lambda offset_vector: ",".join("%s=%g" % (instrument, offset) for instrument, offset in sorted(offset_vector.items()))

		print >>fileobj, "digraph \"Time Slides\" {"
		for nodes in self.generations.values():
			for node in nodes:
				node_name = vectorstring(node.offset_vector)
				print >>fileobj, "\t\"%s\" [shape=box];" % node_name
				for component in node.components:
					component_name = vectorstring(component.offset_vector)
					print >>fileobj, "\t\"%s\" -> \"%s\";" % (component_name, node_name)
		for node in self.head:
			node_name = vectorstring(node.offset_vector)
			for component in node.components:
				component_name = vectorstring(component.offset_vector)
				print >>fileobj, "\t\"%s\" -> \"%s\";" % (component_name, node_name)
		print >>fileobj, "}"


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

	def get_time_slides(self):
		"""
		Return the time_slide_id --> offset vector dictionary.
		"""
		return self.time_slide_table.as_dict()

	def append_coinc(self, process_id, time_slide_id, coinc_def_id, events):
		"""
		Takes a process ID, a time slide ID, and a list of events,
		and adds the events as a new coincidence to the coinc_event
		and coinc_map tables.
		"""
		# so we can iterate over it more than once incase we've
		# been given a generator expression.
		events = tuple(events)

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
	if len([segs for segs in seglistdict.values() if len(segs) == 1]):
		return proc_ids
	elif max(b[0] - a[1] for segs in seglistdict.values() for a, b in zip(segs[:-1], segs[1:])) <= max_segment_gap:
		return proc_ids

	# protract by half the largest coincidence window so as to not miss
	# edge effects
	seglistdict.protract(max_segment_gap / 2)

	# determine what time slides are possible given the instruments in
	# the search summary table
	avail_instruments = set(seglistdict)
	timeslides = [offset_vector for offset_vector in table.get_table(xmldoc, lsctables.TimeSlideTable.tableName).as_dict().values() if set(offset_vector).issubset(avail_instruments)]

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
		if verbose and not (n % 1000):
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
