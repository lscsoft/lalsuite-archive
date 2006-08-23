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

import sys

from glue import segments
from glue.ligolw import lsctables
from pylal import llwapp
from pylal import snglcoinc
from pylal.date import LIGOTimeGPS

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

def sngl_burst_get_start(self):
	# get_start() override to use pylal.date.LIGOTimeGPS instead of
	# glue.lal.LIGOTimeGPS
	return LIGOTimeGPS(self.start_time, self.start_time_ns)

def sngl_burst_get_period(self):
	# get_period() override to use pylal.date.LIGOTimeGPS instead of
	# glue.lal.LIGOTimeGPS
	start = LIGOTimeGPS(self.start_time, self.start_time_ns)
	return segments.segment(start, start + self.duration)


lsctables.SnglBurst.get_start = sngl_burst_get_start
lsctables.SnglBurst.get_period = sngl_burst_get_period


#
# =============================================================================
#
#                           Add Process Information
#
# =============================================================================
#

def append_process(xmldoc, **kwargs):
	process = llwapp.append_process(xmldoc, program = "ligolw_burca", version = __version__, cvs_repository = "lscsoft", cvs_entry_time = __date__, comment = kwargs["comment"])

	params = [("--program", "lstring", kwargs["program"])]
	for key, value in kwargs["window"].iteritems():
		if key[0] < key[1]:
			params += [("--window", "lstring", "%s,%s=%s" % (key[0], key[1], str(value)))]
	llwapp.append_process_params(xmldoc, process, params)

	return process


#
# =============================================================================
#
#                            Event List Management
#
# =============================================================================
#

class ExcessPowerEventList(snglcoinc.EventList):
	"""
	A customization of the EventList class for use with the excess
	power search.
	"""
	def make_index(self):
		"""
		Build a start time look up table and find the maximum of
		the event durations.  The start time look-up table is
		required because Python's bisection search cannot make use
		of an external comparison function.
		"""
		# sort in increasing order by start time
		self.sort(lambda a, b: cmp(a.start_time, b.start_time) or cmp(a.start_time_ns, b.start_time_ns))
		self.start_times = [event.get_start() for event in self]
		self.max_duration = LIGOTimeGPS(max([event.duration for event in self]))

	def __add_offset(self, delta):
		"""
		Add an amount to the start time of each event.
		"""
		for event in self:
			event.set_start(event.get_start() + delta)

	def get_coincs(self, event_a, threshold, comparefunc):
		min_start = event_a.get_start() - threshold - self.max_duration - self.offset
		max_start = event_a.get_start() + event_a.duration + threshold - self.offset
		return [event_b for event_b in self[bisect.bisect_left(self.start_times, min_start) : bisect.bisect_right(self.start_times, max_start)] if not comparefunc(event_a, event_b, threshold)]


class StringEventList(snglcoinc.EventList):
	"""
	A customization of the EventList class for use with the string
	search.
	"""
	def make_index(self):
		"""
		Build a peak time look up table.  The peak time look-up
		table is required because Python's bisection search cannot
		make use of an external comparison function.
		"""
		# sort in increasing order by peak time
		self.sort(lambda a, b: cmp(a.peak_time, b.peak_time) or cmp(a.peak_time_ns, b.peak_time_ns))
		self.peak_times = [event.get_peak() for event in self]

	def __add_offset(self, delta):
		"""
		Add an amount to the start time of each event.
		"""
		for event in self:
			event.set_peak(event.get_peak() + delta)

	def get_coincs(self, event_a, threshold, comparefunc):
		min_peak = event_a.get_peak() - threshold(0) - self.offset
		max_peak = event_a.get_peak() + threshold(0) - self.offset
		return [event_b for event_b in self[bisect.bisect_left(self.peak_times, min_peak) : bisect.bisect_right(self.peak_times, max_peak)] if not comparefunc(event_a, event_b, threshold)]



#
# =============================================================================
#
#                                 Library API
#
# =============================================================================
#

def ligolw_burca(xmldoc, **kwargs):
	# add an entry in the process table
	process = append_process(xmldoc, **kwargs)

	if kwargs["verbose"]:
		print >>sys.stderr, "indexing..."

	# prepare the coincidence table interface
	coinc_tables = snglcoinc.CoincTables(xmldoc, [lsctables.SnglBurstTable.tableName])

	# build the event list accessors, populated with events from those
	# processes that can participate in a coincidence
	# FIXME: correct determination of max_delta_t for string search
	eventlists = snglcoinc.make_eventlists(xmldoc, lsctables.SnglBurstTable.tableName, max(kwargs["window"].itervalues()), kwargs["program"])

	# iterate over time slides
	time_slide_ids = coinc_tables.time_slide_ids()
	for n, time_slide_id in enumerate(time_slide_ids):
		offsetdict = coinc_tables.get_time_slide(time_slide_id)
		if kwargs["verbose"]:
			print >>sys.stderr, "time slide %d/%d: %s" % (n + 1, len(time_slide_ids), str(offsetdict))
		if False in map(eventlists.__contains__, offsetdict.keys()):
			if kwargs["verbose"]:
				print >>sys.stderr, "\twarning: skipping due to insufficient data"
			continue
		if kwargs["verbose"]:
			print >>sys.stderr, "\tapplying time offsets ..."
		eventlists.set_offsetdict(offsetdict)
		if kwargs["verbose"]:
			print >>sys.stderr, "\tsearching ..."
		# search for and record coincidences
		for ntuple in snglcoinc.CoincidentNTuples(eventlists, offsetdict.keys(), kwargs["window"], kwargs["verbose"]):
			coinc_tables.append_coinc(process.process_id, time_slide_id, ntuple)

	# clean up and finish
	eventlists.remove_offsetdict()
	llwapp.set_process_end_time(process)
	return xmldoc
