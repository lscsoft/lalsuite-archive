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

from glue import segments
from glue.ligolw import table
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


def sngl_burst_get_peak(self):
	# get_peak() override to use pylal.date.LIGOTimeGPS instead of
	# glue.lal.LIGOTimeGPS
	return LIGOTimeGPS(self.peak_time, self.peak_time_ns)


def sngl_burst_get_period(self):
	# get_period() override to use pylal.date.LIGOTimeGPS instead of
	# glue.lal.LIGOTimeGPS
	start = LIGOTimeGPS(self.start_time, self.start_time_ns)
	return segments.segment(start, start + self.duration)


def sngl_burst___cmp__(self, other):
	# compare self's peak time to the LIGOTimeGPS instance other
	return cmp(self.peak_time, other.seconds) or cmp(self.peak_time_ns, other.nanoseconds)


lsctables.SnglBurst.get_start = sngl_burst_get_start
lsctables.SnglBurst.get_peak = sngl_burst_get_peak
lsctables.SnglBurst.get_period = sngl_burst_get_period
lsctables.SnglBurst.__cmp__ = sngl_burst___cmp__


#
# =============================================================================
#
#                           Add Process Information
#
# =============================================================================
#


process_program_name = "ligolw_burca"


def append_process(xmldoc, **kwargs):
	process = llwapp.append_process(xmldoc, program = process_program_name, version = __version__, cvs_repository = "lscsoft", cvs_entry_time = __date__, comment = kwargs["comment"])

	params = [
		("--program", "lstring", kwargs["program"]),
		("--coincidence-algorithm", "lstring", kwargs["coincidence_algorithm"])
	]
	if "stringcusp_params" in kwargs:
		params += [("--stringcusp-params", "lstring", kwargs["stringcusp_params"])]
	if "force" in kwargs:
		params += [("--force", "lstring", "")]
	# FIXME: this loop doesn't actually reproduce the original command
	# line arguments correctly.  The information stored in the table is
	# correct and complete, it's just that you have to understand how
	# burca works, internally, to know how to turn the table's contents
	# back into a valid command line.
	for key, value in kwargs["thresholds"].iteritems():
		if key[0] < key[1]:
			params += [("--thresholds", "lstring", "%s,%s=%s" % (key[0], key[1], ",".join(map(str, value))))]
	llwapp.append_process_params(xmldoc, process, params)

	return process


def dbget_thresholds(connection):
	"""
	Extract the --thresholds arguments that had been given to the
	ligolw_burca job recorded in the process table of database to which
	connection points.
	"""
	thresholds = snglcoinc.parse_thresholds(map(str, llwapp.dbget_process_params(connection, process_program_name, "--thresholds")))
	for pair, (dt, df, dhrss) in thresholds.items():
		thresholds[pair] = (float(dt), float(df), float(dhrss))
	return thresholds


#
# =============================================================================
#
#                          CoincTables Customizations
#
# =============================================================================
#


#
# For use with excess power.
#


def make_multi_burst(process_id, coinc_event_id, events):
	multiburst = lsctables.MultiBurst()
	multiburst.process_id = process_id
	multiburst.coinc_event_id = coinc_event_id

	# snr = sum of snrs
	multiburst.snr = sum(event.snr for event in events)

	# duration = snr-weighted average of durations
	multiburst.duration = sum(event.snr * event.duration for event in events) / multiburst.snr

	# central_freq = snr-weighted average of peak frequencies
	multiburst.central_freq = sum(event.snr * event.peak_frequency for event in events) / multiburst.snr

	# bandwidth = snr-weighted average of bandwidths
	multiburst.bandwidth = sum(event.snr * event.bandwidth for event in events) / multiburst.snr

	# confidence = arithmetic mean of confidences
	#multiburst.confidence = sum(event.confidence for event in events) / len(events)
	# confidence = minimum of confidences
	multiburst.confidence = min(event.confidence for event in events)

	# done
	return multiburst


class ExcessPowerCoincTables(snglcoinc.CoincTables):
	def __init__(self, xmldoc, contributor_table_names):
		snglcoinc.CoincTables.__init__(self, xmldoc, contributor_table_names)

		# find the multi_burst table or create one if not found
		try:
			self.multibursttable = table.get_table(xmldoc, lsctables.MultiBurstTable.tableName)
		except ValueError:
			self.multibursttable = lsctables.New(lsctables.MultiBurstTable, ("process_id", "duration", "central_freq", "bandwidth", "snr", "confidence", "coinc_event_id"))
			xmldoc.childNodes[0].appendChild(self.multibursttable)

	def append_coinc(self, process_id, time_slide_id, events):
		coinc = snglcoinc.CoincTables.append_coinc(self, process_id, time_slide_id, events)
		self.multibursttable.append(make_multi_burst(process_id, coinc.coinc_event_id, events))
		return coinc


#
# =============================================================================
#
#                            Event List Management
#
# =============================================================================
#


#
# For use with excess power coincidence test
#


def ExcessPowerMaxSegmentGap(xmldoc, thresholds):
	"""
	Determine the maximum allowed segment gap for use with the excess
	power coincidence test.
	"""
	# longest duration among events
	max_duration = max(table.get_table(xmldoc, lsctables.SnglBurstTable.tableName).getColumnByName("duration"))

	# loosest relative threshold
	max_dt = max([t[0] for t in thresholds.itervalues()])

	# convert to absolute threshold
	max_dt = max_dt * max_duration

	# double for safety
	return 2 * max_dt


class ExcessPowerEventList(snglcoinc.EventList):
	"""
	A customization of the EventList class for use with the excess
	power search.
	"""
	def make_index(self):
		"""
		Sort events by peak time so that a bisection search can
		retrieve them.  Note that the bisection search relies on
		the __cmp__() method of the SnglBurst row class having
		previously been set to compare the event's peak time to a
		LIGOTimeGPS.
		"""
		# sort by peak time
		self.sort(lambda a, b: cmp(a.peak_time, b.peak_time) or cmp(a.peak_time_ns, b.peak_time_ns))
		self.max_duration = LIGOTimeGPS(max([event.duration for event in self]))

	def _add_offset(self, delta):
		"""
		Add an amount to the peak time of each event.
		"""
		for event in self:
			event.set_peak(event.get_peak() + delta)

	def get_coincs(self, event_a, threshold, comparefunc):
		min_peak = event_a.get_peak() - threshold[0] * self.max_duration
		max_peak = event_a.get_peak() + threshold[0] * self.max_duration
		return [event_b for event_b in self[bisect.bisect_left(self, min_peak) : bisect.bisect_right(self, max_peak)] if not comparefunc(event_a, event_b, threshold)]


#
# For use with string coincidence test
#


def StringMaxSegmentGap(xmldoc, thresholds):
	"""
	Determine the maximum allowed segment gap for use with the string
	coincidence test.
	"""
	# loosest threshold
	max_dt = max([t[0] for t in thresholds.itervalues()])

	# double for safety
	return 2 * max_dt


class StringEventList(snglcoinc.EventList):
	"""
	A customization of the EventList class for use with the string
	search.
	"""
	def make_index(self):
		"""
		Sort events by peak time so that a bisection search can
		retrieve them.  Note that the bisection search relies on
		the __cmp__() method of the SnglBurst row class having
		previously been set to compare the event's peak time to a
		LIGOTimeGPS.
		"""
		self.sort(lambda a, b: cmp(a.peak_time, b.peak_time) or cmp(a.peak_time_ns, b.peak_time_ns))

	def _add_offset(self, delta):
		"""
		Add an amount to the peak time of each event.
		"""
		for event in self:
			event.set_peak(event.get_peak() + delta)

	def get_coincs(self, event_a, threshold, comparefunc):
		min_peak = event_a.get_peak() - threshold[0]
		max_peak = event_a.get_peak() + threshold[0]
		return [event_b for event_b in self[bisect.bisect_left(self, min_peak) : bisect.bisect_right(self, max_peak)] if not comparefunc(event_a, event_b, threshold)]


#
# =============================================================================
#
#                              Coincidence Tests
#
# =============================================================================
#


def ExcessPowerCoincCompare(a, b, thresholds):
	"""
	There are three thresholds, dt, df, and dhrss, setting,
	respectively, the maximum allowed fractional difference in peak
	times, the maximum allowed fractional difference in peak
	frequencies, and the maximum allowed fractional difference in
	h_{rss}.

	In the case of the peak times, the difference is taken as a
	fraction of the average of the durations of the two events' most
	significant contributing tile.  For example, dt = 0 means the peak
	times must be exactly equal, while dt = 1 is roughly equivalent to
	requiring the events' most significant contributing tile's time
	intervals to intersect, and dt -> \infty is equivalent to no
	constraint.

	For the peak frequencies, the difference is taken as a fraction of
	the average of the bandwidths of two events' most significant
	contributing tile.  For example, df = 0 means the peak frequencies
	must be exactly equal, while df = 1 is roughly equivalent to
	requiring the events' most significant contributing tile's
	frequency bands to intersect, and df -> \infty is equivalent to no
	constraint (although since the bandwidth of the search is finite,
	and there is a minimum bandwidth for a candidate event, there
	exists a finite df threshold which admits all events)

	For h_{rss}, the difference is taken as a fraction of the average
	of the h_{rss} of the two events' most significant contributing
	tiles.  So dhrss = 0 means the two events' h_{rss}s must be exactly
	equal, while dhrss = 2 is equivalent to no constraint on h_{rss}
	(h_{rss} cannot be negative, so no two values can differ by more
	than twice their average).

	Returns False (a & b are coincident) if the two events match within
	the tresholds.  Retruns non-zero otherwise.
	"""
	# unpack thresholds
	dt, df, dhrss = thresholds

	# convert fractional deltas to absolute deltas
	dt = dt * (a.ms_duration + b.ms_duration) / 2
	df = df * (a.ms_bandwidth + b.ms_bandwidth) / 2
	dhrss = dhrss * (a.ms_hrss + b.ms_hrss) / 2

	# test for coincidence
	coincident = abs(float(a.get_peak() - b.get_peak())) <= dt and abs(a.peak_frequency - b.peak_frequency) <= df and abs(a.ms_hrss - b.ms_hrss) <= dhrss

	# return result
	return not coincident


def StringCoincCompare(a, b, thresholds):
	"""
	Returns False (a & b are coincident) if their peak times agree
	within dt, and in the case of H1+H2 pairs if their amplitudes agree
	according to some kinda test.
	"""
	# unpack thresholds
	dt, kappa, epsilon = thresholds

	# test for time coincidence
	coincident = abs(float(a.get_peak() - b.get_peak())) <= dt

	# for H1+H2, also test for amplitude coincidence
	if a.ifo in ("H1", "H2") and b.ifo in ("H1", "H2"):
		adelta = abs(a.amplitude) * (kappa / a.snr + epsilon)
		bdelta = abs(b.amplitude) * (kappa / b.snr + epsilon)
		coincident = coincident and a.amplitude - adelta <= b.amplitude <= a.amplitude + adelta and b.amplitude - bdelta <= a.amplitude <= b.amplitude + bdelta

	# return result
	return not coincident


#
# =============================================================================
#
#                                 Library API
#
# =============================================================================
#


def ligolw_burca(xmldoc, CoincTables, comparefunc, **kwargs):
	# add an entry in the process table
	process = append_process(xmldoc, **kwargs)

	if kwargs["verbose"]:
		print >>sys.stderr, "indexing ..."

	# prepare the coincidence table interface
	coinc_tables = CoincTables(xmldoc, [lsctables.SnglBurstTable.tableName])

	# cast offsets in time_slide table to LIGOTimeGPS (avoids repeated
	# conversion when applying the time slide to each trigger)
	for row in coinc_tables.time_slide_table:
		row.offset = LIGOTimeGPS(row.offset)

	# build the event list accessors, populated with events from those
	# processes that can participate in a coincidence
	eventlists = snglcoinc.make_eventlists(xmldoc, lsctables.SnglBurstTable.tableName, kwargs["get_max_segment_gap"](xmldoc, kwargs["thresholds"]), kwargs["program"])

	# iterate over time slides
	time_slide_ids = coinc_tables.time_slide_ids()
	for n, time_slide_id in enumerate(time_slide_ids):
		offsetdict = coinc_tables.get_time_slide(time_slide_id)
		if kwargs["verbose"]:
			print >>sys.stderr, "time slide %d/%d: %s" % (n + 1, len(time_slide_ids), ",".join([" %s = %+.16g s" % (i, float(o)) for i, o in offsetdict.iteritems()]))
		if len(offsetdict.keys()) < 2:
			if kwargs["verbose"]:
				print >>sys.stderr, "\tsingle-instrument time slide: skipped"
			continue
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
		for ntuple in snglcoinc.CoincidentNTuples(eventlists, comparefunc, offsetdict.iterkeys(), kwargs["thresholds"], kwargs["verbose"]):
			coinc_tables.append_coinc(process.process_id, time_slide_id, ntuple)

	# clean up and finish
	eventlists.remove_offsetdict()
	llwapp.set_process_end_time(process)
	return xmldoc
