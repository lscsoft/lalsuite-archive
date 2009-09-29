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
Burst injection identification library.  Contains code providing the
capacity to search a list of sngl_burst burst candidates for events
matching entries in a sim_burst list of software injections, recording the
matches as burst <--> injection coincidences using the standard coincidence
infrastructure.  Also, any pre-recorded burst <--> burst coincidences are
checked for cases where all the burst events in a coincidence match an
injection, and these are recorded as coinc <--> injection coincidences,
again using the standard coincidence infrastructure.
"""


import bisect
import sys


from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw.utils import process as ligolw_process
from pylal import ligolw_burca
from pylal import llwapp
from pylal import SimBurstUtils
from pylal.xlal import tools
from pylal.xlal.datatypes.ligotimegps import LIGOTimeGPS


#
# Use a memory-efficient row class written in C for the coinc_event_map
# table
#


lsctables.CoincMapTable.RowType = lsctables.CoincMap = tools.CoincMap


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


lsctables.LIGOTimeGPS = LIGOTimeGPS


def sngl_burst___cmp__(self, other):
	# compare self's peak time to the LIGOTimeGPS instance other
	return cmp(self.peak_time, other.seconds) or cmp(self.peak_time_ns, other.nanoseconds)


lsctables.SnglBurst.__cmp__ = sngl_burst___cmp__


#
# =============================================================================
#
#                              Document Interface
#
# =============================================================================
#


ExcessPowerSBCoincDef = lsctables.CoincDef(search = u"excesspower", search_coinc_type = 1, description = u"sim_burst<-->sngl_burst coincidences")
ExcessPowerSCCoincDef = lsctables.CoincDef(search = u"excesspower", search_coinc_type = 2, description = u"sim_burst<-->coinc_event coincidences (exact)")
ExcessPowerSCNearCoincDef = lsctables.CoincDef(search = u"excesspower", search_coinc_type = 3, description = u"sim_burst<-->coinc_event coincidences (nearby)")


StringCuspSBCoincDef = lsctables.CoincDef(search = u"StringCusp", search_coinc_type = 1, description = u"sim_burst<-->sngl_burst coincidences")
StringCuspSCCoincDef = lsctables.CoincDef(search = u"StringCusp", search_coinc_type = 2, description = u"sim_burst<-->coinc_event coincidences (exact)")
StringCuspSCNearCoincDef = lsctables.CoincDef(search = u"StringCusp", search_coinc_type = 3, description = u"sim_burst<-->coinc_event coincidences (nearby)")


class DocContents(object):
	"""
	A wrapper interface to the XML document.
	"""
	def __init__(self, xmldoc, bbdef, sbdef, scedef, scndef, process):
		#
		# store the process row
		#

		self.process = process

		#
		# locate the sngl_burst and sim_burst tables
		#

		self.snglbursttable = table.get_table(xmldoc, lsctables.SnglBurstTable.tableName)
		self.simbursttable = table.get_table(xmldoc, lsctables.SimBurstTable.tableName)

		#
		# construct the zero-lag time slide needed to cover the
		# instruments listed in all the triggers, then determine
		# its ID (or create it if needed)
		#
		# FIXME:  in the future, the sim_burst table should
		# indicate time slide at which the injection was done
		#

		self.tisi_id = llwapp.get_time_slide_id(xmldoc, {}.fromkeys(self.snglbursttable.getColumnByName("ifo"), 0.0), create_new = process)

		#
		# get coinc_definer row for sim_burst <--> sngl_burst
		# coincs; this creates a coinc_definer table if the
		# document doesn't have one
		#

		self.sb_coinc_def_id = llwapp.get_coinc_def_id(xmldoc, sbdef.search, sbdef.search_coinc_type, create_new = True, description = sbdef.description)

		#
		# get coinc_def_id's for sngl_burst <--> sngl_burst, and
		# both kinds of sim_burst <--> coinc_event coincs.  set all
		# to None if this document does not contain any sngl_burst
		# <--> sngl_burst coincs.
		#

		try:
			bb_coinc_def_id = llwapp.get_coinc_def_id(xmldoc, bbdef.search, bbdef.search_coinc_type, create_new = False)
		except KeyError:
			bb_coinc_def_id = None
			self.sce_coinc_def_id = None
			self.scn_coinc_def_id = None
		else:
			self.sce_coinc_def_id = llwapp.get_coinc_def_id(xmldoc, scedef.search, scedef.search_coinc_type, create_new = True, description = scedef.description)
			self.scn_coinc_def_id = llwapp.get_coinc_def_id(xmldoc, scndef.search, scndef.search_coinc_type, create_new = True, description = scndef.description)

		#
		# get coinc table, create one if needed
		#

		try:
			self.coinctable = table.get_table(xmldoc, lsctables.CoincTable.tableName)
		except ValueError:
			self.coinctable = lsctables.New(lsctables.CoincTable)
			xmldoc.childNodes[0].appendChild(self.coinctable)
		self.coinctable.sync_next_id()

		#
		# get coinc_map table, create one if needed
		#

		try:
			self.coincmaptable = table.get_table(xmldoc, lsctables.CoincMapTable.tableName)
		except ValueError:
			self.coincmaptable = lsctables.New(lsctables.CoincMapTable)
			xmldoc.childNodes[0].appendChild(self.coincmaptable)

		#
		# index the document
		#
		# FIXME:  burst<-->burst coincs should be organized by time
		# slide ID, but since injections are only done at zero lag
		# for now this is ignored.
		#

		# index sngl_burst table
		index = dict((row.event_id, row) for row in self.snglbursttable)
		# find IDs of burst<-->burst coincs
		self.coincs = dict((row.coinc_event_id, []) for row in self.coinctable if row.coinc_def_id == bb_coinc_def_id)
		# construct event list for each burst<-->burst coinc
		for row in self.coincmaptable:
			try:
				self.coincs[row.coinc_event_id].append(index[row.event_id])
			except KeyError:
				continue
		del index
		# sort each event list by peak time and convert to tuples
		# for speed
		for coinc_event_id, events in self.coincs.items():
			events.sort(lambda a, b: cmp(a.peak_time, b.peak_time) or cmp(a.peak_time_ns, b.peak_time_ns))
			self.coincs[coinc_event_id] = tuple(events)
		# convert dictionary to a list
		self.coincs = self.coincs.items()

		#
		# sort sngl_burst table by peak time, and sort the coincs
		# list by the peak time of the first (earliest) event in
		# each coinc (recall that the event tuple for each coinc
		# has been time-ordered)
		#

		self.snglbursttable.sort(lambda a, b: cmp(a.peak_time, b.peak_time) or cmp(a.peak_time_ns, b.peak_time_ns))
		self.coincs.sort(lambda (id_a, a), (id_b, b): cmp(a[0].peak_time, b[0].peak_time) or cmp(a[0].peak_time_ns, b[0].peak_time_ns))

		#
		# set the window for bursts_near_peaktime().  this window
		# is the amount of time such that if an injection's peak
		# time and a burst event's peak time differ by more than
		# this it is *impossible* for them to match one another.
		#

		# the radius of Earth in light seconds.  (the most an
		# injection's peak time column can differ from the time it
		# peaks in an instrument)  6.378140e6 m = mean radius of
		# earth at equator, 299792458 m/s = c, 1.5 = add 50% for
		# good luck.  (constants copied from LALConstants.h)
		self.burst_peak_time_window = 6.378140e6 / 299792458 * 1.5

		# add the duration of the longest burst event (the most an
		# event's peak time could differ from either the start or
		# stop time of the event)
		if len(self.snglbursttable):
			self.burst_peak_time_window += max(self.snglbursttable.getColumnByName("duration"))

		#
		# set the window for identifying coincs near a peak time
		# FIXME:  this is kinda specific to the excess power
		# search.
		#

		self.coinc_peak_time_window = self.burst_peak_time_window + SimBurstUtils.burst_is_near_injection_window

	def bursts_near_peaktime(self, t):
		"""
		Return a list of the burst events whose peak times are
		within self.burst_peak_time_window of t.
		"""
		return self.snglbursttable[bisect.bisect_left(self.snglbursttable, t - self.burst_peak_time_window):bisect.bisect_right(self.snglbursttable, t + self.burst_peak_time_window)]

	def coincs_near_peaktime(self, t):
		"""
		Return a list of the (coinc_event_id, event list) tuples in
		which at least one burst event's peak time is within
		self.coinc_peak_time_window of t.
		"""
		# FIXME:  this test does not consider the time slide
		# offsets that should be applied to the coinc, but for now
		# injections are done at zero lag so this isn't a problem
		# yet
		return [(coinc_event_id, events) for coinc_event_id, events in self.coincs if (t - self.coinc_peak_time_window <= events[-1].get_peak()) and (events[0].get_peak() <= t + self.coinc_peak_time_window)]

	def sort_triggers_by_id(self):
		"""
		Sort the sngl_burst table's rows by ID (tidy-up document
		for output).
		"""
		self.snglbursttable.sort(lambda a, b: cmp(a.event_id, b.event_id))

	def new_coinc(self, coinc_def_id):
		"""
		Construct a new coinc_event row attached to the given
		process, and belonging to the set of coincidences defined
		by the given coinc_def_id.
		"""
		coinc = self.coinctable.RowType()
		coinc.process_id = self.process.process_id
		coinc.coinc_def_id = coinc_def_id
		coinc.coinc_event_id = self.coinctable.get_next_id()
		coinc.time_slide_id = self.tisi_id
		coinc.set_instruments(None)
		coinc.nevents = 0
		coinc.likelihood = None
		self.coinctable.append(coinc)
		return coinc


#
# =============================================================================
#
#                           Add Process Information
#
# =============================================================================
#


process_program_name = "ligolw_binjfind"


def append_process(xmldoc, match_algorithm, comment):
	"""
	Convenience wrapper for adding process metadata to the document.
	"""
	process = llwapp.append_process(xmldoc, program = process_program_name, version = __version__, cvs_repository = u"lscsoft", cvs_entry_time = __date__, comment = comment)

	params = [(u"--match-algorithm", u"lstring", match_algorithm)]
	ligolw_process.append_process_params(xmldoc, process, params)

	return process


#
# =============================================================================
#
#                 Injection <--> Burst Event Comparison Tests
#
# =============================================================================
#


def StringCuspSnglCompare(sim, burst):
	"""
	Return False if the peak time of the injection sim lies within the
	time interval of burst.
	"""
	return SimBurstUtils.time_at_instrument(sim, burst.ifo) not in burst.get_period()


def ExcessPowerSnglCompare(sim, burst):
	"""
	Return False if the peak time and centre frequency of sim lie
	within the time-frequency tile of burst.
	"""
	return StringCuspSnglCompare(sim, burst) or (sim.frequency not in burst.get_band())


def NearCoincCompare(sim, burst):
	"""
	Return False if the peak time of the sim is "near" the burst event.
	"""
	return not SimBurstUtils.burst_is_near_injection(sim, burst.start_time, burst.start_time_ns, burst.duration, burst.ifo)


#
# =============================================================================
#
#                 Build sim_burst <--> sngl_burst Coincidences
#
# =============================================================================
#


def find_sngl_burst_matches(contents, sim, comparefunc):
	"""
	Scan the burst table for triggers matching sim.
	"""
	return [burst for burst in contents.bursts_near_peaktime(sim.get_time_geocent()) if not comparefunc(sim, burst)]


def add_sim_burst_coinc(contents, sim, events):
	"""
	Create a coinc_event in the coinc table, and add arcs in the
	coinc_event_map table linking the sim_burst row and the list of
	sngl_burst rows to the new coinc_event row.
	"""
	coinc = contents.new_coinc(contents.sb_coinc_def_id)
	coinc.set_instruments(event.ifo for event in events)
	coinc.nevents = len(events)

	coincmap = contents.coincmaptable.RowType()
	coincmap.coinc_event_id = coinc.coinc_event_id
	coincmap.table_name = sim.simulation_id.table_name
	coincmap.event_id = sim.simulation_id
	contents.coincmaptable.append(coincmap)

	for event in events:
		coincmap = contents.coincmaptable.RowType()
		coincmap.coinc_event_id = coinc.coinc_event_id
		coincmap.table_name = event.event_id.table_name
		coincmap.event_id = event.event_id
		contents.coincmaptable.append(coincmap)


#
# =============================================================================
#
#                   Build sim_burst <--> coinc Coincidences
#
# =============================================================================
#


def find_exact_coinc_matches(coincs, sim, comparefunc):
	"""
	Return a list of the coinc_event_ids of the burst<-->burst coincs
	in which all burst events match sim.
	"""
	# FIXME:  this test does not consider the time slide offsets that
	# should be applied to the coinc, but for now injections are done
	# at zero lag so this isn't a problem yet
	return [coinc_event_id for coinc_event_id, events in coincs if True not in (bool(comparefunc(sim, event)) for event in events)]


def find_near_coinc_matches(coincs, sim, comparefunc):
	"""
	Return a list of the coinc_event_ids of the burst<-->burst coincs
	in which at least one burst event matches sim.
	"""
	# FIXME:  this test does not consider the time slide offsets that
	# should be applied to the coinc, but for now injections are done
	# at zero lag so this isn't a problem yet
	return [coinc_event_id for coinc_event_id, events in coincs if False in (bool(comparefunc(sim, event)) for event in events)]


def add_sim_coinc_coinc(contents, sim, coinc_event_ids, coinc_def_id):
	"""
	Create a coinc_event in the coinc table, and add arcs in the
	coinc_event_map table linking the sim_burst row and the list of
	coinc_event rows to the new coinc_event row.
	"""
	coinc = contents.new_coinc(coinc_def_id)
	coinc.nevents = len(coinc_event_ids)

	coincmap = contents.coincmaptable.RowType()
	coincmap.coinc_event_id = coinc.coinc_event_id
	coincmap.table_name = sim.simulation_id.table_name
	coincmap.event_id = sim.simulation_id
	contents.coincmaptable.append(coincmap)

	for coinc_event_id in coinc_event_ids:
		coincmap = contents.coincmaptable.RowType()
		coincmap.coinc_event_id = coinc.coinc_event_id
		coincmap.table_name = coinc_event_id.table_name
		coincmap.event_id = coinc_event_id
		contents.coincmaptable.append(coincmap)


#
# =============================================================================
#
#                                 Library API
#
# =============================================================================
#


def ligolw_binjfind(xmldoc, process, search, snglcomparefunc, nearcoinccomparefunc, verbose = False):
	#
	# Analyze the document's contents.
	#

	if verbose:
		print >>sys.stderr, "indexing ..."

	bbdef = {
		"StringCusp": ligolw_burca.StringCuspCoincDef,
		"excesspower": ligolw_burca.ExcessPowerCoincDef
	}[search]
	sbdef = {
		"StringCusp": StringCuspSBCoincDef,
		"excesspower": ExcessPowerSBCoincDef
	}[search]
	scedef = {
		"StringCusp": StringCuspSCCoincDef,
		"excesspower": ExcessPowerSCCoincDef
	}[search]
	scndef = {
		"StringCusp": StringCuspSCNearCoincDef,
		"excesspower": ExcessPowerSCNearCoincDef
	}[search]

	contents = DocContents(xmldoc = xmldoc, bbdef = bbdef, sbdef = sbdef, scedef = scedef, scndef = scndef, process = process)
	N = len(contents.simbursttable)

	#
	# Find sim_burst <--> sngl_burst coincidences.
	#

	if verbose:
		print >>sys.stderr, "constructing %s:" % sbdef.description
	for n, sim in enumerate(contents.simbursttable):
		if verbose:
			print >>sys.stderr, "\t%.1f%%\r" % (100.0 * n / N),
		events = find_sngl_burst_matches(contents, sim, snglcomparefunc)
		if events:
			add_sim_burst_coinc(contents, sim, events)
	if verbose:
		print >>sys.stderr, "\t100.0%"

	#
	# Find sim_burst <--> coinc_event coincidences.
	#

	if contents.sce_coinc_def_id:
		if verbose:
			print >>sys.stderr, "constructing %s and %s:" % (scedef.description, scndef.description)
		for n, sim in enumerate(contents.simbursttable):
			if verbose:
				print >>sys.stderr, "\t%.1f%%\r" % (100.0 * n / N),
			coincs = contents.coincs_near_peaktime(sim.get_time_geocent())
			coinc_event_ids = find_exact_coinc_matches(coincs, sim, snglcomparefunc)
			if coinc_event_ids:
				add_sim_coinc_coinc(contents, sim, coinc_event_ids, contents.sce_coinc_def_id)
			coinc_event_ids = find_near_coinc_matches(coincs, sim, nearcoinccomparefunc)
			if coinc_event_ids:
				add_sim_coinc_coinc(contents, sim, coinc_event_ids, contents.scn_coinc_def_id)
		if verbose:
			print >>sys.stderr, "\t100.0%"

	#
	# Restore the original event order.
	#

	if verbose:
		print >>sys.stderr, "finishing ..."
	contents.sort_triggers_by_id()

	#
	# Done.
	#

	return xmldoc
