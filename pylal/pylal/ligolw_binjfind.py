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
matches as burst<-->injection coincidences using the standard coincidence
infrastructure.
"""


import bisect
import sys

from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw import ilwd
from pylal import llwapp
from pylal import SnglBurstUtils

__author__ = "Kipp Cannon <kipp@gravity.phys.uwm.edu>"
__version__ = "$Revision$"[11:-2]
__date__ = "$Date$"[7:-2]


#
# =============================================================================
#
#                              Document Interface
#
# =============================================================================
#


class DocContents(object):
	"""
	A wrapper interface to the XML document.
	"""
	def __init__(self, xmldoc, process):
		# locate the sngl_burst and sim_burst tables
		self.snglbursttable = table.get_table(xmldoc, lsctables.SnglBurstTable.tableName)
		self.simbursttable = table.get_table(xmldoc, lsctables.SimBurstTable.tableName)

		# construct the zero-lag time slide needed to cover the
		# instruments listed in all the triggers, then determine
		# its ID (or create it if needed)
		time_slide = {}
		for instrument in self.snglbursttable.getColumnByName("ifo"):
			time_slide[instrument] = 0.0
		self.tisi_id = llwapp.get_time_slide_id(xmldoc, time_slide, create_new = process)

		# locate the time_slide table;  get_time_slide_id() above
		# will have created the table if it didn't exist
		self.tisitable = table.get_table(xmldoc, lsctables.TimeSlideTable.tableName)

		# get coinc_def_id for sim_burst <--> sngl_burst coincs
		self.sb_coinc_def_id = llwapp.get_coinc_def_id(xmldoc, [lsctables.SnglBurstTable.tableName, lsctables.SimBurstTable.tableName])

		# get coinc_def_id's for sngl_burst <--> sngl_burst, and
		# sim_burst <--> coinc coincs.
		try:
			self.bb_coinc_def_id = llwapp.get_coinc_def_id(xmldoc, [lsctables.SnglBurstTable.tableName], create_new = False)
			self.sc_coinc_def_id = llwapp.get_coinc_def_id(xmldoc, [lsctables.CoincTable.tableName, lsctables.SimBurstTable.tableName])
		except KeyError:
			self.bb_coinc_def_id = None
			self.sc_coinc_def_id = None

		# get coinc table, create one if needed
		try:
			self.coinctable = table.get_table(xmldoc, lsctables.CoincTable.tableName)
		except ValueError:
			self.coinctable = lsctables.New(lsctables.CoincTable)
			xmldoc.childNodes[0].appendChild(self.coinctable)
		self.coinctable.sync_ids()

		# get coinc_map table, create one if needed
		try:
			self.coincmaptable = table.get_table(xmldoc, lsctables.CoincMapTable.tableName)
		except ValueError:
			self.coincmaptable = lsctables.New(lsctables.CoincMapTable)
			xmldoc.childNodes[0].appendChild(self.coincmaptable)

		# build an index of zero-lag sngl_burst <--> sngl_burst
		# coincs
		burst_event = {}
		for row in self.snglbursttable:
			burst_event[ilwd.ILWDID(row.event_id)] = row
		coinc_maps = {}
		for row in self.coincmaptable:
			try:
				# assumes that all coincs are burst + burst
				burst = burst_event[ilwd.ILWDID(row.event_id)]
			except KeyError:
				continue
			try:
				coinc_maps[ilwd.ILWDID(row.coinc_event_id)].append(burst)
			except KeyError:
				coinc_maps[ilwd.ILWDID(row.coinc_event_id)] = [burst]
		del burst_event
		self.coinc_event = {}
		for coinc in self.coinctable:
			if (coinc.coinc_def_id == self.bb_coinc_def_id) and self.tisitable.is_null(coinc.time_slide_id):
				self.coinc_event[coinc] = coinc_maps[ilwd.ILWDID(coinc.coinc_event_id)]

		# sort triggers by start time truncated to integer seconds,
		# and construct an integer start time look-up table
		self.snglbursttable.sort(lambda a, b: cmp(a.start_time, b.start_time))
		self.starttimes = self.snglbursttable.getColumnByName("start_time")

		# set the window for triggers_near_starttime();  add 2.0 s
		# because of the discreteness in the look-up table
		# FIXME: I think this can be reduced to 1.0 s, check that
		if len(self.starttimes):
			self.starttime_dt = max(self.snglbursttable.getColumnByName("duration")) + 2.0
		else:
			# max() doesn't like empty sequences
			self.starttime_dt = 2.0

	def triggers_near_starttime(self, t):
		"""
		Return a list of the triggers with start times within a
		reasonable window around t.
		"""
		return self.snglbursttable[bisect.bisect_left(self.starttimes, t - self.starttime_dt):bisect.bisect_right(self.starttimes, t + self.starttime_dt)]

	def sort_triggers_by_id(self):
		"""
		Sort the sngl_burst table's rows by ID (tidy-up document
		for output).
		"""
		self.snglbursttable.sort(lambda a, b: cmp(ilwd.ILWDID(a.event_id), ilwd.ILWDID(b.event_id)))

	def new_coinc(self, process, coinc_def_id):
		"""
		Construct a new coinc_event row attached to the given
		process, and belonging to the set of coincidences defined
		by the given coinc_def_id.
		"""
		coinc = lsctables.Coinc()
		coinc.process_id = process.process_id
		coinc.coinc_def_id = coinc_def_id
		coinc.coinc_event_id = self.coinctable.ids.next()
		coinc.time_slide_id = self.tisi_id
		coinc.nevents = 0
		self.coinctable.append(coinc)
		return coinc


#
# =============================================================================
#
#                           Add Process Information
#
# =============================================================================
#


def append_process(xmldoc, **kwargs):
	"""
	Convenience wrapper for adding process metadata to the document.
	"""
	process = llwapp.append_process(xmldoc, program = "ligolw_binjfind", version = __version__, cvs_repository = "lscsoft", cvs_entry_time = __date__, comment = kwargs["comment"])

	params = [("--compare", "lstring", kwargs["compare"])]
	llwapp.append_process_params(xmldoc, process, params)

	return process


#
# =============================================================================
#
#                 Injection <--> Burst Event Comparison Tests
#
# =============================================================================
#


def StringCuspCompare(sim, burst):
	"""
	Return False if the peak time of the injection sim lies within the
	time interval of burst.
	"""
	if sim.coordinates == "ZENITH":
		return sim.get_geocent_peak() not in burst.get_period()
	else:
		return sim.get_peak(burst.ifo) not in burst.get_period()


def ExcessPowerCompare(sim, burst):
	"""
	Return False if the peak time and centre frequency of sim lie
	within the time-frequency tile of burst.
	"""
	return StringCuspCompare(sim, burst) or (sim.freq not in burst.get_band())


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
	return [burst for burst in contents.triggers_near_starttime(sim.geocent_peak_time) if not comparefunc(sim, burst)]


def add_sim_burst_sngl_burst_coinc(contents, process, sim, bursts):
	"""
	Create a coinc_event in the coinc table, and add arcs in the
	coinc_event_map table linking the sim_burst row and the list of
	sngl_burst rows to the new coinc_event row.
	"""
	coinc = contents.new_coinc(process, contents.sb_coinc_def_id)
	coinc.nevents = len(bursts)

	coincmap = lsctables.CoincMap()
	coincmap.coinc_event_id = coinc.coinc_event_id
	coincmap.table_name = ilwd.ILWDTableName(sim.simulation_id)
	coincmap.event_id = sim.simulation_id
	contents.coincmaptable.append(coincmap)

	for event in bursts:
		coincmap = lsctables.CoincMap()
		coincmap.coinc_event_id = coinc.coinc_event_id
		coincmap.table_name = ilwd.ILWDTableName(event.event_id)
		coincmap.event_id = event.event_id
		contents.coincmaptable.append(coincmap)


#
# =============================================================================
#
#                   Build sim_burst <--> coinc Coincidences
#
# =============================================================================
#


def cut_noninjection_coincs(contents, injections):
	"""
	Using the list of burst events known to be injections, delete
	entries from the coinc event list that contain bursts not in this
	list.  Doing this before checking each sim against the coincident
	events greatly speeds up that search.  Infact, after this step the
	only coincs remaining are exactly the coincident injections;  all
	that remains to be done is determine which injection matches which
	coinc.
	"""
	injections.sort()
	non_injection_coincs = []
	for coinc, bursts in contents.coinc_event.iteritems():
		for burst in bursts:
			if not llwapp.bisect_contains(injections, burst):
				non_injection_coincs.append(coinc)
				break
	# can't modify dictionary inside loop
	for coinc in non_injection_coincs:
		del contents.coinc_event[coinc]
	

def find_coinc_matches(contents, sim, comparefunc):
	"""
	Scan the coinc_event table for matching burst coincs coincident
	with sim.
	"""
	matches = []
	for coinc, bursts in contents.coinc_event.iteritems():
		for burst in bursts:
			if comparefunc(sim, burst):
				break
		else:
			matches.append(coinc)
	return matches


def add_sim_burst_coinc_coinc(contents, process, sim, coinc_events):
	"""
	Create a coinc_event in the coinc table, and add arcs in the
	coinc_event_map table linking the sim_burst row and the list of
	coinc_event rows to the new coinc_event row.
	"""
	coinc = contents.new_coinc(process, contents.sc_coinc_def_id)
	coinc.nevents = len(coinc_events)

	coincmap = lsctables.CoincMap()
	coincmap.coinc_event_id = coinc.coinc_event_id
	coincmap.table_name = ilwd.ILWDTableName(sim.simulation_id)
	coincmap.event_id = sim.simulation_id
	contents.coincmaptable.append(coincmap)

	for event in coinc_events:
		coincmap = lsctables.CoincMap()
		coincmap.coinc_event_id = coinc.coinc_event_id
		coincmap.table_name = ilwd.ILWDTableName(event.event_id)
		coincmap.event_id = event.coinc_event_id
		contents.coincmaptable.append(coincmap)


#
# =============================================================================
#
#                                 Library API
#
# =============================================================================
#


def ligolw_binjfind(xmldoc, **kwargs):
	process = append_process(xmldoc, **kwargs)
	if kwargs["verbose"]:
		print >>sys.stderr, "indexing ..."
	contents = DocContents(xmldoc, process)

	if kwargs["verbose"]:
		print >>sys.stderr, "constructing injection-burst coincidences:"
	injections = []	# list of burst events that are injections
	for n, sim in enumerate(contents.simbursttable):
		if kwargs["verbose"] and not (n % (n / 50 or 1)):
			print >>sys.stderr, "\t%.1f%%\r" % (100.0 * n / len(contents.simbursttable)),
		matches = find_sngl_burst_matches(contents, sim, kwargs["comparefunc"])
		if matches:
			add_sim_burst_sngl_burst_coinc(contents, process, sim, matches)
			injections.extend(matches)
	if kwargs["verbose"]:
		print >>sys.stderr, "\t100.0%"

	if contents.sc_coinc_def_id:
		if kwargs["verbose"]:
			print >>sys.stderr, "constructing injection-coinc coincidences:"
		cut_noninjection_coincs(contents, injections)
		for n, sim in enumerate(contents.simbursttable):
			if kwargs["verbose"] and not (n % (n / 50 or 1)):
				print >>sys.stderr, "\t%.1f%%\r" % (100.0 * n / len(contents.simbursttable)),
			matches = find_coinc_matches(contents, sim, kwargs["comparefunc"])
			if matches:
				add_sim_burst_coinc_coinc(contents, process, sim, matches)
		if kwargs["verbose"]:
			print >>sys.stderr, "\t100.0%"
	contents.sort_triggers_by_id()
	llwapp.set_process_end_time(process)
	return xmldoc
