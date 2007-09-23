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
from glue.ligolw import ilwd
from pylal import llwapp
from pylal import SnglBurstUtils
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


def sim_burst_get_geocent_peak(self):
	return LIGOTimeGPS(self.geocent_peak_time, self.geocent_peak_time_ns)


def sngl_burst_get_peak(self):
	return LIGOTimeGPS(self.peak_time, self.peak_time_ns)


def sngl_burst___cmp__(self, other):
	"""
	Compare self's peak time to the LIGOTimeGPS instance other.
	"""
	return cmp(self.peak_time, other.seconds) or cmp(self.peak_time_ns, other.nanoseconds)


lsctables.SimBurst.get_geocent_peak = sim_burst_get_geocent_peak
lsctables.SnglBurst.get_peak = sngl_burst_get_peak
lsctables.SnglBurst.__cmp__ = sngl_burst___cmp__


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

		self.tisi_id = llwapp.get_time_slide_id(xmldoc, {}.fromkeys(self.snglbursttable.getColumnByName("ifo"), 0.0), create_new = process)

		#
		# identify the IDs of zero-lag time slides;
		# get_time_slide_id() above will have created the time
		# slide table if it didn't exist
		#

		zero_lag_time_slides = set(llwapp.get_zero_lag_time_slides(xmldoc).keys())

		#
		# get coinc_def_id for sim_burst <--> sngl_burst coincs;
		# this creates a coinc_definer table if the document
		# doesn't have one
		#

		self.sb_coinc_def_id = llwapp.get_coinc_def_id(xmldoc, [lsctables.SnglBurstTable.tableName, lsctables.SimBurstTable.tableName])

		#
		# get coinc_def_id's for sngl_burst <--> sngl_burst, and
		# sim_burst <--> coinc coincs.  set both to None if this
		# document does not contain any sngl_burst <--> sngl_burst
		# coincs.
		#

		try:
			self.bb_coinc_def_id = llwapp.get_coinc_def_id(xmldoc, [lsctables.SnglBurstTable.tableName], create_new = False)
		except KeyError:
			self.bb_coinc_def_id = None
			self.sc_coinc_def_id = None
		else:
			self.sc_coinc_def_id = llwapp.get_coinc_def_id(xmldoc, [lsctables.CoincTable.tableName, lsctables.SimBurstTable.tableName])

		#
		# get coinc table, create one if needed
		#

		try:
			self.coinctable = table.get_table(xmldoc, lsctables.CoincTable.tableName)
		except ValueError:
			self.coinctable = lsctables.New(lsctables.CoincTable)
			xmldoc.childNodes[0].appendChild(self.coinctable)
		self.coinctable.sync_ids()

		#
		# get coinc_map table, create one if needed
		#

		try:
			self.coincmaptable = table.get_table(xmldoc, lsctables.CoincMapTable.tableName)
		except ValueError:
			self.coincmaptable = lsctables.New(lsctables.CoincMapTable)
			xmldoc.childNodes[0].appendChild(self.coincmaptable)

		#
		# index the document.  after this, coinc_peak_time_window
		# is the time, in seconds, separating the peak times of the
		# most widely-separated event pair of all the burst <-->
		# burst coincs.
		#

		# index sngl_burst table
		index = {}
		for row in self.snglbursttable:
			index[row.event_id] = row
		# find IDs of burst<-->burst coincs
		self.index = {}
		for coinc in self.coinctable:
			if (coinc.coinc_def_id == self.bb_coinc_def_id) and (coinc.time_slide_id in zero_lag_time_slides):
				self.index[coinc.coinc_event_id] = []
		# construct event list for each burst<-->burst coinc
		sngl_burst_table_name = table.StripTableName(lsctables.SnglBurstTable.tableName)
		for row in self.coincmaptable:
			if (row.table_name == sngl_burst_table_name) and (row.coinc_event_id in self.index):
				self.index[row.coinc_event_id].append(index[row.event_id])
		# sort each event list by peak time, convert to tuples for
		# speed, and find peak time window
		self.coinc_peak_time_window = 0
		for id, events in self.index.items():
			events.sort(lambda a, b: cmp(a.peak_time, b.peak_time) or cmp(a.peak_time_ns, b.peak_time_ns))
			self.index[id] = tuple(events)
			self.coinc_peak_time_window = max(self.coinc_peak_time_window, float(events[-1].get_peak() - events[0].get_peak()))

		#
		# sort sngl_burst table by peak time
		#

		self.snglbursttable.sort(lambda a, b: cmp(a.peak_time, b.peak_time) or cmp(a.peak_time_ns, b.peak_time_ns))

		#
		# set the window for bursts_near_peaktime()
		#

		if len(self.snglbursttable):
			self.burst_peak_time_window = max(self.snglbursttable.getColumnByName("duration"))
		else:
			# max() doesn't like empty sequences
			self.burst_peak_time_window = 0.0

		#
		# set the window for identifying coincs near a peak time
		#

		self.coinc_peak_time_window += self.burst_peak_time_window

	def bursts_near_peaktime(self, t):
		"""
		Return a list of the burst events with peak times are
		within self.burst_peak_time_window of t.
		"""
		return self.snglbursttable[bisect.bisect_left(self.snglbursttable, t - self.burst_peak_time_window):bisect.bisect_right(self.snglbursttable, t + self.burst_peak_time_window)]

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
		coinc.likelihood = 1.0
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


def append_process(xmldoc, **kwargs):
	"""
	Convenience wrapper for adding process metadata to the document.
	"""
	process = llwapp.append_process(xmldoc, program = process_program_name, version = __version__, cvs_repository = "lscsoft", cvs_entry_time = __date__, comment = kwargs["comment"])

	params = [("--match-algorithm", "lstring", kwargs["match_algorithm"])]
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
	return [burst for burst in contents.bursts_near_peaktime(sim.get_geocent_peak()) if not comparefunc(sim, burst)]


def add_sim_burst_coinc(contents, process, sim, bursts):
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


def cut_noninjection_coincs(contents, injection_burst_ids):
	"""
	Using the list of burst events known to be injections, delete
	entries from the coinc event list that contain bursts not in this
	list.  Doing this before checking each sim against the coincident
	events greatly speeds up that search.  Infact, after this step the
	only coincs remaining are exactly the coincident injections;  all
	that remains to be done is determine which injection matches which
	coinc.
	"""
	coincs_to_delete = []
	for coinc, bursts in contents.index.iteritems():
		for burst in bursts:
			if id(burst) not in injection_burst_ids:
				coincs_to_delete.append(coinc)
				break
	# can't modify dictionary inside loop
	for coinc in coincs_to_delete:
		del contents.index[coinc]
	

def find_coinc_matches(contents, sim, comparefunc):
	"""
	Return a list of the burst<-->burst coincs in which all burst
	events match sim.
	"""
	# the first half of the conditional is a speed hack to reduce the
	# number of times the real test, which is costly, needs to be run.
	# perhaps imap() could be used instead of map to provide a simpler
	# early bail out
	return [coinc_event_id for coinc_event_id, bursts in contents.index.iteritems() if (abs(float(bursts[0].get_peak() - sim.get_geocent_peak())) <= contents.coinc_peak_time_window) and True not in map(lambda burst: comparefunc(sim, burst), bursts)]


def add_sim_coinc_coinc(contents, process, sim, coinc_event_ids):
	"""
	Create a coinc_event in the coinc table, and add arcs in the
	coinc_event_map table linking the sim_burst row and the list of
	coinc_event rows to the new coinc_event row.
	"""
	coinc = contents.new_coinc(process, contents.sc_coinc_def_id)
	coinc.nevents = len(coinc_event_ids)

	coincmap = lsctables.CoincMap()
	coincmap.coinc_event_id = coinc.coinc_event_id
	coincmap.table_name = ilwd.ILWDTableName(sim.simulation_id)
	coincmap.event_id = sim.simulation_id
	contents.coincmaptable.append(coincmap)

	for coinc_event_id in coinc_event_ids:
		coincmap = lsctables.CoincMap()
		coincmap.coinc_event_id = coinc.coinc_event_id
		coincmap.table_name = ilwd.ILWDTableName(coinc_event_id)
		coincmap.event_id = coinc_event_id
		contents.coincmaptable.append(coincmap)


#
# =============================================================================
#
#                                 Library API
#
# =============================================================================
#


def ligolw_binjfind(xmldoc, **kwargs):
	#
	# Add process metadata to document.
	#

	process = append_process(xmldoc, **kwargs)

	#
	# Analyze the document's contents.
	#

	if kwargs["verbose"]:
		print >>sys.stderr, "indexing ..."
	contents = DocContents(xmldoc, process)
	N = len(contents.simbursttable)

	#
	# Find sim_burst <--> sngl_burst coincidences and record them.
	# injection_burst_ids is a list of the Python "ids" (as returned by
	# the id() builtin) of the burst events that have been identified
	# as being the result of injections.
	#

	if kwargs["verbose"]:
		print >>sys.stderr, "constructing injection--burst coincidences:"
	injection_burst_ids = set()
	for n, sim in enumerate(contents.simbursttable):
		if kwargs["verbose"] and not (n % (N / 50 or 1)):
			print >>sys.stderr, "\t%.1f%%\r" % (100.0 * n / N),
		matches = find_sngl_burst_matches(contents, sim, kwargs["comparefunc"])
		if matches:
			add_sim_burst_coinc(contents, process, sim, matches)
			injection_burst_ids |= set(map(id, matches))
	if kwargs["verbose"]:
		print >>sys.stderr, "\t100.0%"

	#
	# Find sim_burst <--> coinc_event coincidences.
	#

	if contents.sc_coinc_def_id:
		if kwargs["verbose"]:
			print >>sys.stderr, "constructing injection--coinc coincidences:"
		cut_noninjection_coincs(contents, injection_burst_ids)
		for n, sim in enumerate(contents.simbursttable):
			if kwargs["verbose"] and not (n % (N / 50 or 1)):
				print >>sys.stderr, "\t%.1f%%\r" % (100.0 * n / N),
			matches = find_coinc_matches(contents, sim, kwargs["comparefunc"])
			if matches:
				add_sim_coinc_coinc(contents, process, sim, matches)
		if kwargs["verbose"]:
			print >>sys.stderr, "\t100.0%"

	#
	# Restore the original event order, and close out the process
	# metadata.
	#

	contents.sort_triggers_by_id()
	llwapp.set_process_end_time(process)

	#
	# Done.
	#

	return xmldoc
