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


import math
import matplotlib
matplotlib.rcParams.update({
	"font.size": 8.0,
	"axes.titlesize": 10.0,
	"axes.labelsize": 10.0,
	"xtick.labelsize": 8.0,
	"ytick.labelsize": 8.0,
	"legend.fontsize": 8.0,
	"figure.dpi": 300,
	"savefig.dpi": 300,
	"text.usetex": True	# render all text with TeX
})
from matplotlib import figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import re
import sys
# Python 2.3 compatibility
try:
	set
except NameError:
	from sets import Set as set

from glue import segments
from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw import types
from glue.ligolw import ilwd
from pylal import llwapp
from pylal.date import LIGOTimeGPS


#
# =============================================================================
#
#                                   Database
#
# =============================================================================
#


class SnglBurst(lsctables.SnglBurst):
	# use pylal.date.LIGOTimeGPS for speed
	def get_peak(self):
		return LIGOTimeGPS(self.peak_time, self.peak_time_ns)

	def get_start(self):
		return LIGOTimeGPS(self.start_time, self.start_time_ns)

	def get_stop(self):
		return LIGOTimeGPS(self.stop_time, self.stop_time_ns)

	def get_period(self):
		start = LIGOTimeGPS(self.start_time, self.start_time_ns)
		return segments.segment(start, start + self.duration)


class CoincDatabase(object):
	def __init__(self):
		from glue.ligolw import dbtables
		self.connection = dbtables.DBTable_get_connection()
		dbtables.SnglBurstTable.RowType = SnglBurst

	def summarize(self, xmldoc, live_time_program, verbose = False):
		"""
		Compute and record some summary information about the
		database.  Call this after all the data has been inserted,
		and before you want any of this information.
		"""
		cursor = self.connection.cursor()

		# find the tables
		try:
			self.sngl_burst_table = table.get_table(xmldoc, lsctables.SnglBurstTable.tableName)
		except ValueError:
			self.sngl_burst_table = None
		try:
			self.sim_burst_table = table.get_table(xmldoc, lsctables.SimBurstTable.tableName)
		except ValueError:
			self.sim_burst_table = None
		try:
			self.coinc_def_table = table.get_table(xmldoc, lsctables.CoincDefTable.tableName)
			self.coinc_table = table.get_table(xmldoc, lsctables.CoincTable.tableName)
			self.time_slide_table = table.get_table(xmldoc, lsctables.TimeSlideTable.tableName)
		except ValueError:
			self.coinc_def_table = None
			self.coinc_table = None
			self.time_slide_table = None

		# get the segment lists
		self.seglists = llwapp.segmentlistdict_fromsearchsummary(xmldoc, live_time_program)
		self.instruments = set(self.seglists.keys())

		# determine a few coinc_definer IDs
		if self.coinc_def_table is not None:
			try:
				self.bb_definer_id = self.coinc_def_table.get_coinc_def_id([lsctables.SnglBurstTable.tableName], create_new = False)
			except KeyError:
				self.bb_definer_id = None
			try:
				self.sb_definer_id = self.coinc_def_table.get_coinc_def_id([lsctables.SnglBurstTable.tableName, lsctables.SimBurstTable.tableName], create_new = False)
			except KeyError:
				self.sb_definer_id = None
			try:
				self.sc_definer_id = self.coinc_def_table.get_coinc_def_id([lsctables.CoincTable.tableName, lsctables.SimBurstTable.tableName], create_new = False)
			except KeyError:
				self.sc_definer_id = None
		else:
			self.bb_definer_id = None
			self.sb_definer_id = None
			self.sc_definer_id = None

		# verbosity
		if verbose:
			print >>sys.stderr, "database stats:"
			print >>sys.stderr, "\tburst events: %d" % len(self.sngl_burst_table)
			if self.sim_burst_table is not None:
				print >>sys.stderr, "\tinjections: %d" % len(self.sim_burst_table)
			if self.time_slide_table is not None:
				print >>sys.stderr, "\ttime slides: %d" % cursor.execute("SELECT COUNT(DISTINCT(time_slide_id)) FROM time_slide").fetchone()[0]
			if self.bb_definer_id is not None:
				print >>sys.stderr, "\tburst + burst coincidences: %d" % cursor.execute("SELECT COUNT(*) FROM coinc_event WHERE coinc_def_id = ?", (self.bb_definer_id,)).fetchone()[0]
			if self.sb_definer_id is not None:
				print >>sys.stderr, "\tinjection + burst coincidences: %d" % cursor.execute("SELECT COUNT(*) FROM coinc_event WHERE coinc_def_id = ?", (self.sb_definer_id,)).fetchone()[0]
			if self.sc_definer_id is not None:
				print >>sys.stderr, "\tinjection + (burst + burst) coincidences: %d" % cursor.execute("SELECT COUNT(*) FROM coinc_event WHERE coinc_def_id = ?", (self.sc_definer_id,)).fetchone()[0]

		return self

	def coinc_select_by_def_id(self, coinc_def_id):
		for values in self.connection.cursor().execute("SELECT * FROM coinc_event WHERE coinc_def_id == ?", (coinc_def_id,)):
			yield self.coinc_table._row_from_cols(values)

	def coinc_select_by_time_slide_id(self, time_slide_id):
		for values in self.connection.cursor().execute("SELECT * FROM coinc_event WHERE time_slide_id == ?", (time_slide_id,)):
			yield self.coinc_table._row_from_cols(values)

	def coinc_sim_bursts(self, coinc):
		for values in self.connection.cursor().execute("""
SELECT sim_burst.* FROM
	sim_burst
	JOIN coinc_event_map ON (
		sim_burst.simulation_id == coinc_event_map.event_id
		AND coinc_event_map.table_name == 'sim_burst'
	)
WHERE
	coinc_event_map.coinc_event_id == ?
		""", (coinc.coinc_event_id,)):
			yield self.sim_burst_table._row_from_cols(values)

	def coinc_sngl_bursts(self, coinc):
		for values in self.connection.cursor().execute("""
SELECT sngl_burst.* FROM
	sngl_burst
	JOIN coinc_event_map ON (
		sngl_burst.event_id == coinc_event_map.event_id
		AND coinc_event_map.table_name == 'sngl_burst'
	)
WHERE
	coinc_event_map.coinc_event_id == ?
		""", (coinc.coinc_event_id,)):
			yield self.sngl_burst_table._row_from_cols(values)

	def coinc_coincs(self, coinc):
		for values in self.connection.cursor().execute("""
SELECT coinc_event.* FROM
	coinc_event
	JOIN coinc_event_map ON (
		coinc_event.coinc_event_id == coinc_event_map.event_id
		AND coinc_event_map.table_name == 'coinc_event'
	)
WHERE
	coinc_event_map.coinc_event_id == ?
		""", (coinc.coinc_event_id,)):
			yield self.coinc_table._row_from_cols(values)

	def incomplete_injection_coincs(self):
		# determine burst <--> burst coincidences for which at
		# least one burst, *but not all*, was identified as an
		# injection;  these are places in the data where an
		# injection was done, a coincident event was seen, but
		# where, later, the injection was not found to match all
		# events in the coincidence;  these perhaps indicate power
		# leaking from the injection into nearby tiles, or
		# accidental coincidence of a marginal injection with
		# near-by noise, etc, and so although they aren't "bang-on"
		# reconstructions of injections they are nevertheless
		# injections that are found and survive a coincidence cut.
		#
		# the select inside the first select finds a list of the
		# burst event_ids that were marked as coincident with an
		# injection; the outer select finds a list of the
		# burst+burst coinc_event_ids pointing to at least one of
		# those bursts;  the except clause removes coinc_event_ids
		# for which all bursts were marked as injections;  the
		# result is a list of all coinc_event_ids for burst+burst
		# coincidences in which at least 1 but not all burst events
		# was identified as an injection
		for (id,) in self.connection.cursor().execute("""
SELECT DISTINCT coinc_event.coinc_event_id FROM
	coinc_event
	JOIN coinc_event_map ON (
		coinc_event.coinc_event_id == coinc_event_map.coinc_event_id
	)
WHERE
	coinc_event.coinc_def_id == ?
	AND coinc_event_map.table_name == 'sngl_burst'
	AND coinc_event_map.event_id IN (
		SELECT DISTINCT coinc_event_map.event_id FROM
			coinc_event_map
			JOIN coinc_event ON (
				coinc_event_map.coinc_event_id == coinc_event.coinc_event_id
			)
		WHERE
			coinc_event_map.table_name == 'sngl_burst'
			AND coinc_event.coinc_def_id == ?
	)
	-- Get only zero-lag events (? why ?  I can't remember ...)
	AND NOT EXISTS (
		SELECT * FROM
			time_slide
		WHERE
			time_slide_id == coinc_event.time_slide_id
			AND offset != 0.0
	)
EXCEPT SELECT DISTINCT coinc_event_map.event_id FROM
	coinc_event_map
	JOIN coinc_event ON (
		coinc_event_map.coinc_event_id == coinc_event.coinc_event_id
	)
WHERE
	coinc_event_map.table_name == 'coinc_event'
	AND coinc_event.coinc_def_id == ?
		""", (self.bb_definer_id, self.sb_definer_id, self.sc_definer_id)):
			yield self.coinc_table[id]


	def found_injections(self, instrument):
		"""
		Iterate over found injections.
		"""
		for values in self.connection.cursor().execute("""
SELECT DISTINCT sim_burst.* FROM
	sim_burst
	JOIN coinc_event_map AS a ON (
		sim_burst.simulation_id == a.event_id
		AND a.table_name == 'sim_burst'
	)
	JOIN coinc_event_map AS b ON (
		a.coinc_event_id == b.coinc_event_id
		AND b.table_name == 'sngl_burst'
	)
	JOIN sngl_burst ON (
		b.event_id == sngl_burst.event_id
	)
WHERE
	sngl_burst.ifo == ?
		""", (instrument,)):
			yield self.sim_burst_table._row_from_cols(values)

	def missed_injections(self, instrument):
		"""
		Iterate over missed injections.
		"""
		for values in self.connection.cursor().execute("""
SELECT * FROM
	sim_burst
WHERE
	simulation_id NOT IN (
		SELECT sim_burst.simulation_id FROM
			sim_burst
			JOIN coinc_event_map AS a ON (
				a.event_id == sim_burst.simulation_id
				AND a.table_name == 'sim_burst'
			)
			JOIN coinc_event_map AS b ON (
				b.coinc_event_id == a.coinc_event_id
				AND b.table_name == 'sngl_burst'
			)
			JOIN sngl_burst ON (
				sngl_burst.event_id == b.event_id
			)
		WHERE
			sngl_burst.ifo == ?
	)
		""", (instrument,)):
			yield self.sim_burst_table._row_from_cols(values)


#
# =============================================================================
#
#                                 TeX Helpers
#
# =============================================================================
#


def latexnumber(s):
	"""
	Convert a string of the form "d.dddde-dd" to "d.dddd \times
	10^{-dd}"
	"""
	return re.sub(r"([+-]?[.0-9]+)[Ee]?([+-]?[0-9]+)", r"\1 \\times 10^{\2}", s)


#
# =============================================================================
#
#                                    Plots
#
# =============================================================================
#


class BurstPlotError(Exception):
	"""
	Used to relay error messages from plotting routines to
	applications.
	"""
	pass


class BurstPlot(object):
	def __init__(self, x_label, y_label):
		self.nevents = 0
		self.fig = figure.Figure()
		FigureCanvas(self.fig)
		# 6.5" wide, golden ratio high
		self.fig.set_size_inches(6.5, 6.5 / ((1 + math.sqrt(5)) / 2))
		#self.fig.set_size_inches(16, 8)
		self.axes = self.fig.gca()
		self.axes.grid(True)
		self.axes.set_xlabel(x_label)
		self.axes.set_ylabel(y_label)

	def add_contents(self, doc):
		raise NotImplementedError

	def finish(self):
		pass
