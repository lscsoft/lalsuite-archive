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
matplotlib.use("Agg")	# use Agg backend
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
from matplotlib.backends.backend_agg import FigureCanvasAgg
import re
import sys

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


class SnglBurstTable(table.DBTable):
	tableName = lsctables.SnglBurstTable.tableName
	validcolumns = lsctables.SnglBurstTable.validcolumns
	constraints = "PRIMARY KEY (event_id)"

	def __getitem__(self, id):
		self.cursor.execute("SELECT * FROM sngl_burst WHERE event_id == ?", (id,))
		return self._row_from_cols(self.cursor.fetchone())


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


SnglBurstTable.RowType = SnglBurst


class SimBurstTable(table.DBTable):
	tableName = lsctables.SimBurstTable.tableName
	validcolumns = lsctables.SimBurstTable.validcolumns
	constraints = "PRIMARY KEY (simulation_id)"

	def __getitem__(self, id):
		self.cursor.execute("SELECT * FROM sim_burst WHERE simulation_id == ?", (id,))
		return self._row_from_cols(self.cursor.fetchone())


class SimBurst(lsctables.SimBurst):
	# use pylal.date.LIGOTimeGPS for speed
	def get_geocent_peak(self):
		return LIGOTimeGPS(self.geocent_peak_time, self.geocent_peak_time_ns)

	def get_peak(self, instrument):
		observatory = instrument[0]
		if observatory == "H":
			return LIGOTimeGPS(self.h_peak_time, self.h_peak_time_ns)
		if observatory == "L":
			return LIGOTimeGPS(self.l_peak_time, self.l_peak_time_ns)
		raise ValueError, instrument


SimBurstTable.RowType = SimBurst


class TimeSlideTable(table.DBTable):
	# this is a little different because multiple rows share ID
	tableName = lsctables.TimeSlideTable.tableName
	validcolumns = lsctables.TimeSlideTable.validcolumns
	constraints = "PRIMARY KEY (time_slide_id, instrument)"

	def __len__(self):
		return self.cursor.execute("SELECT COUNT(DISTINCT time_slide_id) FROM time_slide").fetchone()[0]

	def __getitem__(self, id):
		offsets = {}
		for instrument, offset in self.cursor.execute("SELECT instrument, offset FROM time_slide WHERE time_slide_id == ?", (id,)):
			offsets[instrument] = offset
		return offsets

	def iterkeys(self):
		for (id,) in self.connection.cursor().execute("SELECT DISTINCT time_slide_id FROM time_slide"):
			yield id

	def is_null(self, id):
		return not self.cursor.execute("SELECT EXISTS (SELECT * FROM time_slide WHERE time_slide_id == ? AND offset != 0.0)", (id,)).fetchone()[0]

	def all_offsets(self):
		return [self[id] for (id,) in self.connection.cursor().execute("SELECT DISTINCT time_slide_id FROM time_slide")]


class CoincDefTable(table.DBTable):
	# this is a little different because multiple rows share an ID
	tableName = lsctables.CoincDefTable.tableName
	validcolumns = lsctables.CoincDefTable.validcolumns

	def __len__(self):
		return self.cursor.execute("SELECT COUNT(DISTINCT coinc_def_id) FROM coinc_definer").fetchone()[0]

	def get_table_names(self, id):
		"""
		From an ID, return a sorted list of table names or raise
		KeyError if no matching ID is found.
		"""
		l = [table_name for (table_name,) in self.cursor.execute("SELECT table_name FROM coinc_definer WHERE coinc_def_id == ?", (id,))]
		if not l:
			raise KeyError, id
		l.sort()
		return l

	def get_id(self, table_names):
		"""
		From a list of table names, return an ID or raise KeyError
		if no matching ID is found.
		"""
		table_names = list(table_names)	# so we can modify it
		table_names.sort()
		for (id,) in self.connection.cursor().execute("SELECT DISTINCT coinc_def_id FROM coinc_definer"):
			if self.get_table_names(id) == table_names:
				return id
		raise KeyError, table_names


class CoincTable(table.DBTable):
	tableName = lsctables.CoincTable.tableName
	validcolumns = lsctables.CoincTable.validcolumns
	constraints = "PRIMARY KEY (coinc_event_id)"

	def __getitem__(self, id):
		self.cursor.execute("SELECT * FROM coinc_event WHERE coinc_event_id == ?", (id,))
		return self._row_from_cols(self.cursor.fetchone())

	def selectByDefID(self, coinc_def_id):
		for values in self.connection.cursor().execute("SELECT * FROM coinc_event WHERE coinc_def_id == ?", (coinc_def_id,)):
			yield self._row_from_cols(values)

	def selectByTimeSlideID(self, time_slide_id):
		for values in self.connection.cursor().execute("SELECT * FROM coinc_event WHERE time_slide_id == ?", (time_slide_id,)):
			yield self._row_from_cols(values)


class Coinc(lsctables.Coinc):
	def get_time_slide(self):
		offsets = {}
		for instrument, offset in CoincTable.connection.cursor().execute("SELECT instrument, offset FROM time_slide WHERE time_slide_id == ?", (self.time_slide_id,)):
			offsets[instrument] = offset
		return offsets

	def is_zero_lag(self):
		return not CoincTable.connection.cursor().execute("SELECT EXISTS (SELECT * FROM time_slide WHERE time_slide_id == ? AND offset != 0.0)", (self.time_slide_id,)).fetchone()[0]


CoincTable.RowType = Coinc


class CoincMapTable(table.DBTable):
	tableName = lsctables.CoincMapTable.tableName
	validcolumns = lsctables.CoincMapTable.validcolumns

	def _end_of_rows(self):
		table.DBTable._end_of_rows(self)
		self.cursor.execute("CREATE INDEX coinc_event_id_index ON coinc_event_map (table_name, coinc_event_id)")


class CoincDatabase(object):
	def __init__(self, connection):
		self.connection = connection
		table.DBTable.connection = connection
		lsctables.TableByName.update({
			table.StripTableName(SnglBurstTable.tableName): SnglBurstTable,
			table.StripTableName(SimBurstTable.tableName): SimBurstTable,
			table.StripTableName(TimeSlideTable.tableName): TimeSlideTable,
			table.StripTableName(CoincDefTable.tableName): CoincDefTable,
			table.StripTableName(CoincTable.tableName): CoincTable,
			table.StripTableName(CoincMapTable.tableName): CoincMapTable
		})

	def summarize(self, xmldoc, live_time_program, verbose = False):
		"""
		Compute and record some summary information about the
		database.  Call this after all the data has been inserted,
		and before you want any of this information.
		"""
		cursor = self.connection.cursor()

		# find the tables
		self.sngl_burst_table = table.get_table(xmldoc, lsctables.SnglBurstTable.tableName)
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
		self.instruments = self.seglists.keys()

		# determine a few coinc_definer IDs
		try:
			self.bb_definer_id = self.coinc_def_table.get_id([lsctables.SnglBurstTable.tableName])
		except KeyError:
			self.bb_definer_id = None
		try:
			self.sb_definer_id = self.coinc_def_table.get_id([lsctables.SnglBurstTable.tableName, lsctables.SimBurstTable.tableName])
		except KeyError:
			self.sb_definer_id = None
		try:
			self.sc_definer_id = self.coinc_def_table.get_id([lsctables.CoincTable.tableName, lsctables.SimBurstTable.tableName])
		except KeyError:
			self.sc_definer_id = None

		# compute the missed injections by instrument;  first
		# generate a copy of all injection IDs for each instrument,
		# then remove the ones that each instrument didn't find
		if self.sb_definer_id:
			self.missed_injections = {}
			for instrument in self.instruments:
				self.missed_injections[instrument] = [id for (id,) in cursor.execute(
					"""
SELECT simulation_id FROM
	sim_burst
WHERE
	simulation_id NOT IN (
		SELECT DISTINCT sim_burst.simulation_id FROM
			sim_burst
			JOIN coinc_event_map AS a ON (
				sim_burst.simulation_id == a.event_id
				AND a.table_name == 'sim_burst'
			)
			JOIN coinc_event_map AS b ON (
				b.table_name == 'sngl_burst'
				AND b.coinc_event_id == a.coinc_event_id
			)
			JOIN sngl_burst ON (
				b.event_id == sngl_burst.event_id
				AND sngl_burst.ifo == ?
			)
	)
					""",
					(instrument,))]
		else:
			self.missed_injections = {}
			for instrument in self.instruments:
				self.missed_injections[instrument] = []

		# determine burst <--> burst coincidences for which at
		# least one burst, but not all, was identified as an
		# injection;  these are places in the data where an
		# injection was done, a coincident event was seen, but
		# where, later, the injection was not found to match all
		# events in the coincidence;  these perhaps indicate power
		# leaking from the injection into nearby tiles, or
		# accidental coincidence with near-by noise, etc, and so
		# although they aren't "bang-on" reconstructions of
		# injections they are nevertheless injections that are
		# found and survive a coincidence cut.
		if self.sim_burst_table:
			# the select inside the first select finds a list
			# of the burst event_ids that were marked as
			# coincident with an injection; the outer select
			# finds a list of the burst+burst coinc_event_ids
			# pointing to at least one of those bursts;  the
			# except clause removes coinc_event_ids for which
			# all bursts were marked as injections;  the result
			# is a list of all coinc_event_ids for burst+burst
			# coincidences in which at least 1 but not all
			# burst events was identified as an injection
			self.incomplete_injection_coinc_ids = [coinc_event_id for (coinc_event_id,) in cursor.execute(
				"""
SELECT DISTINCT coinc_event.coinc_event_id FROM
	coinc_event
	JOIN coinc_event_map ON (
		coinc_event.coinc_event_id == coinc_event_map.coinc_event_id
	)
WHERE
	coinc_def_id == ?
	AND table_name == 'sngl_burst'
	AND event_id IN (
		SELECT DISTINCT event_id FROM
			coinc_event_map
			JOIN coinc_event ON (
				coinc_event_map.coinc_event_id == coinc_event.coinc_event_id
			)
			WHERE
				coinc_event_map.table_name == 'sngl_burst'
				AND coinc_event.coinc_def_id == ?
	)
EXCEPT
	SELECT DISTINCT event_id FROM
		coinc_event_map
		JOIN coinc_event ON (
			coinc_event_map.coinc_event_id == coinc_event.coinc_event_id
		)
		WHERE
			table_name == 'coinc_event'
			AND coinc_def_id == ?
				""",
				(self.bb_definer_id, self.sb_definer_id, self.sc_definer_id))]

			# remove coinc_event_ids for coincs that are
			# not at zero-lag
			self.incomplete_injection_coinc_ids = filter(lambda id: self.coinc_table[id].is_zero_lag(), self.incomplete_injection_coinc_ids)
		else:
			self.incomplete_injection_coinc_ids = []

		# verbosity
		if verbose:
			print >>sys.stderr, "database stats:"
			print >>sys.stderr, "\tburst events: %d" % len(self.sngl_burst_table)
			if self.sim_burst_table:
				print >>sys.stderr, "\tinjections: %d" % len(self.sim_burst_table)
			print >>sys.stderr, "\ttime slides: %d" % len(self.time_slide_table)
			print >>sys.stderr, "\tburst + burst coincidences: %d" % cursor.execute("SELECT COUNT(*) FROM coinc_event WHERE coinc_def_id = ?", (self.bb_definer_id,)).fetchone()[0]
			if self.sim_burst_table:
				print >>sys.stderr, "\tinjection + burst coincidences: %d" % cursor.execute("SELECT COUNT(*) FROM coinc_event WHERE coinc_def_id = ?", (self.sb_definer_id,)).fetchone()[0]
				print >>sys.stderr, "\tinjection + (burst + burst) coincidences: %d" % cursor.execute("SELECT COUNT(*) FROM coinc_event WHERE coinc_def_id = ?", (self.sc_definer_id,)).fetchone()[0]
				print >>sys.stderr, "\tburst + burst coincidences involving at least one injection: %d" % len(self.incomplete_injection_coinc_ids)

	def coinc_sim_bursts(self, coinc):
		for values in self.connection.cursor().execute("SELECT sim_burst.* FROM sim_burst JOIN coinc_event_map ON sim_burst.simulation_id == coinc_event_map.event_id WHERE coinc_event_map.table_name == 'sim_burst' AND coinc_event_map.coinc_event_id == ?", (coinc.coinc_event_id,)):
			yield self.sim_burst_table._row_from_cols(values)

	def coinc_sngl_bursts(self, coinc):
		for values in self.connection.cursor().execute("SELECT sngl_burst.* FROM sngl_burst JOIN coinc_event_map ON sngl_burst.event_id == coinc_event_map.event_id WHERE coinc_event_map.table_name == 'sngl_burst' AND coinc_event_map.coinc_event_id == ?", (coinc.coinc_event_id,)):
			yield self.sngl_burst_table._row_from_cols(values)

	def coinc_coincs(self, coinc):
		for values in self.connection.cursor().execute("SELECT coinc_event.* FROM coinc_event JOIN coinc_event_map ON coinc_event.coinc_event_id == coinc_event_map.event_id WHERE coinc_event_map.table_name == 'coinc_event' AND coinc_event_map.coinc_event_id == ?", (self.coinc_event_id,)):
			yield self.coinc_table._row_from_cols(values)


#
# =============================================================================
#
#                                  Clustering
#
# =============================================================================
#


def cmp_segs(a, b):
	"""
	Returns 1 if a covers an interval above b's interval, -1 if a
	covers an interval below b's, and 0 if the two intervals overlap
	(including if their edges touch).
	"""
	if a[0] > b[1]:
		return 1
	if a[1] < b[0]:
		return -1
	return 0


def CompareSnglBurstByPeakTime(a, b):
	"""
	Orders a and b by peak time.
	"""
	return cmp(a.get_peak(), b.get_peak())


def CompareSnglBurstByPeakTimeAndFreq(a, b):
	"""
	Orders a and b by peak time, then by frequency band.  Returns 0 if
	a and b have the same peak time, and their frequency bands
	intersect.
	"""
	return cmp(a.get_peak(), b.get_peak()) or cmp_segs(a.get_band(), b.get_band())


def CompareSnglBurst(a, b, twindow = LIGOTimeGPS(0)):
	"""
	Orders a and b by time interval, then by frequency band.  Returns 0
	if a and b's time-frequency tiles intersect.  A time window can be
	optionally applied, and the time-frequency tiles will continue to
	compare as equal if they do not overlap by as much as the window
	amount.
	"""
	return cmp_segs(a.get_period().protract(twindow), b.get_period()) or cmp_segs(a.get_band(), b.get_band())


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
		FigureCanvasAgg(self.fig)
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
