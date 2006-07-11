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

import matplotlib
matplotlib.use("Agg")	# use Agg backend
matplotlib.rcParams["text.usetex"] = True	# render all text with TeX
from matplotlib import figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
import sys

from glue.ligolw import table
from glue.ligolw import lsctables
from pylal import llwapp
from pylal.date import LIGOTimeGPS


#
# =============================================================================
#
#                                   Database
#
# =============================================================================
#

#
# This is a FUGLY work-in-progress;  please be kind.
#

class SnglBurstTable(table.Table):
	tableName = lsctables.SnglBurstTable.tableName
	validcolumns = lsctables.SnglBurstTable.validcolumns
	connection = None

	def __init__(self, *attrs):
		table.Table.__init__(self, *attrs)
		self.cursor = self.connection.cursor()
		self.cursor.execute("CREATE TABLE sngl_burst (process_id INTEGER, ifo TEXT, search TEXT, channel TEXT, start_time INTEGER, start_time_ns INTEGER, peak_time INTEGER, peak_time_ns INTEGER, duration REAL, central_freq REAL, bandwidth REAL, amplitude REAL, snr REAL, confidence REAL, tfvolume REAL, event_id INTEGER UNIQUE PRIMARY KEY)")

	def append(self, row):
		self.cursor.execute("INSERT INTO sngl_burst VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", (lsctables.ILWDID(row.process_id), row.ifo, row.search, row.channel, row.start_time, row.start_time_ns, row.peak_time, row.peak_time_ns, row.duration, row.central_freq, row.bandwidth, row.amplitude, row.snr, row.confidence, row.tfvolume, lsctables.ILWDID(row.event_id)))

	def __len__(self):
		return self.cursor.execute("SELECT COUNT(*) FROM sngl_burst").fetchone()[0]

	def _row_from_cols(cls, values):
		row = cls.RowType()
		row.process_id, row.ifo, row.search, row.channel, row.start_time, row.start_time_ns, row.peak_time, row.peak_time_ns, row.duration, row.central_freq, row.bandwidth, row.amplitude, row.snr, row.confidence, row.tfvolume, row.event_id = values
		return row
	_row_from_cols = classmethod(_row_from_cols)

	def __getitem__(self, id):
		self.cursor.execute("SELECT * FROM sngl_burst WHERE event_id == ?", (id,))
		return self._row_from_cols(self.cursor.fetchone())

	def select(self):
		for values in self.connection.cursor().execute("SELECT * FROM sngl_burst"):
			yield self._row_from_cols(values)

	def unlink(self):
		table.Table.unlink(self)
		self.cursor.execute("DROP TABLE sngl_burst")


class SnglBurst(lsctables.SnglBurst):
	def get_peak(self):
		return LIGOTimeGPS(self.peak_time, self.peak_time_ns)

SnglBurstTable.RowType = SnglBurst


class SimBurstTable(table.Table):
	tableName = lsctables.SimBurstTable.tableName
	validcolumns = lsctables.SimBurstTable.validcolumns
	connection = None

	def __init__(self, *attrs):
		table.Table.__init__(self, *attrs)
		self.cursor = self.connection.cursor()
		self.cursor.execute("CREATE TABLE sim_burst (process_id INTEGER, waveform TEXT, geocent_peak_time INTEGER, geocent_peak_time_ns INTEGER, h_peak_time INTEGER, h_peak_time_ns INTEGER, l_peak_time INTEGER, l_peak_time_ns INTEGER, peak_time_gmst REAL, dtminus REAL, dtplus REAL, longitude REAL, latitude REAL, coordinates TEXT, polarization REAL, hrss REAL, hpeak REAL, distance REAL, freq REAL, tau REAL, zm_number INTEGER, simulation_id INTEGER UNIQUE PRIMARY KEY)")

	def append(self, row):
		self.cursor.execute("INSERT INTO sim_burst VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", (lsctables.ILWDID(row.process_id), row.waveform, row.geocent_peak_time, row.geocent_peak_time_ns, row.h_peak_time, row.h_peak_time_ns, row.l_peak_time, row.l_peak_time_ns, row.peak_time_gmst, row.dtminus, row.dtplus, row.longitude, row.latitude, row.coordinates, row.polarization, row.hrss, row.hpeak, row.distance, row.freq, row.tau, row.zm_number, lsctables.ILWDID(row.simulation_id)))

	def __len__(self):
		return self.cursor.execute("SELECT COUNT(*) FROM sim_burst").fetchone()[0]

	def _row_from_cols(cls, values):
		row = cls.RowType()
		row.process_id, row.waveform, row.geocent_peak_time, row.geocent_peak_time_ns, row.h_peak_time, row.h_peak_time_ns, row.l_peak_time, row.l_peak_time_ns, row.peak_time_gmst, row.dtminus, row.dtplus, row.longitude, row.latitude, row.coordinates, row.polarization, row.hrss, row.hpeak, row.distance, row.freq, row.tau, row.zm_number, row.simulation_id = values
		return row
	_row_from_cols = classmethod(_row_from_cols)

	def __getitem__(self, id):
		self.cursor.execute("SELECT * FROM sim_burst WHERE simulation_id == ?", (id,))
		return self._row_from_cols(self.cursor.fetchone())

	def select(self):
		for values in self.connection.cursor().execute("SELECT * FROM sim_burst"):
			yield self._row_from_cols(values)

	def unlink(self):
		table.Table.unlink(self)
		self.cursor.execute("DROP TABLE sim_burst")


class SimBurst(lsctables.SimBurst):
	def get_geocent_peak(self):
		return LIGOTimeGPS(self.peak_time, self.peak_time_ns)

SimBurstTable.RowType = SimBurst


class TimeSlideTable(table.Table):
	tableName = lsctables.TimeSlideTable.tableName
	validcolumns = lsctables.TimeSlideTable.validcolumns
	connection = None

	def __init__(self, *attrs):
		table.Table.__init__(self, *attrs)
		self.cursor = self.connection.cursor()
		self.cursor.execute("CREATE TABLE time_slide (process_id INTEGER, time_slide_id INTEGER, instrument TEXT, offset REAL)")

	def append(self, row):
		self.cursor.execute("INSERT INTO time_slide VALUES (?, ?, ?, ?)", (lsctables.ILWDID(row.process_id), lsctables.ILWDID(row.time_slide_id), row.instrument, row.offset))

	def __len__(self):
		return self.cursor.execute("SELECT COUNT(DISTINCT time_slide_id) FROM time_slide").fetchone()[0]

	def __getitem__(self, id):
		offsets = {}
		for instrument, offset in self.cursor.execute("SELECT instrument, offset FROM time_slide WHERE time_slide_id == ?", (id,)):
			offsets[instrument] = offset
		return offsets

	def is_null(self, id):
		return not self.cursor.execute("SELECT COUNT(time_slide_id) FROM time_slide WHERE time_slide_id == ? AND offset != 0 LIMIT 1", (id,)).fetchone()[0]

	def all_offsets(self):
		return [self[id] for (id,) in self.connection.cursor().execute("SELECT DISTINCT time_slide_id FROM time_slide")]

	def unlink(self):
		table.Table.unlink(self)
		self.cursor.execute("DROP TABLE time_slide")


class CoincDefTable(table.Table):
	tableName = lsctables.CoincDefTable.tableName
	validcolumns = lsctables.CoincDefTable.validcolumns
	connection = None

	def __init__(self, *attrs):
		table.Table.__init__(self, *attrs)
		self.cursor = self.connection.cursor()
		self.cursor.execute("CREATE TABLE coinc_definer (coinc_def_id INTEGER, table_name TEXT)")

	def append(self, row):
		self.cursor.execute("INSERT INTO coinc_definer VALUES (?, ?)", (lsctables.ILWDID(row.coinc_def_id), row.table_name))

	def __len__(self):
		return self.cursor.execute("SELECT COUNT(DISTINCT coinc_def_id) FROM coinc_definer").fetchone()[0]

	def get_table_names(self, id):
		"""
		From a numeric ID, return a sorted list of table names or
		raise KeyError if no matching ID is found.
		"""
		l = [table_name for (table_name,) in self.cursor.execute("SELECT table_name FROM coinc_definer WHERE coinc_def_id == ?", (id,))]
		if not l:
			raise KeyError, id
		l.sort()
		return l

	def get_id(self, table_names):
		"""
		From a list of table names, return a numeric ID or raise
		KeyError if no matching ID is found.
		"""
		table_names = list(table_names)	# so we can modify it
		table_names.sort()
		for id in [id for (id,) in self.cursor.execute("SELECT DISTINCT coinc_def_id FROM coinc_definer")]:
			if self.get_table_names(id) == table_names:
				return id
		raise KeyError, table_names

	def unlink(self):
		table.Table.unlink(self)
		self.cursor.execute("DROP TABLE coinc_definer")


class CoincTable(table.Table):
	tableName = lsctables.CoincTable.tableName
	validcolumns = lsctables.CoincTable.validcolumns
	connection = None

	def __init__(self, *attrs):
		table.Table.__init__(self, *attrs)
		self.cursor = self.connection.cursor()
		self.cursor.execute("CREATE TABLE coinc_event (process_id INTEGER, coinc_def_id INTEGER, time_slide_id INTEGER, nevents INTEGER, coinc_event_id INTEGER UNIQUE PRIMARY KEY)")

	def append(self, row):
		self.cursor.execute("INSERT INTO coinc_event VALUES (?, ?, ?, ?, ?)", (lsctables.ILWDID(row.process_id), lsctables.ILWDID(row.coinc_def_id), lsctables.ILWDID(row.time_slide_id), row.nevents, lsctables.ILWDID(row.coinc_event_id)))

	def _row_from_cols(cls, values):
		row = cls.RowType()
		row.process_id, row.coinc_def_id, row.time_slide_id, row.nevents, row.coinc_event_id = values
		return row
	_row_from_cols = classmethod(_row_from_cols)

	def __getitem__(self, id):
		self.cursor.execute("SELECT * FROM coinc_event WHERE coinc_event_id == ?", (id,))
		return self._row_from_cols(self.cursor.fetchone())

	def __len__(self):
		return self.cursor.execute("SELECT COUNT(*) FROM coinc_event").fetchone()[0]

	def selectByDefID(self, coinc_def_id):
		for values in self.connection.cursor().execute("SELECT * FROM coinc_event WHERE coinc_def_id == ?", (coinc_def_id,)):
			yield self._row_from_cols(values)

	def unlink(self):
		table.Table.unlink(self)
		self.cursor.execute("DROP TABLE coinc_event")


class Coinc(lsctables.Coinc):
	def get_time_slide(self):
		offsets = {}
		for instrument, offset in CoincTable.connection.cursor().execute("SELECT instrument, offset FROM time_slide WHERE time_slide_id == ?", (self.time_slide_id,)):
			offsets[instrument] = offset
		return offsets

	def is_zero_lag(self):
		return not CoincTable.connection.cursor().execute("SELECT COUNT(offset) FROM time_slide WHERE time_slide_id == ? AND offset != 0", (self.time_slide_id,)).fetchone()[0]

	def sim_bursts(self):
		for values in CoincTable.connection.cursor().execute("SELECT * FROM sim_burst WHERE simulation_id IN (SELECT event_id FROM coinc_event_map WHERE table_name == 'sim_burst' AND coinc_event_id == ?)", (self.coinc_event_id,)):
			yield SimBurstTable._row_from_cols(values)

	def sngl_bursts(self):
		for values in CoincTable.connection.cursor().execute("SELECT * FROM sngl_burst WHERE event_id in (SELECT event_id FROM coinc_event_map WHERE table_name == 'sngl_burst' AND coinc_event_id == ?)", (self.coinc_event_id,)):
			yield SnglBurstTable._row_from_cols(values)

	def coincs(self):
		for values in CoincTable.connection.cursor().execute("SELECT * FROM coinc_event WHERE coinc_event_id IN (SELECT event_id FROM coinc_event_map WHERE table_name == 'coinc_event' AND coinc_event_id == ?)", (self.coinc_event_id,)):
			yield CoincTable._row_from_cols(values)

CoincTable.RowType = Coinc


class CoincMapTable(table.Table):
	tableName = lsctables.CoincMapTable.tableName
	validcolumns = lsctables.CoincMapTable.validcolumns
	connection = None

	def __init__(self, *attrs):
		table.Table.__init__(self, *attrs)
		self.cursor = self.connection.cursor()
		self.cursor.execute("CREATE TABLE coinc_event_map (coinc_event_id INTEGER, table_name TEXT, event_id INTEGER)")

	def append(self, row):
		self.cursor.execute("INSERT INTO coinc_event_map VALUES (?, ?, ?)", (lsctables.ILWDID(row.coinc_event_id), lsctables.ILWDTableName(row.event_id), lsctables.ILWDID(row.event_id)))

	def unlink(self):
		table.Table.unlink(self)
		self.cursor.execute("DROP TABLE coinc_event_map")


class CoincDatabase(object):
	def __init__(self, connection):
		self.connection = connection
		SnglBurstTable.connection = connection
		SimBurstTable.connection = connection
		TimeSlideTable.connection = connection
		CoincDefTable.connection = connection
		CoincTable.connection = connection
		CoincMapTable.connection = connection
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
		self.connection.commit()
		self.connection.execute("CREATE INDEX time_slide_id_index ON time_slide (time_slide_id)")
		self.connection.execute("CREATE INDEX coinc_event_id_index ON coinc_event_map (coinc_event_id)")

		# find the tables
		self.sngl_burst_table = llwapp.get_table(xmldoc, lsctables.SnglBurstTable.tableName)
		self.sim_burst_table = llwapp.get_table(xmldoc, lsctables.SimBurstTable.tableName)
		self.coinc_def_table = llwapp.get_table(xmldoc, lsctables.CoincDefTable.tableName)
		self.coinc_table = llwapp.get_table(xmldoc, lsctables.CoincTable.tableName)
		self.time_slide_table = llwapp.get_table(xmldoc, lsctables.TimeSlideTable.tableName)

		# get the segment lists
		self.seglists = llwapp.segmentlistdict_fromsearchsummary(xmldoc, live_time_program)
		self.instruments = self.seglists.keys()

		# determine a few coinc_definer IDs
		self.bb_definer_id = self.coinc_def_table.get_id([lsctables.SnglBurstTable.tableName])
		self.sb_definer_id = self.coinc_def_table.get_id([lsctables.SnglBurstTable.tableName, lsctables.SimBurstTable.tableName])
		self.sc_definer_id = self.coinc_def_table.get_id([lsctables.CoincTable.tableName, lsctables.SimBurstTable.tableName])

		# compute the missed injections
		self.missed_injections = [simulation_id for (simulation_id,) in self.sim_burst_table.cursor.execute("SELECT simulation_id FROM sim_burst")]
		for coinc in self.coinc_table.selectByDefID(self.sb_definer_id):
			for sim in coinc.sim_bursts():
				try:
					self.missed_injections.remove(sim.simulation_id)
				except ValueError:
					# already removed
					pass

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

		# the select inside the outer select finds a list of the
		# burst event_ids that were marked as coincident with an
		# injection; the outer select finds a list of the
		# burst+burst coinc_event_ids pointing to at least one of
		# those bursts
		self.incomplete_injection_coinc_ids = [coinc_event_id for (coinc_event_id,) in self.coinc_table.cursor.execute("SELECT DISTINCT coinc_event.coinc_event_id FROM coinc_event JOIN coinc_event_map ON (coinc_event.coinc_event_id == coinc_event_map.coinc_event_id) WHERE coinc_def_id == ? AND table_name == 'sngl_burst' AND event_id IN (SELECT DISTINCT event_id FROM coinc_event_map JOIN coinc_event ON (coinc_event_map.coinc_event_id == coinc_event.coinc_event_id) WHERE coinc_event_map.table_name == 'sngl_burst' AND coinc_event.coinc_def_id == ?)", (self.bb_definer_id, self.sb_definer_id))]

		# now remove the coinc_event_ids for which all bursts were
		# marked as injections
		map(self.incomplete_injection_coinc_ids.remove, [coinc_event_id for (coinc_event_id,) in self.coinc_table.cursor.execute("SELECT DISTINCT event_id FROM coinc_event_map JOIN coinc_event ON (coinc_event_map.coinc_event_id == coinc_event.coinc_event_id) WHERE table_name == 'coinc_event' AND coinc_def_id == ?", (self.sc_definer_id,))])

		# now remove coinc_event_ids for coincs that are not at
		# zero-lag
		self.incomplete_injection_coinc_ids = filter(lambda id: self.coinc_table[id].is_zero_lag(), self.incomplete_injection_coinc_ids)

		# verbosity
		if verbose:
			print >>sys.stderr, "database stats:"
			print >>sys.stderr, "\tburst events: %d" % len(self.sngl_burst_table)
			print >>sys.stderr, "\tinjections: %d" % len(self.sim_burst_table)
			print >>sys.stderr, "\ttime slides: %d" % len(self.time_slide_table)
			print >>sys.stderr, "\tburst + burst coincidences: %d" % self.coinc_table.cursor.execute("SELECT COUNT(*) FROM coinc_event WHERE coinc_def_id = ?", (self.bb_definer_id,)).fetchone()[0]
			print >>sys.stderr, "\tinjection + burst coincidences: %d" % self.coinc_table.cursor.execute("SELECT COUNT(*) FROM coinc_event WHERE coinc_def_id = ?", (self.sb_definer_id,)).fetchone()[0]
			print >>sys.stderr, "\tinjection + (burst + burst) coincidences: %d" % self.coinc_table.cursor.execute("SELECT COUNT(*) FROM coinc_event WHERE coinc_def_id = ?", (self.sc_definer_id,)).fetchone()[0]
			print >>sys.stderr, "\tburst + burst coincidences involving at least one injection: %d" % len(self.incomplete_injection_coinc_ids)


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
#                              Injection Related
#
# =============================================================================
#

def CompareSimBurstAndSnglBurstByTime(sim, burst):
	"""
	Return True if the peak time of the injection sim lies within the
	time interval of burst.
	"""
	if sim.coordinates == "ZENITH":
		return sim.get_geocent_peak() in burst.get_period()
	else:
		return sim.get_peak(burst.ifo) in burst.get_period()

def CompareSimBurstAndSnglBurstByTimeandFreq(sim, burst):
	"""
	Return True if the peak time and centre frequency of sim lie within
	the time-frequency tile of burst.
	"""
	return CompareSimBurstAndSnglBurstByTime(sim, burst) and (sim.freq in burst.get_band())


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
		self.fig.set_figsize_inches(16, 8)
		self.axes = self.fig.gca()
		self.axes.grid(True)
		self.axes.set_xlabel(x_label)
		self.axes.set_ylabel(y_label)

	def add_contents(self, doc):
		raise NotImplementedError

	def finish(self):
		pass
