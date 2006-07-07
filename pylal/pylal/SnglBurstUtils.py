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
import sqlobject
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

def CoincDatabaseConnection():
	return sqlobject.connectionForURI("sqlite:/:memory:")


class CoincDatabase(object):
	#
	# Table classes, and append overrides for glue.ligolw
	#

	class SnglBurst(sqlobject.SQLObject):
		class sqlmeta:
			table = "sngl_burst"
			idName = "event_id"
		ifo = sqlobject.StringCol(length = 2)
		central_freq = sqlobject.FloatCol()
		confidence = sqlobject.FloatCol()
		peak_time = sqlobject.IntCol()
		peak_time_ns = sqlobject.IntCol()

		def get_peak(self):
			return LIGOTimeGPS(self.peak_time, self.peak_time_ns)


	def sngl_burst_append(self, row):
		CoincDatabase.SnglBurst(
			id = lsctables.ILWDID(row.event_id),
			ifo = row.ifo,
			central_freq = row.central_freq,
			confidence = row.confidence,
			peak_time = row.peak_time,
			peak_time_ns = row.peak_time_ns
		)
	sngl_burst_append = staticmethod(sngl_burst_append)


	class SimBurst(sqlobject.SQLObject):
		class sqlmeta:
			table = "sim_burst"
			idName = "simulation_id"
		geocent_peak_time = sqlobject.IntCol()
		geocent_peak_time_ns = sqlobject.IntCol()
		h_peak_time = sqlobject.IntCol()
		h_peak_time_ns = sqlobject.IntCol()
		l_peak_time = sqlobject.IntCol()
		l_peak_time_ns = sqlobject.IntCol()
		freq = sqlobject.FloatCol()
		hrss = sqlobject.FloatCol()

		def get_geocent_peak(self):
			return LIGOTimeGPS(self.peak_time, self.peak_time_ns)


	def sim_burst_append(self, row):
		CoincDatabase.SimBurst(
			id = lsctables.ILWDID(row.simulation_id),
			geocent_peak_time = row.geocent_peak_time,
			geocent_peak_time_ns = row.geocent_peak_time_ns,
			h_peak_time = row.h_peak_time,
			h_peak_time_ns = row.h_peak_time_ns,
			l_peak_time = row.l_peak_time,
			l_peak_time_ns = row.l_peak_time_ns,
			freq = row.freq,
			hrss = row.hrss
		)
	sim_burst_append = staticmethod(sim_burst_append)


	class TimeSlide(sqlobject.SQLObject):
		class sqlmeta:
			table = "time_slide"
		process_id = sqlobject.IntCol()
		time_slide_id = sqlobject.IntCol()
		instrument = sqlobject.StringCol(length = 2)
		offset = sqlobject.FloatCol()

		def get_offsets(cls, id):
			offsets = {}
			for instrument, offset in cls._connection.queryAll("SELECT instrument, offset FROM time_slide WHERE time_slide_id == %d" % id):
				offsets[instrument] = LIGOTimeGPS(offset)
			return offsets
		get_offsets = classmethod(get_offsets)

		def is_null(cls, id):
			return not cls._connection.queryOne("SELECT COUNT(offset) FROM time_slide WHERE time_slide_id == %d AND offset != 0" % id)[0]
		is_null = classmethod(is_null)

		def all_offsets(cls):
			return [cls.get_offsets(id) for (id, ) in cls._connection.queryAll("SELECT DISTINCT time_slide_id FROM time_slide")]
		all_offsets = classmethod(all_offsets)


	def time_slide_append(self, row):
		CoincDatabase.TimeSlide(
			process_id = lsctables.ILWDID(row.process_id),
			time_slide_id = lsctables.ILWDID(row.time_slide_id),
			instrument = row.instrument,
			offset = row.offset
		)
	time_slide_append = staticmethod(time_slide_append)


	class CoincDef(sqlobject.SQLObject):
		class sqlmeta:
			table = "coinc_definer"
		coinc_def_id = sqlobject.IntCol()
		table_name = sqlobject.StringCol()

		def get_table_names(cls, id):
			"""
			From a numeric ID, return a sorted list of table
			names or raise KeyError if no matching ID is found.
			"""
			l = [table_name for (table_name, ) in cls._connection.queryAll("SELECT table_name FROM coinc_definer WHERE coinc_def_id = %d" % id)]
			if not l:
				raise KeyError, id
			l.sort()
			return l
		get_table_names = classmethod(get_table_names)

		def get_id(cls, table_names):
			"""
			From a list of table names, return a numeric ID or
			raise KeyError if no matching ID is found.
			"""
			table_names = list(table_names)	# so we can modify it
			table_names.sort()
			for id in [id for (id, ) in cls._connection.queryAll("SELECT DISTINCT coinc_def_id FROM coinc_definer")]:
				if cls.get_table_names(id) == table_names:
					return id
			raise KeyError, table_names
		get_id = classmethod(get_id)


	def coinc_def_append(self, row):
		CoincDatabase.CoincDef(
			coinc_def_id = lsctables.ILWDID(row.coinc_def_id),
			table_name = row.table_name
		)
	coinc_def_append = staticmethod(coinc_def_append)


	class Coinc(sqlobject.SQLObject):
		class sqlmeta:
			table = "coinc_event"
			idName = "coinc_event_id"
		process_id = sqlobject.IntCol()
		coinc_def_id = sqlobject.IntCol()
		sngl_bursts = sqlobject.RelatedJoin("SnglBurst")
		sim_bursts = sqlobject.RelatedJoin("SimBurst")
		time_slide_id = sqlobject.IntCol()
		nevents = sqlobject.IntCol(default = 0)

		def get_time_slide(self):
			return CoincDatabase.TimeSlide.get_offsets(self.time_slide_id)

		def is_zero_lag(self):
			return CoincDatabase.TimeSlide.is_null(self.time_slide_id)


	def coinc_append(self, row):
		CoincDatabase.Coinc(
			id = lsctables.ILWDID(row.coinc_event_id),
			process_id = lsctables.ILWDID(row.process_id),
			coinc_def_id = lsctables.ILWDID(row.coinc_def_id),
			time_slide_id = lsctables.ILWDID(row.time_slide_id),
			nevents = row.nevents
		)
	coinc_append = staticmethod(coinc_append)


	def coinc_map_append(self, row):
		coinc_event_id = lsctables.ILWDID(row.coinc_event_id)
		event_table = lsctables.ILWDTableName(row.event_id)
		event_id = lsctables.ILWDID(row.event_id)
		if not table.CompareTableNames(event_table, lsctables.SnglBurstTable.tableName):
			CoincDatabase.Coinc.get(coinc_event_id).addSnglBurst(CoincDatabase.SnglBurst.get(event_id))
		elif not table.CompareTableNames(event_table, lsctables.SimBurstTable.tableName):
			CoincDatabase.Coinc.get(coinc_event_id).addSimBurst(CoincDatabase.SimBurst.get(event_id))
		elif not table.CompareTableNames(event_table, lsctables.CoincTable.tableName):
			# resolve sim/coincs to the burst events in the
			# coinc FIXME:  we're assuming here that (a) this
			# is a sim/coinc coinc, and (b) that the coinc
			# target is a burst/burst coinc;  add checks of
			# coinc_def_id's to confirm?
			map(CoincDatabase.Coinc.get(coinc_event_id).addSnglBurst, CoincDatabase.Coinc.get(event_id).sngl_bursts)
		else:
			raise TypeError, row.event_id
	coinc_map_append = staticmethod(coinc_map_append)


	#
	# Database class methods
	#

	def __init__(self, connection):
		for cls in [self.SnglBurst, self.SimBurst, self.TimeSlide, self.CoincDef, self.Coinc]:
			cls._connection = connection
			cls.createTable()
		self.TimeSlide._connection.query("CREATE INDEX time_slide_id_index ON time_slide (time_slide_id)")

		lsctables.SnglBurstTable.append = self.sngl_burst_append
		lsctables.SimBurstTable.append = self.sim_burst_append
		lsctables.TimeSlideTable.append = self.time_slide_append
		lsctables.CoincDefTable.append = self.coinc_def_append
		lsctables.CoincTable.append = self.coinc_append
		lsctables.CoincMapTable.append = self.coinc_map_append

		self.bb_definer_id = None
		self.sb_definer_id = None
		self.sc_definer_id = None
		self.missed_injections = []


	def summarize(self, xmldoc, live_time_program, verbose = False):
		"""
		Compute and record some summary information about the
		database.  Call this after all the data has been inserted,
		and before you want any of this information.
		"""
		# get the segment lists
		self.seglists = llwapp.segmentlistdict_fromsearchsummary(xmldoc, live_time_program)
		self.instruments = self.seglists.keys()

		# determine a few coinc_definer IDs
		self.bb_definer_id = self.CoincDef.get_id([lsctables.SnglBurstTable.tableName])
		self.sb_definer_id = self.CoincDef.get_id([lsctables.SnglBurstTable.tableName, lsctables.SimBurstTable.tableName])
		self.sc_definer_id = self.CoincDef.get_id([lsctables.CoincTable.tableName, lsctables.SimBurstTable.tableName])

		# compute the missed injections
		self.missed_injections = [id for (id, ) in self.SimBurst._connection.queryAll("SELECT DISTINCT sim_burst.simulation_id FROM sim_burst")]
		for coinc in self.Coinc.selectBy(coinc_def_id = self.sb_definer_id):
			for sim in coinc.sim_bursts:
				try:
					self.missed_injections.remove(sim.id)
				except ValueError:
					# already removed
					pass

		# determine burst <--> burst coincidences for which at
		# least one burst, but not all, was identified as an
		# injection;  these are places in the data where an
		# injection was done, a coincident event was seen, but
		# where the injection was not found to match all events in
		# the coincidence;  these perhaps indicate power leaking
		# from the injection into nearby tiles, or accidental
		# coincidence with near-by noise, etc, and so although they
		# aren't "bang-on" reconstructions of injections they are
		# nevertheless injections that are found and survive a
		# coincidence cut
		# FIXME: I don't know how to do this

		# verbosity
		if verbose:
			print >>sys.stderr, "database stats:"
			print >>sys.stderr, "\tburst events: %d" % self.SnglBurst.select().count()
			print >>sys.stderr, "\tinjections: %d" % self.SimBurst.select().count()
			print >>sys.stderr, "\ttime slides: %d" % len(self.TimeSlide.all_offsets())
			print >>sys.stderr, "\tburst + burst coincidences: %d" % self.Coinc.selectBy(coinc_def_id = self.bb_definer_id).count()
			print >>sys.stderr, "\tinjection + burst coincidences: %d" % self.Coinc.selectBy(coinc_def_id = self.sb_definer_id).count()
			print >>sys.stderr, "\tinjection + (burst + burst) coincidences: %d" % self.Coinc.selectBy(coinc_def_id = self.sc_definer_id).count()


	def clear(self):
		"""
		Empty the contents of the database
		"""
		for cls in [self.SnglBurst, self.SimBurst, self.TimeSlide, self.CoincDef, self.Coinc]:
			cls.clearTable()

		self.bb_definer_id = None
		self.sb_definer_id = None
		self.sc_definer_id = None
		self.missed_injections = []


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
