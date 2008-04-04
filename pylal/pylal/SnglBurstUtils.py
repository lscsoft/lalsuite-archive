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
#from matplotlib.backends.backend_pdf import FigureCanvasPdf as FigureCanvas
import re
import sys
# Python 2.3 compatibility
try:
	set
except NameError:
	from sets import Set as set


from glue.ligolw import table
from glue.ligolw import lsctables
from pylal import llwapp


#
# =============================================================================
#
#                                   Database
#
# =============================================================================
#


class CoincDatabase(object):
	def __init__(self, live_time_program, search = "excesspower", verbose = False):
		"""
		Compute and record some summary information about the
		database.  Call this after all the data has been inserted,
		and before you want any of this information.
		"""

		from glue.ligolw import dbtables
		self.connection = dbtables.DBTable_get_connection()
		self.xmldoc = dbtables.DBTable_get_xml()

		cursor = self.connection.cursor()

		# find the tables
		try:
			self.sngl_burst_table = table.get_table(self.xmldoc, lsctables.SnglBurstTable.tableName)
		except ValueError:
			self.sngl_burst_table = None
		try:
			self.sim_burst_table = table.get_table(self.xmldoc, lsctables.SimBurstTable.tableName)
		except ValueError:
			self.sim_burst_table = None
		try:
			self.coinc_def_table = table.get_table(self.xmldoc, lsctables.CoincDefTable.tableName)
			self.coinc_table = table.get_table(self.xmldoc, lsctables.CoincTable.tableName)
			self.time_slide_table = table.get_table(self.xmldoc, lsctables.TimeSlideTable.tableName)
		except ValueError:
			self.coinc_def_table = None
			self.coinc_table = None
			self.time_slide_table = None
		try:
			self.multi_burst_table = table.get_table(self.xmldoc, lsctables.MultiBurstTable.tableName)
		except ValueError:
			self.multi_burst_table = None

		# get the segment lists
		self.seglists = llwapp.segmentlistdict_fromsearchsummary(self.xmldoc, live_time_program)
		self.instruments = set(self.seglists.keys())

		# determine a few coinc_definer IDs
		# FIXME:  don't hard-code the numbers
		if self.coinc_def_table is not None:
			try:
				self.bb_definer_id = self.coinc_def_table.get_coinc_def_id(search, 0, create_new = False)
			except KeyError:
				self.bb_definer_id = None
			try:
				self.sb_definer_id = self.coinc_def_table.get_coinc_def_id(search, 1, create_new = False)
			except KeyError:
				self.sb_definer_id = None
			try:
				self.sce_definer_id = self.coinc_def_table.get_coinc_def_id(search, 2, create_new = False)
			except KeyError:
				self.sce_definer_id = None
			try:
				self.scn_definer_id = self.coinc_def_table.get_coinc_def_id(search, 3, create_new = False)
			except KeyError:
				self.scn_definer_id = None
		else:
			self.bb_definer_id = None
			self.sb_definer_id = None
			self.sce_definer_id = None
			self.scn_definer_id = None

		# verbosity
		if verbose:
			print >>sys.stderr, "database stats:"
			if self.sngl_burst_table is not None:
				print >>sys.stderr, "\tburst events: %d" % len(self.sngl_burst_table)
			if self.sim_burst_table is not None:
				print >>sys.stderr, "\tinjections: %d" % len(self.sim_burst_table)
			if self.time_slide_table is not None:
				print >>sys.stderr, "\ttime slides: %d" % cursor.execute("SELECT COUNT(DISTINCT(time_slide_id)) FROM time_slide").fetchone()[0]
			if self.coinc_def_table is not None:
				for description, n in self.connection.cursor().execute("SELECT description, COUNT(*) FROM coinc_definer NATURAL JOIN coinc_event GROUP BY coinc_def_id"):
					print >>sys.stderr, "\t%s: %d" % (description, n)


def coinc_sngl_bursts(contents, coinc_event_id):
	for values in contents.connection.cursor().execute("""
SELECT sngl_burst.* FROM
	sngl_burst
	JOIN coinc_event_map ON (
		sngl_burst.event_id == coinc_event_map.event_id
		AND coinc_event_map.table_name == 'sngl_burst'
	)
WHERE
	coinc_event_map.coinc_event_id == ?
	""", (coinc_event_id,)):
		yield contents.sngl_burst_table._row_from_cols(values)


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
