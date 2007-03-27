# $Id$
#
# Copyright (C) 2007  Kipp C. Cannon
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


import numpy
from scipy.stats import stats
from xml import sax


from glue import segments
from glue.ligolw import ligolw
from glue.ligolw import array
from glue.ligolw import param
from glue.ligolw import lsctables
from pylal import ligolw_burca
from pylal import llwapp
from pylal import rate
from pylal.date import LIGOTimeGPS


__author__ = "Kipp Cannon <kipp@gravity.phys.uwm.edu>"
__version__ = "$Revision$"[11:-2]
__date__ = "$Date$"[7:-2]


#
# =============================================================================
#
#                             Thresholds Recovery
#
# =============================================================================
#


def dbget_thresholds(connection):
	"""
	Extract ligolw_burca's --thresholds arguments from the
	process_params table, and munge into the form desired by the rest
	of the code in this module.
	"""

	#
	# Retrieve the --thresholds arguments from the process_params
	# table.
	#

	thresholds = ligolw_burca.dbget_thresholds(connection)

	#
	# Convert to symmetric intervals
	#

	for (inst1, inst2), (dt, df, dh) in thresholds.items():
		if inst1 > inst2:
			continue
		dt_other, df_other, dh_other = thresholds[(inst2, inst1)]
		dt = segments.segment(-dt, +dt) | segments.segment(-dt_other, +dt_other)
		df = segments.segment(-df, +df) | segments.segment(-df_other, +df_other)
		dh = segments.segment(-dh, +dh) | segments.segment(-dh_other, +dh_other)
		thresholds[(inst1, inst2)] = thresholds[(inst2, inst1)] = (dt, df, dh)

	#
	# Remove duplicates.
	#

	for pair in thresholds.keys():
		if pair[0] > pair[1]:
			del thresholds[pair]

	#
	# Done.
	#

	return thresholds


#
# =============================================================================
#
#                                 Bookkeeping
#
# =============================================================================
#


#
# Track the distributions of the parameters associated with background and
# injection coincidences.
#


class Delta_Distributions(object):
	def __init__(self, thresholds):
		self.thresholds = thresholds

		# initialize the binnings
		self.bak_dt = {}
		self.inj_dt = {}
		self.bak_df = {}
		self.inj_df = {}
		self.bak_dh = {}
		self.inj_dh = {}
		for pair, (dtinterval, dfinterval, dhinterval) in thresholds.items():
			self.bak_dt[pair] = rate.Rate(dtinterval, abs(dtinterval) / 75.0)
			self.inj_dt[pair] = rate.Rate(dtinterval, abs(dtinterval) / 75.0)
			self.bak_df[pair] = rate.Rate(dfinterval, abs(dfinterval) / 75.0)
			self.inj_df[pair] = rate.Rate(dfinterval, abs(dfinterval) / 75.0)
			self.bak_dh[pair] = rate.Rate(dhinterval, abs(dhinterval) / 75.0)
			self.inj_dh[pair] = rate.Rate(dhinterval, abs(dhinterval) / 75.0)

	def add_background(self, pair, dt, df, dh):
		# IndexError == not within thresholds
		try:
			self.bak_dt[pair][dt] = 1.0
		except IndexError:
			pass
		try:
			self.bak_df[pair][df] = 1.0
		except IndexError:
			pass
		try:
			self.bak_dh[pair][dh] = 1.0
		except IndexError:
			pass

	def add_injections(self, pair, dt, df, dh):
		# IndexError == not within thresholds
		try:
			self.inj_dt[pair][dt] = 1.0
		except IndexError:
			pass
		try:
			self.inj_df[pair][df] = 1.0
		except IndexError:
			pass
		try:
			self.inj_dh[pair][dh] = 1.0
		except IndexError:
			pass

	def normalize(self):
		# normalize the distributions
		for pair in self.thresholds.keys():
			self.bak_dt[pair].array /= numpy.sum(self.bak_dt[pair].array)
			self.inj_dt[pair].array /= numpy.sum(self.inj_dt[pair].array)
			self.bak_df[pair].array /= numpy.sum(self.bak_df[pair].array)
			self.inj_df[pair].array /= numpy.sum(self.inj_df[pair].array)
			self.bak_dh[pair].array /= numpy.sum(self.bak_dh[pair].array)
			self.inj_dh[pair].array /= numpy.sum(self.inj_dh[pair].array)


#
# Scatter plot data
#


class Scatter(object):
	def __init__(self):
		self.bak_x = []
		self.bak_y = []
		self.inj_x = []
		self.inj_y = []

	def add_background(self, pair, x, y):
		self.bak_x.append(x)
		self.bak_y.append(y)

	def add_injections(self, pair, x, y):
		self.inj_x.append(x)
		self.inj_y.append(y)

	def finish(self):
		pass


#
# Covariance matrix
#


def covariance_normalize(c):
	"""
	Normalize a covariance matrix so that the variances (diagonal
	elements) are 1.
	"""
	std_dev = numpy.sqrt(numpy.diagonal(c))
	return c / numpy.outer(std_dev, std_dev)


class Covariance(object):
	def __init__(self):
		self.bak_observations = []
		self.inj_observations = []

	def add_background(self, *args):
		self.bak_observations.append(args)

	def add_injections(self, *args):
		self.inj_observations.append(args)

	def finish(self):
		self.bak_observations = numpy.array(self.bak_observations)
		self.inj_observations = numpy.array(self.inj_observations)
		self.bak_cov = covariance_normalize(stats.cov(self.bak_observations))
		self.inj_cov = covariance_normalize(stats.cov(self.inj_observations))


#
# =============================================================================
#
#                                  Interface
#
# =============================================================================
#


class Stats(object):
	def __init__(self, thresholds):
		self.thresholds = thresholds
		self.n_time_slides = None
		self.n_background_events = 0

		self.deltas = Delta_Distributions(thresholds)
		self.scatter = Scatter()
		self.covariance = Covariance()


	def add_background(self, database):
		# count the number of time slides (assume all input files
		# list the exact same time slides)
		if self.n_time_slides is None:
			self.n_time_slides = database.connection.cursor().execute("""SELECT COUNT(DISTINCT time_slide_id) FROM time_slide""").fetchone()[0]

		# iterate over non-zero-lag burst+burst coincidences
		# involving the two desired instruments
		for pair in self.thresholds.keys():
			for b1_confidence, b1_peak_time, b1_peak_time_ns, b1_duration, b1_peak_frequency, b1_bandwidth, b1_hrss, b2_confidence, b2_peak_time, b2_peak_time_ns, b2_duration, b2_peak_frequency, b2_bandwidth, b2_hrss in database.connection.cursor().execute("""
SELECT b1.confidence, b1.peak_time + t1.offset, b1.peak_time_ns, b1.ms_duration, b1.peak_frequency, b1.ms_bandwidth, b1.ms_hrss, b2.confidence, b2.peak_time + t2.offset, b2.peak_time_ns, b2.ms_duration, b2.peak_frequency, b2.ms_bandwidth, b2.ms_hrss FROM
	sngl_burst AS b1
	JOIN coinc_event_map AS a ON (
		a.event_id == b1.event_id
		AND a.table_name == 'sngl_burst'
	)
	JOIN coinc_event_map AS b ON (
		b.coinc_event_id == a.coinc_event_id
	)
	JOIN sngl_burst AS b2 ON (
		b.event_id == b2.event_id
		AND b.table_name == 'sngl_burst'
	)
	JOIN coinc_event ON (
		coinc_event.coinc_event_id == a.coinc_event_id
	)
	JOIN time_slide AS t1 ON (
		coinc_event.time_slide_id == t1.time_slide_id
		AND b1.ifo == t1.instrument
	)
	JOIN time_slide AS t2 ON (
		coinc_event.time_slide_id == t2.time_slide_id
		AND b2.ifo == t2.instrument
	)
WHERE
	coinc_event.coinc_def_id == ?
	AND EXISTS (
		SELECT * FROM
			time_slide
		WHERE
			time_slide.time_slide_id == coinc_event.time_slide_id
			AND time_slide.offset != 0
	)
	AND b1.ifo == ?
	AND b2.ifo == ?
			""", (database.bb_definer_id, pair[0], pair[1])):
				self.n_background_events += 1

				dt = float(LIGOTimeGPS(b1_peak_time, b1_peak_time_ns) - LIGOTimeGPS(b2_peak_time, b2_peak_time_ns)) / ((b1_duration + b2_duration) / 2)
				df = (b1_peak_frequency - b2_peak_frequency) / ((b1_bandwidth + b2_bandwidth) / 2)
				dh = (b1_hrss - b2_hrss) / ((b1_hrss + b2_hrss) / 2)

				self.deltas.add_background(pair, dt, df, dh)
				self.scatter.add_background(pair, dt, df)
				self.covariance.add_background(dt, df, dh)


	def add_injections(self, database):
		# iterate over injections recovered in both of the two
		# desired instruments
		for pair in self.thresholds.keys():
			for b1_confidence, b1_peak_time, b1_peak_time_ns, b1_duration, b1_peak_frequency, b1_bandwidth, b1_hrss, b2_confidence, b2_peak_time, b2_peak_time_ns, b2_duration, b2_peak_frequency, b2_bandwidth, b2_hrss in database.connection.cursor().execute("""
SELECT b1.confidence, b1.peak_time + t1.offset, b1.peak_time_ns, b1.ms_duration, b1.peak_frequency, b1.ms_bandwidth, b1.ms_hrss, b2.confidence, b2.peak_time + t2.offset, b2.peak_time_ns, b2.ms_duration, b2.peak_frequency, b2.ms_bandwidth, b2.ms_hrss FROM
	sngl_burst AS b1
	JOIN coinc_event_map AS a ON (
		a.event_id == b1.event_id
		AND a.table_name == 'sngl_burst'
	)
	JOIN coinc_event_map AS b ON (
		b.coinc_event_id == a.coinc_event_id
	)
	JOIN sngl_burst AS b2 ON (
		b.event_id == b2.event_id
		AND b.table_name == 'sngl_burst'
	)
	JOIN coinc_event ON (
		coinc_event.coinc_event_id == a.coinc_event_id
	)
	JOIN time_slide AS t1 ON (
		coinc_event.time_slide_id == t1.time_slide_id
		AND b1.ifo == t1.instrument
	)
	JOIN time_slide AS t2 ON (
		coinc_event.time_slide_id == t2.time_slide_id
		AND b2.ifo == t2.instrument
	)
WHERE
	coinc_event.coinc_def_id == ?
	AND b1.ifo == ?
	AND b2.ifo == ?
			""", (database.sb_definer_id, pair[0], pair[1])):
				dt = float(LIGOTimeGPS(b1_peak_time, b1_peak_time_ns) - LIGOTimeGPS(b2_peak_time, b2_peak_time_ns)) / ((b1_duration + b2_duration) / 2)
				df = (b1_peak_frequency - b2_peak_frequency) / ((b1_bandwidth + b2_bandwidth) / 2)
				dh = (b1_hrss - b2_hrss) / ((b1_hrss + b2_hrss) / 2)

				self.deltas.add_injections(pair, dt, df, dh)
				self.scatter.add_injections(pair, dt, df)
				self.covariance.add_injections(dt, df, dh)

	def finish(self):
		self.deltas.normalize()
		self.scatter.finish()
		self.covariance.finish()


#
# =============================================================================
#
#                             Process Information
#
# =============================================================================
#


process_program_name = "ligolw_burca_tailor"


def append_process(xmldoc, **kwargs):
	process = llwapp.append_process(xmldoc, program = process_program_name, version = __version__, cvs_repository = "lscsoft", cvs_entry_time = __date__, comment = kwargs["comment"])

	params = []
	for (inst1, inst2), value in kwargs["thresholds"].iteritems():
		params += [(u"--thresholds", u"lstring", u"%s,%s=%s" % (inst1, inst2, ",".join(map(unicode, value))))]
	llwapp.append_process_params(xmldoc, process, params)

	return process


#
# =============================================================================
#
#                           Likelihood Control File
#
# =============================================================================
#


#
# Construct LIGO Light Weight likelihood distributions document.
#


def gen_likelihood_control(deltas):
	xmldoc = ligolw.Document()
	node = xmldoc.appendChild(ligolw.LIGO_LW())

	node.appendChild(lsctables.New(lsctables.ProcessTable))
	node.appendChild(lsctables.New(lsctables.ProcessParamsTable))
	process = append_process(xmldoc, thresholds = deltas.thresholds, comment = u"")

	node = node.appendChild(ligolw.LIGO_LW(sax.xmlreader.AttributesImpl({u"Name": process_program_name})))
	node.appendChild(param.new_param(u"process_id", u"ilwd:char", process.process_id))

	node.appendChild(llwapp.pickle_to_param(deltas.thresholds, u"thresholds"))
	for pair in deltas.thresholds.keys():
		node.appendChild(array.from_array(u"%s_%s_dt" % pair, numpy.array([deltas.inj_dt[pair].xvals(), deltas.inj_dt[pair].array, deltas.bak_dt[pair].array]), (u"dt", u"dt,P_inj,P_bak")))
		node.appendChild(array.from_array(u"%s_%s_df" % pair, numpy.array([deltas.inj_df[pair].xvals(), deltas.inj_df[pair].array, deltas.bak_df[pair].array]), (u"df", u"df,P_inj,P_bak")))
		node.appendChild(array.from_array(u"%s_%s_dh" % pair, numpy.array([deltas.inj_dh[pair].xvals(), deltas.inj_dh[pair].array, deltas.bak_dh[pair].array]), (u"dh", u"dh,P_inj,P_bak")))

	return xmldoc

