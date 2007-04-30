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
from glue.ligolw import ilwd
from glue.ligolw import array
from glue.ligolw import param
from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw import utils
from pylal import itertools
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
		# break order degeneracy by requiring instrument pairs to
		# be in alphabetical order
		if inst1 > inst2:
			continue

		# now retrieve the thresholds for the other order
		dt_other, df_other, dh_other = thresholds[(inst2, inst1)]

		# determine the smallest intervals enclosing both intervals
		dt = segments.segment(-dt, +dt) | segments.segment(-dt_other, +dt_other)
		df = segments.segment(-df, +df) | segments.segment(-df_other, +df_other)
		dh = segments.segment(-dh, +dh) | segments.segment(-dh_other, +dh_other)

		# record the result
		thresholds[(inst1, inst2)] = thresholds[(inst2, inst1)] = (dt, df, dh)

	#
	# Remove duplicates.  Leave only instrument pairs in which the
	# instruments are listed in alphabetical order.
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
#                                 Book-keeping
#
# =============================================================================
#


#
# How to compute coincidence parameters from burst events
#


# FIXME:  should probably introduct a CoincParams class, initialized from a
# list of events, rather than mucking with dictionaries.


def coinc_params(events, offsetdict):
	events.sort(lambda a, b: cmp(a.ifo, b.ifo))
	params = {}
	for event1, event2 in itertools.choices(events, 2):
		if event1.ifo == event2.ifo:
			# a coincidence is parameterized only by
			# inter-instrument deltas
			continue

		prefix = "%s_%s_" % (event1.ifo, event2.ifo)

		# in each of the following, if the list of events contains
		# more than one event from a given instrument, the smallest
		# deltas are recorded
		dt = float(event1.get_peak() + offsetdict[event1.ifo] - event2.get_peak() - offsetdict[event2.ifo]) / ((event1.ms_duration + event2.ms_duration) / 2)
		if prefix + "dt" not in params or abs(params[prefix + "dt"]) > abs(dt):
			params[prefix + "dt"] = dt

		df = (event1.peak_frequency - event2.peak_frequency) / ((event1.ms_bandwidth + event2.ms_bandwidth) / 2)
		if prefix + "df" not in params or abs(params[prefix + "df"]) > abs(df):
			params[prefix + "df"] = df

		dh = (event1.ms_hrss - event2.ms_hrss) / ((event1.ms_hrss + event2.ms_hrss) / 2)
		if prefix + "dh" not in params or abs(params[prefix + "dh"]) > abs(dh):
			params[prefix + "dh"] = dh

		dband = (event1.ms_bandwidth - event2.ms_bandwidth) / ((event1.ms_bandwidth + event2.ms_bandwidth) / 2)
		if prefix + "dband" not in params or abs(params[prefix + "dband"]) > abs(dband):
			params[prefix + "dband"] = dband

	return params


#
# A class for measuring 1-D parameter distributions
#


class CoincParamsDistributions(object):
	def __init__(self, rate_args):
		self.background_rates = {}
		self.injection_rates = {}
		for name, args in rate_args.iteritems():
			self.background_rates[name] = rate.Rate(*args)
			self.injection_rates[name] = rate.Rate(*args)

	def __iadd__(self, other):
		if type(other) != type(self):
			raise TypeError, other
		for name, rate in other.background_rates.iteritems():
			if name in self.background_rates:
				self.background_rates[name] += rate
			else:
				self.background_rates[name] = rate
		for name, rate in other.injection_rates.iteritems():
			if name in self.injection_rates:
				self.injection_rates[name] += rate
			else:
				self.injection_rates[name] = rate
		return self

	def add_background(self, param_func, events, offsetdict):
		for name, value in param_func(events, offsetdict).iteritems():
			rate = self.background_rates[name]
			try:
				rate[value] = 1.0
			except IndexError:
				# param value out of range
				pass

	def add_injection(self, param_func, events, offsetdict):
		for name, value in param_func(events, offsetdict).iteritems():
			rate = self.injection_rates[name]
			try:
				rate[value] = 1.0
			except IndexError:
				# param value out of range
				pass

	def set_filter(self, name, filterwidth, windowfunc):
		self.background_rates[name].set_filter(filterwidth, windowfunc)
		self.injection_rates[name].set_filter(filterwidth, windowfunc)

	def finish(self):
		for rate in self.background_rates.itervalues():
			rate.array /= numpy.sum(rate.array)
			rate.filter()
		for rate in self.injection_rates.itervalues():
			rate.array /= numpy.sum(rate.array)
			rate.filter()
		return self


#
# Scatter plot data
#


class Scatter(object):
	def __init__(self):
		self.bak_x = []
		self.bak_y = []
		self.inj_x = []
		self.inj_y = []

	def add_background(self, x, y):
		self.bak_x.append(x)
		self.bak_y.append(y)

	def add_injection(self, x, y):
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

	def add_injection(self, *args):
		self.inj_observations.append(args)

	def finish(self):
		self.bak_cov = covariance_normalize(stats.cov(self.bak_observations))
		self.inj_cov = covariance_normalize(stats.cov(self.inj_observations))


#
# =============================================================================
#
#                                  Interface
#
# =============================================================================
#


#
# Base class used to hook the database contents into a statistics analyzer.
#


class Stats(object):
	def __init__(self, thresholds):
		self.thresholds = thresholds
		self.n_time_slides = None
		self.n_background_events = 0


	def _add_background(self, param_func, events, offsetdict):
		pass


	def _add_injections(self, param_func, events, offsetdict):
		pass


	def add_background(self, database):
		# count the number of time slides (assume all input files
		# list the exact same time slides, so only do this once)
		if self.n_time_slides is None:
			self.n_time_slides = database.connection.cursor().execute("""SELECT COUNT(DISTINCT time_slide_id) FROM time_slide""").fetchone()[0]

		# iterate over non-zero-lag burst<-->burst coincs
		for (coinc_event_id,) in database.connection.cursor().execute("""
SELECT coinc_event_id FROM
	coinc_event
WHERE
	coinc_def_id == ?
	AND EXISTS (
		SELECT * FROM
			time_slide
		WHERE
			time_slide.time_slide_id == coinc_event.time_slide_id
			AND time_slide.offset != 0
	)
		""", (database.bb_definer_id,)):
			# add 1 to the count of background coincs
			self.n_background_events += 1

			# retrieve the list of the sngl_bursts in this
			# coinc, and their time slide dictionary
			events = []
			offsetdict = {}
			for values in database.connection.cursor().execute("""
SELECT sngl_burst.*, time_slide.offset FROM
	sngl_burst
	JOIN coinc_event_map ON (
		coinc_event_map.event_id == sngl_burst.event_id
		AND coinc_event_map.table_name == 'sngl_burst'
	)
	JOIN coinc_event ON (
		coinc_event.coinc_event_id == coinc_event_map.coinc_event_id
	)
	JOIN time_slide ON (
		coinc_event.time_slide_id == time_slide.time_slide_id
		AND sngl_burst.ifo == time_slide.instrument
	)
WHERE
	coinc_event.coinc_event_id == ?
			""", (coinc_event_id,)):
				# reconstruct the event
				event = database.sngl_burst_table._row_from_cols(values[:-1])

				# add to list
				events.append(event)

				# store the time slide offset
				offsetdict[event.ifo] = values[-1]

			self._add_background(coinc_params, events, offsetdict)


	def add_injections(self, database):
		# iterate over burst<-->injection coincs
		for (coinc_event_id,) in database.connection.cursor().execute("""
SELECT coinc_event_id FROM
	coinc_event
WHERE
	coinc_def_id == ?
		""", (database.sb_definer_id,)):
			# retrieve the list of the sngl_bursts in this
			# coinc, and their time slide dictionary
			events = []
			offsetdict = {}
			for values in database.connection.cursor().execute("""
SELECT sngl_burst.*, time_slide.offset FROM
	sngl_burst
	JOIN coinc_event_map ON (
		coinc_event_map.event_id == sngl_burst.event_id
		AND coinc_event_map.table_name == 'sngl_burst'
	)
	JOIN coinc_event ON (
		coinc_event.coinc_event_id == coinc_event_map.coinc_event_id
	)
	JOIN time_slide ON (
		coinc_event.time_slide_id == time_slide.time_slide_id
		AND sngl_burst.ifo == time_slide.instrument
	)
WHERE
	coinc_event.coinc_event_id == ?
			""", (coinc_event_id,)):
				# reconstruct the burst events
				event = database.sngl_burst_table._row_from_cols(values[:-1])

				# add to list
				events.append(event)

				# store the time slide offset
				offsetdict[event.ifo] = values[-1]

			# pass the events to whatever wants them
			self._add_injections(coinc_params, events, offsetdict)

	def finish(self):
		pass


#
# Scatter plot data
#


class ScatterStats(Stats):
	def __init__(self, thresholds):
		Stats.__init__(self, thresholds)
		self.scatter = Scatter()

	def _add_background(self, param_func, events, offsetdict):
		params = param_func(events, offsetdict)
		for event1, event2 in itertools.choices(events, 2):
			if event1.ifo == event2.ifo:
				continue
			prefix = "%s_%s_" % (event1.ifo, event2.ifo)
			self.scatter.add_background(params[prefix + "dt"], params[prefix + "df"])

	def _add_injections(self, param_func, events, offsetdict):
		params = param_func(events, offsetdict)
		for event1, event2 in itertools.choices(events, 2):
			if event1.ifo == event2.ifo:
				continue
			prefix = "%s_%s_" % (event1.ifo, event2.ifo)
			self.scatter.add_injection(params[prefix + "dt"], params[prefix + "df"])

	def finish(self):
		self.scatter.finish()


#
# Covariance data
#


class CovarianceStats(Stats):
	def __init__(self, thresholds):
		Stats.__init__(self, thresholds)
		self.covariance = Covariance()

	def _add_background(self, param_func, events, offsetdict):
		params = param_func(events, offsetdict)
		items = params.items()
		items.sort()
		self.covariance.add_background([value for name, value in items])

	def _add_injections(self, param_func, events, offsetdict):
		params = param_func(events, offsetdict)
		items = params.items()
		items.sort()
		self.covariance.add_injection([value for name, value in items])

	def finish(self):
		self.covariance.finish()


class DistributionsStats(Stats):
	filter_widths = {
		"H1_H2_dband": 1.0 / 25,
		"H1_L1_dband": 1.0 / 25,
		"H2_L1_dband": 1.0 / 25,
		"H1_H2_df": 1.0 / 60,
		"H1_L1_df": 1.0 / 60,
		"H2_L1_df": 1.0 / 60,
		"H1_H2_dh": 1.0 / 150,
		"H1_L1_dh": 1.0 / 9,
		"H2_L1_dh": 1.0 / 9,
		"H1_H2_dt": 1.0 / 300,
		"H1_L1_dt": 1.0 / 75,
		"H2_L1_dt": 1.0 / 75
	}

	def __init__(self, thresholds):
		Stats.__init__(self, thresholds)
		# careful, the intervals have to be unpacked in the order
		# in which they were packed by dbget_thresholds(); also
		# each instrument pair must list the instruments in
		# alphabetical order, or the parameters returned by
		# coinc_params() won't match Rate instances
		rate_args = {}
		for pair, (dtinterval, dfinterval, dhinterval) in thresholds.items():
			name = "%s_%s_dt" % pair
			rate_args[name] = (dtinterval, self.filter_widths[name])
			name = "%s_%s_df" % pair
			rate_args[name] = (dfinterval, self.filter_widths[name])
			name = "%s_%s_dh" % pair
			rate_args[name] = (dhinterval, self.filter_widths[name])
			name = "%s_%s_dband" % pair
			rate_args[name] = (segments.segment(-2.0, 2.0), self.filter_widths[name])
		self.distributions = CoincParamsDistributions(rate_args)


	def _add_background(self, param_func, events, offsetdict):
		self.distributions.add_background(param_func, events, offsetdict)

	def _add_injections(self, param_func, events, offsetdict):
		self.distributions.add_injection(param_func, events, offsetdict)

	def finish(self):
		self.distributions.finish()


#
# =============================================================================
#
#                                     I/O
#
# =============================================================================
#


def coinc_params_distributions_to_xml(coinc_params_distributions, name):
	xml = ligolw.LIGO_LW({u"Name": u"%s:pylal_ligolw_burca_tailor_coincparamsdistributions" % name})
	for name, rateobj in coinc_params_distributions.background_rates.iteritems():
		xml.appendChild(rate.rate_to_xml(rateobj, "background:%s" % name))
	for name, rateobj in coinc_params_distributions.injection_rates.iteritems():
		xml.appendChild(rate.rate_to_xml(rateobj, "injection:%s" % name))
	return xml


def coinc_params_distributions_from_xml(xml, name):
	xml, = [elem for elem in xml.getElementsByTagName(ligolw.LIGO_LW.tagName) if elem.hasAttribute(u"Name") and elem.getAttribute(u"Name") == u"%s:pylal_ligolw_burca_tailor_coincparamsdistributions" % name]
	names = [elem.getAttribute("Name").split(":")[1] for elem in xml.childNodes if elem.getAttribute("Name")[:11] == "background:"]
	c = CoincParamsDistributions({})
	for name in names:
		c.background_rates[name] = rate.rate_from_xml(xml, "background:%s" % name)
		c.injection_rates[name] = rate.rate_from_xml(xml, "injection:%s" % name)
	return c


def coinc_params_distributions_from_filenames(filenames, name, verbose = False):
	c = CoincParamsDistributions({})
	for filename in filenames:
		c += coinc_params_distributions_from_xml(utils.load_filename(filename, verbose = verbose, gz = (filename or "stdin")[-3:] == ".gz"), name)
	return c


#
# =============================================================================
#
#                             Process Information
#
# =============================================================================
#


process_program_name = "ligolw_burca_tailor"


def append_process(xmldoc, **kwargs):
	return llwapp.append_process(xmldoc, program = process_program_name, version = __version__, cvs_repository = "lscsoft", cvs_entry_time = __date__, comment = kwargs["comment"])


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


def gen_likelihood_control(coinc_params_distributions):
	xmldoc = ligolw.Document()
	node = xmldoc.appendChild(ligolw.LIGO_LW())

	node.appendChild(lsctables.New(lsctables.ProcessTable))
	node.appendChild(lsctables.New(lsctables.ProcessParamsTable))
	process = append_process(xmldoc, comment = u"")

	node.appendChild(coinc_params_distributions_to_xml(coinc_params_distributions, u"ligolw_burca_tailor"))

	return xmldoc


#
# =============================================================================
#
#                           param_dist_definer:table
#
# =============================================================================
#


class ParamDistDefinerIDs(ilwd.ILWD):
	def __init__(self, n = 0):
		ilwd.ILWD.__init__(self, "param_dist_definer", "param_dist_def_id", n)


class ParamDistDefinerTable(table.Table):
	tableName = "param_dist_definer:table"
	validcolumns = {
		"process_id": "ilwd:char",
		"param_dist_def_id": "ilwd:char",
		"search": "lstring",
		"distribution_name": "lstring",
		"start_time": "int_4s",
		"start_time_ns": "int_4s",
		"end_time": "int_4s",
		"end_time_ns": "int_4s",
		"array_name": "lstring"
	}
	constraints = "PRIMARY KEY (param_dist_def_id)"
	ids = ParamDistDefinerIDs()


class ParamDistDefiner(object):
	__slots__ = ParamDistDefinerTable.validcolumns.keys()


ParamDistDefinerTable.RowType = ParamDistDefiner
