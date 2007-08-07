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
from pylal import inject
from pylal import itertools
from pylal import ligolw_burca
from pylal import llwapp
from pylal import rate
from pylal import SimBurstUtils
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


def coinc_params(events, offsetdict):
	params = {}
	events.sort(lambda a, b: cmp(a.ifo, b.ifo))
	for event1, event2 in itertools.choices(events, 2):
		if event1.ifo == event2.ifo:
			# a coincidence is parameterized only by
			# inter-instrument deltas
			continue

		prefix = "%s_%s_" % (event1.ifo, event2.ifo)

		# in each of the following, if the list of events contains
		# more than one event from a given instrument, the smallest
		# deltas are recorded

		dt = float(event1.get_peak() + offsetdict[event1.ifo] - event2.get_peak() - offsetdict[event2.ifo])
		name = prefix + "dt"
		if name not in params or abs(params[name]) > abs(dt):
			params[name] = dt

		df = (event1.peak_frequency - event2.peak_frequency) / ((event1.peak_frequency + event2.peak_frequency) / 2)
		name = prefix + "df"
		if name not in params or abs(params[name]) > abs(df):
			params[name] = df

		dh = (event1.ms_hrss - event2.ms_hrss) / ((event1.ms_hrss + event2.ms_hrss) / 2)
		name = prefix + "dh"
		if name not in params or abs(params[name]) > abs(dh):
			params[name] = dh

		dband = (event1.ms_bandwidth - event2.ms_bandwidth) / ((event1.ms_bandwidth + event2.ms_bandwidth) / 2)
		name = prefix + "dband"
		if name not in params or abs(params[name]) > abs(dband):
			params[name] = dband

		ddur = (event1.ms_duration - event2.ms_duration) / ((event1.ms_duration + event2.ms_duration) / 2)
		name = prefix + "ddur"
		if name not in params or abs(params[name]) > abs(ddur):
			params[name] = ddur

	return params


#
# A class for measuring 1-D parameter distributions
#


class CoincParamsDistributions(object):
	def __init__(self, **kwargs):
		self.background_rates = {}
		self.injection_rates = {}
		for param, rateargs in kwargs.iteritems():
			self.background_rates[param] = rate.Rate(*rateargs)
			self.injection_rates[param] = rate.Rate(*rateargs)

	def __iadd__(self, other):
		if type(other) != type(self):
			raise TypeError, other
		for param, rate in other.background_rates.iteritems():
			if param in self.background_rates:
				self.background_rates[param] += rate
			else:
				self.background_rates[param] = rate
		for param, rate in other.injection_rates.iteritems():
			if param in self.injection_rates:
				self.injection_rates[param] += rate
			else:
				self.injection_rates[param] = rate
		return self

	def add_background(self, param_func, events, timeslide):
		for param, value in param_func(events, timeslide).iteritems():
			rate = self.background_rates[param]
			try:
				rate[value] += 1.0
			except IndexError:
				# param value out of range
				pass

	def add_injection(self, param_func, events, timeslide):
		for param, value in param_func(events, timeslide).iteritems():
			rate = self.injection_rates[param]
			try:
				rate[value] += 1.0
			except IndexError:
				# param value out of range
				pass

	def set_filter(self, param, filterwidth, windowfunc):
		self.background_rates[param].set_filter(filterwidth, windowfunc)
		self.injection_rates[param].set_filter(filterwidth, windowfunc)

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
#                               Injection Filter
#
# =============================================================================
#


def good_injection_matches(sim, events, max_hrss_ratio, max_frequency_ratio):
	"""
	Return a list of the sngl_burst events from the events list that
	are "good" matches for the sim_burst.  binjfind will any old thing
	that happens to be at the same time and frequency of an injection,
	but we want to teach the Bayesian coincidence filter what a "good"
	found injection looks like.  So we remove poor reconstructions from
	the event list to show it what we're really looking for.
	"""
	return [sngl_burst for sngl_burst in events if
		(1.0 / max_hrss_ratio <= sngl_burst.ms_hrss / SimBurstUtils.hrss_in_instrument(sim, sngl_burst.ifo) <= max_hrss_ratio)
		and
		(1.0 / max_frequency_ratio <= sngl_burst.peak_frequency / sim.freq <= max_frequency_ratio)
	]


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


	def _add_background(self, param_func, events, timeslide):
		pass


	def _add_injections(self, param_func, sim, events, timeslide):
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
		coinc_event_map.table_name == 'sngl_burst'
		AND coinc_event_map.event_id == sngl_burst.event_id
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
ORDER BY
	sngl_burst.ifo
			""", (coinc_event_id,)):
				# reconstruct the event
				event = database.sngl_burst_table._row_from_cols(values)

				# add to list
				events.append(event)

				# store the time slide offset
				offsetdict[event.ifo] = values[-1]

			self._add_background(coinc_params, events, offsetdict)


	def add_injections(self, database):
		# iterate over burst<-->burst coincs in which at least one
		# burst was identified as being the result of an injection
		for values in database.connection.cursor().execute("""
SELECT
	sim_burst.*,
	burst_coinc_event.coinc_event_id
FROM
	sim_burst
	JOIN coinc_event AS burst_coinc_event
	JOIN coinc_event AS sim_coinc_event ON (
		sim_coinc_event.time_slide_id == burst_coinc_event.time_slide_id
	)
	JOIN coinc_event_map AS c ON (
		-- Each injection can result in at most 1 associated entry
		-- in the coinc_event table, so doing this join in the
		-- outer select will not result in any injections being
		-- counted more than once (and it hugely increases the
		-- speed of the query).
		c.coinc_event_id == sim_coinc_event.coinc_event_id
		AND c.table_name == 'sim_burst'
		AND c.event_id == sim_burst.simulation_id
	)
WHERE
	burst_coinc_event.coinc_def_id == ?
	AND sim_coinc_event.coinc_def_id == ?
	AND EXISTS (
		-- Find a link from the sim coinc to the burst coinc
		-- through the coinc_event_map table
		SELECT
			*
		FROM
			coinc_event_map AS a
			JOIN coinc_event_map AS b ON (
				a.table_name == 'sngl_burst'
				AND b.table_name == 'sngl_burst'
				AND b.event_id == a.event_id
			)
		WHERE
			a.coinc_event_id == burst_coinc_event.coinc_event_id
			AND b.coinc_event_id == sim_coinc_event.coinc_event_id
	)
		""", (database.bb_definer_id, database.sb_definer_id)):
			# retrieve the injection and the coinc_event_id
			sim = database.sim_burst_table._row_from_cols(values)
			coinc_event_id = values[-1]

			# retrieve the list of the sngl_bursts in this
			# coinc, and their time slide dictionary
			events = []
			offsetdict = {}
			for values in database.connection.cursor().execute("""
SELECT sngl_burst.*, time_slide.offset FROM
	sngl_burst
	JOIN coinc_event_map ON (
		coinc_event_map.table_name == 'sngl_burst'
		AND coinc_event_map.event_id == sngl_burst.event_id
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
ORDER BY
	sngl_burst.ifo
			""", (coinc_event_id,)):
				# reconstruct the burst events
				event = database.sngl_burst_table._row_from_cols(values)

				# add to list
				events.append(event)

				# store the time slide offset
				offsetdict[event.ifo] = values[-1]

			# pass the events to whatever wants them
			self._add_injections(coinc_params, sim, events, offsetdict)

	def finish(self):
		pass


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

	def _add_injections(self, param_func, sim, events, offsetdict):
		params = param_func(events, offsetdict)
		items = params.items()
		items.sort()
		self.covariance.add_injection([value for name, value in items])

	def finish(self):
		self.covariance.finish()


class DistributionsStats(Stats):
	filter_widths = {
		"H1_H2_dband": 1.0 / 400,
		"H1_L1_dband": 1.0 / 400,
		"H2_L1_dband": 1.0 / 400,
		"H1_H2_ddur": 1.0 / 400,
		"H1_L1_ddur": 1.0 / 400,
		"H2_L1_ddur": 1.0 / 400,
		"H1_H2_df": 1.0 / 400,
		"H1_L1_df": 1.0 / 400,
		"H2_L1_df": 1.0 / 400,
		"H1_H2_dh": 1.0 / 200,
		"H1_L1_dh": 1.0 / 200,
		"H2_L1_dh": 1.0 / 200,
		"H1_H2_dt": 1.0 / 7000,
		"H1_L1_dt": 1.0 / 7000,
		"H2_L1_dt": 1.0 / 7000
	}

	def __init__(self, max_hrss_ratio, max_frequency_ratio, thresholds):
		Stats.__init__(self, thresholds)
		self.max_hrss_ratio = max_hrss_ratio
		self.max_frequency_ratio = max_frequency_ratio
		# careful, the intervals have to be unpacked in the order
		# in which they were packed by dbget_thresholds(); also
		# each instrument pair must list the instruments in
		# alphabetical order, or the parameters returned by
		# coinc_params() won't match Rate instances
		rate_args = {}
		for pair, (dtinterval, dfinterval, dhinterval) in thresholds.items():
			name = "%s_%s_dt" % pair
			rate_args[name] = (dtinterval.protract(inject.light_travel_time(*pair)), self.filter_widths[name])
			name = "%s_%s_df" % pair
			rate_args[name] = (dfinterval, self.filter_widths[name])
			name = "%s_%s_dh" % pair
			rate_args[name] = (dhinterval, self.filter_widths[name])
			name = "%s_%s_dband" % pair
			rate_args[name] = (segments.segment(-2.0, +2.0), self.filter_widths[name])
			name = "%s_%s_ddur" % pair
			rate_args[name] = (segments.segment(-2.0, +2.0), self.filter_widths[name])
		self.distributions = CoincParamsDistributions(**rate_args)

	def _add_background(self, param_func, events, offsetdict):
		self.distributions.add_background(param_func, events, offsetdict)

	def _add_injections(self, param_func, sim, events, offsetdict):
		# remove events whose h_rss differs from the correct value
		# by more than a factor of max_hrss_ratio, and whose
		# frequency differs from the correct value by more than a
		# factor of max_frequency_ratio.
		events = good_injection_matches(sim, events, self.max_hrss_ratio, self.max_frequency_ratio)
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
	c = CoincParamsDistributions()
	for name in names:
		c.background_rates[name] = rate.rate_from_xml(xml, "background:%s" % name)
		c.injection_rates[name] = rate.rate_from_xml(xml, "injection:%s" % name)
	return c


def coinc_params_distributions_from_filenames(filenames, name, verbose = False):
	c = CoincParamsDistributions()
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
