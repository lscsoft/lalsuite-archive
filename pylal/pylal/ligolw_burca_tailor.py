# Copyright (C) 2007-2010  Kipp Cannon
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


import itertools
import math
import numpy
from scipy.stats import stats
import sys
import threading


from glue import iterutils
try:
	any, all
except NameError:
	# Python <2.5
	from glue.iterutils import any, all
from glue.ligolw import ligolw
from glue.ligolw import ilwd
from glue.ligolw import param
from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw import utils
from pylal import date
from pylal import git_version
from pylal import inject
from pylal import llwapp
from pylal import rate
from pylal.xlal.datatypes.ligotimegps import LIGOTimeGPS


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date


#
# =============================================================================
#
#                                 Speed Hacks
#
# =============================================================================
#


lsctables.LIGOTimeGPS = LIGOTimeGPS


#
# =============================================================================
#
#             Generating Coincidence Parameters from Burst Events
#
# =============================================================================
#


#
# All sky version.
#


def coinc_params(events, offsetvector):
	params = {}

	if events:
		# the "time" is the ms_snr squared weighted average of the
		# peak times neglecting light-travel times.  because
		# LIGOTimeGPS objects have overflow problems in this sort
		# of a calculation, the first event's peak time is used as
		# an epoch and the calculations are done w.r.t. that time.

		# FIXME: this time is available as the peak_time in the
		# multi_burst table, and it should be retrieved from that
		# table instead of being recomputed

		t = events[0].get_peak()
		t += sum(float(event.get_peak() - t) * event.ms_snr**2.0 for event in events) / sum(event.ms_snr**2.0 for event in events)
		gmst = date.XLALGreenwichMeanSiderealTime(t) % (2 * math.pi)

	for event1, event2 in iterutils.choices(sorted(events, lambda a, b: cmp(a.ifo, b.ifo)), 2):
		if event1.ifo == event2.ifo:
			# a coincidence is parameterized only by
			# inter-instrument deltas
			continue

		prefix = "%s_%s_" % (event1.ifo, event2.ifo)

		# in each of the following, if the list of events contains
		# more than one event from a given instrument, the smallest
		# deltas are recorded

		dt = float(event1.get_peak() + offsetvector[event1.ifo] - event2.get_peak() - offsetvector[event2.ifo])
		name = "%sdt" % prefix
		if name not in params or abs(params[name][0]) > abs(dt):
			#params[name] = (dt,)
			params[name] = (dt, gmst)

		df = (event1.peak_frequency - event2.peak_frequency) / ((event1.peak_frequency + event2.peak_frequency) / 2)
		name = "%sdf" % prefix
		if name not in params or abs(params[name][0]) > abs(df):
			#params[name] = (df,)
			params[name] = (df, gmst)

		dh = (event1.ms_hrss - event2.ms_hrss) / ((event1.ms_hrss + event2.ms_hrss) / 2)
		name = "%sdh" % prefix
		if name not in params or abs(params[name][0]) > abs(dh):
			#params[name] = (dh,)
			params[name] = (dh, gmst)

		dband = (event1.ms_bandwidth - event2.ms_bandwidth) / ((event1.ms_bandwidth + event2.ms_bandwidth) / 2)
		name = "%sdband" % prefix
		if name not in params or abs(params[name][0]) > abs(dband):
			#params[name] = (dband,)
			params[name] = (dband, gmst)

		ddur = (event1.ms_duration - event2.ms_duration) / ((event1.ms_duration + event2.ms_duration) / 2)
		name = "%sddur" % prefix
		if name not in params or abs(params[name][0]) > abs(ddur):
			#params[name] = (ddur,)
			params[name] = (ddur, gmst)

	return params


#
# Galactic core version.
#


def delay_and_amplitude_correct(event, ra, dec):
	# retrieve station metadata

	detector = inject.cached_detector[inject.prefix_to_name[event.ifo]]

	# delay-correct the event to the geocentre

	peak = event.get_peak()
	delay = date.XLALTimeDelayFromEarthCenter(detector.location, ra, dec, peak)
	event.set_peak(peak - delay)
	event.set_start(event.get_start() - delay)
	event.set_ms_start(event.get_ms_start() - delay)

	# amplitude-correct the event using the polarization-averaged
	# antenna response

	fp, fc = inject.XLALComputeDetAMResponse(detector.response, ra, dec, 0, date.XLALGreenwichMeanSiderealTime(peak))
	mean_response = math.sqrt(fp**2 + fc**2)
	event.amplitude /= mean_response
	event.ms_hrss /= mean_response

	# done

	return event


def targeted_coinc_params(events, offsetvector, ra, dec):
	return coinc_params((delay_and_amplitude_correct(event, ra, dec) for event in events), offsetvector)


#
# =============================================================================
#
#                                 Book-keeping
#
# =============================================================================
#


#
# A class for measuring parameter distributions
#


class FilterThread(threading.Thread):
	# allow at most 5 threads
	cpu = threading.Semaphore(5)
	# allow at most one to write to stderr
	stderr = threading.Semaphore(1)

	def __init__(self, binnedarray, filter, verbose = False, name = None):
		threading.Thread.__init__(self, name = name)
		self.binnedarray = binnedarray
		self.filter = filter
		self.verbose = verbose

	def run(self):
		self.cpu.acquire()
		if self.verbose:
			self.stderr.acquire()
			print >>sys.stderr, "\tstarting %s" % self.getName()
			self.stderr.release()

		self.binnedarray.array /= numpy.sum(self.binnedarray.array)
		rate.to_moving_mean_density(self.binnedarray, self.filter)

		if self.verbose:
			self.stderr.acquire()
			print >>sys.stderr, "\tcompleted %s" % self.getName()
			self.stderr.release()
		self.cpu.release()


class CoincParamsDistributions(object):
	def __init__(self, **kwargs):
		self.zero_lag_rates = {}
		self.background_rates = {}
		self.injection_rates = {}
		for param, binning in kwargs.items():
			self.zero_lag_rates[param] = rate.BinnedArray(binning)
			self.background_rates[param] = rate.BinnedArray(binning)
			self.injection_rates[param] = rate.BinnedArray(binning)

	def __iadd__(self, other):
		if type(other) != type(self):
			raise TypeError, other
		for param, rate in other.zero_lag_rates.items():
			if param in self.zero_lag_rates:
				self.zero_lag_rates[param] += rate
			else:
				self.zero_lag_rates[param] = rate
		for param, rate in other.background_rates.items():
			if param in self.background_rates:
				self.background_rates[param] += rate
			else:
				self.background_rates[param] = rate
		for param, rate in other.injection_rates.items():
			if param in self.injection_rates:
				self.injection_rates[param] += rate
			else:
				self.injection_rates[param] = rate
		return self

	@classmethod
	def copy(cls, other):
		new = cls(**dict((param, other.zero_lag_rates[param].bins) for param in other.zero_lag_rates))
		new += other
		return new

	def add_zero_lag(self, param_dict, weight = 1.0):
		for param, value in (param_dict or {}).items():
			rate = self.zero_lag_rates[param]
			try:
				rate[value] += weight
			except IndexError:
				# param value out of range
				pass

	def add_background(self, param_dict, weight = 1.0):
		for param, value in (param_dict or {}).items():
			rate = self.background_rates[param]
			try:
				rate[value] += weight
			except IndexError:
				# param value out of range
				pass

	def add_injection(self, param_dict, weight = 1.0):
		for param, value in (param_dict or {}).items():
			rate = self.injection_rates[param]
			try:
				rate[value] += weight
			except IndexError:
				# param value out of range
				pass

	def finish(self, filters = {}, verbose = False):
		default_filter = rate.gaussian_window(21)
		# normalizing each array so that its sum is 1 has the
		# effect of making the integral of P(x) dx equal 1 after
		# the array is transformed to an array of densities (which
		# is done by dividing each bin by dx).
		N = len(self.zero_lag_rates) + len(self.background_rates) + len(self.injection_rates)
		n = 0
		threads = []
		for group, (name, binnedarray) in itertools.chain(zip(["zero lag"] * len(self.zero_lag_rates), self.zero_lag_rates.items()), zip(["background"] * len(self.background_rates), self.background_rates.items()), zip(["injections"] * len(self.injection_rates), self.injection_rates.items())):
			n += 1
			threads.append(FilterThread(binnedarray, filters.get(name, default_filter), verbose = verbose, name = "%d / %d: %s \"%s\"" % (n, N, group, name)))
			threads[-1].start()
		for thread in threads:
			thread.join()
		return self

	@classmethod
	def from_xml(cls, xml, name):
		xml, = [elem for elem in xml.getElementsByTagName(ligolw.LIGO_LW.tagName) if elem.hasAttribute(u"Name") and elem.getAttribute(u"Name") == u"%s:pylal_ligolw_burca_tailor_coincparamsdistributions" % name]
		process_id = param.get_pyvalue(xml, u"process_id")
		names = [elem.getAttribute("Name").split(":")[1] for elem in xml.childNodes if elem.getAttribute("Name").startswith("background:")]
		self = cls()
		for name in names:
			self.zero_lag_rates[str(name)] = rate.binned_array_from_xml(xml, "zero_lag:%s" % name)
			self.background_rates[str(name)] = rate.binned_array_from_xml(xml, "background:%s" % name)
			self.injection_rates[str(name)] = rate.binned_array_from_xml(xml, "injection:%s" % name)
		return self, process_id

	def to_xml(self, process, name):
		xml = ligolw.LIGO_LW({u"Name": u"%s:pylal_ligolw_burca_tailor_coincparamsdistributions" % name})
		xml.appendChild(param.new_param(u"process_id", u"ilwd:char", process.process_id))
		for name, binnedarray in self.zero_lag_rates.items():
			xml.appendChild(rate.binned_array_to_xml(binnedarray, u"zero_lag:%s" % name))
		for name, binnedarray in self.background_rates.items():
			xml.appendChild(rate.binned_array_to_xml(binnedarray, u"background:%s" % name))
		for name, binnedarray in self.injection_rates.items():
			xml.appendChild(rate.binned_array_to_xml(binnedarray, u"injection:%s" % name))
		return xml


#
# =============================================================================
#
#                                  Interface
#
# =============================================================================
#


def get_noninjections(contents):
	"""
	Generator function to return

		is_background, event_list, offsetvector

	tuples by querying the coinc_event and sngl_burst tables in the
	database described by contents.  Only coincs corresponding to
	sngl_burst<-->sngl_burst coincs will be retrieved.
	"""
	cursor = contents.connection.cursor()
	for coinc_event_id, time_slide_id in contents.connection.cursor().execute("""
SELECT
	coinc_event_id,
	time_slide_id
FROM
	coinc_event
WHERE
	coinc_def_id == ?
	""", (contents.bb_definer_id,)):
		rows = [(contents.sngl_burst_table.row_from_cols(row), row[-1]) for row in cursor.execute("""
SELECT
	sngl_burst.*,
	time_slide.offset
FROM
	coinc_event_map
	JOIN sngl_burst ON (
		coinc_event_map.table_name == 'sngl_burst'
		AND sngl_burst.event_id == coinc_event_map.event_id
	)
	JOIN time_slide ON (
		time_slide.instrument == sngl_burst.ifo
	)
WHERE
	coinc_event_map.coinc_event_id == ?
	AND time_slide.time_slide_id == ?
		""", (coinc_event_id, time_slide_id))]
		offsetvector = dict((event.ifo, offset) for event, offset in rows)
		if any(offsetvector.values()):
			yield True, [event for event, offset in rows], offsetvector
		else:
			yield False, [event for event, offset in rows], offsetvector
	cursor.close()


def get_injections(contents):
	"""
	Generator function to return

		sim, event_list, offsetvector

	tuples by querying the sim_burst, coinc_event and sngl_burst tables
	in the database described by contents.  Only coincs corresponding
	to "exact" sim_burst<-->coinc_event coincs will be retrieved.
	"""
	cursor = contents.connection.cursor()
	for values in contents.connection.cursor().execute("""
SELECT
	sim_burst.*,
	burst_coinc_event_map.event_id
FROM
	sim_burst
	JOIN coinc_event_map AS sim_coinc_event_map ON (
		sim_coinc_event_map.table_name == 'sim_burst'
		AND sim_coinc_event_map.event_id == sim_burst.simulation_id
	)
	JOIN coinc_event AS sim_coinc_event ON (
		sim_coinc_event.coinc_event_id == sim_coinc_event_map.coinc_event_id
	)
	JOIN coinc_event_map AS burst_coinc_event_map ON (
		burst_coinc_event_map.coinc_event_id == sim_coinc_event_map.coinc_event_id
		AND burst_coinc_event_map.table_name == 'coinc_event'
	)
WHERE
	sim_coinc_event.coinc_def_id == ?
	""", (contents.sce_definer_id,)):
		# retrieve the injection and the coinc_event_id
		sim = contents.sim_burst_table.row_from_cols(values)
		coinc_event_id = values[-1]

		# retrieve the list of the sngl_bursts in this
		# coinc, and their time slide dictionary
		rows = [(contents.sngl_burst_table.row_from_cols(row), row[-1]) for row in cursor.execute("""
SELECT
	sngl_burst.*,
	time_slide.offset
FROM
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
		AND time_slide.instrument == sngl_burst.ifo
	)
WHERE
	coinc_event.coinc_event_id == ?
		""", (coinc_event_id,))]
		# pass the events to whatever wants them
		yield sim, [event for event, offset in rows], dict((event.ifo, offset) for event, offset in rows)
	cursor.close()


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

	def add_noninjections(self, param_func, database):
		# iterate over burst<-->burst coincs
		for is_background, events, offsetvector in get_noninjections(database):
			if is_background:
				self.bak_observations.append(tuple(value for name, value in sorted(param_func(events, offsetvector).items())))
			else:
				# zero-lag not used
				pass

	def add_injections(self, param_func, database):
		# iterate over burst<-->burst coincs matching injections
		# "exactly"
		for sim, events, offsetvector in get_injections(database):
			self.inj_observations.append(tuple(value for name, value in sorted(param_func(events, offsetvector).items())))

	def finish(self):
		self.bak_cov = covariance_normalize(stats.cov(self.bak_observations))
		self.inj_cov = covariance_normalize(stats.cov(self.inj_observations))


#
# Parameter distributions
#


def dt_binning(instrument1, instrument2):
	# FIXME:  hard-coded for directional search
	#dt = 0.02 + inject.light_travel_time(instrument1, instrument2)
	dt = 0.02
	return rate.NDBins((rate.ATanBins(-dt, +dt, 12001), rate.LinearBins(0.0, 2 * math.pi, 61)))


class DistributionsStats(object):
	"""
	A class used to populate a CoincParamsDistribution instance with
	the data from the outputs of ligolw_burca and ligolw_binjfind.
	"""

	binnings = {
		"H1_H2_dband": rate.NDBins((rate.LinearBins(-2.0, +2.0, 12001), rate.LinearBins(0.0, 2 * math.pi, 61))),
		"H1_L1_dband": rate.NDBins((rate.LinearBins(-2.0, +2.0, 12001), rate.LinearBins(0.0, 2 * math.pi, 61))),
		"H2_L1_dband": rate.NDBins((rate.LinearBins(-2.0, +2.0, 12001), rate.LinearBins(0.0, 2 * math.pi, 61))),
		"H1_H2_ddur": rate.NDBins((rate.LinearBins(-2.0, +2.0, 12001), rate.LinearBins(0.0, 2 * math.pi, 61))),
		"H1_L1_ddur": rate.NDBins((rate.LinearBins(-2.0, +2.0, 12001), rate.LinearBins(0.0, 2 * math.pi, 61))),
		"H2_L1_ddur": rate.NDBins((rate.LinearBins(-2.0, +2.0, 12001), rate.LinearBins(0.0, 2 * math.pi, 61))),
		"H1_H2_df": rate.NDBins((rate.LinearBins(-2.0, +2.0, 12001), rate.LinearBins(0.0, 2 * math.pi, 61))),
		"H1_L1_df": rate.NDBins((rate.LinearBins(-2.0, +2.0, 12001), rate.LinearBins(0.0, 2 * math.pi, 61))),
		"H2_L1_df": rate.NDBins((rate.LinearBins(-2.0, +2.0, 12001), rate.LinearBins(0.0, 2 * math.pi, 61))),
		"H1_H2_dh": rate.NDBins((rate.LinearBins(-2.0, +2.0, 12001), rate.LinearBins(0.0, 2 * math.pi, 61))),
		"H1_L1_dh": rate.NDBins((rate.LinearBins(-2.0, +2.0, 12001), rate.LinearBins(0.0, 2 * math.pi, 61))),
		"H2_L1_dh": rate.NDBins((rate.LinearBins(-2.0, +2.0, 12001), rate.LinearBins(0.0, 2 * math.pi, 61))),
		"H1_H2_dt": dt_binning("H1", "H2"),
		"H1_L1_dt": dt_binning("H1", "L1"),
		"H2_L1_dt": dt_binning("H2", "L1")
	}

	filters = {
		"H1_H2_dband": rate.gaussian_window2d(11, 5),
		"H1_L1_dband": rate.gaussian_window2d(11, 5),
		"H2_L1_dband": rate.gaussian_window2d(11, 5),
		"H1_H2_ddur": rate.gaussian_window2d(11, 5),
		"H1_L1_ddur": rate.gaussian_window2d(11, 5),
		"H2_L1_ddur": rate.gaussian_window2d(11, 5),
		"H1_H2_df": rate.gaussian_window2d(11, 5),
		"H1_L1_df": rate.gaussian_window2d(11, 5),
		"H2_L1_df": rate.gaussian_window2d(11, 5),
		"H1_H2_dh": rate.gaussian_window2d(11, 5),
		"H1_L1_dh": rate.gaussian_window2d(11, 5),
		"H2_L1_dh": rate.gaussian_window2d(11, 5),
		"H1_H2_dt": rate.gaussian_window2d(11, 5),
		"H1_L1_dt": rate.gaussian_window2d(11, 5),
		"H2_L1_dt": rate.gaussian_window2d(11, 5)
	}

	def __init__(self):
		self.distributions = CoincParamsDistributions(**self.binnings)

	def add_noninjections(self, param_func, database, *args):
		# iterate over burst<-->burst coincs
		for is_background, events, offsetvector in get_noninjections(database):
			if is_background:
				self.distributions.add_background(param_func(events, offsetvector, *args))
			else:
				self.distributions.add_zero_lag(param_func(events, offsetvector, *args))

	def add_injections(self, param_func, database, *args):
		# iterate over burst<-->burst coincs matching injections
		# "exactly"
		for sim, events, offsetvector in get_injections(database):
			self.distributions.add_injection(param_func(events, offsetvector, *args))

	def finish(self):
		self.distributions.finish(filters = self.filters)


#
# =============================================================================
#
#                                     I/O
#
# =============================================================================
#


#
# XML construction and parsing
#


def get_coincparamsdistributions(xmldoc, name):
	coincparamsdistributions, process_id = CoincParamsDistributions.from_xml(xmldoc, name)
	seglists = lsctables.table.get_table(xmldoc, lsctables.SearchSummaryTable.tableName).get_out_segmentlistdict(set([process_id])).coalesce()
	return coincparamsdistributions, seglists


#
# Construct LIGO Light Weight likelihood distributions document
#


def gen_likelihood_control(coinc_params_distributions, seglists, name = u"ligolw_burca_tailor", comment = u""):
	xmldoc = ligolw.Document()
	node = xmldoc.appendChild(ligolw.LIGO_LW())

	node.appendChild(lsctables.New(lsctables.ProcessTable))
	node.appendChild(lsctables.New(lsctables.ProcessParamsTable))
	node.appendChild(lsctables.New(lsctables.SearchSummaryTable))
	process = append_process(xmldoc, comment = comment)
	llwapp.append_search_summary(xmldoc, process, ifos = seglists.keys(), inseg = seglists.extent_all(), outseg = seglists.extent_all())

	node.appendChild(coinc_params_distributions.to_xml(process, name))

	llwapp.set_process_end_time(process)

	return xmldoc


#
# I/O
#


def load_likelihood_data(filenames, name, verbose = False):
	coincparamsdistributions = None
	for n, filename in enumerate(filenames):
		if verbose:
			print >>sys.stderr, "%d/%d:" % (n + 1, len(filenames)),
		xmldoc = utils.load_filename(filename, verbose = verbose)
		if coincparamsdistributions is None:
			coincparamsdistributions, seglists = get_coincparamsdistributions(xmldoc, name)
		else:
			a, b = get_coincparamsdistributions(xmldoc, name)
			coincparamsdistributions += a
			seglists |= b
			del a, b
		xmldoc.unlink()
	return coincparamsdistributions, seglists


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
#                           param_dist_definer:table
#
# =============================================================================
#


ParamDistDefinerID = ilwd.get_ilwdchar_class(u"param_dist_definer", u"param_dist_def_id")


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
	next_id = ParamDistDefinerID(0)


class ParamDistDefiner(object):
	__slots__ = ParamDistDefinerTable.validcolumns.keys()


ParamDistDefinerTable.RowType = ParamDistDefiner
