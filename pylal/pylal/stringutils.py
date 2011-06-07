# Copyright (C) 2009  Kipp Cannon
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
import scipy.stats
import sys


from glue import iterutils
from glue import segmentsUtils
from glue.ligolw import lsctables
from glue.ligolw import utils
from glue.ligolw.utils import process as ligolw_process
from pylal import ligolw_burca_tailor
from pylal import git_version
from pylal import inject
from pylal import rate
from pylal import snglcoinc


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date


#
# =============================================================================
#
#                             Likelihood Machinery
#
# =============================================================================
#


#
# Make a look-up table of time-of-arrival triangulators
#


def triangulators(timing_uncertainties):
	"""
	Return a dictionary of snglcoinc.TOATriangulator objects
	initialized for a variety of instrument combinations.
	timing_uncertainties is a dictionary of instrument->\Delta t pairs.
	The return value is a dictionary of (instrument
	tuple)->TOATrangulator mappings.  The instrument names in each
	tuple are sorted in alphabetical order, and the triangulators are
	constructed with the instruments in that order (the the
	documentation for snglcoinc.TOATriangulator for more information).

	Example:

	>>> x = triangulators({"H1": 0.005, "L1": 0.005, "V1": 0.005})

	constructs a dictionary of triangulators for every combination of
	two or more instruments that can be constructed from those three.

	The program lalapps_string_plot_binj can be used to measure the
	timing uncertainties for the instruments in a search.
	"""
	allinstruments = sorted(timing_uncertainties.keys())

	triangulators = {}
	for n in range(2, len(allinstruments) + 1):
		for instruments in iterutils.choices(allinstruments, n):
			triangulators[instruments] = snglcoinc.TOATriangulator([inject.cached_detector[inject.prefix_to_name[instrument]].location for instrument in instruments], [timing_uncertainties[instrument] for instrument in instruments])

	return triangulators


#
# A look-up table used to convert instrument names to powers of 2.  Why?
# To create a bidirectional mapping between combinations of instrument
# names and integers so we can use a pylal.rate style binning for the
# instrument combinations
#


instrument_to_factor = dict((instrument, int(2**n)) for n, instrument in enumerate(("G1", "H1", "H2", "H+", "H-", "L1", "V1")))


def instruments_to_category(instruments):
	return sum(instrument_to_factor[instrument] for instrument in instruments)


def category_to_instruments(category):
	return set(instrument for instrument, factor in instrument_to_factor.items() if category & factor)


#
# Coinc params function
#


def coinc_params_func(events, offsetvector, triangulators):
	#
	# check for coincs that have been vetoed entirely
	#

	if len(events) < 2:
		return None

	#
	# Initialize the parameter dictionary, sort the events by
	# instrument name (the multi-instrument parameters are defined for
	# the instruments in this order and the triangulators are
	# constructed this way too), and retrieve the sorted instrument
	# names
	#

	params = {}
	events = tuple(sorted(events, key = lambda event: event.ifo))
	instruments = tuple(event.ifo for event in events)

	#
	# zero-instrument parameters
	#

	ignored, ignored, ignored, rss_timing_residual = triangulators[instruments](tuple(event.get_peak() + offsetvector[event.ifo] for event in events))
	# FIXME:  rss_timing_residual is forced to 0 to disable this
	# feature.  all the code to compute it properly is still here and
	# given suitable initializations, the distribution data is still
	# two-dimensional and has a suitable filter applied to it, but all
	# events are forced into the RSS_{\Delta t} = 0 bin, in effect
	# removing that dimension from the data.  We can look at this again
	# sometime in the future if we're curious why it didn't help.  Just
	# delete the next line and you're back in business.
	rss_timing_residual = 0.0
	params["instrumentgroup,rss_timing_residual"] = (instruments_to_category(instruments), rss_timing_residual)

	#
	# one-instrument parameters
	#

	for event in events:
		prefix = "%s_" % event.ifo

		params["%ssnr2_chi2" % prefix] = (event.snr**2.0, event.chisq / event.chisq_dof)

	#
	# two-instrument parameters.  note that events are sorted by
	# instrument
	#

	for event1, event2 in iterutils.choices(events, 2):
		assert event1.ifo != event2.ifo

		prefix = "%s_%s_" % (event1.ifo, event2.ifo)

		dt = float((event1.get_peak() + offsetvector[event1.ifo]) - (event2.get_peak() + offsetvector[event2.ifo]))
		params["%sdt" % prefix] = (dt,)

		dA = math.log10(abs(event1.amplitude / event2.amplitude))
		params["%sdA" % prefix] = (dA,)

		# f_cut = central_freq + bandwidth/2
		f_cut1 = event1.central_freq + event1.bandwidth / 2
		f_cut2 = event2.central_freq + event2.bandwidth / 2
		df = float((math.log10(f_cut1) - math.log10(f_cut2)) / (math.log10(f_cut1) + math.log10(f_cut2)))
		params["%sdf" % prefix] = (df,)

	#
	# done
	#

	return params


#
# Parameter distributions
#


def dt_binning(instrument1, instrument2):
	dt = 0.005 + inject.light_travel_time(instrument1, instrument2)	# seconds
	return rate.NDBins((rate.ATanBins(-dt, +dt, 801),))


class DistributionsStats(object):
	"""
	A class used to populate a CoincParamsDistribution instance with
	the data from the outputs of ligolw_burca and ligolw_binjfind.
	"""

	binnings = {
		"H1_snr2_chi2": rate.NDBins((rate.ATanLogarithmicBins(10, 1e7, 801), rate.ATanLogarithmicBins(.1, 1e4, 801))),
		"H2_snr2_chi2": rate.NDBins((rate.ATanLogarithmicBins(10, 1e7, 801), rate.ATanLogarithmicBins(.1, 1e4, 801))),
		"L1_snr2_chi2": rate.NDBins((rate.ATanLogarithmicBins(10, 1e7, 801), rate.ATanLogarithmicBins(.1, 1e4, 801))),
		"V1_snr2_chi2": rate.NDBins((rate.ATanLogarithmicBins(10, 1e7, 801), rate.ATanLogarithmicBins(.1, 1e4, 801))),
		"H1_H2_dt": dt_binning("H1", "H2"),
		"H1_L1_dt": dt_binning("H1", "L1"),
		"H1_V1_dt": dt_binning("H1", "V1"),
		"H2_L1_dt": dt_binning("H2", "L1"),
		"H2_V1_dt": dt_binning("H2", "V1"),
		"L1_V1_dt": dt_binning("L1", "V1"),
		"H1_H2_dA": rate.NDBins((rate.ATanBins(-0.5, +0.5, 801),)),
		"H1_L1_dA": rate.NDBins((rate.ATanBins(-0.5, +0.5, 801),)),
		"H1_V1_dA": rate.NDBins((rate.ATanBins(-0.5, +0.5, 801),)),
		"H2_L1_dA": rate.NDBins((rate.ATanBins(-0.5, +0.5, 801),)),
		"H2_V1_dA": rate.NDBins((rate.ATanBins(-0.5, +0.5, 801),)),
		"L1_V1_dA": rate.NDBins((rate.ATanBins(-0.5, +0.5, 801),)),
		"H1_H2_df": rate.NDBins((rate.ATanBins(-0.2, +0.2, 501),)),
		"H1_L1_df": rate.NDBins((rate.ATanBins(-0.2, +0.2, 501),)),
		"H1_V1_df": rate.NDBins((rate.ATanBins(-0.2, +0.2, 501),)),
		"H2_L1_df": rate.NDBins((rate.ATanBins(-0.2, +0.2, 501),)),
		"H2_V1_df": rate.NDBins((rate.ATanBins(-0.2, +0.2, 501),)),
		"L1_V1_df": rate.NDBins((rate.ATanBins(-0.2, +0.2, 501),)),
		# instrument group bin centres are at 1, 2, 3, ...;  only
		# non-negative rss timing residual bins will be used but we
		# want a binning that's linear at the origin so instead of
		# inventing a new one we just use atan bins that are
		# symmetric about 0
		"instrumentgroup,rss_timing_residual": rate.NDBins((rate.LinearBins(0.5, sum(instrument_to_factor.values()) + 0.5, sum(instrument_to_factor.values())), rate.ATanBins(-0.02, +0.02, 1001)))
	}

	filters = {
		"H1_snr2_chi2": rate.gaussian_window2d(11, 11, sigma = 20),
		"H2_snr2_chi2": rate.gaussian_window2d(11, 11, sigma = 20),
		"L1_snr2_chi2": rate.gaussian_window2d(11, 11, sigma = 20),
		"V1_snr2_chi2": rate.gaussian_window2d(11, 11, sigma = 20),
		"H1_H2_dt": rate.gaussian_window(11, sigma = 20),
		"H1_L1_dt": rate.gaussian_window(11, sigma = 20),
		"H1_V1_dt": rate.gaussian_window(11, sigma = 20),
		"H2_L1_dt": rate.gaussian_window(11, sigma = 20),
		"H2_V1_dt": rate.gaussian_window(11, sigma = 20),
		"L1_V1_dt": rate.gaussian_window(11, sigma = 20),
		"H1_H2_dA": rate.gaussian_window(11, sigma = 20),
		"H1_L1_dA": rate.gaussian_window(11, sigma = 20),
		"H1_V1_dA": rate.gaussian_window(11, sigma = 20),
		"H2_L1_dA": rate.gaussian_window(11, sigma = 20),
		"H2_V1_dA": rate.gaussian_window(11, sigma = 20),
		"L1_V1_dA": rate.gaussian_window(11, sigma = 20),
		"H1_H2_df": rate.gaussian_window(11, sigma = 20),
		"H1_L1_df": rate.gaussian_window(11, sigma = 20),
		"H1_V1_df": rate.gaussian_window(11, sigma = 20),
		"H2_L1_df": rate.gaussian_window(11, sigma = 20),
		"H2_V1_df": rate.gaussian_window(11, sigma = 20),
		"L1_V1_df": rate.gaussian_window(11, sigma = 20),
		# instrument group filter is a no-op, should produce a
		# 1-bin top-hat window.
		"instrumentgroup,rss_timing_residual": rate.gaussian_window2d(1e-100, 11, sigma = 20)
	}

	def __init__(self):
		self.distributions = ligolw_burca_tailor.CoincParamsDistributions(**self.binnings)

	def add_noninjections(self, param_func, database, param_func_args = ()):
		# iterate over burst<-->burst coincs
		for is_background, events, offsetvector in ligolw_burca_tailor.get_noninjections(database):
			events = [event for event in events if event.ifo not in database.vetoseglists or event.get_peak() not in database.vetoseglists[event.ifo]]
			if is_background:
				self.distributions.add_background(param_func(events, offsetvector, *param_func_args))
			else:
				self.distributions.add_zero_lag(param_func(events, offsetvector, *param_func_args))

	def add_slidelessbackground(self, param_func, database, experiments, param_func_args = ()):
		# FIXME:  this needs to be taught how to not slide H1 and
		# H2 with respect to each other

		# segment lists
		seglists = database.seglists - database.vetoseglists

		# construct the event list dictionary.  remove vetoed
		# events from the lists and save event peak times so they
		# can be restored later
		eventlists = {}
		orig_peak_times = {}
		for event in database.sngl_burst_table:
			if event.get_peak() in seglists[event.ifo]:
				try:
					eventlists[event.ifo].append(event)
				except KeyError:
					eventlists[event.ifo] = [event]
				orig_peak_times[event] = event.get_peak()

		# parse the --thresholds H1,L1=... command-line options from burca
		delta_t = [float(threshold.split("=")[-1]) for threshold in ligolw_process.get_process_params(database.xmldoc, "ligolw_burca", "--thresholds")]
		if not all(delta_t[0] == threshold for threshold in delta_t[1:]):
			raise ValueError, "\Delta t is not unique in ligolw_burca arguments"
		delta_t = delta_t.pop()

		# construct the coinc generator.  note that H1+H2-only
		# coincs are forbidden, which is affected here by removing
		# that instrument combination from mu_conic
		mu, tau = snglcoinc.slideless_coinc_generator_mu_tau(eventlists, segmentlists, delta_t)
		zero_lag_offset_vector = dict((instrument, 0.0) for instrument in mu)
		mu_coinc = snglcoinc.slideless_coinc_generator_rates(mu, tau)
		if frozenset(("H1", "H2")) in mu_coinc:
			del mu_coinc[frozenset(("H1", "H2"))]
		coinc_generator = snglcoinc.slideless_coinc_generator(eventlists, mu_coinc, tau, lsctables.SnglBurst.get_peak)
		toa_generator = dict((instruments, snglcoinc.slideless_coinc_generator_plausible_toas(instruments, tau)) for instruments in mu_coinc.keys())

		# how many coincs?  the expected number is obtained by
		# multiplying the total zero-lag time for which at least
		# two instruments were on by the sum of the rates for all
		# coincs to get the mean number of coincs per zero-lag
		# observation time, and multiplying that by the number of
		# experiments the background should simulate to get the
		# mean number of background events to simulate.  the actual
		# number simulated is a Poisson-distributed RV with that
		# mean.
		n_coincs, = scipy.stats.poisson.rvs(float(abs(segmentsUtils.vote(seglists.values(), 2))) * sum(mu_coinc.values()) * experiments)

		# generate synthetic background coincs
		for n, events in enumerate(coinc_generator):
			# assign fake peak times
			toas = toa_generator[frozenset(event.ifo for event in events)].next()
			for event in events:
				event.set_peak(toas[event.ifo])
			# compute coincidence parameters
			self.distributions.add_background(param_func(events, zero_lag_offset_vector, *param_func_args))
			if n > n_coincs:
				break

		# restore original peak times
		for event, peak_time in orig_peak_times.iteritems():
			event.set_peak(peak_time)

	def add_injections(self, param_func, database, weight_func = lambda sim: 1.0, param_func_args = ()):
		# iterate over burst<-->burst coincs matching injections
		# "exactly"
		for sim, events, offsetvector in ligolw_burca_tailor.get_injections(database):
			events = [event for event in events if event.ifo not in database.vetoseglists or event.get_peak() not in database.vetoseglists[event.ifo]]
			self.distributions.add_injection(param_func(events, offsetvector, *param_func_args), weight = weight_func(sim))

	def finish(self):
		self.distributions.finish(filters = self.filters)


#
# Livetime
#


def time_slides_livetime(seglists, time_slides, min_instruments, verbose = False, clip = None):
	"""
	seglists is a segmentlistdict of times when each of a set of
	instruments were on, time_slides is a sequence of
	instrument-->offset dictionaries, each vector of offsets in the
	sequence is applied to the segmentlists and the total time during
	which at least min_instruments were on is summed and returned.  If
	clip is not None, after each offset vector is applied to seglists
	the result is intersected with clip before computing the livetime.
	If verbose is True then progress reports are printed to stderr.
	"""
	livetime = 0.0
	seglists = seglists.copy()	# don't modify original
	N = len(time_slides)
	if verbose:
		print >>sys.stderr, "computing the live time for %d time slides:" % N
	for n, time_slide in enumerate(time_slides):
		if verbose:
			print >>sys.stderr, "\t%.1f%%\r" % (100.0 * n / N),
		seglists.offsets.update(time_slide)
		if clip is None:
			livetime += float(abs(segmentsUtils.vote(seglists.values(), min_instruments)))
		else:
			livetime += float(abs(segmentsUtils.vote((seglists & clip).values(), min_instruments)))
	if verbose:
		print >>sys.stderr, "\t100.0%"
	return livetime


def time_slides_livetime_for_instrument_combo(seglists, time_slides, instruments, verbose = False, clip = None):
	"""
	like time_slides_livetime() except computes the time for which
	exactly the instruments given by the sequence instruments were on
	(and nothing else).
	"""
	livetime = 0.0
	# segments for instruments that must be on
	onseglists = seglists.copy(keys = instruments)
	# segments for instruments that must be off
	offseglists = seglists.copy(keys = set(seglists) - set(instruments))
	N = len(time_slides)
	if verbose:
		print >>sys.stderr, "computing the live time for %s in %d time slides:" % (", ".join(instruments), N)
	for n, time_slide in enumerate(time_slides):
		if verbose:
			print >>sys.stderr, "\t%.1f%%\r" % (100.0 * n / N),
		onseglists.offsets.update(time_slide)
		offseglists.offsets.update(time_slide)
		if clip is None:
			livetime += float(abs(onseglists.intersection(onseglists.keys()) - offseglists.union(offseglists.keys())))
		else:
			livetime += float(abs((onseglists & clip).intersection(onseglists.keys()) - offseglists.union(offseglists.keys())))
	if verbose:
		print >>sys.stderr, "\t100.0%"
	return livetime


#
# I/O
#


def get_coincparamsdistributions(xmldoc):
	coincparamsdistributions, process_id = ligolw_burca_tailor.coinc_params_distributions_from_xml(xmldoc, u"string_cusp_likelihood")
	seglists = lsctables.table.get_table(xmldoc, lsctables.SearchSummaryTable.tableName).get_out_segmentlistdict(set([process_id])).coalesce()
	return coincparamsdistributions, seglists


def load_likelihood_data(filenames, verbose = False):
	coincparamsdistributions = None
	for n, filename in enumerate(filenames):
		if verbose:
			print >>sys.stderr, "%d/%d:" % (n + 1, len(filenames)),
		xmldoc = utils.load_filename(filename, verbose = verbose)
		if coincparamsdistributions is None:
			coincparamsdistributions, seglists = get_coincparamsdistributions(xmldoc)
		else:
			a, b = get_coincparamsdistributions(xmldoc)
			coincparamsdistributions += a
			seglists |= b
			del a, b
		xmldoc.unlink()
	return coincparamsdistributions, seglists


def write_likelihood_data(filename, coincparamsdistributions, seglists, verbose = False):
	utils.write_filename(ligolw_burca_tailor.gen_likelihood_control(coincparamsdistributions, seglists, name = u"string_cusp_likelihood"), filename, verbose = verbose, gz = (filename or "stdout").endswith(".gz"))


#
# =============================================================================
#
#                              Database Utilities
#
# =============================================================================
#


def create_recovered_likelihood_table(connection, bb_coinc_def_id):
	"""
	Create a temporary table named "recovered_likelihood" containing
	two columns:  "simulation_id", the simulation_id of an injection,
	and "likelihood", the highest likelihood ratio at which that
	injection was recovered by a coincidence of type bb_coinc_def_id.
	"""
	cursor = connection.cursor()
	cursor.execute("""
CREATE TEMPORARY TABLE recovered_likelihood (simulation_id TEXT PRIMARY KEY, likelihood REAL)
	""")
	cursor.execute("""
INSERT OR REPLACE INTO
	recovered_likelihood
SELECT
	sim_burst.simulation_id AS simulation_id,
	MAX(coinc_event.likelihood) AS likelihood
FROM
	sim_burst
	JOIN coinc_event_map AS a ON (
		a.table_name == "sim_burst"
		AND a.event_id == sim_burst.simulation_id
	)
	JOIN coinc_event_map AS b ON (
		b.coinc_event_id == a.coinc_event_id
	)
	JOIN coinc_event ON (
		b.table_name == "coinc_event"
		AND b.event_id == coinc_event.coinc_event_id
	)
WHERE
	coinc_event.coinc_def_id == ?
GROUP BY
	sim_burst.simulation_id
	""", (bb_coinc_def_id,))
	cursor.close()
