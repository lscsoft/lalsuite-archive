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


import bisect
import math
import numpy
from scipy.interpolate import interpolate
import sys


from glue.ligolw import lsctables


__author__ = "Kipp Cannon <kipp@gravity.phys.uwm.edu>"
__version__ = "$Revision$"[11:-2]
__date__ = "$Date$"[7:-2]


#
# =============================================================================
#
#                                  Likelihood
#
# =============================================================================
#


#
# starting from Bayes' theorem:
#
# P(coinc is a g.w. | its parameters)
#     P(those parameters | a coinc known to be a g.w.) * P(coinc is g.w.)
#   = -------------------------------------------------------------------
#                                P(parameters)
#
#     P(those parameters | a coinc known to be a g.w.) * P(coinc is g.w.)
#   = -------------------------------------------------------------------
#     P(noise params) * P(coinc is not g.w.) + P(inj params) * P(coinc is g.w.)
#
#                       P(inj params) * P(coinc is g.w.)
#   = -------------------------------------------------------------------
#     P(noise params) * [1 - P(coinc is g.w.)] + P(inj params) * P(coinc is g.w.)
#
#                        P(inj params) * P(coinc is g.w.)
#   = ----------------------------------------------------------------------
#     P(noise params) + [P(inj params) - P(noise params)] * P(coinc is g.w.)
#


#
# Interpolator wrappers.
#


class interp1d(interpolate.interp1d):
	def __init__(self, x, y, fill_value = 0.0):
		# Extrapolate x and y arrays by half a bin at each end.
		# This is done because the BinnedArray class in pylal.rate
		# returns x co-ordinates as the bin centres, which is
		# correct but it means that an event can have a set of
		# parameter values that lie beyond the end of the x
		# co-ordinate array.  The parameters are still in the last
		# bin, but in the outer half of it.

		x = numpy.hstack((x[0] + (x[0] - x[1]) / 2.0, x, x[-1] + (x[-1] - x[-2]) / 2.0))
		y = numpy.hstack((y[0] + (y[0] - y[1]) / 2.0, y, y[-1] + (y[-1] - y[-2]) / 2.0))

		# Because the y values are assumed to be probability
		# densities they cannot be negative.  This constraint is
		# imposed by clipping the y array to 0, and using a linear
		# interpolator below.  Assuming the input is valid, the
		# only reason the y array could have negative numbers in it
		# is the extrapolation that has just been done.

		y = numpy.clip(y, 0.0, float("inf"))

		# Build the interpolator.  Note the use of fill_value as
		# the return value for x co-ordinates outside the domain of
		# the interpolator.

		interpolate.interp1d.__init__(self, x, y, kind = "linear", bounds_error = False, fill_value = fill_value)


class interp2d(interpolate.interp2d):
	# unlike the real interp2d, this version requires the x and y
	# arrays to be just the co-ordinates of the columns and rows of a
	# regular mesh.
	def __init__(self, x, y, z, fill_value = 0.0):
		# Extrapolate x, y and z arrays by half a bin at each end.
		# See interp1d for an explanation.

		x = numpy.hstack((x[0] + (x[0] - x[1]) / 2.0, x, x[-1] + (x[-1] - x[-2]) / 2.0))
		y = numpy.hstack((y[0] + (y[0] - y[1]) / 2.0, y, y[-1] + (y[-1] - y[-2]) / 2.0))
		z = numpy.vstack((z[0] + (z[0] - z[1]) / 2.0, z, z[-1] + (z[-1] - z[-2]) / 2.0))
		z = numpy.transpose(z)
		z = numpy.vstack((z[0] + (z[0] - z[1]) / 2.0, z, z[-1] + (z[-1] - z[-2]) / 2.0))
		z = numpy.transpose(z)

		# Clip the z array to 0.  See interp1d for an explanation.

		z = numpy.clip(z, 0.0, float("inf"))

		# Build the interpolator.  Note the use of fill_value as
		# the return value for co-ordinates outside the domain of
		# the interpolator.

		# FIXME:  scipy's 2D interpolator is busted, put this back
		# when it's fixed.
		#interpolate.interp2d.__init__(self, x, y, z, kind = "linear", bounds_error = False, fill_value = fill_value)
		self.x = x
		self.y = y
		self.z = z

	# FIXME:  scipy's 2D interpolator is busted.  remove this when it's
	# fixed.
	def __call__(self, x, y):
		return (self.z[bisect.bisect(self.x, x), bisect.bisect(self.y, y)],)


#
# Class for computing foreground likelihoods from the measurements in a
# CoincParamsDistributions instance.
#


class Likelihood(object):
	def __init__(self, coinc_param_distributions):
		# check input
		if set(coinc_param_distributions.background_rates.keys()) != set(coinc_param_distributions.injection_rates.keys()):
			raise ValueError, "distribution density name mismatch"
		for name, binnedarray in coinc_param_distributions.background_rates.items():
			if len(binnedarray.array.shape) != len(coinc_param_distributions.injection_rates[name].array.shape):
				raise ValueError, "distribution density dimension mismatch"

		# construct interpolators from the distribution data
		self.background_rates = {}
		self.injection_rates = {}
		for name, binnedarray in coinc_param_distributions.background_rates.items():
			if len(binnedarray.array.shape) == 1:
				x, = binnedarray.centres()
				self.background_rates[name] = interp1d(x, binnedarray.array)
			elif len(binnedarray.array.shape) == 2:
				x, y = binnedarray.centres()
				self.background_rates[name] = interp2d(x, y, binnedarray.array)
			else:
				raise ValueError, "distribution densities with >2 dimensions not supported"
		for name, binnedarray in coinc_param_distributions.injection_rates.items():
			if len(binnedarray.array.shape) == 1:
				x, = binnedarray.centres()
				self.injection_rates[name] = interp1d(x, binnedarray.array)
			elif len(binnedarray.array.shape) == 2:
				x, y = binnedarray.centres()
				self.injection_rates[name] = interp2d(x, y, binnedarray.array)
			else:
				raise ValueError, "distribution densities with >2 dimensions not supported"

	def set_P_gw(self, P):
		self.P_gw = P

	def P(self, params_func, events, offsetdict, *params_func_extra_args):
		P_bak = 1.0
		P_inj = 1.0
		for name, value in params_func(events, offsetdict, *params_func_extra_args).items():
			P_bak *= self.background_rates[name](*value)[0]
			P_inj *= self.injection_rates[name](*value)[0]
		return P_bak, P_inj

	def __call__(self, params_func, events, offsetdict, *params_func_extra_args):
		"""
		Compute the likelihood that the coincident n-tuple of
		events is the result of a gravitational wave:  the
		probability that the hypothesis "the events are a
		gravitational wave" is correct, in the context of the
		measured background and foreground distributions, and the
		intrinsic rate of gravitational wave coincidences.  offsets
		is a dictionary of instrument --> offset mappings to be
		used to time shift the events before comparison.
		"""
		P_bak, P_inj = self.P(params_func, events, offsetdict, *params_func_extra_args)
		return (P_inj * self.P_gw) / (P_bak + (P_inj - P_bak) * self.P_gw)


class Confidence(Likelihood):
	def __call__(self, params_func, events, offsetdict, *params_func_extra_args):
		"""
		Compute the confidence that the list of events are the
		result of a gravitational wave:  -ln[1 - P(gw)], where
		P(gw) is the likelihood returned by the Likelihood class.
		A set of events very much like gravitational waves will
		have a likelihood of being a gravitational wave very close
		to 1, so 1 - P is a small positive number, and so -ln of
		that is a large positive number.
		"""
		P_bak, P_inj = self.P(params_func, events, offsetdict, *params_func_extra_args)
		return  math.log(P_bak + (P_inj - P_bak) * self.P_gw) - math.log(P_inj) - math.log(self.P_gw)


class LikelihoodRatio(Likelihood):
	def set_P_gw(self, P):
		"""
		Raises NotImplementedError.  The likelihood ratio is
		computed without using this parameter.
		"""
		raise NotImplementedError

	def __call__(self, params_func, events, offsetdict, *params_func_extra_args):
		"""
		Compute the likelihood ratio for the hypothesis that the
		list of events are the result of a gravitational wave.  The
		likelihood ratio is the ratio P(inj params) / P(noise
		params).  The probability that the events are the result of
		a gravitiational wave is a monotonically increasing
		function of the likelihood ratio, so ranking events from
		"most like a gravitational wave" to "least like a
		gravitational wave" can be performed by calculating the
		likelihood ratios, which has the advantage of not requiring
		a prior probability to be provided.
		"""
		P_bak, P_inj = self.P(params_func, events, offsetdict, *params_func_extra_args)
		if P_bak == 0.0 and P_inj == 0.0:
			# this can happen.  "correct" answer is 0, not NaN,
			# because if a tuple of events has been found in a
			# region of parameter space where the probability
			# of an injection occuring is 0 then there is no
			# way this is an injection.  there is also,
			# aparently, no way it's a noise event, but that's
			# irrelevant because we are supposed to be
			# computing something that is a monotonically
			# increasing function of the probability that an
			# event tuple is a gravitational wave, which is 0
			# in this part of the parameter space.
			return 0.0
		return  P_inj / P_bak


#
# =============================================================================
#
#                              Library Interface
#
# =============================================================================
#


#
# Main routine
#


def ligolw_burca2(database, likelihood_ratio, params_func, verbose = False, params_func_extra_args = ()):
	"""
	Assigns likelihood ratio values to excess power coincidences.
	database is pylal.SnglBurstUtils.CoincDatabase instance, and
	likelihood_ratio is a LikelihoodRatio class instance.
	"""
	#
	# Find document parts.
	#

	time_slides = database.time_slide_table.as_dict()

	#
	# Iterate over all coincs, assigning likelihood ratios to
	# burst+burst coincs if the document contains them.
	#

	if verbose:
		print >>sys.stderr, "computing likelihood ratios ..."
		n_coincs = len(database.coinc_table)

	cursor = database.connection.cursor()
	cursor.execute("""
CREATE TEMPORARY VIEW
	coinc_burst_map
AS
	SELECT
		*
	FROM
		sngl_burst
		JOIN coinc_event_map ON (
			coinc_event_map.table_name == 'sngl_burst'
			AND coinc_event_map.event_id == sngl_burst.event_id
		)
	""")
	for n, (coinc_event_id, time_slide_id) in enumerate(database.connection.cursor().execute("SELECT coinc_event_id, time_slide_id FROM coinc_event WHERE coinc_def_id == ?", (database.bb_definer_id,))):
		if verbose and not n % 200:
			print >>sys.stderr, "\t%.1f%%\r" % (100.0 * n / n_coincs),

		# retrieve sngl_burst events
		events = map(database.sngl_burst_table.row_from_cols, cursor.execute("""SELECT * FROM coinc_burst_map WHERE coinc_event_id == ?""", (coinc_event_id,)))

		# compute and store likelihood ratio
		cursor.execute("""
UPDATE
	coinc_event
SET
	likelihood = ?
WHERE
	coinc_event_id == ?
		""", (likelihood_ratio(params_func, events, time_slides[time_slide_id], *params_func_extra_args), coinc_event_id))
	if verbose:
		print >>sys.stderr, "\t100.0%"

	cursor.execute("""DROP VIEW coinc_burst_map""")

	#
	# Done
	#

	return database
