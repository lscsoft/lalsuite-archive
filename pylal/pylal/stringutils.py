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


from glue import iterutils
from pylal import ligolw_burca_tailor
from pylal import git_version
from pylal import inject
from pylal import rate


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
# Coinc params function
#


def coinc_params_func(events, offsetdict):
	params = {}

	#
	# one-instrument parameters
	#

	for event in events:
		prefix = "%s_" % event.ifo

		params["%ssnr2_chi2" % prefix] = (event.snr**2.0, event.chisq / event.chisq_dof)

	#
	# two-instrument parameters
	#

	for event1, event2 in iterutils.choices(sorted(events, lambda a, b: cmp(a.ifo, b.ifo)), 2):
		if event1.ifo == event2.ifo:
			# shouldn't happen, but might as well check for it
			continue

		prefix = "%s_%s_" % (event1.ifo, event2.ifo)

		dt = float((event1.get_peak() + offsetdict[event1.ifo]) - (event2.get_peak() + offsetdict[event2.ifo]))
		params["%sdt" % prefix] = (dt,)

		dA = math.log10(abs(event1.amplitude / event2.amplitude))
		params["%sdA" % prefix] = (dA,)

		df = float((event1.central_freq + 0.5*event1.bandwidth - event2.central_freq - 0.5*event2.bandwidth)/(event1.central_freq + 0.5*event1.bandwidth + event2.central_freq + 0.5*event2.bandwidth))
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
	return rate.NDBins((rate.ATanBins(-dt, +dt, 3001),))


class DistributionsStats(ligolw_burca_tailor.Stats):
	"""
	A subclass of the Stats class used to populate a
	CoincParamsDistribution instance with the data from the outputs of
	ligolw_burca and ligolw_binjfind.
	"""

	binnings = {
		"H1_snr2_chi2": rate.NDBins((rate.ATanLogarithmicBins(10, 1e7, 1201), rate.ATanLogarithmicBins(.1, 1e4, 1201))),
		"H2_snr2_chi2": rate.NDBins((rate.ATanLogarithmicBins(10, 1e7, 1201), rate.ATanLogarithmicBins(.1, 1e4, 1201))),
		"L1_snr2_chi2": rate.NDBins((rate.ATanLogarithmicBins(10, 1e7, 1201), rate.ATanLogarithmicBins(.1, 1e4, 1201))),
		"V1_snr2_chi2": rate.NDBins((rate.ATanLogarithmicBins(10, 1e7, 1201), rate.ATanLogarithmicBins(.1, 1e4, 1201))),
		"H1_H2_dt": dt_binning("H1", "H2"),
		"H1_L1_dt": dt_binning("H1", "L1"),
		"H1_V1_dt": dt_binning("H1", "V1"),
		"H2_L1_dt": dt_binning("H2", "L1"),
		"H2_V1_dt": dt_binning("H2", "V1"),
		"L1_V1_dt": dt_binning("L1", "V1"),
		"H1_H2_dA": rate.NDBins((rate.ATanBins(-0.5, +0.5, 6001),)),
		"H1_L1_dA": rate.NDBins((rate.ATanBins(-0.5, +0.5, 6001),)),
		"H1_V1_dA": rate.NDBins((rate.ATanBins(-0.5, +0.5, 6001),)),
		"H2_L1_dA": rate.NDBins((rate.ATanBins(-0.5, +0.5, 6001),)),
		"H2_V1_dA": rate.NDBins((rate.ATanBins(-0.5, +0.5, 6001),)),
		"L1_V1_dA": rate.NDBins((rate.ATanBins(-0.5, +0.5, 6001),)),
		"H1_H2_df": rate.NDBins((rate.ATanBins(-0.5, +0.5, 6001),)),
		"H1_L1_df": rate.NDBins((rate.ATanBins(-0.5, +0.5, 6001),)),
		"H1_V1_df": rate.NDBins((rate.ATanBins(-0.5, +0.5, 6001),)),
		"H2_L1_df": rate.NDBins((rate.ATanBins(-0.5, +0.5, 6001),)),
		"H2_V1_df": rate.NDBins((rate.ATanBins(-0.5, +0.5, 6001),)),
		"L1_V1_df": rate.NDBins((rate.ATanBins(-0.5, +0.5, 6001),))
	}

	filters = {
		"H1_snr2_chi2": rate.gaussian_window2d(11, 11),
		"H2_snr2_chi2": rate.gaussian_window2d(11, 11),
		"L1_snr2_chi2": rate.gaussian_window2d(11, 11),
		"V1_snr2_chi2": rate.gaussian_window2d(11, 11),
		"H1_H2_dt": rate.gaussian_window(11),
		"H1_L1_dt": rate.gaussian_window(11),
		"H1_V1_dt": rate.gaussian_window(11),
		"H2_L1_dt": rate.gaussian_window(11),
		"H2_V1_dt": rate.gaussian_window(11),
		"L1_V1_dt": rate.gaussian_window(11),
		"H1_H2_dA": rate.gaussian_window(11),
		"H1_L1_dA": rate.gaussian_window(11),
		"H1_V1_dA": rate.gaussian_window(11),
		"H2_L1_dA": rate.gaussian_window(11),
		"H2_V1_dA": rate.gaussian_window(11),
		"L1_V1_dA": rate.gaussian_window(11),
		"H1_H2_df": rate.gaussian_window(11),
		"H1_L1_df": rate.gaussian_window(11),
		"H1_V1_df": rate.gaussian_window(11),
		"H2_L1_df": rate.gaussian_window(11),
		"H2_V1_df": rate.gaussian_window(11),
		"L1_V1_df": rate.gaussian_window(11)
	}

	def __init__(self):
		ligolw_burca_tailor.Stats.__init__(self)
		self.distributions = ligolw_burca_tailor.CoincParamsDistributions(**self.binnings)

	def _add_zero_lag(self, param_func, events, offsetdict, *args):
		self.distributions.add_zero_lag(param_func, events, offsetdict, *args)

	def _add_background(self, param_func, events, offsetdict, *args):
		self.distributions.add_background(param_func, events, offsetdict, *args)

	def _add_injections(self, param_func, sim, events, offsetdict, *args):
		self.distributions.add_injection(param_func, events, offsetdict, *args)

	def finish(self):
		self.distributions.finish(filters = self.filters)


#
# I/O
#


def get_coincparamsdistributions(xmldoc):
	coincparamsdistributions, process_id = ligolw_burca_tailor.coinc_params_distributions_from_xml(xmldoc, u"string_cusp_likelihood")
	return coincparamsdistributions
