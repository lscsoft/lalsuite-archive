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


from pylal import date
from pylal import inject
from pylal import rate
from pylal import SnglBurstUtils


__author__ = "Kipp Cannon <kipp@gravity.phys.uwm.edu>"
__version__ = "$Revision$"[11:-2]
__date__ = "$Date$"[7:-2]


#
# =============================================================================
#
#                                    Misc.
#
# =============================================================================
#


def injection_was_made(sim, seglist, instruments):
	"""
	Return True if the peak times for sim for all of the instruments
	lie within the segment list.
	"""
	return False not in map(lambda i: sim.get_peak(i) in seglist, instruments)


#
# =============================================================================
#
#                            h_{rss} in Instrument
#
# =============================================================================
#


def hrss_in_instrument(sim, instrument):
	"""
	Given an injection and an instrument, compute and return the h_rss
	of the injection as should be observed in the instrument.  That is,
	project the waveform onto the instrument, and return the integrated
	strain squared.
	"""
	if sim.coordinates != "EQUATORIAL":
		raise ValueError, sim.coordinates
	detector = inject.cached_detector[{
		"H1": "LHO_4k",
		"H2": "LHO_2k",
		"L1": "LLO_4k"
	}[instrument]]
	gps = {
		"H1": date.LIGOTimeGPS(sim.h_peak_time, sim.h_peak_time_ns),
		"H2": date.LIGOTimeGPS(sim.h_peak_time, sim.h_peak_time_ns),
		"L1": date.LIGOTimeGPS(sim.l_peak_time, sim.l_peak_time_ns)
	}[instrument]
	fplus, fcross = inject.XLALComputeDetAMResponse(
		detector.response,
		sim.longitude,	# ra
		sim.latitude,	# dec
		sim.polarization,
		date.XLALGreenwichMeanSiderealTime(gps)
	)
	return abs(fplus) * sim.hrss


#
# =============================================================================
#
#                        Iterate Over Found Injections
#
# =============================================================================
#


def found_injections(contents, instrument):
	for values in contents.connection.cursor().execute("""
SELECT DISTINCT sim_burst.* FROM
	sim_burst
	JOIN coinc_event_map AS a ON (
		sim_burst.simulation_id == a.event_id
		AND a.table_name == 'sim_burst'
	)
	JOIN coinc_event_map AS b ON (
		a.coinc_event_id == b.coinc_event_id
		AND b.table_name == 'sngl_burst'
	)
	JOIN sngl_burst ON (
		b.event_id == sngl_burst.event_id
	)
WHERE
	sngl_burst.ifo == ?
	""", (instrument,)):
		yield self.sim_burst_table._row_from_cols(values)


#
# =============================================================================
#
#                             Efficiency Contours
#
# =============================================================================
#


#
# h_{rss} vs. peak frequency
#


class Efficiency_hrss_vs_freq(object):
	def __init__(self, instrument, hrss_func, error):
		self.instrument = instrument
		self.hrss_func = hrss_func
		self.error = error
		self.num_injections = 0
		self.injected_x = []
		self.injected_y = []
		self.found_x = []
		self.found_y = []

	def add_contents(self, contents):
		self.num_injections += len(contents.sim_burst_table)
		for sim in contents.sim_burst_table:
			if injection_was_made(sim, contents.seglists[self.instrument], [self.instrument]):
				self.injected_x.append(sim.freq)
				self.injected_y.append(self.hrss_func(sim, self.instrument))
		for sim in found_injections(contents, self.instrument):
			self.found_x.append(sim.freq)
			self.found_y.append(self.hrss_func(sim, self.instrument))

	def finish(self):
		self.efficiency = rate.BinnedRatios(rate.Bins(min(self.injected_x), max(self.injected_x), 256, min(self.injected_y), max(self.injected_y), 256, spacing = ["log", "log"]))
		map(self.efficiency.incdenominator, zip(self.injected_x, self.injected_y))
		map(self.efficiency.incnumerator, zip(self.found_x, self.found_y))

		# 1 / error^2 is the number of injections that need to be
		# within the window in order for the uncertainty in that
		# number to be = error.  multiplying by bins_per_inj tells
		# us how many bins the window needs to cover, and taking
		# the square root translates that into the window's length
		# on a side in bins.
		bins_per_inj = self.efficiency.used() / float(self.num_injections)
		window_size = (bins_per_inj / self.error**2)**0.5
		if window_size > 100:
			raise ValueError, "smoothing filter too large (not enough injections)"

		# smooth the efficiency data.
		rate.filter_binned_ratios(self.efficiency, rate.gaussian_window2d(window_size, window_size))

		# regularize to prevent divide-by-zero errors
		self.efficiency.regularize()


def plot_Efficiency_hrss_vs_freq(efficiency):
	"""
	Generate a plot from an Efficiency_hrss_vs_freq instance.
	"""
	plot = SnglBurstUtils.BurstPlot("Frequency (Hz)", r"$h_{\mathrm{rss}}$")
	plot.axes.loglog()

	xcoords, ycoords = efficiency.efficiency.centres()
	zvals = efficiency.efficiency.ratio()
	cset = plot.axes.contour(xcoords, ycoords, numpy.transpose(zvals), [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
	plot.axes.set_title(r"Injection Detection Efficiency (%d Injections, Contours at 10\%% Intervals, %g\%% Uncertainty)" % (efficiency.num_injections, 100 * efficiency.error))
	return plot.fig
