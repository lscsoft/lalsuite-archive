# Copyright (C) 2007  Kipp Cannon
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
import numpy
import sys


from pylal import date
from pylal import git_version
from pylal import inject
from pylal import rate
from pylal import SnglBurstUtils


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date


#
# =============================================================================
#
#                                    Misc.
#
# =============================================================================
#


def time_at_instrument(sim, instrument, offsetvector):
	"""
	Return the "time" of the injection, delay corrected for the
	displacement from the geocentre to the given instrument.

	NOTE:  this function does not account for the rotation of the Earth
	that occurs during the transit of the plane wave from the detector
	to the geocentre.  That is, it is assumed the Earth is in the same
	orientation with respect to the celestial sphere when the wave
	passes through the detector as when it passes through the
	geocentre.  The Earth rotates by about 1.5 urad during the 21 ms it
	takes light to travel the radius of the Earth, which corresponds to
	10 m of displacement at the equator, or 33 ns in units of light
	travel time.  Therefore, the failure to do a proper retarded time
	calculation here results in errors no larger than 33 ns, which
	should be insignificant.
	"""
	# the offset is subtracted from the time of the injection.
	# injections are done this way so that when the triggers that
	# result from an injection have the offset vector added to their
	# times the triggers will form a coinc
	t_geocent = sim.get_time_geocent() - offsetvector[instrument]
	ra, dec = sim.get_ra_dec()
	return t_geocent + date.XLALTimeDelayFromEarthCenter(inject.cached_detector[inject.prefix_to_name[instrument]].location, ra, dec, t_geocent)


def on_instruments(sim, seglists, offsetvector):
	"""
	Return a set of the names of the instruments that were on at the
	time of the injection.
	"""
	return set(instrument for instrument, seglist in seglists.items() if time_at_instrument(sim, instrument, offsetvector) in seglist)


def create_sim_burst_best_string_sngl_map(connection, coinc_def_id):
	"""
	Construct a sim_burst --> best matching coinc_event mapping.
	"""
	connection.cursor().execute("""
CREATE TEMPORARY TABLE
	sim_burst_best_string_sngl_map
AS
	SELECT
		sim_burst.simulation_id AS simulation_id,
		(
			SELECT
				sngl_burst.event_id
			FROM
				coinc_event_map AS a
				JOIN coinc_event_map AS b ON (
					b.coinc_event_id == a.coinc_event_id
				)
				JOIN coinc_event ON (
					coinc_event.coinc_event_id == a.coinc_event_id
				)
				JOIN sngl_burst ON (
					b.table_name == 'sngl_burst'
					AND b.event_id == sngl_burst.event_id
				)
			WHERE
				a.table_name == 'sim_burst'
				AND a.event_id == sim_burst.simulation_id
				AND coinc_event.coinc_def_id == ?
			ORDER BY
				(sngl_burst.chisq / sngl_burst.chisq_dof) / (sngl_burst.snr * sngl_burst.snr)
			LIMIT 1
		) AS event_id
	FROM
		sim_burst
	WHERE
		event_id IS NOT NULL
	""", (coinc_def_id,))


def create_sim_burst_best_string_coinc_map(connection, coinc_def_id):
	"""
	Construct a sim_burst --> best matching coinc_event mapping for
	string cusp injections and coincs.
	"""
	# FIXME:  this hasn't finished being ported from the inspiral code
	connection.cursor().execute("""
CREATE TEMPORARY TABLE
	sim_burst_best_string_coinc_map
AS
	SELECT
		sim_burst.simulation_id AS simulation_id,
		(
			SELECT
				coinc_inspiral.coinc_event_id
			FROM
				coinc_event_map AS a
				JOIN coinc_event_map AS b ON (
					b.coinc_event_id == a.coinc_event_id
				)
				JOIN coinc_inspiral ON (
					b.table_name == 'coinc_event'
					AND b.event_id == coinc_inspiral.coinc_event_id
				)
			WHERE
				a.table_name == 'sim_burst'
				AND a.event_id == sim_burst.simulation_id
				AND coinc_event.coinc_def_id == ?
			ORDER BY
				(sngl_burst.chisq / sngl_burst.chisq_dof) / (sngl_burst.snr * sngl_burst.snr)
			LIMIT 1
		) AS coinc_event_id
	FROM
		sim_burst
	WHERE
		coinc_event_id IS NOT NULL
	""", (coinc_def_id,))


#
# =============================================================================
#
#                           Amplitudes in Instrument
#
# =============================================================================
#


def hrss_in_instrument(sim, instrument, offsetvector):
	"""
	Given an injection and an instrument, compute and return the h_rss
	of the injection as should be observed in the instrument.  That is,
	project the waveform onto the instrument, and return the root
	integrated strain squared.
	"""
	# FIXME:  this function is really only correct for sine-Gaussian
	# injections.  that's OK because I only quote sensitivities in
	# units of hrss when discussing sine-Gaussians.
	#
	# the problem is the following.  first,
	#
	#	h = F+ h+ + Fx hx
	#
	# so
	#
	#	h^{2} = F+^2 h+^2 + Fx^2 hx^2 + 2 F+ Fx h+ hx
	#
	# which means to calculate the hrss in the instrument you need to
	# know:  mean-square h in the + polarization, mean-square h in the
	# x polarization, and the correlation between the polarizations <h+
	# hx>.  these could be recorded in the sim_burst table, but they
	# aren't at present.

	# semimajor and semiminor axes of polarization ellipse

	a = 1.0 / math.sqrt(2.0 - sim.pol_ellipse_e**2)
	b = a * math.sqrt(1.0 - sim.pol_ellipse_e**2)

	# hrss in plus and cross polarizations

	hplusrss  = sim.hrss * (a * math.cos(sim.pol_ellipse_angle) - b * math.sin(sim.pol_ellipse_angle))
	hcrossrss = sim.hrss * (b * math.cos(sim.pol_ellipse_angle) + a * math.sin(sim.pol_ellipse_angle))

	# antenna response factors

	fplus, fcross = inject.XLALComputeDetAMResponse(
		inject.cached_detector[inject.prefix_to_name[instrument]].response,
		sim.ra,
		sim.dec,
		sim.psi,
		date.XLALGreenwichMeanSiderealTime(time_at_instrument(sim, instrument, offsetvector))
	)

	# hrss in detector

	return math.sqrt((fplus * hplusrss)**2 + (fcross * hcrossrss)**2)


def string_amplitude_in_instrument(sim, instrument, offsetvector):
	"""
	Given a string cusp injection and an instrument, compute and return
	the amplitude of the injection as should be observed in the
	instrument.
	"""
	assert sim.waveform == "StringCusp"

	# antenna response factors

	fplus, fcross = inject.XLALComputeDetAMResponse(
		inject.cached_detector[inject.prefix_to_name[instrument]].response,
		sim.ra,
		sim.dec,
		sim.psi,
		date.XLALGreenwichMeanSiderealTime(time_at_instrument(sim, instrument, offsetvector))
	)

	# amplitude in detector

	return fplus * sim.amplitude


#
# =============================================================================
#
#                             Efficiency Contours
#
# =============================================================================
#


#
# typical width of a string cusp signal's autocorrelation function when the
# signal is normalized to the interferometer noise (whitened).  this is
# used to provide a window around each injection for the purpose of
# identifying burst events that match the injection.  this number is
# O(1/f_{bucket}).
#


stringcusp_autocorrelation_width = .016	# seconds


#
# used to find burst events near injections for the purpose of producing a
# short list of coinc events for use in more costly comparisons
#


burst_is_near_injection_window = 2.0	# seconds


#
# h_{rss} vs. peak frequency
#


class Efficiency_hrss_vs_freq(object):
	def __init__(self, instruments, amplitude_func, amplitude_lbl, error):
		self.instruments = set(instruments)
		self.amplitude_func = amplitude_func
		self.amplitude_lbl = amplitude_lbl
		self.error = error
		self.injected_x = []
		self.injected_y = []
		self.found_x = []
		self.found_y = []


	def add_contents(self, contents):
		# FIXME:  the first left outer join can yield multiple
		# rows.
		cursor = contents.connection.cursor()
		for values in contents.connection.cursor().execute("""
SELECT
	sim_burst.*,
	coinc_event.coinc_event_id
FROM
	sim_burst
	-- The rest of this join can yield at most 1 row for each sim_burst
	-- row
	LEFT OUTER JOIN coinc_event_map ON (
		coinc_event_map.table_name == 'sim_burst'
		AND coinc_event_map.event_id == sim_burst.simulation_id
	)
	LEFT OUTER JOIN coinc_event ON (
		coinc_event.coinc_event_id == coinc_event_map.coinc_event_id
	)
WHERE
	coinc_event.coinc_def_id == ?
		""", (contents.sb_definer_id,)):
			sim = contents.sim_burst_table.row_from_cols(values)
			coinc_event_id = values[-1]
			instruments = set(cursor.execute("""
SELECT
	sngl_burst.ifo
FROM
	coinc_event_map
	JOIN sngl_burst ON (
		coinc_event_map.table_name == 'sngl_burst'
		AND coinc_event_map.event_id == sngl_burst.event_id
	)
WHERE
	coinc_event_map.coinc_event_id == ?
			""", (coinc_event_id,)))
			found = self.instruments.issubset(instruments)
			# FIXME:  this following assumes all injections are
			# done at zero lag (which is correct, for now, but
			# watch out for this)
			if injection_was_made(sim, contents.seglists, self.instruments):
				for instrument in self.instruments:
					amplitude = self.amplitude_func(sim, instrument)
					self.injected_x.append(sim.frequency)
					self.injected_y.append(amplitude)
					if found:
						self.found_x.append(sim.frequency)
						self.found_y.append(amplitude)
			elif found:
				print >>sys.stderr, "odd, injection %s was found in %s but not injected..." % (sim.simulation_id, "+".join(self.instruments))

	def _bin_events(self, binning = None):
		# called internally by finish()
		if binning is None:
			minx, maxx = min(self.injected_x), max(self.injected_x)
			miny, maxy = min(self.injected_y), max(self.injected_y)
			binning = rate.NDBins((rate.LogarithmicBins(minx, maxx, 256), rate.LogarithmicBins(miny, maxy, 256)))

		self.efficiency = rate.BinnedRatios(binning)

		for xy in zip(self.injected_x, self.injected_y):
			self.efficiency.incdenominator(xy)
		for xy in zip(self.found_x, self.found_y):
			self.efficiency.incnumerator(xy)

		# 1 / error^2 is the number of injections that need to be
		# within the window in order for the fractional uncertainty
		# in that number to be = error.  multiplying by
		# bins_per_inj tells us how many bins the window needs to
		# cover, and taking the square root translates that into
		# the window's length on a side in bins.  because the
		# contours tend to run parallel to the x axis, the window
		# is dilated in that direction to improve resolution.

		bins_per_inj = self.efficiency.used() / float(len(self.injected_x))
		self.window_size_x = self.window_size_y = math.sqrt(bins_per_inj / self.error**2)
		self.window_size_x *= math.sqrt(2)
		self.window_size_y /= math.sqrt(2)
		if self.window_size_x > 100 or self.window_size_y > 100:
			# program will take too long to run
			raise ValueError, "smoothing filter too large (not enough injections)"

		print >>sys.stderr, "The smoothing window for %s is %g x %g bins" % ("+".join(self.instruments), self.window_size_x, self.window_size_y),
		print >>sys.stderr, "which is %g%% x %g%% of the binning" % (100.0 * self.window_size_x / binning[0].n, 100.0 * self.window_size_y / binning[1].n)

	def finish(self, binning = None):
		# compute the binning if needed, and set the injections
		# into the numerator and denominator bins.  also compute
		# the smoothing window's parameters.

		self._bin_events(binning)

		# smooth the efficiency data.
		print >>sys.stderr, "Sum of numerator bins before smoothing = %g" % self.efficiency.numerator.array.sum()
		print >>sys.stderr, "Sum of denominator bins before smoothing = %g" % self.efficiency.denominator.array.sum()
		rate.filter_binned_ratios(self.efficiency, rate.gaussian_window2d(self.window_size_x, self.window_size_y))
		print >>sys.stderr, "Sum of numerator bins after smoothing = %g" % self.efficiency.numerator.array.sum()
		print >>sys.stderr, "Sum of denominator bins after smoothing = %g" % self.efficiency.denominator.array.sum()

		# regularize to prevent divide-by-zero errors
		self.efficiency.regularize()


def plot_Efficiency_hrss_vs_freq(efficiency):
	"""
	Generate a plot from an Efficiency_hrss_vs_freq instance.
	"""
	fig, axes = SnglBurstUtils.make_burst_plot("Frequency (Hz)", efficiency.amplitude_lbl)
	axes.loglog()

	xcoords, ycoords = efficiency.efficiency.centres()
	zvals = efficiency.efficiency.ratio()
	cset = axes.contour(xcoords, ycoords, numpy.transpose(zvals), (.1, .2, .3, .4, .5, .6, .7, .8, .9))
	cset.clabel(inline = True, fontsize = 5, fmt = r"$%%g \pm %g$" % efficiency.error, colors = "k")
	axes.set_title(r"%s Injection Detection Efficiency (%d of %d Found)" % ("+".join(sorted(efficiency.instruments)), len(efficiency.found_x), len(efficiency.injected_x)))
	return fig


#
# =============================================================================
#
#                             Source Co-ordinates
#
# =============================================================================
#


#
# Location of the galactic core
#	ra = 27940.04 s = 7 h 45 m 40.04 s
#	dec = -29o 00' 28.1"
#


MW_CENTER_J2000_RA_RAD = 2.0318570464121519
MW_CENTER_J2000_DEC_RAD = -0.50628171572274738

