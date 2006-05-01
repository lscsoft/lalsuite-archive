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
# Preamble
#
# =============================================================================
#

"""
This module provides code to convolve impulsive events with a window
function to produce a smooth data set.  The convolution is
integral-preserving, and so is particularly-well suited for use in
eliminating binning artifacts from histograms, and other rate-like data
sets.
"""

import math
import numarray
from numarray import convolve

__author__ = "Kipp Cannon <kipp@gravity.phys.uwm.edu>"
__version__ = "$Revision$"[11:-2]
__date__ = "$Date$"[7:-2]


#
# =============================================================================
#
#                                   Windows
#
# =============================================================================
#

def gaussian_window(halfwidth):
	"""
	Generate a normalized (integral = 1) Gaussian window.
	"""
	bins_per_unit = 10.0 / halfwidth
	return numarray.exp(-numarray.arrayrange(-10.0 * halfwidth, +10.0 * halfwidth, 1.0/bins_per_unit)**2.0 / (2.0 * halfwidth**2.0)) / math.sqrt(2.0 * math.pi) / halfwidth


def tophat_window(halfwidth):
	"""
	Generate a normalized (integral = 1) top-hat window.
	"""
	bins_per_unit = 10.0 / halfwidth
	return numarray.ones(2.0 * halfwidth / bins_per_unit) / (2.0 * halfwidth)


#
# =============================================================================
#
#                                     1-D
#
# =============================================================================
#

class Rate1D(object):
	"""
	An object for binning and smoothing impulsive data.
	"""
	def __init__(self, segment, halfwidth, windowfunc = gaussian_window):
		"""
		Initialize the bins for the given segment and width.
		"""
		self.halfwidth = halfwidth
		self.bins_per_unit = 10.0 / halfwidth
		self.start = segment[0]
		self.xvals = numarray.arrayrange(0.0, float(segment.duration()), 1.0/self.bins_per_unit) + float(segment[0])
		self.yvals = numarray.zeros(len(self.xvals), "Float32")
		self.set_window(windowfunc)

	def bin(self, x):
		"""
		Return the index for the bin corresponding to x.
		"""
		return int((x - self.start) * self.bins_per_unit)

	def __getitem__(self, x):
		"""
		Retrieve the weight in bin corresponding to x.
		"""
		return self.yvals[self.bin(x)]

	def __setitem__(self, x, weight):
		"""
		Add weight to the bin corresponding to x.
		"""
		self.yvals[self.bin(x)] += weight

	def set_window(self, windowfunc):
		"""
		Set the window function.
		"""
		self.window = windowfunc

	def convolve(self):
		"""
		Convolve the binned weights with the window to smooth the
		data set.
		"""
		self.yvals = convolve.convolve(self.yvals, self.window(self.halfwidth), mode=convolve.SAME)
		return self


#
# =============================================================================
#
#                                  Utilities
#
# =============================================================================
#

def smooth(impulses, segment, width, weights = None):
	"""
	Given a list of the co-ordinates of impulsive events, an optional
	matching list of weights for those events, a segment and a
	characteristic length, this function returns two arrays:  an array
	of x co-ordinates, and an array of corresponding y co-ordinates
	describing the rate of impulse events.

	The y values returned carry units of "weight per unit of x".  The
	default weight of each impulse is 1.0 if weights are not provided.

	Example:
		# generate regular impulses at a rate of 2 Hz from 10 s to
		# 90 s.
		impulses = numarray.arrayrange(10, 90, 0.5)

		# using the impulses, generate rate data in the interval
		# 0 s to 100 s with a 1 s Gaussian window.
		x, y = smooth(impulses, segment(0, 100), 1)

		# plot
		plot(x, y)
	"""
	rate = Rate1D(segment, width / 2.0)
	if weights != None:
		for n, x in enumerate(impulses):
			if segment[0] <= x < segment[1]:
				rate[x] = weights[n]
	else:
		for x in impulses:
			if segment[0] <= x < segment[1]:
				rate[x] = 1.0
	rate.convolve()
	return rate.xvals, rate.yvals

