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
#                                   Preamble
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
#                                     Bins
#
# =============================================================================
#

class _Bins(object):
	def __init__(self, min, max, n):
		if type(n) != int:
			raise TypeError, n
		if max <= min:
			raise ValueError, (min, max)
		self.min = min
		self.max = max
		self.n = n
		self._set_delta(min, max, n)

	def _set_delta(self, min, max, n):
		raise NotImplementedError

	def __getitem__(self, x):
		raise NotImplementedError


class _LinBins(_Bins):
	def _set_delta(self, min, max, n):
		self.delta = float(max - min) / (n - 1)

	def __getitem__(self, x):
		if self.min <= x <= self.max:
			return int((x - self.min) / self.delta + 0.5)
		raise IndexError, x

	def centres(self):
		return self.min + self.delta * numarray.arange(0, self.n, 1, "Float64")


class _LogBins(_Bins):
	def _set_delta(self, min, max, n):
		self.delta = math.log(float(max / min) ** (1.0 / (n - 1)))

	def __getitem__(self, x):
		if self.min <= x <= self.max:
			return int(math.log(x / self.min) / self.delta + 0.5)
		raise IndexError, x

	def centres(self):
		return self.min * numarray.exp(self.delta * numarray.arange(0, self.n, 1, "Float64"))


class Bins(object):
	def __init__(self, *args, **kwargs):
		if len(args) % 3:
			raise TypeError, "arguments must be min,max,n[,min,max,n]..."
		spacing = kwargs.get("spacing", ["lin"] * (len(args) / 3))
		if len(spacing) != len(args) / 3:
			raise ValueError, spacing
		self.__bins = []
		try:
			it = iter(args)
			spacing = iter(spacing)
			while True:
				s = spacing.next()
				if s == "lin":
					self.__bins.append(_LinBins(it.next(), it.next(), it.next()))
				elif s == "log":
					self.__bins.append(_LogBins(it.next(), it.next(), it.next()))
				else:
					raise ValueError, s
		except StopIteration:
			pass
		self.shape = tuple([b.n for b in self.__bins])

	def __getitem__(self, coords):
		return tuple([self.__bins[i][coords[i]] for i in xrange(len(self.__bins))])

	def centres(self):
		return tuple([self.__bins[i].centres() for i in xrange(len(self.__bins))])


#
# =============================================================================
#
#                                 Binned Array
#
# =============================================================================
#

class BinnedArray(object):
	# Argh.  Ideally, we would subclass numarray.NumArray and construct
	# a variant that accepts arbitrary co-ordinates in the
	# __getitem__() and __setitem__() methods.  Unforunately, numarray
	# is somewhat poorly designed, and so subclasses that do not export
	# *exactly* the same interface as that exported by
	# numarray.NumArray break the library.  If you can't construct a
	# subclass that modifies the interface in any way whatsoever, that
	# is if all subclasses must behave in exactly the same way as the
	# parent class, then what would ever be the point?
	def __init__(self, bins):
		self.bins = bins
		self.array = numarray.zeros(bins.shape, "Float64")

	def __getitem__(self, coords):
		return self.array[self.bins[coords]]

	def __setitem__(self, coords, val):
		self.array[self.bins[coords]] = val

	def centres(self):
		return self.bins.centres()


class BinnedRatios(object):
	def __init__(self, bins):
		self.bins = bins
		self.numerator = numarray.zeros(bins.shape, "Float64")
		self.denominator = numarray.zeros(bins.shape, "Float64")

	def incnumerator(self, coords, weight = 1.0):
		self.numerator[self.bins[coords]] += weight

	def incdenominator(self, coords, weight = 1.0):
		self.denominator[self.bins[coords]] += weight

	def regularize(self):
		self.denominator = numarray.where(self.denominator > 0, self.denominator, 1.0)

	def centres(self):
		return self.bins.centres()


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
	def __init__(self, segment, width, windowfunc = gaussian_window):
		"""
		Initialize the bins for the given segment and width.
		"""
		self.halfwidth = width / 2.0
		self.xvals = numarray.arrayrange(0.0, float(segment.duration()) + self.halfwidth / 10.0, self.halfwidth / 10.0) + float(segment[0])
		self.__yvals = BinnedArray(Bins(segment[0], segment[1], len(self.xvals)))
		self.yvals = self.__yvals.array
		self.set_window(windowfunc)

	def __getitem__(self, x):
		"""
		Retrieve the weight in bin corresponding to x.
		"""
		return self.__yvals[x,]

	def __setitem__(self, x, weight):
		"""
		Add weight to the bin corresponding to x.
		"""
		self.__yvals[x,] += weight

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
	rate = Rate1D(segment, width)
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

