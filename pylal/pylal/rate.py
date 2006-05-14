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
This module provides facilities for handling impulsive events.  A number
multi-dimensional binning functions are provided, as well as code to
convolve binned data with integral-preserving window functions to produce
smoothed representations of data sets.  This is particularly well suited
for use in computing moving-average rate data from collections of impulsive
events, elliminating binning artifacts from histograms, and improving the
appearance of contour plots.
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
#                                     Bins
#
# =============================================================================
#

class _Bins(object):
	"""
	Parent class for 1-D bins.  For internal use only.
	"""
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
	"""
	Linearly-spaced 1-D bins.  For internal use only.
	"""
	def _set_delta(self, min, max, n):
		self.delta = float(max - min) / (n - 1)

	def __getitem__(self, x):
		if self.min <= x <= self.max:
			return int((x - self.min) / self.delta + 0.5)
		raise IndexError, x

	def centres(self):
		return self.min + self.delta * numarray.arange(0, self.n, 1, "Float64")


class _LogBins(_Bins):
	"""
	Logarithmically-spaced 1-D bins.  For internal use only.
	"""
	def _set_delta(self, min, max, n):
		self.delta = math.log(float(max / min) ** (1.0 / (n - 1)))

	def __getitem__(self, x):
		if self.min <= x <= self.max:
			return int(math.log(x / self.min) / self.delta + 0.5)
		raise IndexError, x

	def centres(self):
		return self.min * numarray.exp(self.delta * numarray.arange(0, self.n, 1, "Float64"))


class Bins(object):
	"""
	Multi-dimensional co-ordinate binning.  An instance of this object
	is used to convert a tuple of co-ordinates into a tuple of array
	indeces, thereby allowing the contents of an array object to be
	accessed with real-valued coordinates.  When creating a Bins
	object, the arguments describe the range along each co-ordinate
	direction, and the number of bins.  The arguments come in groups of
	three, one for each direction.  An optional keyword argument is
	used to set the spacing of the bins to linear or logarithmic in
	each of the directions independantly.

	Example:
	
	>>> b = Bins(1, 25, 3, 1, 25, 3, spacing = ["lin", "log"])
	>>> b[1, 1]
	(0, 0)
	>>> b[1.5, 1]
	(0, 0)
	>>> b[10, 1]
	(1, 0)
	>>> b[1, 5]
	(0, 1)
	>>> b.centres()
	(array([  1.,  13.,  25.]), array([  1.,   5.,  25.]))
	"""
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
		"""
		Return a tuple of arrays containing the bin centres for
		each dimension.
		"""
		return tuple([self.__bins[i].centres() for i in xrange(len(self.__bins))])


#
# =============================================================================
#
#                                 Binned Array
#
# =============================================================================
#

class BinnedArray(object):
	"""
	A convenience wrapper, using the Bins class to provide access to
	the elements of an array object.  Technical reasons preclude
	providing a subclass of the array object, so the array data is made
	available as the "array" attribute of this class.

	Example:

	>>> x = BinnedArray(Bins(0, 10, 5))
	>>> x.array
	array([ 0.,  0.,  0.,  0.,  0.])
	>>> x[0,] += 1
	>>> x[0.5,] += 1
	>>> x.array
	array([ 2.,  0.,  0.,  0.,  0.])

	Note that even for 1 dimensional array the index but be a tuple.
	"""
	def __init__(self, bins):
		self.bins = bins
		self.array = numarray.zeros(bins.shape, "Float64")

	def __getitem__(self, coords):
		return self.array[self.bins[coords]]

	def __setitem__(self, coords, val):
		self.array[self.bins[coords]] = val

	def centres(self):
		"""
		Return a tuple of arrays containing the bin centres for
		each dimension.
		"""
		return self.bins.centres()


class BinnedRatios(object):
	"""
	Like BinnedArray, but provides a numerator array and a denominator
	array.  The incnumerator() method increments a bin in the numerator
	by the given weight, and the incdenominator() method increments a
	bin in the denominator by the given weight.  There are no methods
	provided for setting or decrementing either, but the arrays can be
	accessed directly.
	"""
	def __init__(self, bins):
		self.bins = bins
		self.numerator = numarray.zeros(bins.shape, "Float64")
		self.denominator = numarray.zeros(bins.shape, "Float64")

	def incnumerator(self, coords, weight = 1.0):
		"""
		Add weight to the numerator bin at coords.
		"""
		self.numerator[self.bins[coords]] += weight

	def incdenominator(self, coords, weight = 1.0):
		"""
		Add weight to the denominator bin at coords.
		"""
		self.denominator[self.bins[coords]] += weight

	def ratio(self):
		"""
		Compute and return the array of ratios.
		"""
		return self.numerator / self.denominator

	def regularize(self):
		"""
		Find bins in the denominator that are 0, and set them to 1.
		Presumably the corresponding bin in the numerator is also
		0, so this has the effect of allowing the ratio array to be
		evaluated without error, returning zeros in those bins that
		have had no weight added to them.
		"""
		numarray.where(self.denominator > 0, self.denominator, 1.0, out = self.denominator)

	def logregularize(self, epsilon = 2**-1074):
		"""
		Find bins in the denominator that are 0, and set them to 1,
		while setting the corresponding bin in the numerator to
		float epsilon.  This has the effect of allowing the
		logarithm of the ratio array to be evaluated without error,
		return log(epsilon) to be returned.
		"""
		numarray.where(self.denominator > 0, self.numerator, epsilon, out = self.numerator)
		numarray.where(self.denominator > 0, self.denominator, 1.0, out = self.denominator)

	def centres(self):
		"""
		Return a tuple of arrays containing the bin centres for
		each dimension.
		"""
		return self.bins.centres()


#
# =============================================================================
#
#                                   Windows
#
# =============================================================================
#

def gaussian_window(bins):
	"""
	Generate a normalized (integral = 1) Gaussian window in 1
	dimension.  bins sets the width of the window in bin count.
	"""
	bins /= 2.0	# half-width
	return numarray.exp(-numarray.arrayrange(-6 * int(bins), 6 * int(bins) + 1, 1, "Float64")**2.0 / (2.0 * bins**2.0)) / math.sqrt(2.0 * math.pi) / bins


def gaussian_window2d(bins_x, bins_y):
	"""
	Generate a normalized (integral = 1) Gaussian window in 2
	dimensions.  bins_x and bins_y set the widths of the window in bin
	counts.
	"""
	return numarray.outerproduct(gaussian_window(bins_x), gaussian_window(bins_y))


def tophat_window(bins):
	"""
	Generate a normalized (integral = 1) top-hat window in 1 dimension.
	bins sets the width of the window in bin counts.
	"""
	return numarray.ones(bins, "Float64") / bins


def tophat_window2d(bins_x, bins_y):
	"""
	Generate a normalized (integral = 1) top-hat window in 2
	dimensions.  bins_x and bins_y set the widths of the window in bin
	counts.
	"""
	return numarray.outerproduct(tophat_window(bins_x), tophat_window(bins_y))


#
# =============================================================================
#
#                                  Filtering
#
# =============================================================================
#

def filter(binned, window):
	"""
	Filter the binned data using the window function.  The
	transformation is done in place.
	"""
	dims = len(binned.array.shape)
	if dims != len(window.shape):
		raise ValueError, "binned array and window dimensions mismatch"
	if 1 in map(int.__cmp__, window.shape, binned.array.shape):
		raise ValueError, "window data excedes size of array data"
	if dims == 1:
		binned.array = convolve.convolve(binned.array, window, mode = convolve.SAME)
	elif dims == 2:
		convolve.convolve2d(binned.array, window, output = binned.array, mode = "constant")
	else:
		raise ValueError, "can only filter 1 and 2 dimensional arrays"


def filter_ratios(binned, window):
	"""
	Filter the binned ratio data using the window function.  The
	transformation is done in place.
	"""
	dims = len(binned.numerator.shape)
	if dims != len(window.shape):
		raise ValueError, "binned array and window dimensions mismatch"
	if 1 in map(int.__cmp__, window.shape, binned.array.shape):
		raise ValueError, "window data excedes size of array data"
	if dims == 1:
		binned.numerator = convolve.convolve2d(binned.numerator, window, mode = "constant")
		binned.denominator = convolve.convolve2d(binned.denominator, window, mode = "constant")
	elif dims == 2:
		binned.numerator = convolve.convolve2d(binned.numerator, window, mode = "constant")
		binned.denominator = convolve.convolve2d(binned.denominator, window, mode = "constant")
	else:
		raise ValueError, "can only filter 1 and 2 dimensional arrays"


#
# =============================================================================
#
#                                    Rates
#
# =============================================================================
#

class Rate(BinnedArray):
	"""
	An object for binning and smoothing impulsive data in 1 dimension,
	normalized so as to measure events (or event weight) per filter
	width.
	"""
	bins_per_filterwidth = 20
	def __init__(self, segment, filterwidth, windowfunc = gaussian_window):
		"""
		Initialize the bins for the given segment and filter width.
		"""
		self.filterwidth = filterwidth
		BinnedArray.__init__(self, Bins(segment[0], segment[1], int(segment.duration() / filterwidth * self.bins_per_filterwidth + 1) ) )
		self.set_window(windowfunc)

	def __setitem__(self, x, weight):
		"""
		Add weight to the bin corresponding to x.
		"""
		self.array[self.bins[x,]] += weight

	def set_window(self, windowfunc):
		"""
		Set the window function.
		"""
		self.window = windowfunc(self.bins_per_filterwidth) / self.filterwidth * self.bins_per_filterwidth

	def filter(self):
		"""
		Convolve the binned weights with the window to smooth the
		data set.
		"""
		filter(self, self.window)
		return self

	def xvals(self):
		return self.centres()[0]

	def yvals(self):
		return self.array
