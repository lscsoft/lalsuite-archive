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
This module provides facilities for studying impulsive events.  A number of
multi-dimensional binning functions are provided, as well as code to
convolve binned data with integral- and phase-preserving window functions
to produce smoothed representations of data sets.  This is particularly
well suited for use in computing moving-average rate data from collections
of impulsive events, elliminating binning artifacts from histograms, and
smoothing contour plots.
"""


import bisect
import math
import numpy
from scipy.signal import signaltools


from glue import iterutils
from glue import segments
from pylal import window


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


class Bins(object):
	"""
	Parent class for 1-dimensional binnings.  This class is not
	intended to be used directly, but to be subclassed for use in real
	bins classes.
	"""
	def __init__(self, min, max, n):
		"""
		Initialize a Bins instance.  The three arguments are the
		minimum and maximum of the values spanned by the bins, and
		the number of bins to place between them.  Subclasses may
		require additional arguments, or different arguments
		altogether.
		"""
		# convenience code to do some common initialization and
		# input checking
		if not isinstance(n, int):
			raise TypeError, n
		if max <= min:
			raise ValueError, (min, max)
		self.min = min
		self.max = max
		self.n = n

	def __len__(self):
		return self.n

	def __cmp__(self, other):
		"""
		Two binnings are the same if they are instances of the same
		class, have the same lower and upper bounds, and the same
		count of bins.
		"""
		if not isinstance(other, type(self)):
			return -1
		return cmp((type(self), self.min, self.max, len(self)), (type(other), other.min, other.max, len(other)))

	def __getitem__(self, x):
		"""
		Convert a co-ordinate into a bin index.  The co-ordinate
		can be a single number, or a Python slice instance.  If a
		single number is given, it is mapped to the bin in which it
		falls.  If a slice is given, it is converted to a slice
		whose upper and lower bounds are the bins in which the
		input slice's upper and lower bounds fall.
		"""
		raise NotImplementedError

	def lower(self):
		"""
		Return an array containing the locations of the lower
		boundaries of the bins.
		"""
		raise NotImplementedError

	def centres(self):
		"""
		Return an array containing the locations of the bin
		centres.
		"""
		raise NotImplementedError

	def upper(self):
		"""
		Return an array containing the locations of the upper
		boundaries of the bins.
		"""
		raise NotImplementedError


class LinearBins(Bins):
	"""
	Linearly-spaced bins.  There are n bins of equal size, the first
	bin starts on the lower bound and the last bin ends on the upper
	bound inclusively.

	Example:

	>>> x = LinearBins(1.0, 25.0, 3)
	>>> x[1]
	0
	>>> x[1.5]
	0
	>>> x[10]
	1
	>>> x[25]
	2
	"""
	def __init__(self, min, max, n):
		Bins.__init__(self, min, max, n)
		self.delta = float(max - min) / n

	def __getitem__(self, x):
		if isinstance(x, slice):
			if x.step is not None:
				raise NotImplementedError, x
			if x.start is None:
				start = 0
			else:
				start = self[x.start]
			if x.stop is None:
				stop = len(self)
			else:
				stop = self[x.stop]
			return slice(start, stop)
		if self.min <= x < self.max:
			return int(math.floor((x - self.min) / self.delta))
		if x == self.max:
			# special "measure zero" corner case
			return len(self) - 1
		raise IndexError, x

	def lower(self):
		return self.min + self.delta * numpy.arange(len(self))

	def centres(self):
		return self.min + self.delta * (numpy.arange(len(self)) + 0.5)

	def upper(self):
		return self.min + self.delta * (numpy.arange(len(self)) + 1)


class LogarithmicBins(Bins):
	"""
	Logarithmically-spaced bins.  There are n bins, each of whose upper
	and lower bounds differ by the same factor.  The first bin starts
	on the lower bound, and the last bin ends on the upper bound
	inclusively.

	Example:

	>>> x = LogarithmicBins(1.0, 25.0, 3)
	>>> x[1]
	0
	>>> x[5]
	1
	>>> x[25]
	2
	"""
	def __init__(self, min, max, n):
		Bins.__init__(self, min, max, n)
		self.delta = math.log(float(max / min)) / n

	def __getitem__(self, x):
		if isinstance(x, slice):
			if x.step is not None:
				raise NotImplementedError, x
			if x.start is None:
				start = 0
			else:
				start = self[x.start]
			if x.stop is None:
				stop = len(self)
			else:
				stop = self[x.stop]
			return slice(start, stop)
		if self.min <= x < self.max:
			return int(math.floor(math.log(x / self.min) / self.delta))
		if x == self.max:
			# special "measure zero" corner case
			return len(self) - 1
		raise IndexError, x

	def lower(self):
		return self.min * numpy.exp(self.delta * numpy.arange(len(self)))

	def centres(self):
		return self.min * numpy.exp(self.delta * (numpy.arange(len(self)) + 0.5))

	def upper(self):
		return self.min * numpy.exp(self.delta * (numpy.arange(len(self)) + 1))


class ATanBins(Bins):
	"""
	Bins spaced uniformly in tan^-1 x.  Provides approximately linear
	binning in the middle portion, with the bin density dropping
	asymptotically to 0 as x goes to +/- \infty.  The min and max
	parameters set the bounds of the region of approximately
	uniformly-spaced bins.  In a sense, these are where the roll-over
	from uniformly-spaced bins to asymptotically diminishing bin
	density occurs.  There is a total of n bins.

	Example:

	>>> x = ATanBins(-1.0, +1.0, 11)
	>>> x[float("-inf")]
	0
	>>> x[0]
	5
	>>> x[float("+inf")]
	10
	>>> x.centres()
	array([-4.42778777, -1.39400285, -0.73469838, -0.40913068, -0.18692843,
	        0.        ,  0.18692843,  0.40913068,  0.73469838,  1.39400285,
                4.42778777])
	"""
	def __init__(self, min, max, n):
		Bins.__init__(self, min, max, n)
		self.mid = (min + max) / 2.0
		self.scale = math.pi / float(max - min)
		self.delta = 1.0 / n

	def __getitem__(self, x):
		if isinstance(x, slice):
			if x.step is not None:
				raise NotImplementedError, x
			if x.start is None:
				start = 0
			else:
				start = self[x.start]
			if x.stop is None:
				stop = len(self)
			else:
				stop = self[x.stop]
			return slice(start, stop)
		# map to the domain [0, 1]
		x = math.atan(float(x - self.mid) * self.scale) / math.pi + 0.5
		if x < 1:
			return int(math.floor(x / self.delta))
		# x == 1, special "measure zero" corner case
		return len(self) - 1

	def lower(self):
		x = numpy.tan(-math.pi / 2 + math.pi * self.delta * numpy.arange(len(self))) / self.scale + self.mid
		x[0] = float("-inf")
		return x

	def centres(self):
		return numpy.tan(-math.pi / 2 + math.pi * self.delta * (numpy.arange(len(self)) + 0.5)) / self.scale + self.mid

	def upper(self):
		x = numpy.tan(-math.pi / 2 + math.pi * self.delta * (numpy.arange(len(self)) + 1)) / self.scale + self.mid
		x[-1] = float("+inf")
		return x

class IrregularBins(Bins):
	"""
	Bins with arbitrary, irregular spacing.  We only require strict
	monotonicity of the bin boundaries.  N boundaries define N-1 bins.

	Example:

	>>> x = IrregularBins([0.0, 11.0, 15.0, numpy.inf])
	>>> len(x)
	3
	>>> x[1]
	0
	>>> x[1.5]
	0
	>>> x[13]
	1
	>>> x[25]
	2
	"""
	def __init__(self, boundaries):
		"""
		Initialize a set of custom bins with the bin boundaries.
		This includes all left edges plus the right edge.  The
		boundaries must be monotonic and there must be at least two
		elements.
		"""
		# check pre-conditions
		if len(boundaries) < 2:
			raise ValueError, "less than two boundaries provided"
		boundaries = numpy.array(boundaries)
		if (boundaries[:-1] > boundaries[1:]).any():
			raise ValueError, "non-monotonic boundaries provided"

		self.boundaries = boundaries
		self.n = len(boundaries) - 1
		self.min = boundaries[0]
		self.max = boundaries[-1]

	def __cmp__(self, other):
		"""
		Two binnings are the same if they are instances of the same
		class, and have the same boundaries.
		"""
		if not isinstance(other, type(self)):
			return -1
		return cmp(len(self), len(other)) or (self.boundaries != other.boundaries).any()

	def __getitem__(self, x):
		if isinstance(x, slice):
			if x.step is not None:
				raise NotImplementedError, x
			if x.start is None:
				start = 0
			else:
				start = self[x.start]
			if x.stop is None:
				stop = len(self)
			else:
				stop = self[x.stop]
			return slice(start, stop)
		if x < self.min or x > self.max:
			raise IndexError, x
		# special measure-zero edge case
		if x == self.max:
			return len(self.boundaries) - 2
		return bisect.bisect_right(self.boundaries, x) - 1

	def lower(self):
		return self.boundaries[:-1]

	def upper(self):
		return self.boundaries[1:]

	def centres(self):
		return (self.lower() + self.upper()) / 2

class NDBins(tuple):
	"""
	Multi-dimensional co-ordinate binning.  An instance of this object
	is used to convert a tuple of co-ordinates into a tuple of bin
	indices.  This can be used to allow the contents of an array object
	to be accessed with real-valued coordinates.

	NDBins is a subclass of the tuple builtin, and is initialized with
	an iterable of instances of subclasses of Bins.  Each Bins subclass
	instance describes the binning to apply in the corresponding
	co-ordinate direction, and the number of them sets the dimensions
	of the binning.

	Example:
	
	>>> x = NDBins((LinearBins(1, 25, 3), LogarithmicBins(1, 25, 3)))
	>>> x[1, 1]
	(0, 0)
	>>> x[1.5, 1]
	(0, 0)
	>>> x[10, 1]
	(1, 0)
	>>> x[1, 5]
	(0, 1)
	>>> x[1, 1:5]
	(0, slice(0, 1, None))
	>>> x.centres()
	(array([  5.,  13.,  21.]), array([  1.70997595,   5.,  14.62008869]))

	Note that the co-ordinates to be converted must be a tuple, even if
	it is only a 1-dimensional co-ordinate.
	"""
	def __init__(self, *args):
		tuple.__init__(self, *args)
		self.min = tuple([b.min for b in self])
		self.max = tuple([b.max for b in self])
		self.shape = tuple([len(b) for b in self])

	def __getitem__(self, coords):
		"""
		When coords is a tuple, it is interpreted as an
		N-dimensional co-ordinate which is converted to an N-tuple
		of bin indices by the Bins instances in this object.
		Otherwise coords is interpeted as an index into the tuple,
		and the corresponding Bins instance is returned.

		Example:

		>>> x[1, 1]
		(0, 0)
		>>> x[1]
		<pylal.rate.LinearBins object at 0xb5cfa9ac>

		When used to convert co-ordinates to bin indices, each
		co-ordinate can be anything the corresponding Bins instance
		will accept.  Note that the co-ordinates to be converted
		must be a tuple, even if it is only a 1-dimensional
		co-ordinate.
		"""
		if isinstance(coords, tuple):
			if len(coords) != len(self):
				raise ValueError, "dimension mismatch"
			return tuple(map(lambda b, c: b[c], self, coords))
		else:
			return tuple.__getitem__(self, coords)

	def lower(self):
		"""
		Return a tuple of arrays, where each array contains the
		locations of the lower boundaries of the bins in the
		corresponding dimension.
		"""
		return tuple([b.lower() for b in self])

	def centres(self):
		"""
		Return a tuple of arrays, where each array contains the
		locations of the bin centres for the corresponding
		dimension.
		"""
		return tuple([b.centres() for b in self])

	def upper(self):
		"""
		Return a tuple of arrays, where each array contains the
		locations of the upper boundaries of the bins in the
		corresponding dimension.
		"""
		return tuple([b.upper() for b in self])

	def volumes(self):
		"""
		Return an n-dimensional array of the bin volumes.
		"""
		# FIXME:  works for 1 and 2 dimensions, but goes funny
		# after that.  always returns a 2 dimensional array no
		# matter what...
		bounds = iter(zip(self.upper(), self.lower()))
		u, l = bounds.next()
		volumes = u - l
		for u, l in bounds:
			volumes = numpy.outer(volumes, u - l)
		return volumes


#
# =============================================================================
#
#                              Segments and Bins
#
# =============================================================================
#


def bins_spanned(bins, seglist, dtype = "double"):
	"""
	Input is a Bins subclass instance and a glue.segments.segmentlist
	instance.  The output is an array object the length of the binning,
	which each element in the array set to the interval in the
	corresponding bin spanned by the segment list.

	Example:

	>>> from glue.segments import *
	>>> s = segmentlist([segment(1.5, 10.333), segment(15.8, 24)])
	>>> b = LinearBins(0, 30, 100)
	>>> bins_spanned(b, s)
	array([ 0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.3  ,  0.3  ,  0.3  ,
	        0.3  ,  0.3  ,  0.3  ,  0.3  ,  0.3  ,  0.3  ,  0.3  ,  0.3  ,
	        0.3  ,  0.3  ,  0.3  ,  0.3  ,  0.3  ,  0.3  ,  0.3  ,  0.3  ,
	        0.3  ,  0.3  ,  0.3  ,  0.3  ,  0.3  ,  0.3  ,  0.3  ,  0.3  ,
	        0.3  ,  0.3  ,  0.133,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,
	        0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,
	        0.   ,  0.   ,  0.   ,  0.   ,  0.1  ,  0.3  ,  0.3  ,  0.3  ,
	        0.3  ,  0.3  ,  0.3  ,  0.3  ,  0.3  ,  0.3  ,  0.3  ,  0.3  ,
	        0.3  ,  0.3  ,  0.3  ,  0.3  ,  0.3  ,  0.3  ,  0.3  ,  0.3  ,
	        0.3  ,  0.3  ,  0.3  ,  0.3  ,  0.3  ,  0.3  ,  0.3  ,  0.3  ,
	        0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,
	        0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,
	        0.   ,  0.   ,  0.   ,  0.   ])
	"""
	array = numpy.zeros((len(bins),), dtype = dtype)
	for i, (a, b) in enumerate(zip(bins.lower(), bins.upper())):
		array[i] = abs(seglist & segments.segmentlist([segments.segment(a, b)]))
	return array


#
# =============================================================================
#
#                                 Binned Array
#
# =============================================================================
#


class BinnedArray(object):
	"""
	A convenience wrapper, using the NDBins class to provide access to
	the elements of an array object.  Technical reasons preclude
	providing a subclass of the array object, so the array data is made
	available as the "array" attribute of this class.

	Example:

	>>> x = BinnedArray(NDBins((LinearBins(0, 10, 5),)))
	>>> x.array
	array([ 0.,  0.,  0.,  0.,  0.])
	>>> x[0,] += 1
	>>> x[0.5,] += 1
	>>> x.array
	array([ 2.,  0.,  0.,  0.,  0.])

	Note that even for 1 dimensional arrays the index must be a tuple.
	"""
	def __init__(self, bins, dtype = "double"):
		self.bins = bins
		self.array = numpy.zeros(bins.shape, dtype = dtype)

	def __getitem__(self, coords):
		return self.array[self.bins[coords]]

	def __setitem__(self, coords, val):
		self.array[self.bins[coords]] = val

	def __len__(self):
		return len(self.array)

	def __iadd__(self, other):
		"""
		Add the contents of another BinnedArray object to this one.
		It is not necessary for the binnings to be identical, but
		an integer number of the bins in other must fit into each
		bin in self.
		"""
		# identical binning? (fast path)
		if not cmp(self.bins, other.bins):
			self.array += other.array
			return self
		# can other's bins be put into ours?
		if self.bins.min != other.bins.min or self.bins.max != other.bins.max or False in map(lambda a, b: (b % a) == 0, self.bins.shape, other.bins.shape):
			raise TypeError, "incompatible binning: %s" % repr(other)
		for coords in iterutils.MultiIter(*other.bins.centres()):
			self[coords] += other[coords]
		return self

	def centres(self):
		"""
		Return a tuple of arrays containing the bin centres for
		each dimension.
		"""
		return self.bins.centres()

	def to_density(self):
		"""
		Divide each bin's value by the volume of the bin.
		"""
		self.array /= self.bins.volumes()

	def logregularize(self, epsilon = 2**-1074):
		"""
		Find bins <= 0, and set them to epsilon, This has the
		effect of allowing the logarithm of the array to be
		evaluated without error.
		"""
		self.array[self.array <= 0] = epsilon
		return self


class BinnedRatios(object):
	"""
	Like BinnedArray, but provides a numerator array and a denominator
	array.  The incnumerator() method increments a bin in the numerator
	by the given weight, and the incdenominator() method increments a
	bin in the denominator by the given weight.  There are no methods
	provided for setting or decrementing either, but the they are
	accessible as the numerator and denominator attributes, which are
	both BinnedArray objects.
	"""
	def __init__(self, bins, dtype = "double"):
		self.numerator = BinnedArray(bins, dtype = dtype)
		self.denominator = BinnedArray(bins, dtype = dtype)

	def bins(self):
		return self.numerator.bins

	def __iadd__(self, other):
		"""
		Add the weights from another BinnedRatios object's
		numerator and denominator to the numerator and denominator
		of this one.  Note that this is not the same as adding the
		ratios.  It is not necessary for the binnings to be
		identical, but an integer number of the bins in other must
		fit into each bin in self.
		"""
		try:
			self.numerator += other.numerator
			self.denominator += other.denominator
		except TypeError:
			raise TypeError, "incompatible binning: %s" % repr(other)
		return self

	def incnumerator(self, coords, weight = 1):
		"""
		Add weight to the numerator bin at coords.
		"""
		self.numerator[coords] += weight

	def incdenominator(self, coords, weight = 1):
		"""
		Add weight to the denominator bin at coords.
		"""
		self.denominator[coords] += weight

	def ratio(self):
		"""
		Compute and return the array of ratios.
		"""
		return self.numerator.array / self.denominator.array

	def regularize(self):
		"""
		Find bins in the denominator that are 0, and set them to 1.
		Presumably the corresponding bin in the numerator is also
		0, so this has the effect of allowing the ratio array to be
		evaluated without error, returning zeros in those bins that
		have had no weight added to them.
		"""
		self.denominator.array[self.denominator.array == 0] = 1
		return self

	def logregularize(self, epsilon = 2**-1074):
		"""
		Find bins in the denominator that are 0, and set them to 1,
		while setting the corresponding bin in the numerator to
		float epsilon.  This has the effect of allowing the
		logarithm of the ratio array to be evaluated without error.
		"""
		self.numerator.array[self.denominator.array == 0] = epsilon
		self.denominator.array[self.denominator.array == 0] = 1
		return self

	def centres(self):
		"""
		Return a tuple of arrays containing the bin centres for
		each dimension.
		"""
		return self.numerator.bins.centres()

	def used(self):
		"""
		Return the number of bins with non-zero denominator.
		"""
		return numpy.sum(self.denominator.array != 0)


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
	dimension.  bins sets the width of the window in bin count.  The
	vector is created with an odd number of samples, and the Gaussian
	is peaked on the middle sample.
	"""
	if bins <= 0:
		raise ValueError, bins
	l = 10 * int(bins / 2.0)
	w = window.XLALCreateGaussREAL8Window(l + 1, l / float(bins))
	return w.data / w.sum


def gaussian_window2d(bins_x, bins_y):
	"""
	Generate a normalized (integral = 1) Gaussian window in 2
	dimensions.  bins_x and bins_y set the widths of the window in bin
	counts.
	"""
	return numpy.outer(gaussian_window(bins_x), gaussian_window(bins_y))


def tophat_window(bins):
	"""
	Generate a normalized (integral = 1) top-hat window in 1 dimension.
	bins sets the width of the window in bin counts, which is rounded
	up to the nearest odd integer.
	"""
	if bins <= 0:
		raise ValueError, bins
	w = window.XLALCreateRectangularREAL8Window(int(bins / 2.0) * 2 + 1)
	return w.data / w.sum


def tophat_window2d(bins_x, bins_y):
	"""
	Generate a normalized (integral = 1) top-hat window in 2
	dimensions.  bins_x and bins_y set the widths of the window in bin
	counts, which are both rounded up to the nearest odd integer.  The
	result is a rectangular array, with an elliptical pattern of
	elements set to a constant value centred on the array's mid-point,
	and all other elements set to 0.
	"""
	if bins_x <= 0:
		raise ValueError, bins_x
	if bins_y <= 0:
		raise ValueError, bins_y

	# This might appear to be using a screwy, slow, algorithm but it's
	# the only way I have found to get a window with the correct bins
	# set and cleared as appropriate.  I'd love this to be replaced by
	# something that's easier to know is correct.

	# fill rectangle with ones, making the number of bins odd in each
	# direction
	window = numpy.ones((int(bins_x / 2.0) * 2 + 1, int(bins_y / 2.0) * 2 + 1), "Float64")

	# zero the bins outside the window
	for x, y in iterutils.MultiIter(*map(range, window.shape)):
		if ((x - window.shape[0] / 2) / float(bins_x) * 2.0)**2 + ((y - window.shape[1] / 2) / float(bins_y) * 2.0)**2 > 1.0:
			window[x, y] = 0.0

	# normalize
	window /= numpy.sum(window)

	return window


#
# =============================================================================
#
#                                  Filtering
#
# =============================================================================
#


def filter_array(a, window, cyclic = False):
	"""
	Filter an array using the window function.  The transformation is
	done in place.  If cyclic = True, then the data are assumed to be
	periodic in each dimension, otherwise the data are assumed to be 0
	outside of their domain of definition.  The window function must
	have an odd number of samples in each dimension;  this is done so
	that it is always clear which sample is at the window's centre,
	which helps prevent phase errors.  If the window function's size
	exceeds that of the data in one or more dimensions, the largest
	allowed central portion of the window function in the affected
	dimensions will be used.  This is done silently;  to determine if
	window function truncation will occur, check for yourself that your
	window function is smaller than your data in all directions.
	"""
	# check that the window and the data have the same number of
	# dimensions
	dims = len(a.shape)
	if dims != len(window.shape):
		raise ValueError, "array and window dimensions mismatch"
	# check that all of the window's dimensions have an odd size
	if 0 in map((1).__and__, window.shape):
		raise ValueError, "window size is not an odd integer in at least 1 dimension"
	# determine how much of the window function can be used
	window_slices = []
	for d in xrange(dims):
		if window.shape[d] > a.shape[d]:
			# largest odd integer <= size of a
			n = ((a.shape[d] + 1) / 2) * 2 - 1
			first = (window.shape[d] - n) / 2
			window_slices.append(slice(first, first + n))
		else:
			window_slices.append(slice(0, window.shape[d]))
	if dims == 1:
		if cyclic:
			# FIXME: use fftconvolve for increase in speed when
			# cyclic boundaries are wanted.  but look into
			# accuracy.
			a[:] = signaltools.convolve(numpy.concatenate((a, a, a)), window[window_slices], mode = "same")[len(a) : 2 * len(a)]
		else:
			a[:] = signaltools.convolve(a, window[window_slices], mode = "same")
	elif dims == 2:
		if cyclic:
			a[:,:] = signaltools.convolve2d(a, window[window_slices], mode = "same", boundary = "wrap")
		else:
			a[:,:] = signaltools.convolve2d(a, window[window_slices], mode = "same")
	else:
		raise ValueError, "can only filter 1 and 2 dimensional arrays"
	return a


def filter_binned_ratios(ratios, window, cyclic = False):
	"""
	Convolve the numerator and denominator of a BinnedRatios instance
	each with the same window function.  This has the effect of
	interpolating the ratio of the two between bins where it has been
	measured, weighting bins by the number of measurements made in
	each.  For example, consider a 1-dimensional binning, with zeros in
	the denominator and numerator bins everywhere except in one bin
	where both are set to 1.0.  The ratio is 1.0 in that bin, and
	undefined everywhere else, where it has not been measured.
	Convolving both numerator and denominator with a Gaussian window
	will replace the "delta function" in each with a smooth hill
	spanning some number of bins.  Since the same smooth hill will be
	seen in both the numerator and the denominator bins, the ratio of
	the two is now 1.0 --- the ratio from the bin where a measurement
	was made --- everywhere the window function had support.  Contrast
	this to the result of convolving the ratio with a window function.

	Convolving the numerator and denominator bins separately preserves
	the integral of each.  In other words the total number of events in
	each of the denominator and numerator is conserved, only their
	locations are shuffled about in order to smooth out irregularities
	in their distributions.  Convolving, instead, the ratios with a
	window function would preserve the integral of the efficiency,
	which is probably meaningless.

	Note that you should be using the window functions defined in this
	module, which are carefully designed to be norm preserving (the
	integrals of the numerator and denominator bins are preserved), and
	phase preserving.

	Note, also, that you should apply this function *before* using
	either of the regularize() methods of the BinnedRatios object.
	"""
	filter_array(ratios.numerator.array, window, cyclic = cyclic)
	filter_array(ratios.denominator.array, window, cyclic = cyclic)


#
# =============================================================================
#
#                                    Rates
#
# =============================================================================
#


def to_moving_mean_density(binned_array, filterdata, cyclic = False):
	"""
	Convolve a BinnedArray with a filter function, then divide all bins
	by their volumes.  The result is the density function smoothed by
	the filter.  The default is to assume 0 values beyond the ends of
	the array when convolving with the filter function.  Set the
	optional cyclic parameter to True for periodic boundaries.

	Example:

	>>> x = BinnedArray(NDBins((LinearBins(0, 10, 5),)))
	>>> x[5.0,] = 1
	>>> x.array
	array([ 0.,  0.,  1.,  0.,  0.])
	>>> to_moving_mean_density(x, tophat_window(3))
	>>> x.array
	array([ 0.        ,  0.16666667,  0.16666667,  0.16666667,  0.
	])

	Explanation.  There are five bins spanning the interval [0, 10],
	making each bin 2 "units" in size.  A single count is placed at
	5.0, which is bin number 2.  The averaging filter is a top-hat
	window 3 bins wide.  The single count in bin #2, when averaged over
	the three bins around it, is equivalent to a mean density of 1/6
	events / unit.

	Example:

	>>> x = BinnedArray(NDBins((LinearBins(0, 10, 5),)))
	>>> x[1,] = 1
	>>> x[3,] = 1
	>>> x[5,] = 1
	>>> x[7,] = 1
	>>> x[9,] = 1
	>>> x.array
	array([ 1.,  1.,  1.,  1.,  1.])
	>>> to_moving_mean_density(x, tophat_window(3))
	>>> x.array
	array([ 0.33333333,  0.5       ,  0.5       ,  0.5       ,  0.33333333])

	We have uniformly distributed events at 2 unit intervals (the first
	is at 1, the second at 3, etc.).  The event density is 0.5 events /
	unit, except at the edges where the smoothing window has picked up
	zero values from beyond the ends of the array.
	"""
	filter_array(binned_array.array, filterdata, cyclic = cyclic)
	binned_array.to_density()


#
# =============================================================================
#
#                                     I/O
#
# =============================================================================
#


from glue.ligolw import ligolw
from glue.ligolw import array
from glue.ligolw import param
from glue.ligolw import table
from glue.ligolw import lsctables


class BinsTable(table.Table):
	"""
	LIGO Light Weight XML table defining a binning.
	"""
	tableName = "pylal_rate_bins:table"
	validcolumns = {
		"order": "int_4u",
		"type": "lstring",
		"min": "real_8",
		"max": "real_8",
		"n": "int_4u"
	}


def bins_to_xml(bins):
	"""
	Construct a LIGO Light Weight XML table representation of the
	NDBins instance bins.
	"""
	xml = lsctables.New(BinsTable)
	for order, bin in enumerate(bins):
		row = xml.RowType()
		row.order = order
		row.type = {
			LinearBins: "lin",
			LogarithmicBins: "log",
			ATanBins: "atan"
		}[bin.__class__]
		row.min = bin.min
		row.max = bin.max
		row.n = len(bin)
		xml.append(row)
	return xml


def bins_from_xml(xml):
	"""
	From the XML document tree rooted at xml, retrieve the table
	describing a binning, and construct and return a rate.NDBins object
	from it.
	"""
	xml = table.get_table(xml, BinsTable.tableName)
	binnings = [None] * (len(xml) and (max(xml.getColumnByName("order")) + 1))
	for row in xml:
		if binnings[row.order] is not None:
			raise ValueError, "duplicate binning for dimension %d" % row.order
		binnings[row.order] = {
			"lin": LinearBins,
			"log": LogarithmicBins,
			"atan": ATanBins
		}[row.type](row.min, row.max, row.n)
	if None in binnings:
		raise ValueError, "no binning for dimension %d" % binnings.find(None)
	return NDBins(binnings)


def binned_array_to_xml(binnedarray, name):
	"""
	Retrun an XML document tree describing a rate.BinnedArray object.
	"""
	xml = ligolw.LIGO_LW({u"Name": u"%s:pylal_rate_binnedarray" % name})
	xml.appendChild(bins_to_xml(binnedarray.bins))
	xml.appendChild(array.from_array(u"array", binnedarray.array))
	return xml


def binned_array_from_xml(xml, name):
	"""
	Search for the description of a rate.BinnedArray object named
	"name" in the XML document tree rooted at xml, and construct and
	return a new rate.BinnedArray object from the data contained
	therein.
	"""
	xml = [elem for elem in xml.getElementsByTagName(ligolw.LIGO_LW.tagName) if elem.getAttribute(u"Name") == u"%s:pylal_rate_binnedarray" % name]
	try:
		xml, = xml
	except ValueError:
		raise ValueError, "document must contain exactly 1 BinnedArray named '%s'" % name
	# an empty binning is used for the initial object creation instead
	# of using the real binning to avoid the creation of a (possibly
	# large) array that would otherwise accompany this step
	binnedarray = BinnedArray(NDBins())
	binnedarray.bins = bins_from_xml(xml)
	binnedarray.array = array.get_array(xml, u"array").array
	return binnedarray


def binned_ratios_to_xml(ratios, name):
	"""
	Return an XML document tree describing a rate.BinnedRatios object.
	"""
	xml = ligolw.LIGO_LW({u"Name": u"%s:pylal_rate_binnedratios" % name})
	xml.appendChild(binned_array_to_xml(ratios.numerator, u"numerator"))
	xml.appendChild(binned_array_to_xml(ratios.denominator, u"denominator"))
	return xml


def binned_ratios_from_xml(xml, name):
	"""
	Search for the description of a rate.BinnedRatios object named
	"name" in the XML document tree rooted at xml, and construct and
	return a new rate.BinnedRatios object from the data contained
	therein.
	"""
	xml, = [elem for elem in xml.getElementsByTagName(ligolw.LIGO_LW.tagName) if elem.getAttribute(u"Name") == u"%s:pylal_rate_binnedratios" % name]
	ratios = BinnedRatios(NDBins())
	ratios.numerator = binned_array_from_xml(xml, u"numerator")
	ratios.denominator = binned_array_from_xml(xml, u"denominator")
	# normally they share a single NDBins instance
	ratios.denominator.bins = ratios.numerator.bins
	return ratios
