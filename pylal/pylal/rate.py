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


import math
import numpy
from scipy.signal import signaltools


from glue import segments
from pylal import itertools
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

	def __cmp__(self, other):
		return cmp((type(self), self.min, self.max, self.n), (type(other), other.min, other.max, other.n))

	def _set_delta(self, min, max, n):
		raise NotImplementedError

	def __getitem__(self, x):
		raise NotImplementedError


class _LinBins(_Bins):
	"""
	Linearly-spaced 1-D bins.  For internal use only.
	"""
	def _set_delta(self, min, max, n):
		self.delta = float(max - min) / n

	def __getitem__(self, x):
		if isinstance(x, segments.segment):
			return slice(self[x[0]], self[x[1]] + 1)
		if isinstance(x, slice):
			if x.step is not None:
				raise NotImplementedError, "slices with steps not yet supported"
			return slice(self[x.start], self[x.stop])
		if self.min <= x < self.max:
			return int((x - self.min) / self.delta)
		if x == self.max:
			# special "measure zero" corner case
			return self.n - 1
		raise IndexError, x

	def centres(self):
		return self.min + self.delta * (numpy.arange(0, self.n) + 0.5)


class _LogBins(_Bins):
	"""
	Logarithmically-spaced 1-D bins.  For internal use only.
	"""
	def _set_delta(self, min, max, n):
		self.delta = math.log(float(max / min)) / n

	def __getitem__(self, x):
		if isinstance(x, segments.segment):
			return slice(self[x[0]], self[x[1]] + 1)
		if isinstance(x, slice):
			if x.step is not None:
				raise NotImplementedError, x
			return slice(self[x.start], self[x.stop])
		if self.min <= x < self.max:
			return int(math.log(x / self.min) / self.delta)
		if x == self.max:
			# special "measure zero" corner case
			return self.n - 1
		raise IndexError, x

	def centres(self):
		return self.min * numpy.exp(self.delta * (numpy.arange(0, self.n) + 0.5))


class Bins(object):
	"""
	Multi-dimensional co-ordinate binning.  An instance of this object
	is used to convert a tuple of co-ordinates into a tuple of array
	indices, thereby allowing the contents of an array object to be
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
	>>> b[1, 1:5]
	(0, slice(0, 1, None))
	>>> b[1, segment(1, 5)]
	(0, slice(0, 2, None))
	>>> b.centres
	(array([  5.,  13.,  21.]), array([  1.70997595,   5.,  14.62008869]))
	"""
	def __init__(self, *args, **kwargs):
		if len(args) % 3:
			raise TypeError, "arguments must be min,max,n[,min,max,n]..."
		spacing = kwargs.get("spacing", ["lin"] * (len(args) / 3))
		if len(spacing) != len(args) / 3:
			raise ValueError, spacing
		it = iter(args)
		# FIXME: don't construct the intermediate lists when we
		# no longer care about Python 2.3 compatibility.
		self.bins = tuple([{
			"lin": _LinBins,
			"log": _LogBins
		}[s](it.next(), it.next(), it.next()) for s in spacing])
		self.min = tuple([b.min for b in self.bins])
		self.max = tuple([b.max for b in self.bins])
		self.shape = tuple([b.n for b in self.bins])
		self.centres = tuple([b.centres() for b in self.bins])

	def __call__(self, *args):
		return self[args]

	def __getitem__(self, coords):
		"""
		Return the indices corresponding to the tuple of
		co-ordinates, coords.  Note that a the co-ordinates must be
		a tuple even if there is only 1 dimension.  Each
		co-ordinate can be a single number, a Python slice object,
		or a glue.segments.segment object.  The difference between
		co-ordinate ranges given as slices and ranges given as
		segments is that a slice is exclusive of the upper bound
		while a segment is inclusive of the upper bound.
		"""
		if len(coords) != len(self.bins):
			raise ValueError, "dimension mismatch"
		return tuple(map(lambda b, c: b[c], self.bins, coords))

	def __cmp__(self, other):
		"""
		Return 0 if the Bins objects are "compatible", meaning they
		have the same upper and lower bounds, the same number of
		bins, and the same bin spacing (linear, logarithmic, etc.).
		Return non-zero otherwise.
		"""
		return reduce(int.__or__, map(cmp, self.bins, other.bins))


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

	Note that even for 1 dimensional arrays the index must be a tuple.
	"""
	def __init__(self, bins, dtype = numpy.float64):
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
		for coords in itertools.MultiIter(*other.bins.centres):
			self[coords] += other[coords]
		return self

	def centres(self):
		"""
		Return a tuple of arrays containing the bin centres for
		each dimension.
		"""
		return self.bins.centres

	def logregularize(self, epsilon = 2**-1074):
		"""
		Find bins <= 0, and set them to epsilon, This has the
		effect of allowing the logarithm of the array to be
		evaluated without error.
		"""
		# FIXME:  assign to array's contents instead
		self.array = numpy.where(self.array > 0, self.array, epsilon)
		return self

	def used(self):
		"""
		Return the number of bins that are non-zero.
		"""
		return len(numpy.nonzero(self.array)[0])


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
	def __init__(self, bins, dtype = numpy.float64):
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
		# FIXME: assign to denominator's contents instead
		self.denominator.array = numpy.where(self.denominator.array, self.denominator.array, 1)
		return self

	def logregularize(self, epsilon = 2**-1074):
		"""
		Find bins in the denominator that are 0, and set them to 1,
		while setting the corresponding bin in the numerator to
		float epsilon.  This has the effect of allowing the
		logarithm of the ratio array to be evaluated without error.
		"""
		# FIXME: assign to contents instead
		self.numerator.array = numpy.where(self.denominator.array, self.numerator.array, epsilon)
		self.denominator.array = numpy.where(self.denominator.array, self.denominator.array, 1)
		return self

	def centres(self):
		"""
		Return a tuple of arrays containing the bin centres for
		each dimension.
		"""
		return self.numerator.bins.centres

	def used(self):
		"""
		Return the number of bins with non-zero denominator.
		"""
		return numpy.sum(numpy.where(self.denominator.array, 1, 0))


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
	if bins <= 0:
		raise ValueError, bins
	l = 20 * int(bins / 2.0)
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
	bins sets the width of the window in bin counts.
	"""
	if bins <= 0:
		raise ValueError, bins
	w = window.XLALCreateRectangularREAL8Window(int(bins / 2.0) * 2 + 1)
	return w.data / w.sum


def tophat_window2d(bins_x, bins_y):
	"""
	Generate a normalized (integral = 1) top-hat window in 2
	dimensions.  bins_x and bins_y set the widths of the window in bin
	counts.  The result is a rectangular array, with an elliptical
	pattern of elements set to a constant value centred on the array's
	mid-point, and all other elements set to 0.
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
	window = numpy.ones((int(bins_x / 2) * 2 + 1, int(bins_y / 2) * 2 + 1), "Float64")

	# zero the bins outside the window
	for x in xrange(window.shape[0]):
		for y in xrange(window.shape[1]):
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
	done in place.  If cyclic = True, then the data is assumed to be
	periodic in each dimension, otherwise the data is assumed to be 0
	outside of its domain of definition.  The window function must have
	an odd number of samples in each dimension;  this is done so that
	it is always clear which sample is at the window's centre, which
	helps prevent phase errors.  If the window function's size exceeds
	that of the data in one or more dimensions, the largest allowed
	central portion of the window function in the affected dimensions
	will be used.  This is done silently;  to determine if window
	function truncation will occur, check for yourself that your window
	function is smaller than your data in all directions.
	"""
	# check that the window and the data have the same number of
	# dimensions
	dims = len(a.shape)
	if dims != len(window.shape):
		raise ValueError, "array and window dimensions mismatch"
	# check that all of the window's dimensions have an odd size
	if 0 in map(int(1).__and__, window.shape):
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


class Rate(BinnedArray):
	"""
	An object for binning and smoothing impulsive data in 1 dimension,
	normalized so as to measure events (or event weight) per filter
	width.
	"""
	def __init__(self, segment, filterwidth, windowfunc = gaussian_window):
		"""
		Initialize the bins for the given segment and filter width.
		"""
		#
		# bin size is 1/20th of the filter's width, but adjusted so
		# that there is an integer number of bins in the interval
		#

		BinnedArray.__init__(self, Bins(segment[0], segment[1], int(abs(segment) / (filterwidth / 20.0)) + 1))

		#
		# determine the true bin size from the final integer bin
		# count and the interval's length
		#

		self.binsize = float(abs(segment)) / self.bins.shape[0]

		#
		# generate the filter data
		#

		self.set_filter(filterwidth, windowfunc)

	def __getitem__(self, x):
		return self.array[self.bins[x,]]

	def __setitem__(self, x, value):
		self.array[self.bins[x,]] = value

	def xvals(self):
		return self.centres()[0]

	def set_filter(self, filterwidth, windowfunc):
		"""
		Update the filter function.  Note that this is less ideal
		than specifying the correct filter function at object
		creation time because the binning may not be appropriate
		for this new filter (too many bins or too few).  Returns
		self.
		"""
		#
		# safety-check filter width
		#

		if filterwidth / self.binsize < 3:
			raise ValueError, "filter too narrow (less than 3 bins)"
		self.filterwidth = filterwidth

		#
		# Construct and store the window function
		#

		self.filterdata = windowfunc(self.filterwidth / self.binsize) / self.binsize

		#
		# Done
		#

		return self

	def filter(self, cyclic = False):
		"""
		Convolve the binned weights with the window to generate the
		rate data.
		"""
		return filter_array(self.array, self.filterdata, cyclic = cyclic)


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
	rate.Bins instance bins.
	"""
	xml = lsctables.New(BinsTable)
	for order, bin in enumerate(bins.bins):
		row = xml.RowType()
		row.order = order
		row.type = {
			_LinBins: "lin",
			_LogBins: "log"
		}[bin.__class__]
		row.min = bin.min
		row.max = bin.max
		row.n = bin.n
		xml.append(row)
	return xml


def bins_from_xml(xml):
	"""
	From the XML document tree rooted at xml, retrieve the table
	describing a binning, and construct and return a rate.Bins object
	from it.
	"""
	xml = table.get_table(xml, BinsTable.tableName)
	xml.sort(lambda a, b: cmp(a.order, b.order))
	args = [None] * 3 * (len(xml) and (xml[-1].order + 1))
	kwargs = {"spacing" : [None] * (len(xml) and (xml[-1].order + 1))}
	for row in xml:
		args[row.order * 3 : row.order * 3 + 3] = row.min, row.max, row.n
		kwargs["spacing"][row.order] = row.type
	if None in args:
		raise ValueError, "incomplete bin spec: %s" % str(args)
	return Bins(*args, **kwargs)


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
	xml, = [elem for elem in xml.getElementsByTagName(ligolw.LIGO_LW.tagName) if elem.getAttribute(u"Name") == u"%s:pylal_rate_binnedarray" % name]
	binnedarray = BinnedArray(Bins())
	binnedarray.bins = bins_from_xml(xml)
	binnedarray.array = array.get_array(xml, u"array").array
	return binnedarray


def rate_to_xml(rate, name):
	"""
	Retrun an XML document tree describing a rate.BinnedArray object.
	"""
	xml = binned_array_to_xml(rate, name)
	xml.appendChild(param.from_pyvalue(u"binsize", rate.binsize))
	xml.appendChild(param.from_pyvalue(u"filterwidth", rate.filterwidth))
	xml.appendChild(array.from_array(u"filterdata", rate.filterdata))
	return xml


def rate_from_xml(xml, name):
	"""
	Search for the description of a rate.Rate object named "name" in
	the XML document tree rooted at xml, and construct and return a new
	rate.Rate instance initialized from the data contained therein.
	"""
	# FIXME:  figure out how to not duplicate code from
	# binned_array_from_xml()
	xml, = [elem for elem in xml.getElementsByTagName(ligolw.LIGO_LW.tagName) if elem.getAttribute(u"Name") == u"%s:pylal_rate_binnedarray" % name]
	rate = Rate(segments.segment(0, 1), 1)
	rate.bins = bins_from_xml(xml)
	rate.array = array.get_array(xml, u"array").array
	rate.binsize = param.get_pyvalue(xml, u"binsize")
	rate.filterwidth = param.get_pyvalue(xml, u"filterwidth")
	rate.filterdata = array.get_array(xml, u"filterdata").array
	return rate


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
	ratios = BinnedRatios(Bins())
	ratios.numerator = binned_array_from_xml(xml, u"numerator")
	ratios.denominator = binned_array_from_xml(xml, u"denominator")
	# normally they share a single Bins instance
	ratios.denominator.bins = ratios.numerator.bins
	return ratios
