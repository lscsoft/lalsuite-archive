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


#
# NOTE:  the logic in this code is unintuitively complicated.  Small,
# apparently irrelevant, changes to conditionals can have subtly unexpected
# consequences to the behaviour of the class methods.  ALWAYS make sure that
# the test suite returns OK on ALL tests after any changes you make.
#


"""
This module defines the segment and segmentlist objects, as well as the
infinity object used to define semi-infinite and infinite segments.
"""

from bisect import bisect_left, bisect_right
from copy import copy as shallowcopy


__author__ = "Kipp Cannon <kipp@gravity.phys.uwm.edu>"
__date__ = "$Date$"[7:-2]
__version__ = "$Revision$"[11:-2]


#
# =============================================================================
#
#                                   infinity
#
# =============================================================================
#


class infinity(object):
	"""
	The infinity object possess the algebraic properties necessary for
	use as a bound on semi-infinite and infinite segments.

	This class uses comparison-by-identity rather than
	comparison-by-value.  What this means, is there are only ever two
	instances of this class, representing positive and negative
	infinity respectively.  All other "instances" of this class are
	infact simply references to one of these two, and comparisons are
	done by checking which one you've got.  This improves speed and
	reduces memory use, and is similar in implementation to Python's
	boolean True and False objects.

	The normal way to obtain references to positive or negative
	infinity is to do infinity() or -infinity() respectively.  It is
	also possible to select the sign by passing a single numeric
	argument to the constructor.  The sign of the argument causes a
	reference to either positive or negative infinity to be returned,
	respectively.  For example infinity(-1) is equivalent to
	-infinity().  However, this feature is a little slower and not
	recommended for normal use;  it is provided only to simplify the
	pickling and unpickling of instances of the class.

	Example:

	>>> x = infinity()
	>>> x > 0
	True
	>>> x + 10 == x
	True
	>>> segment(-10, 10) - segment(-x, 0)
	segment(0, 10)
	"""
	__slots__ = []

	def __new__(cls, *args):
		if args:
			# pickle support
			(sign,) = args
			if sign > 0:
				return PosInfinity
			if sign < 0:
				return NegInfinity
			raise ValueError, sign
		return PosInfinity

	def __repr__(self):
		"""
		Returns a string.
		"""
		if self is PosInfinity:
			return "infinity"
		return "-infinity"

	__str__ = __repr__

	# pickle support

	def __reduce__(self):
		"""
		Pickle support.
		"""
		if self is PosInfinity:
			return (infinity, (1,))
		# self is NegInfinity
		return (infinity, (-1,))

	# tests

	def __cmp__(self, other):
		"""
		Positive infinity compares as greater than everything
		except itself, negative infinity compares as less than
		everything except itself.
		"""
		if self is other:
			return 0
		if self is PosInfinity:
			return 1
		# self is NegInfinity
		return -1

	def __nonzero__(self):
		"""
		Returns True.
		"""
		return True

	# arithmetic

	def __neg__(self):
		"""
		Returns -self.
		"""
		if self is PosInfinity:
			return NegInfinity
		# self is NegInfinity
		return PosInfinity

	def __pos__(self):
		"""
		Returns self.
		"""
		return self

	def __add__(self, other):
		"""
		Returns self.
		"""
		return self

	def __radd__(self, other):
		"""
		Returns self.
		"""
		return self

	def __sub__(self, other):
		"""
		Returns self.
		"""
		return self

	def __rsub__(self, other):
		"""
		Returns -self.
		"""
		if self is PosInfinity:
			return NegInfinity
		# self is NegInfinity
		return PosInfinity


PosInfinity = object.__new__(infinity)
NegInfinity = object.__new__(infinity)


#
# =============================================================================
#
#                                   segment
#
# =============================================================================
#


class segment(tuple):
	"""
	The segment class defines objects that represent a range of values.
	A segment has a start and an end, and is taken to represent the
	range of values in the semi-open interval [start, end).  Some
	limited arithmetic operations are possible with segments, but
	because the set of (single) segments is not closed under the
	sensible definitions of the standard arithmetic operations, the
	behaviour of the arithmetic operators on segments may not be as you
	would expect.  For general arithmetic on segments, use segmentlist
	objects.  The methods for this class exist mostly for purpose of
	simplifying the implementation of the segmentlist class.

	Example:

	>>> segment(0, 10) & segment(5, 15)
	segment(5, 10)
	>>> segment(0, 10) | segment(5, 15)
	segment(0, 15)
	>>> segment(0, 10) - segment(5, 15)
	segment(0, 5)
	>>> segment(0, 10) < segment(5, 15)
	True
	>>> segment(1, 2) in segment(0, 10)
	True
	>>> bool(segment(0, 0))
	False
	>>> segment("AAA Towing", "York University") & segment("Pool", "Zoo")
	segment('Pool', 'York University')
	>>> x = [0, 1]
	>>> segment(x)
	segment(0, 1)
	"""

	# basic class methods

	def __new__(cls, *args):
		if len(args) == 1:
			args = args[0]
		if len(args) != 2:
			raise TypeError, "__new__() takes 2 arguments, or 1 argument when it is a sequence of length 2"
		if args[0] <= args[1]:
			return tuple.__new__(cls, args)
		else:
			return tuple.__new__(cls, (args[1], args[0]))

	def __repr__(self):
		return "segment(" + repr(self[0]) + ", " + repr(self[1]) + ")"

	def __str__(self):
		return "[" + str(self[0]) + " ... " + str(self[1]) + ")"

	# accessors

	def __abs__(self):
		"""
		Returns the length of the interval represented by the
		segment.
		"""
		return self[1] - self[0]

	duration = __abs__

	# comparisons

	def __nonzero__(self):
		"""
		Test for segment having non-zero duration.
		"""
		return self[0] != self[1]

	def not_continuous(self, other):
		"""
		Returns >0 if self covers an interval above other's
		interval, <0 if self covers an interval below other's, or 0
		if the two intervals are not disjoint (intersect or touch).
		A return value of 0 indicates the two segments would
		coalesce.
		"""
		if self[0] > other[1]:
			return 1
		if self[1] < other[0]:
			return -1
		return 0

	def order(self, other):
		"""
		Equivalent to cmp(self, other) except that a result of 0
		indicates that other is contained in self rather than being
		identically equal to self.
		"""
		# FIXME: is this used?  It is not used anywhere in glue or
		# pylal.  I want to deprecate it.
		if other in self:
			return 0
		return cmp(self, other)

	# some arithmetic operations that (mostly) make sense for segments

	def __and__(self, other):
		"""
		Return the segment that is the intersection of the given
		segments.  Raises ValueError if the result cannot be
		presented as a single segment.
		"""
		if (self[1] <= other[0]) or (self[0] >= other[1]):
			# self and other don't intersect
			raise ValueError, other
		return tuple.__new__(segment, (max(self[0], other[0]), min(self[1], other[1])))

	def __or__(self, other):
		"""
		Return the segment that is the union of the given segments.
		Raises ValueError if the result cannot be represented as a
		single segment.
		"""
		if (self[1] < other[0]) or (self[0] > other[1]):
			# self and other are disjoint
			raise ValueError, other
		return tuple.__new__(segment, (min(self[0], other[0]), max(self[1], other[1])))

	# addition is union
	__add__ = __or__

	def __sub__(self, other):
		"""
		Return the segment that is that part of self which is not
		contained in other.  Raises ValueError if the result cannot
		be represented as a single segment.
		"""
		if (self[1] <= other[0]) or (self[0] >= other[1]):
			# self and other do not intersect
			return self
		if (self in other) or ((self[0] < other[0]) and (self[1] > other[1])):
			# result is not exactly 1 segment
			raise ValueError, other
		if self[0] < other[0]:
			return tuple.__new__(segment, (self[0], other[0]))
		return tuple.__new__(segment, (other[1], self[1]))

	# check for proper intersection, containment, and continuity

	def intersects(self, other):
		"""
		Return True if the intersection of self and other is not a
		null segment.
		"""
		return (self[1] > other[0]) and (self[0] < other[1])

	def __contains__(self, other):
		"""
		Return True if other is wholly contained in self.  other
		can be another segment or an object of the same type as the
		bounds of self.
		"""
		if type(other) == segment:
			return (self[0] <= other[0]) and (self[1] >= other[1])
		else:
			return self[0] <= other < self[1]

	def continuous(self, other):
		"""
		Return True if self and other are not disjoint.
		"""
		# FIXME: deprecate this in favour of not_continuous()
		# above.
		return not self.not_continuous(other)

	# protraction and contraction and shifting

	def protract(self, x):
		"""
		Move both the start and the end of the segment a distance x
		away from the other.
		"""
		return segment(self[0] - x, self[1] + x)

	def contract(self, x):
		"""
		Move both the start and the end of the segment a distance x
		towards the the other.
		"""
		return segment(self[0] + x, self[1] - x)

	def shift(self, x):
		"""
		Return a new segment by adding x to the upper and lower
		bounds of this segment.
		"""
		return tuple.__new__(segment, (self[0] + x, self[1] + x))


#
# =============================================================================
#
#                                 segmentlist
#
# =============================================================================
#


class segmentlist(list):
	"""
	The segmentlist class defines a list of segments, and is an
	extension of the built-in list class.  This class provides
	addtional methods that assist in the manipulation of lists of
	segments.  In particular, arithmetic operations such as union and
	intersection are provided.  Unlike the segment class, the
	segmentlist class is closed under all supported arithmetic
	operations.

	All standard Python sequence-like operations are supported, like
	slicing, iteration and so on, but the arithmetic and other methods
	in this class generally expect the segmentlist to be in what is
	refered to as a "coalesced" state --- consisting solely of disjoint
	segments listed in ascending order.  Using the standard Python
	sequence-like operations, a segmentlist can be easily constructed
	that is not in this state;  for example by simply appending a
	segment to the end of the list that overlaps some other segment
	already in the list.  The class provides a coalesce() method that
	can be called to put it in the coalesced state.  Following
	application of the coalesce method, all arithmetic operations will
	function reliably.  All arithmetic methods themselves return
	coalesced results, so there is never a need to call the coalesce
	method when manipulating segmentlists exclusively via the
	arithmetic operators.

	Example:

	>>> x = segmentlist([segment(-10, 10)])
	>>> x |= segmentlist([segment(20, 30)])
	>>> x -= segmentlist([segment(-5, 5)])
	>>> print x
	[segment(-10, -5), segment(5, 10), segment(20, 30)]
	>>> print ~x
	[segment(-infinity, -10), segment(-5, 5), segment(10, 20), segment(30, infinity)]
	"""

	# container method over-rides.

	def __contains__(self, item):
		"""
		Returns True if the given object is wholly contained within
		one of the segments in self.  Does not require the
		segmentlist to be coalesced.

		Note the difference between this operator, and the standard
		Python "in" operator for sequence-like objects:  in the
		case of standard sequence-like objects the in operator
		checks for an exact match between the given item and one of
		the contents of the list; for segmentlists, the in operator
		checks if the given item is contained within any of the
		segments in the segmentlist.
		"""
		for seg in self:
			if item in seg:
				return True
		return False

	# supplementary accessors

	def __abs__(self):
		"""
		Return the sum of the durations of all segments in self.
		Does not require the segmentlist to be coalesced.
		"""
		d = 0
		for seg in self:
			d += abs(seg)
		return d

	duration = __abs__

	def extent(self):
		"""
		Return the segment whose end-points denote the maximum and
		minimum extent of the segmentlist.  Does not require the
		segmentlist to be coalesced.
		"""
		if not len(self):
			raise ValueError, "empty list"
		min, max = self[0]
		for seg in self:
			if min > seg[0]:
				min = seg[0]
			if max < seg[1]:
				max = seg[1]
		return segment(min, max)

	def find(self, item):
		"""
		Return the smallest i such that i is the index of an
		element that wholly contains the given segment.  Raises
		ValueError if no such element exists.  Does not require the
		segmentlist to be coalesced.
		"""
		for i, seg in enumerate(self):
			if item in seg:
				return i
		raise ValueError, item

	# arithmetic operations that are sensible with segment lists

	def __iand__(self, other):
		"""
		Replace the segmentlist with the intersection of itself and
		another.  This operation is O(n).
		"""
		return self.__isub__(~other)

	def __and__(self, other):
		"""
		Return the intersection of the segmentlist and another.
		This operation is O(n).
		"""
		return segmentlist(self[:]).__iand__(other)

	def __ior__(self, other):
		"""
		Replace the segmentlist with the union of itself and
		another.  If the two lists have numbers of elements m and n
		respectively, then this algorithm is O(n log m), which
		means it is optimized for the case when self is large and
		other is small.  In practice, this has been found to be the
		more common scenario.  If you have two large lists of
		comparable size, n, then it is faster to do

		>>> list1.extend(list2)
		>>> list1.coalesce()

		This is still O(n log n), but with a smaller leading
		coefficient.
		"""
		i = 0
		for seg in other:
			i = j = bisect_right(self, seg, i)
			if i and not self[i - 1].not_continuous(seg):
				i -= 1
				seg |= self[i]
			n = len(self)
			while j < n and not seg.not_continuous(self[j]):
				j += 1
			if j > i:
				self[i : j] = [seg | self[j - 1]]
			else:
				self.insert(i, seg)
			i += 1
		return self

	def __or__(self, other):
		"""
		Return the union of the segment list and another.  The
		comment in the segmentlist.__ior__() method about
		performance optimization when computing the union of large
		lists of similar size applies here as well.
		"""
		if len(self) >= len(other):
			return segmentlist(self[:]).__ior__(other)
		return segmentlist(other[:]).__ior__(self)

	def __xor__(self, other):
		"""
		Return the segmentlist that is the list of all intervals
		contained in exactly one of this and another list.  This
		operation is O(n log n).
		"""
		return (self - other) | (other - self)

	# addition is union
	__iadd__ = __ior__
	__add__ = __or__

	def __isub__(self, other):
		"""
		Replace the segmentlist with the difference between itself
		and another.  This operation is O(n).
		"""
		if not other:
			return self
		i = j = 0
		otherseg = other[j]
		while i < len(self):
			seg = self[i]
			while (not otherseg) or otherseg[1] <= seg[0]:
				j += 1
				if j >= len(other):
					return self
				otherseg = other[j]
			if seg[1] <= otherseg[0]:
				i += 1
			elif otherseg[0] <= seg[0]:
				if otherseg[1] >= seg[1]:
					del self[i]
				else:
					self[i] = segment(otherseg[1], seg[1])
			else:
				self[i] = segment(seg[0], otherseg[0])
				i += 1
				if otherseg[1] < seg[1]:
					self.insert(i, segment(otherseg[1], seg[1]))
		return self

	def __sub__(self, other):
		"""
		Return the difference between the segmentlist and another.
		This operation is O(n).
		"""
		return segmentlist(self[:]).__isub__(other)

	def __invert__(self):
		"""
		Return the segmentlist that is the inversion of the given
		list.  This operation is O(n).
		"""
		if len(self) == 0:
			return segmentlist([segment(NegInfinity, PosInfinity)])
		l = segmentlist()
		if self[0][0] > NegInfinity:
			l.append(segment(NegInfinity, self[0][0]))
		last = self[0][1]
		for j in xrange(1, len(self)):
			l.append(segment(last, self[j][0]))
			last = self[j][1]
		if last < PosInfinity:
			l.append(segment(last, PosInfinity))
		return l

	# other operations

	def split(self, value):
		"""
		Break all segments that stradle the given value at that
		value.  Does not require the segmentlist to be coalesced,
		and the result is not coalesced by definition.  This
		operation is O(n).
		"""
		for i, seg in enumerate(self):
			if value in seg:
				self[i:i+1] = [segment(seg[0], value), segment(value, seg[1])]

	def intersects_segment(self, other):
		"""
		Returns True if the intersection of self and the segment
		other is not the null set, otherwise returns False.  The
		algorithm is O(log n).  Requires the list to be coalesced.
		"""
		i = bisect_left(self, other)
		return ((i != 0) and (other[0] < self[i-1][1])) or ((i != len(self)) and (other[1] > self[i][0]))

	def intersects(self, other):
		"""
		Returns True if the intersection of self and the
		segmentlist other is not the null set, otherwise returns
		False.  The algorithm is O(n), but faster than explicit
		calculation of the intersection, i.e. by testing len(self &
		other).  Requires both lists to be coalesced.
		"""
		if not (self and other):
			return False
		i = j = 0
		seg = self[0]
		otherseg = other[0]
		while True:
			if seg[1] <= otherseg[0]:
				i += 1
				if i >= len(self):
					return False
				seg = self[i]
			elif otherseg[1] <= seg[0]:
				j += 1
				if j >= len(other):
					return False
				otherseg = other[j]
			else:
				return True

	def coalesce(self):
		"""
		Sort the elements of a list into ascending order, and merge
		continuous segments into single segments.  This operation
		is O(n log n).
		"""
		self.sort()
		i = j = 0
		n = len(self)
		while j < n:
			seg = self[j]
			j += 1
			while j < n and not seg.not_continuous(self[j]):
				seg |= self[j]
				j += 1
			self[i] = seg
			i += 1
		del self[i : ]
		return self

	def protract(self, x):
		"""
		For each segment in the list, move both the start and the
		end a distance x away from the other.  Coalesce the result.
		"""
		for i in xrange(len(self)):
			self[i] = self[i].protract(x)
		return self.coalesce()

	def contract(self, x):
		"""
		For each segment in the list, move both the start and the
		end a distance x towards the other.  Coalesce the result.
		"""
		for i in xrange(len(self)):
			self[i] = self[i].contract(x)
		return self.coalesce()

	def shift(self, x):
		"""
		Shift the segmentlist by adding x to the upper and lower
		bounds of all segments.  The algorithm is O(n) and does not
		require the list to be coalesced.
		"""
		for i in xrange(len(self)):
			self[i] = self[i].shift(x)
		return self


#
# =============================================================================
#
#                               segmentlistdict
#
# =============================================================================
#


class _offsets(dict):
	"""
	Implements the segmentlist offset book-keeping in the
	segmentlistdict class.  Not intended for use outside of the
	segmentlistdict class.
	"""
	def __new__(cls, parent):
		return dict.__new__(cls)

	def __init__(self, parent):
		self.__parent = parent
		for key in self.__parent:
			dict.__setitem__(self, key, 0.0)

	def __setitem__(self, key, value):
		"""
		Set an offset.  If the new offset is identical to the
		current offset this is a no-op, otherwise the corresponding
		segmentlist object is shifted.
		"""
		if key not in self:
			raise KeyError, key
		delta = value - self[key]
		if delta:
			self.__parent[key].shift(delta)
			dict.__setitem__(self, key, self[key] + delta)

	def update(self, d):
		"""
		From a dictionary of offsets, apply each offset to the
		corresponding segmentlist.  NOTE:  it is acceptable for the
		offset dictionary to contain entries for which there is no
		matching segmentlist; no error will be raised, but the
		offset will be ignored.  This simplifies the case of
		updating several segmentlistdict objects from a common
		offset dictionary, when one or more of the segmentlistdicts
		contains only a subset of the keys.
		"""
		for key, value in d.iteritems():
			try:
				self[key] = value
			except KeyError:
				pass

	def clear(self):
		"""
		Remove the offsets from all segmentlists.
		"""
		for key in self:
			self[key] = 0.0

	# stubs to prevent bugs
	def __delitem__(*args):
		raise NotImplementedError
	def fromkeys(*args):
		raise NotImplementedError
	def pop(*args):
		raise NotImplementedError
	def popitem(*args):
		raise NotImplementedError


class segmentlistdict(dict):
	"""
	A dictionary associating a unique label and numeric offset with
	each of a set of segmentlist objects.

	This class implements a standard mapping interface, with additional
	features added to assist with the manipulation of a collection of
	segmentlist objects.  In particular, methods for taking unions and
	intersections of the lists in the dictionary are available, as well
	as the ability to record and apply numeric offsets to the
	boundaries of the segments in each list.

	The numeric offsets are stored in the "offsets" attribute, which
	itself is a dictionary, associating a number with each key in the
	main dictionary.  Assigning to one of the entries of the offsets
	attribute has the effect of shifting the corresponding segmentlist
	from its original position (not its current position) by the given
	amount.

	Example:

	>>> x = segmentlistdict()
	>>> x["H1"] = segmentlist([segment(0, 10)])
	>>> print x
	{'H1': [segment(0, 10)]}
	>>> x.offsets["H1"] = 6
	>>> print x
	{'H1': [segment(6.0, 16.0)]}
	>>> x.offsets.clear()
	>>> print x
	{'H1': [segment(0.0, 10.0)]}
	>>> x["H2"] = segmentlist([segment(5, 15)])
	>>> x.intersection(["H1", "H2"])
	[segment(5, 10.0)]
	>>> x.offsets["H1"] = 6
	>>> x.intersection(["H1", "H2"])
	[segment(11.0, 15)]
	>>> c = x.extract_common(["H1", "H2"])
	>>> c.offsets.clear()
	>>> c
	{'H2': [segment(11.0, 15)], 'H1': [segment(5.0, 9.0)]}
	"""
	def __init__(self, *args):
		dict.__init__(self, *args)
		self.offsets = _offsets(self)
		if args and isinstance(args[0], segmentlistdict):
			dict.update(self.offsets, args[0].offsets)

	def copy(self):
		"""
		Return a copy of the segmentlistdict object.  The return
		value is a new object with a new offsets attribute, with
		references to the original keys, and shallow copies of the
		segment lists.  In summary, modifications made to the
		offset dictionary or segmentlists in the object returned by
		this method will not affect the original, but without using
		much memory until such modifications are made.
		"""
		new = self.__class__()
		for key, value in self.iteritems():
			new[key] = shallowcopy(value)
		dict.update(new.offsets, self.offsets)
		return new

	def __setitem__(self, key, value):
		"""
		Set the segmentlist associated with a key.  If key is not
		already in the dictionary, the corresponding offset is
		initialized to 0.0, otherwise it is left unchanged.
		"""
		dict.__setitem__(self, key, value)
		if key not in self.offsets:
			dict.__setitem__(self.offsets, key, 0.0)

	def __delitem__(self, key):
		dict.__delitem__(self, key)
		dict.__delitem__(self.offsets, key)

	# supplementary accessors

	def map(self, func):
		"""
		Return a dictionary of the results of func applied to each
		of the segmentlist objects in self.
		"""
		d = {}
		for key, value in self.iteritems():
			d[key] = func(value)
		return d

	def __abs__(self):
		"""
		Return a dictionary of the results of running abs() on
		each of the segmentlists.
		"""
		return self.map(abs)

	duration = __abs__

	def extent(self):
		"""
		Return a dictionary of the results of running extent() on
		each of the segmentlists.
		"""
		return self.map(segmentlist.extent)

	def extent_all(self):
		"""
		Return the result of running extent on the union of all
		lists in the dictionary.
		"""
		return self.union(self.iterkeys()).extent()

	def find(self, seg):
		"""
		Return a dictionary of the results of running find() on
		each of the segmentlists.
		"""
		return self.map(lambda l: segmentlist.find(l, seg))

	# list-by-list arithmetic

	def __iand__(self, other):
		for key, value in other.iteritems():
			if key in self:
				self[key] &= value
		return self

	def __and__(self, other):
		return self.copy().__iand__(other)

	def __ior__(self, other):
		for key, value in other.iteritems():
			if key in self:
				self[key] |= value
			else:
				self[key] = shallowcopy(value)
		return self

	def __or__(self, other):
		return self.copy().__ior__(other)

	__iadd__ = __ior__
	__add__ = __or__

	def __isub__(self, other):
		for key, value in other.iteritems():
			if key in self:
				self[key] -= value
		return self

	def __sub__(self, other):
		return self.copy().__isub__(other)

	def __invert__(self):
		new = self.copy()
		for key, value in new.iteritems():
			dict.__setitem__(new, key, ~value)
		return new

	# other list-by-list operations

	def intersects_segment(self, seg):
		"""
		Returns True if any segmentlist in self intersects the
		segment, otherwise returns False.
		"""
		for value in self.itervalues():
			if value.intersects_segment(seg):
				return True
		return False

	def intersects_all(self, other):
		"""
		Returns True if each segmentlist in other intersects the
		corresponding segmentlist in self;  returns False
		otherwise.
		"""
		for key, value in other.iteritems():
			if key not in self or not self[key].intersects(value):
				return False
		return True

	def all_intersects(self, other):
		"""
		Returns True if each segmentlist in self intersects the
		corresponding segmentlist in other;  returns False
		otherwise.
		"""
		for key, value in self.iteritems():
			if key not in other or not other[key].intersects(value):
				return False
		return True

	# FIXME: deprecate this
	intersects = intersects_all

	def coalesce(self):
		"""
		Run coalesce() on all segmentlists.
		"""
		for value in self.itervalues():
			value.coalesce()
		return self

	def contract(self, x):
		"""
		Run contract(x) on all segmentlists.
		"""
		for value in self.itervalues():
			value.contract(x)
		return self

	def protract(self, x):
		"""
		Run protract(x) on all segmentlists.
		"""
		for value in self.itervalues():
			value.protract(x)
		return self

	def extract_common(self, keys):
		"""
		Return a new segmentlistdict containing only those
		segmentlists associated with the keys in keys, with each
		set to their mutual intersection.  The offsets are
		preserved.
		"""
		new = segmentlistdict()
		intersection = self.intersection(keys)
		for key in keys:
			dict.__setitem__(new, key, shallowcopy(intersection))
			dict.__setitem__(new.offsets, key, self.offsets[key])
		return new

	# multi-list operations

	def is_coincident(self, other, keys = None):
		"""
		Return True if any segment in any list in self intersects
		any segment in any list in other.  If the optional keys
		argument is not None, then only segment lists for the given
		keys are considered.  Keys not represented in both segment
		lists are ignored.  If keys is None (the default) then all
		segment lists are considered.
		"""
		keys1 = set(self.iterkeys())
		keys2 = set(other.iterkeys())
		if keys is not None:
			keys = set(keys)
			keys1 &= keys
			keys2 &= keys
		for key1 in keys1:
			l1 = self[key1]
			for key2 in keys2:
				if l1.intersects(other[key2]):
					return True
		return False

	def intersection(self, keys):
		"""
		Return the intersection of the segmentlists associated with
		the keys in keys.
		"""
		if not keys:
			return segmentlist()
		seglist = ~segmentlist()
		for key in keys:
			seglist &= self[key]
		return seglist

	def union(self, keys):
		"""
		Return the union of the segmentlists associated with the
		keys in keys.
		"""
		seglist = segmentlist()
		for key in keys:
			seglist |= self[key]
		return seglist

