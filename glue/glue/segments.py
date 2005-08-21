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

__author__ = "Kipp Cannon <kipp@gravity.phys.uwm.edu>"
__date__ = "$Date$"
__version__ = "$Revision$"

import bisect


#
# Infinity for use with segments
#

class infinity:
	"""
	The infinity object possess the algebraic properties necessary for use
	as a bound on semi-infinite and infinite segments.

	Example usage:
		x = infinity()
		x > 0
		x + 10
		segment(-10, 10) - segment(-x, 0)
	"""

	def __init__(self):
		self.__sign = 1
	
	def __repr__(self):
		if self.__sign > 0:
			return "infinity"
		return "-infinity"
	
	__str__ = __repr__
	
	def __nonzero__(self):
		return True
	
	def __cmp__(self, other):
		if type(self) != type(other):
			return self.__sign
		return cmp(self.__sign, other.__sign)
	
	def __add__(self, other):
		return self
	def __radd__(self, other):
		return self
	
	def __sub__(self, other):
		if type(self) != type(other):
			return self
		if self.__sign != other.__sign:
			return self
		return None
	def __rsub__(self, other):
		if type(self) != type(other):
			return -self
		if self.__sign != other.__sign:
			return -self
		return None

	def __neg__(self):
		x = infinity()
		x.__sign = -self.__sign
		return x

	def __pos__(self):
		return self


#
# The segment class
#

class segment(tuple):
	"""
	The segment class defines objects that represent a range of values.  A
	segment has a start and an end, and is taken to represent the range of
	values from the start to the end inclusively.  Some limited arithmetic
	operations are possible with segments, but because the set of (single)
	segments is not closed under the sensible definitions of the standard
	arithmetic operations, the behaviour of the arithmetic operators on
	segments may not be as you would expect.  For general arithmetic on
	segments, use segmentlist objects.  The methods for this class exist
	mostly for purpose of simplifying the implementation of the segmentlist
	class.

	Example use:
		segment(0, 10) & segment(5, 15)
		segment(0, 10) | segment(5, 15)
		segment(0, 10) - segment(5, 15)
		segment(0, 10) < segment(5, 15)
		segment(1, 2) in segment(0, 10)
		bool(segment(0, 0))

	Notes:
	It is also possible to cast 2-element tuples, lists, and other
	container types to segments.  For example segment([0, 1])
	"""

	# basic class methods

	def __new__(cls, *args):
		if len(args) == 1:
			args = args[0]
		if len(args) != 2:
			raise TypeError, "__new__() takes 3 arguments or 2 arguments when the second is a 2-element container type"
		if args[0] <= args[1]:
			return tuple.__new__(cls, args)
		else:
			return tuple.__new__(cls, (args[1], args[0]))

	def __repr__(self):
		return "segment(" + str(self[0]) + ", " + str(self[1]) + ")"

	def __str__(self):
		return str(self[0]) + "_" + str(self[1])

	# accessors

	def duration(self):
		"""
		Returns the length of the interval represented by the segment.
		"""
		return self[1] - self[0]

	# comparisons

	def __nonzero__(self):
		"""
		Test for segment having non-zero duration.
		"""
		return self[0] != self[1]

	def order(self, other):
		"""
		Equivalent to cmp(self, other) except that a result of 0
		indicates that other is contained in self rather than being
		identically equal to self.
		"""
		if other in self:
			return 0
		return cmp(self, other)

	# some arithmetic operations that (mostly) make sense for segments

	def __and__(self, other):
		"""
		Return the segment that is the intersection of the given
		segments, or None if the segments do not intersect.
		"""
		if not self.intersects(other):
			return None
		return segment(max(self[0], other[0]), min(self[1], other[1]))

	def __or__(self, other):
		"""
		Return the segment that is the union of the given segments, or
		None if the result cannot be represented as a single segment.
		"""
		if not self.continuous(other):
			return None
		return segment(min(self[0], other[0]), max(self[1], other[1]))

	# addition is defined to be the union operation
	__add__ = __or__

	def __sub__(self, other):
		"""
		Return the segment that is that part of self which is not
		contained in other, or None if the result cannot be represented
		as a single segment.
		"""
		if not self.intersects(other):
			return self
		if (self in other) or ((self[0] < other[0]) and (self[1] > other[1])):
			return None
		if self[0] < other[0]:
			return segment(self[0], other[0])
		return segment(other[1], self[1])

	# check for proper intersection, containment, and continuity

	def intersects(self, other):
		"""
		Return True if the intersection of self and other is not a null
		segment.
		"""
		return (self[1] > other[0]) and (self[0] < other[1])

	def __contains__(self, other):
		"""
		Return True if other is wholly contained in self.
		"""
		return (self[0] <= other[0]) and (self[1] >= other[1])

	def continuous(self, other):
		"""
		Return True if self and other are not disjoint.
		"""
		return (self[1] >= other[0]) and (self[0] <= other[1])

	# protraction and contraction

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


#
# A segment list class derived from the builtin list class.
#

class segmentlist(list):
	"""
	The segmentlist class defines a list of segments, and is an extension
	of the built-in list class.  This class provides addtional methods that
	assist in the manipulation of lists of segments.  In particular,
	arithmetic operations such as union and intersection are provided.
	Unlike the segment class, the segmentlist class is closed under all
	supported arithmetic operations.

	All standard Python sequence-like operations are supported, like
	slicing, iteration and so on, but the arithmetic and other methods in
	this class generally expect the segmentlist to be in what is refered to
	as a "coalesced" state --- consisting solely of disjoint segments
	listed in ascending order.  Using the standard Python sequence-like
	operations, a segmentlist can be easily constructed that is not in this
	state;  for example by simply appending a segment to the end of the
	list that overlaps some other segment already in the list.  The class
	provides a coalesce() method that can be called to put it in the
	coalesced state.  Following application of the coalesce method, all
	arithmetic operations will function reliably.  All arithmetic methods
	themselves return coalesced results, so there is never a need to call
	the coalesce method when manipulating segmentlists exclusively via the
	arithmetic operators.

	Example use:
		x = segmentlist([segment(-10, 10)])
		x |= segmentlist([segment(20, 30)])
		x -= segmentlist([segment(-5, 5)])
		print x
		print ~x
	"""

	# container method over-rides.

	def __contains__(self, item):
		"""
		Returns True if the given segment is wholly contained within
		one of the segments in self.  Does not require the segmentlist
		to be coalesced.

		Note the difference between this operator, and the standard
		Python "in" operator for sequence-like objects:  in the case of
		standard sequence-like objects the in operator checks for an
		exact match between the given item and one of the contents of
		the list; for segmentlists, the in operator checks if the given
		item is contained within any of the segments in the
		segmentlist.
		"""
		for seg in self:
			if item in seg:
				return True
		return False

	# suplementary accessors

	def duration(self):
		"""
		Return the sum of the durations of all segments in self.  Does
		not require the segmentlist to be coalesced.
		"""
		d = 0
		for seg in self:
			d += seg.duration()
		return d

	def extent(self):
		"""
		Return the segment whose end-points denote the maximum and
		minimum extent of the segmentlist.  Does not require the
		segmentlist to be coalesced.
		"""
		if not len(self):
			raise ValueError, "segmentlist.extent(): empty list"
		(min, max) = self[0]
		for seg in self:
			if min > seg[0]:
				min = seg[0]
			if max < seg[1]:
				max = seg[1]
		return segment(min, max)

	def find(self, item):
		"""
		Return the smallest i such that i is the index of an element
		that wholly contains the given segment.  Raises ValueError if
		no such element exists.  Does not require the segmentlist to be
		coalesced.
		"""
		for i, seg in enumerate(self):
			if item in seg:
				return i
		raise ValueError, "segmentlist.find(x): x not contained in segmentlist"

	# arithmetic operations that are sensible with segment lists

	def __iand__(self, other):
		"""
		Replace the segmentlist with the intersection of itself and
		another.  This operation is O(n).
		"""
		self -= self - other
		return self

	def __and__(self, other):
		"""
		Return the intersection of the segmentlist and another.  This
		operation is O(n).
		"""
		x = segmentlist(self[:])
		x &= other
		return x

	def __ior__(self, other):
		"""
		Replace the segmentlist with the union of itself and another.
		If the two lists have numbers of elements m and n respectively,
		then this algorithm is O(n log m), which means it is optimized
		for the case when the latter list contains a small number of
		segments.  If you have two large lists of n elements each, then
		it is faster to do

			list1.extend(list2)
			list1.coalesce()

		which is also O(n log n), but with a smaller (by about 25%)
		leading coefficient.
		"""
		low = 0
		for seg in other:
			low = bisect.bisect_right(self, seg, low)
			self.insert(low, seg)
			if low:
				if self[low-1].continuous(self[low]):
					self[low-1:low+1] = [ self[low-1] | self[low] ]
					low -= 1
			try:
				while self[low].continuous(self[low+1]):
					self[low:low+2] = [ self[low] | self[low+1] ]
			except IndexError:
				pass
			low += 1
		return self

	def __or__(self, other):
		"""
		Return the union of the segment list and another.  The comment
		in the segmentlist.__ior__() method about performance
		optimization when computing the union of large lists of similar
		size applies here as well.
		"""
		if len(self) >= len(other):
			x = segmentlist(self[:])
			x |= other
		else:
			x = segmentlist(other[:])
			x |= self
		return x

	def __xor__(self, other):
		"""
		Return the segmentlist that is the list of all intervals
		contained in exactly one of this and another list.  This
		operation is O(n log n).
		"""
		return (self - other) | (other - self)

	# addition is defined to be the union operation
	__iadd__ = __ior__
	__add__ = __or__

	def __isub__(self, other):
		"""
		Replace the segmentlist with the difference between itself and
		another.  This operation is O(n).
		"""
		try:
			other_it = iter(other)
			other_seg = other_it.next()
			i = 0
			while 1:
				seg = self[i]
				while (seg[0] >= other_seg[1]) or not bool(other_seg):
					other_seg = other_it.next()
				if seg[1] <= other_seg[0]:
					i += 1
				elif seg in other_seg:
					self[i:i+1] = []
				elif (other_seg[0] > seg[0]) and (other_seg[1] < seg[1]):
					self[i:i+1] = [segment(seg[0], other_seg[0]), segment(other_seg[1], seg[1])]
					i += 1
				else:
					self[i] -= other_seg
		except (StopIteration, IndexError):
			pass
		return self

	def __sub__(self, other):
		"""
		Return the difference between the segmentlist and another.
		This operation is O(n).
		"""
		x = segmentlist(self[:])
		x -= other
		return x

	def __invert__(self):
		"""
		Return the segmentlist that is the inversion of the given list.
		This operation is O(n).
		"""
		return segmentlist([segment(-infinity(), infinity())]) - self

	# other operations

	def split(self, value):
		"""
		Break all segments that stradle the given value at that value.
		Does not require the segmentlist to be coalesced, and the
		result is not coalesced by definition.  This operation is O(n).
		"""
		for i, seg in enumerate(self):
			if (seg[0] < value) and (seg[1] > value):
				self[i:i+1] = [segment(seg[0], value), segment(value, seg[1])]

	def coalesce(self):
		"""
		Sort the elements of a list into ascending order, and merge
		continuous segments into single segments.  This operation is
		O(n log n), and is dominated by the sort.
		"""
		self.sort()
		try:
			for i in range(len(self) - 1):
				while self[i].continuous(self[i+1]):
					self[i:i+2] = [ self[i] | self[i+1] ]
		except IndexError:
			pass
		return self

	def protract(self, x):
		"""
		For each segment in the list, move both the start and the end a
		distance x away from the other.  Coalesce the result.
		"""
		return segmentlist([seg.protract(x) for seg in self]).coalesce()

	def contract(self, x):
		"""
		For each segment in the list, move both the start and the end a
		distance x towards the other.  Coalesce the result.
		"""
		return segmentlist([seg.contract(x) for seg in self]).coalesce()
