#
# NOTE:  the logic in this code is unintuitively complicated.  Small,
# apparently irrelevant, changes to conditionals can have subtly unexpected
# consequences to the behaviour of the class methods.  ALWAYS make sure that
# the module_test() function returns pass on ALL tests after any changes you
# make.
#

"""
This module defines the segment and segmentlist objects, as well as the
infinity object used to define semi-infinite and infinite segments.
"""

__author__ = "Kipp Cannon <kipp@gravity.phys.uwm.edu>"
__date__ = "$Date:"
__version__ = "$Revision:"


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


#
# The segment class
#

class segment:
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
	"""

	# basic class methods

	def __init__(self, start, end):
		if start <= end:
			self.start, self.end = start, end
		else:
			self.start, self.end = end, start

	def __repr__(self):
		return "segment(" + str(self.start) + ", " + str(self.end) + ")"

	def __str__(self):
		return str(self.start) + "_" + str(self.end)

	# accessors

	def get_start(self):
		return self.start

	def get_end(end):
		return self.end

	def duration(self):
		"""
		Returns the length of the interval represented by the segment.
		"""
		return self.end - self.start

	# comparisons

	def __cmp__(self, other):
		"""
		Compare one segment to another using the rules for comparing
		two two-element tuples
		"""
		if type(self) != type(other):
			return -1
		return cmp((self.start, self.end), (other.start, other.end))
	
	def __nonzero__(self):
		"""
		Test for segment having non-zero duration.
		"""
		return self.start != self.end

	# some arithmetic operations that (mostly) make sense for segments

	def __and__(self, other):
		"""
		Return the segment that is the intersection of the given
		segments, or None if the segments do not intersect.
		"""
		if not self.intersects(other):
			return None
		return segment(max(self.start, other.start), min(self.end, other.end))

	def __or__(self, other):
		"""
		Return the segment that is the union of the given segments, or
		None if the result cannot be represented as a single segment.
		"""
		if not self.continuous(other):
			return None
		return segment(min(self.start, other.start), max(self.end, other.end))

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
		if (self in other) or ((self.start < other.start) and (self.end > other.end)):
			return None
		if self.start < other.start:
			return segment(self.start, other.start)
		return segment(other.end, self.end)

	# check for proper intersection, containment, and continuity

	def intersects(self, other):
		"""
		Return True if the intersection of self and other is not a null
		segment.
		"""
		return (self.end > other.start) and (self.start < other.end)

	def __contains__(self, other):
		"""
		Return True if other is wholly contained in self.
		"""
		return (self.start <= other.start) and (self.end >= other.end)

	def continuous(self, other):
		"""
		Return True if self and other are not disjoint.
		"""
		return (self.end >= other.start) and (self.start <= other.end)


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

	Note the difference between the standard Python "in" operator for
	sequence-like objects, and the same operator for segmentlists:  in the
	case of standard sequence-like objects the in operator checks for an
	exact match between the given item and one of the contents of the list;
	for segmentlists, the in operator checks if the given item (a segment)
	is contained within any of the segments in the segmentlist.

	The union operation, and the coalesce method which it uses, are O(n log
	n).  The subtraction operation is O(n^2), as are all other arithmetic
	operations (as they are implemented on top of the subtraction
	operation).

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
		Returns True if item (a segment) is wholly contained within one
		of the segments in self.
		"""
		for seg in self:
			if item in seg:
				return True
		return False

	# suplementary accessors

	def duration(self):
		"""
		Return the sum of the durations of all segments in self.
		"""
		return sum([seg.duration() for seg in self])
	
	def extent(self):
		"""
		Return the segment whose end-points denote the maximum and
		minimum extent of the segmentlist self.
		"""
		return segment(min([seg.start for seg in self]), max([seg.end for seg in self]))

	# arithmetic operations that are sensible with segment lists

	def __iand__(self, other):
		"""
		Replace the segmentlist with the intersection of itself and
		another.
		"""
		self -= self - other
		return self

	def __and__(self, other):
		"""
		Return the intersection of the segmentlist and another.
		"""
		x = segmentlist(self[:])
		x &= other
		return x

	def __ior__(self, other):
		"""
		Replace the segmentlist with the union of itself and another.
		"""
		self.extend(other)
		self.coalesce()
		return self

	def __or__(self, other):
		"""
		Return the union of the segment list and another.
		"""
		x = segmentlist(self[:])
		x |= other
		return x

	def __xor__(self, other):
		"""
		Return the segmentlist that is the list of all intervals
		contained in exactly one of this and another list.
		"""
		return (self | other) - (self & other)

	# addition is defined to be the union operation
	__iadd__ = __ior__
	__add__ = __or__

	def __isub__(self, other):
		"""
		Replace the segmentlist with the difference between itself and
		another.
		"""
		for b in other:
			self.split(b.start)
		try:
			i = 0
			for b in other:
				while self[i].end <= b.start:
					i += 1
				while self[i] in b:
					self[i:i+1] = []
				while self[i].start < b.end:
					self[i] -= b
					i += 1
		except IndexError:
			pass
		return self

	def __sub__(self, other):
		"""
		Return the difference between the segmentlist and another.
		"""
		x = segmentlist(self[:])
		x -= other
		return x

	def __invert__(self):
		"""
		Return the segmentlist that is the inversion of the given list.
		"""
		return segmentlist([segment(-infinity(), infinity())]) - self

	def split(self, value):
		"""
		Break all segments that stradle the given value at that value.
		Note that the result is, by definition, not coalesced.
		"""
		try:
			i = 0
			while 1:
				if self[i].intersects(segment(value,value)):
					self[i:i+1] = [segment(self[i].start, value), segment(value, self[i].end)]
					i += 1
				i += 1
		except IndexError:
			pass
		return

	def coalesce(self):
		"""
		Sort the elements of a list into ascending order, and merge
		continuous segments into single segments.
		"""
		self.sort()
		try:
			for i in range(len(self) - 1):
				while self[i].continuous(self[i+1]):
					self[i:i+2] = [ self[i] | self[i+1] ]
		except IndexError:
			pass


#
# segments module verification tests
#

def module_verify():
	print "=== test infinity class"
	def test_inf_cmp(r,a,b):
		if cmp(a, b) == r:
			s = "pass:  "
		else:
			s = "FAIL:  "
		print s + "cmp(" + str(a) + ", " + str(b) + ") = " + str(cmp(a, b))
	def test_inf_add(r,a,b):
		if a + b == r:
			s = "pass:  "
		else:
			s = "FAIL:  "
		print s + str(a) + " + " + str(b) + " = " + str(a + b)
	def test_inf_sub(r,a,b):
		if a - b == r:
			s = "pass:  "
		else:
			s = "FAIL:  "
		print s + str(a) + " - " + str(b) + " = " + str(a - b)
	a = infinity()
	print str(a) + " " + str(-a)
	print str(a) + " " + str(-a)
	test_inf_cmp(-1, 0, a)
	test_inf_cmp(1, 0, -a)
	test_inf_cmp(1, a, 0)
	test_inf_cmp(-1, -a, 0)
	test_inf_cmp(1, a, -a)
	test_inf_cmp(-1, -a, a)
	test_inf_cmp(0, a, a)
	test_inf_cmp(0, -a, -a)
	test_inf_add(infinity(), a, 10)
	test_inf_add(infinity(), a, -10)
	test_inf_add(-infinity(), -a, 10)
	test_inf_add(-infinity(), -a, -10)
	test_inf_add(infinity(), 10, a)
	test_inf_add(infinity(), -10, a)
	test_inf_add(-infinity(), 10, -a)
	test_inf_add(-infinity(), -10, -a)
	test_inf_add(a, a, a)
	test_inf_add(-a, -a, -a)
	test_inf_sub(a, a, 10)
	test_inf_sub(a, a, -10)
	test_inf_sub(-a, -a, 10)
	test_inf_sub(-a, -a, -10)
	test_inf_sub(-a, 10, a)
	test_inf_sub(-a, -10, a)
	test_inf_sub(a, 10, -a)
	test_inf_sub(a, -10, -a)
	test_inf_sub(None, a, a)
	test_inf_sub(a, a, -a)
	test_inf_sub(-a, -a, a)
	test_inf_sub(None, -a, -a)
	
	print "=== test segment definition"
	print "segment(-2, 2) = " + str(segment(-2, 2))
	print "segment(2, -2) = " + str(segment(2, -2))
	print "segment(-infinity, 2) = " + str(segment(-infinity(), 2))
	print "segment(2, -infinity) = " + str(segment(2, -infinity()))
	print "segment(infinity, 2) = " + str(segment(infinity(), 2))
	print "segment(2, infinity) = " + str(segment(2, infinity()))
	print "segment(-inf, inf) = " + str(segment(-infinity(), infinity()))

	print "segment(-2,2).duration() = " + str(segment(-2, 2).duration())
	print "segment(0, infinity).duration() = " + str(segment(0, infinity()).duration())
	print "segment(-infinity, 0).duration() = " + str(segment(-infinity(), 0).duration())
	print "segment(-infinity, infinity).duration() = " + str(segment(-infinity(), infinity()).duration())

	print "=== test segment comparisons"
	def test_intersects(r,a,b):
		if a.intersects(b) == r:
			s = "pass:  "
		else:
			s = "FAIL:  "
		print s + str(a) + " intersects with " + str(b) + " = " + str(a.intersects(b))
	test_intersects(False, segment(-2, 2), segment(-4, -3))
	test_intersects(False, segment(-2, 2), segment(-4, -2))
	test_intersects(True, segment(-2, 2), segment(-4, 0))
	test_intersects(True, segment(-2, 2), segment(-4, 2))
	test_intersects(True, segment(-2, 2), segment(-4, 4))
	test_intersects(True, segment(-2, 2), segment(-2, 4))
	test_intersects(True, segment(-2, 2), segment(0, 4))
	test_intersects(False, segment(-2, 2), segment(2, 4))
	test_intersects(False, segment(-2, 2), segment(3, 4))
	test_intersects(True, segment(-2, 2), segment(-2, 2))
	test_intersects(True, segment(-2, 2), segment(-1, 1))

	def test_contains(r,a,b):
		if (b in a) == r:
			s = "pass:  "
		else:
			s = "fail:  "
		print s + str(b) + " in " + str(a) + " = " + str(b in a)
	test_contains(False, segment(-2, 2), segment(-4, -3))
	test_contains(False, segment(-2, 2), segment(-4, -2))
	test_contains(False, segment(-2, 2), segment(-4, 0))
	test_contains(False, segment(-2, 2), segment(-4, 2))
	test_contains(False, segment(-2, 2), segment(-4, 4))
	test_contains(False, segment(-2, 2), segment(-2, 4))
	test_contains(False, segment(-2, 2), segment(0, 4))
	test_contains(False, segment(-2, 2), segment(2, 4))
	test_contains(False, segment(-2, 2), segment(3, 4))
	test_contains(True, segment(-2, 2), segment(-2, 2))
	test_contains(True, segment(-2, 2), segment(-1, 1))

	def test_continuous(r,a,b):
		if a.continuous(b) == r:
			s = "pass:  "
		else:
			s = "fail:  "
		print s + str(a) + " continuous " + str(b) + " = " + str(a.continuous(b))
	test_continuous(False, segment(-2, 2), segment(-4, -3))
	test_continuous(True, segment(-2, 2), segment(-4, -2))
	test_continuous(True, segment(-2, 2), segment(-4, 0))
	test_continuous(True, segment(-2, 2), segment(-4, 2))
	test_continuous(True, segment(-2, 2), segment(-4, 4))
	test_continuous(True, segment(-2, 2), segment(-2, 4))
	test_continuous(True, segment(-2, 2), segment(0, 4))
	test_continuous(True, segment(-2, 2), segment(2, 4))
	test_continuous(False, segment(-2, 2), segment(3, 4))
	test_continuous(True, segment(-2, 2), segment(-2, 2))
	test_continuous(True, segment(-2, 2), segment(-1, 1))

	print "=== test segment list excision"
	def test_excision(r,a,b):
		if a - b == r:
			s = "pass:  "
		else:
			s = "fail:  "
		print s + str(a) + " - " + str(b) + " = " + str(a - b)
	test_excision(segmentlist([]), segmentlist([]), segmentlist([]))
	test_excision(segmentlist([]), segmentlist([]), segmentlist([segment(-1,1)]))
	test_excision(segmentlist([segment(-1,1)]), segmentlist([segment(-1,1)]), segmentlist([]))
	test_excision(segmentlist([]), segmentlist([segment(-1,1)]), segmentlist([segment(-1,1)]))

	test_excision(segmentlist([segment(0,1)]), segmentlist([segment(0,1)]), segmentlist([segment(2,3)]))
	test_excision(segmentlist([segment(0,1)]), segmentlist([segment(0,1)]), segmentlist([segment(2,3), segment(4,5)]))
	test_excision(segmentlist([segment(0,1)]), segmentlist([segment(0,1), segment(2,3)]), segmentlist([segment(2,3)]))
	test_excision(segmentlist([segment(2,3)]), segmentlist([segment(0,1), segment(2,3)]), segmentlist([segment(0,1)]))
	test_excision(segmentlist([segment(0,1), segment(4,5)]), segmentlist([segment(0,1), segment(2,3), segment(4,5)]), segmentlist([segment(2,3)]))


	test_excision(segmentlist([segment(-5, 10)]), segmentlist([segment(-10,10)]), segmentlist([segment(-15,-5)]))
	test_excision(segmentlist([segment(-10, -5), segment(5, 10)]), segmentlist([segment(-10,10)]), segmentlist([segment(-5,5)]))
	test_excision(segmentlist([segment(-10, 5)]), segmentlist([segment(-10,10)]), segmentlist([segment(5,15)]))

	test_excision(segmentlist([segment(0,5), segment(45,50)]), segmentlist([segment(0,10), segment(20,30), segment(40,50)]), segmentlist([segment(5, 45)]))

	print "=== test segment list inversion"
	def test_inversion(r,a):
		if ~a == r:
			s = "pass:  "
		else:
			s = "fail:  "
		print s + " ~" + str(a) + " = "  + str(~a)
	test_inversion(segmentlist([segment(-infinity(), infinity())]), segmentlist([]))
	test_inversion(segmentlist([]), segmentlist([segment(-infinity(), infinity())]))
	test_inversion(segmentlist([segment(-infinity(), -5), segment(5, infinity())]), segmentlist([segment(-5,5)]))
