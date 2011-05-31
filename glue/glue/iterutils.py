# Copyright (C) 2007  Kipp Cannon, Nickolas Fotopoulos
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
A collection of iteration utilities.
"""


from glue import git_version


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date


#
# =============================================================================
#
#                               Iteration Tools
#
# =============================================================================
#


def MultiIter(*sequences):
	"""
	A generator for iterating over the elements of multiple sequences
	simultaneously.  With N sequences given as input, the generator
	yields all possible distinct N-tuples that contain one element from
	each of the input sequences.

	Example:

	>>> x = MultiIter([0, 1, 2], [10, 11])
	>>> list(x)
	[(0, 10), (1, 10), (2, 10), (0, 11), (1, 11), (2, 11)]

	The elements in each output tuple are in the order of the input
	sequences, and the left-most input sequence is iterated over first.

	The input sequences are each iterated over only once, so it is safe
	to pass generators as arguments.  Also, this generator is
	significantly faster if the longest input sequence is given as the
	first argument.  For example, this code

	>>> lengths = range(1, 12)
	>>> for x in MultiIter(*map(range, lengths)):
	...	pass
	...

	runs approximately 5 times faster if the lengths list is reversed.
	"""
	if len(sequences) > 1:
		# FIXME:  this loop is about 5% faster if done the other
		# way around, if the last list is iterated over in the
		# inner loop.  but there is code, like snglcoinc.py in
		# pylal, that has been optimized for the current order and
		# would need to be reoptimized if this function were to be
		# reversed.
		head = tuple((x,) for x in sequences[0])
		for t in MultiIter(*sequences[1:]):
			for h in head:
				yield h + t
	elif sequences:
		for t in sequences[0]:
			yield (t,)


def choices(vals, n):
	"""
	A generator for iterating over all choices of n elements from the
	input sequence vals.  In each result returned, the original order
	of the values is preserved.

	Example:

	>>> x = choices(["a", "b", "c"], 2)
	>>> list(x)
	[('a', 'b'), ('a', 'c'), ('b', 'c')]

	The order of combinations in the output sequence is always the
	same, so if choices() is called twice with two different sequences
	of the same length the first combination in each of the two output
	sequences will contain elements from the same positions in the two
	different input sequences, and so on for each subsequent pair of
	output combinations.

	Example:

	>>> x = choices(["a", "b", "c"], 2)
	>>> y = choices(["1", "2", "3"], 2)
	>>> zip(x, y)
	[(('a', 'b'), ('1', '2')), (('a', 'c'), ('1', '3')), (('b', 'c'), ('2', '3'))]

	Furthermore, the order of combinations in the output sequence is
	such that if the input list has n elements, and one constructs the
	combinations choices(input, m), then each combination in
	choices(input, n-m).reverse() contains the elements discarded in
	forming the corresponding combination in the former.

	Example:

	>>> x = ["a", "b", "c", "d", "e"]
	>>> X = list(choices(x, 2))
	>>> Y = list(choices(x, len(x) - 2))
	>>> Y.reverse()
	>>> zip(X, Y)
	[(('a', 'b'), ('c', 'd', 'e')), (('a', 'c'), ('b', 'd', 'e')), (('a', 'd'), ('b', 'c', 'e')), (('a', 'e'), ('b', 'c', 'd')), (('b', 'c'), ('a', 'd', 'e')), (('b', 'd'), ('a', 'c', 'e')), (('b', 'e'), ('a', 'c', 'd')), (('c', 'd'), ('a', 'b', 'e')), (('c', 'e'), ('a', 'b', 'd')), (('d', 'e'), ('a', 'b', 'c'))]
	"""
	if n == len(vals):
		yield tuple(vals)
	elif n > 1:
		n -= 1
		for i, v in enumerate(vals[:-n]):
			v = (v,)
			for c in choices(vals[i+1:], n):
				yield v + c
	elif n == 1:
		for v in vals:
			yield (v,)
	elif n == 0:
		yield ()
	else:
		# n < 0
		raise ValueError, n


def uniq(iterable):
	"""
	Yield the unique items of an iterable, preserving order.
	http://mail.python.org/pipermail/tutor/2002-March/012930.html

	Example:

	>>> x = uniq([0, 0, 2, 6, 2, 0, 5])
	>>> list(x)
	[0, 2, 6, 5]
	"""
	temp_dict = {}
	for e in iterable:
		if e not in temp_dict:
			yield temp_dict.setdefault(e, e)

def flatten(sequence, levels = 1):
	"""
	Example:
	>>> nested = [[1,2], [[3]]]
	>>> list(flatten(nested))
	[1, 2, [3]]
	"""
	if levels == 0:
		for x in sequence:
			yield x
	else:
		for x in sequence:
			for y in flatten(x, levels - 1):
				yield y

#
# =============================================================================
#
#    any() and all() are built-ins in Python 2.5, but I don't want to wait.
#
# =============================================================================
#


try:
	any = any
	all = all
except NameError:
	# These short-circuit, returning as soon as the return value can be
	# determined.  These are a factor of a few slower than Python 2.5's
	# implementation.
	def any(S):
		"""
		any(iterable) -> bool

		Return True if bool(x) is True for any x in the iterable.
		"""
		for x in S:
			if x: return True
		return False
	def all(S):
		"""
		all(iterable) -> bool

		Return True if bool(x) is True for all values x in the iterable.
		"""
		for x in S:
			if not x: return False
		return True


#
# =============================================================================
#
#    itertools.groupby() was added in Python 2.4, but I don't want to wait.
#
# =============================================================================
#

try:
	from itertools import groupby
except ImportError:
	class groupby(object):
		"""
		Python 2.3 compatibility: reimplement itertools.groupby.
		Taken directly from the itertools documentation.
		"""
		def __init__(self, iterable, key=None):
			if key is None:
				key = lambda x: x
			self.keyfunc = key
			self.it = iter(iterable)
			self.tgtkey = self.currkey = self.currvalue = xrange(0)
		def __iter__(self):
			return self
		def next(self):
			while self.currkey == self.tgtkey:
				self.currvalue = self.it.next() # Exit on StopIteration
				self.currkey = self.keyfunc(self.currvalue)
			self.tgtkey = self.currkey
			return (self.currkey, self._grouper(self.tgtkey))
		def _grouper(self, tgtkey):
			while self.currkey == tgtkey:
				yield self.currvalue
				self.currvalue = self.it.next() # Exit on StopIteration
				self.currkey = self.keyfunc(self.currvalue)


#
# =============================================================================
#
#                              In-Place filter()
#
# =============================================================================
#


def inplace_filter(func, sequence):
	"""
	Like Python's filter() builtin, but modifies the sequence in place.

	Example:

	>>> l = range(10)
	>>> inplace_filter(lambda x: x > 5, l)
	>>> l
	[6, 7, 8, 9]
	"""
	target = 0
	for source in xrange(len(sequence)):
		if func(sequence[source]):
			sequence[target] = sequence[source]
			target += 1
	del sequence[target:]


#
# =============================================================================
#
#          Return the Values from Several Ordered Iterables in Order
#
# =============================================================================
#


def inorder(*iterables, **kwargs):
	"""
	A generator that yields the values from several ordered iterables
	in order.

	Example:

	>>> x = [0, 1, 2, 3]
	>>> y = [1.5, 2.5, 3.5, 4.5]
	>>> z = [1.75, 2.25, 3.75, 4.25]
	>>> list(inorder(x, y, z))
	[0, 1, 1.5, 1.75, 2, 2.25, 2.5, 3, 3.5, 3.75, 4.25, 4.5]

	>>> x = [3, 2, 1, 0]
	>>> y = [4.5, 3.5, 2.5, 1.5]
	>>> z = [4.25, 3.75, 2.25, 1.75]
	>>> list(inorder(x, y, z, reverse = True))
	[4.5, 4.25, 3.75, 3.5, 3, 2.5, 2.25, 2, 1.75, 1.5, 1, 0]

	NOTE:  this function will never reverse the order of elements in
	the input iterables.  If the reverse keyword argument is False (the
	default) then the input sequences must yield elements in increasing
	order, likewise if the keyword argument is True then the input
	sequences must yield elements in decreasing order.  Failure to
	adhere to this yields undefined results, and for performance
	reasons no check is performed to validate the element order in the
	input sequences.
	"""
	reverse = False
	if kwargs:
		try:
			reverse = kwargs.pop("reverse")
		except KeyError:
			pass
		if kwargs:
			raise TypeError, "invalid keyword argument '%s'" % kwargs.keys()[0]
	nextvals = {}
	for iterable in iterables:
		next = iter(iterable).next
		try:
			nextvals[next] = next(), next
		except StopIteration:
			pass
	if not nextvals:
		# all sequences are empty
		return
	if reverse:
		select = max
	else:
		select = min
	values = nextvals.itervalues
	if len(nextvals) > 1:
		while True:
			val, next = select(values())
			yield val
			try:
				nextvals[next] = next(), next
			except StopIteration:
				del nextvals[next]
				if len(nextvals) < 2:
					break
	# exactly one sequence remains, short circuit and drain it
	(val, next), = values()
	yield val
	while True:
		yield next()


#
# =============================================================================
#
#           Thing for keeping only the highest values in a sequence
#
# =============================================================================
#


class Highest(list):
	"""
	A class for use when you need to collect the largest in a very long
	sequence of things, a sequence too long to hold in memory all at
	once and sort.  This class behaves like a list, in fact it is a
	Python list, but one that stores only some fraction of all items
	that have been added to it.  The list is always sorted in
	decreasing order, so the 0th element is the highest-valued item
	added to the list.

	For this class to function correctly the list must remain sorted,
	so the insert() and __setitem__() methods are disabled.  Only the
	append() and extend() methods can be used to add elements to the
	list.  The __len__() method returns the total number of items that
	have been added to the list ever, not the number that are actually
	stored in it at any given time.  To retrieve the number of items
	actually stored in memory, use the list class' __len__() method.
	See the example below.

	Example:

	>>> l = Highest(max = 3)
	>>> for i in range(10000):
	...	l.append(i)
	...
	>>> l
	[9999, 9998, 9997]
	>>> len(l)
	10000
	>>> list.__len__(l)
	3
	"""
	def __init__(self, sequence = tuple(), max = None):
		list.__init__(self, sequence)
		self.n = list.__len__(self)
		self.max = int(max)
		list.sort(self, reverse = True)
		del self[self.max:]

	def __len__(self):
		return self.n

	def __reduce__(self):
		return (type(self), ((), self.max), {"n": self.n}, iter(self))

	def append(self, value):
		# can't use bisect module because list is sorted in reverse
		# order
		hi = list.__len__(self)
		lo = 0
		while lo < hi:
			mid = (lo + hi) // 2
			if value > self[mid]:
				hi = mid
			else:
				lo = mid + 1
		list.insert(self, lo, value)
		self.n += 1
		del self[self.max:]

	def __iadd__(self, other):
		"""
		In-place addition.  Add the contents of another Highest
		object to this one.  This method only works correctly with
		other Highest objects.

		Example:

		>>> x = Highest(max = 10)
		>>> x.extend(range(50))
		>>> x
		[49, 48, 47, 46, 45, 44, 43, 42, 41, 40]
		>>> y = Highest(max = 5)
		>>> y.append(100)
		>>> y.append(49)
		>>> y
		[100, 49]
		>>> x += y
		>>> x
		[100, 49, 49, 48, 47]
		>>> len(x)
		52
		"""
		assert isinstance(other, type(self))
		self.n += other.n
		# we can only claim to provide the minimum of the two max's
		# worth of the highest valued objects
		self.max = min(self.max, other.max)
		# FIXME:  since the lists are sorted, lots could be done to
		# speed this up
		i = j = 0
		while i + j < self.max:
			if i >= list.__len__(self):
				j = min(self.max - i, list.__len__(other))
				break
			if j >= list.__len__(other):
				i = min(self.max - j, list.__len__(self))
				break
			if self[i] > other[j]:
				i += 1
			elif self[i] < other[j]:
				j += 1
			else:
				# self[i] == other[j]
				i += 1
				j += 1
		self[i:] = other[:j]
		list.sort(self, reverse = True)
		del self[self.max:]
		return self

	def extend(self, sequence):
		# this method will not work correctly with other Highest
		# objects
		assert not isinstance(sequence, type(self))
		# FIXME:  could be implemented with heapq.merge() when 2.6
		# is generally available
		# n is updated in a way that allows sequence to be a
		# generator or other object without a __len__() method
		before = list.__len__(self)
		list.extend(self, sequence)
		self.n += list.__len__(self) - before
		list.sort(self, reverse = True)
		del self[self.max:]

	#
	# Stubs to prevent bugs
	#

	def __setitem__(*args, **kwargs):
		raise NotImplementedError

	def reverse(*args, **kwargs):
		raise NotImplementedError

	def remove(*args, **kwargs):
		raise NotImplementedError

	def pop(*args, **kwargs):
		raise NotImplementedError

	def insert(*args, **kwargs):
		raise NotImplementedError

	def index(*args, **kwargs):
		raise NotImplementedError

	def count(*args, **kwargs):
		raise NotImplementedError

	def sort(*args, **kwargs):
		raise NotImplementedError
