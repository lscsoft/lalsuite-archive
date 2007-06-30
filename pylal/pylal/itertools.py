# $Id$
#
# Copyright (C) 2007  Kipp C. Cannon, Nickolas Fotopoulos
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


import math


__author__ = "Kipp Cannon <kipp@gravity.phys.uwm.edu>"
__version__ = "$Revision$"[11:-2]
__date__ = "$Date$"[7:-2]


#
# =============================================================================
#
#                               Iteration Tools
#
# =============================================================================
#


def MultiIter(*lists):
	"""
	A generator for iterating over the elements of multiple lists
	simultaneously.  With N lists given as input, the output sequence
	consists of all possible N-element lists containing one element
	from each of the input lists.

	Example:

	>>> x = MultiIter([0, 1, 2], [10, 11])
	>>> list(x)
	[[0, 10], [1, 10], [2, 10], [0, 11], [1, 11], [2, 11]]
	"""
	if len(lists) > 0:
		head = lists[0]
		for t in MultiIter(*lists[1:]):
			for h in head:
				yield [h] + t
	else:
		yield []


def choices(vals, n):
	"""
	Iterate over all choices of n elements from the list vals.  In each
	result returned, the original order of the values is preserved.

	Example:

	>>> x = choices(["a", "b", "c"], 2)
	>>> list(x)
	[['a', 'b'], ['a', 'c'], ['b', 'c']]
	"""
	if n == len(vals):
		yield vals
	elif n > 1:
		for i, v in enumerate(vals[:1 - n]):
			v = [v]
			for c in choices(vals[i+1:], n - 1):
				yield v + c
	elif n == 1:
		for v in vals:
			yield [v]
	else:
		# n < 1
		raise ValueError, n


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
	sequence of things, too long a sequence to hold in memory all at
	once and sort.  This class behaves like a list, in fact it is a
	Python list, but one that stores only some fixed fraction of all
	items that have been added to it.  The list is always ordered, so
	the insert() and __setitem__() methods are not supported, only
	append() and extend().

	Example:

	>>> import random
	>>> l = Highest(fraction = 0.0002)
	>>> for i in range(10000):
	...	l.append(random.random())
	...
	>>> l
	[0.99999998092197417, 0.99955052111790088]
	>>> len(l)
	10000

	Notes:

	- Because the list contains a fixed fraction of the total number of
	  elements, but the length must be an integer, there are times when
	  appending one additional element causes the number of elements
	  retained in the list to increase by 1.  When this occurs, the new
	  element is always the one just added, even if it is smaller than
	  elements that have previously been discarded.  To mitigate this
	  effect, a "guard" element is retained, so the actual number of
	  elements retained is equal to the desired fraction of the total
	  plus 1.  However, pathological input sets can always be
	  constructed that defeats the mechanism, for example

	>>> l = Highest(fraction = 0.5)
	>>> l.append(1)
	>>> l.append(1)
	>>> l.append(.125)
	>>> l.append(.125)
	>>> l.append(.125)
	>>> l.append(.125)
	>>> l.append(.125)
	>>> l
	[1, 0.125, 0.125, 0.125]

	  These are not the four highest elements in the sequence, nor are
	  the first three even the three highest elements (discarding the
	  guard element).

	- What is true is that the first N elements are guaranteed to be
	  correct where N is the smallest number of elements that have ever
	  been in the list.  In the example above, after the first append()
	  l has just 1 element in it, and so only ever the first element
	  can be guaranteed to be correct (which is the case).  If you will
	  require to know, without error, the 100 highest values in the
	  sequence, then you must be sure to initialize the list with at
	  least 100/fraction elements from the sequence.
	"""
	def __init__(self, sequence = [], fraction = 1.0):
		list.__init__(self, sequence)
		self.n = list.__len__(self)
		self.fraction = fraction
		list.sort(self, reverse = True)
		del self[int(math.ceil(self.fraction * self.n)) + 1:]

	def __len__(self):
		return self.n

	def append(self, value):
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
		del self[int(math.ceil(self.fraction * self.n)) + 1:]

	def extend(self, sequence):
		before = list.__len__(self)
		list.extend(self, sequence)
		self.n += list.__len__(self) - before
		list.sort(self, reverse = True)
		del self[int(math.ceil(self.fraction * self.n)) + 1:]

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

