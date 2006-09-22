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
A collection of iteration utilities.
"""

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

class MultiIter(object):
	"""
	An iterator class for iterating over the elements of multiple lists
	simultaneously.  An instance of the class is initialized with a
	list of lists.  A call to next() returns a list of elements, one
	from each of the lists.  Subsequent calls to next() iterate over
	all combinations of elements from the lists.

	Example:

	>>> x = MultiIter([[0, 1, 2], [10, 11, 12]])
	>>> list(x)
	[[0, 10], [1, 10], [2, 10], [0, 11], [1, 11], [2, 11], [0, 12], [1,
	12], [2, 12]]
	"""
	def __init__(self, lists):
		self.lists = tuple(lists)
		self.index = [0] * len(lists)
		self.length = tuple(map(len, lists))
		self.stop = 0 in self.length

	def __len__(self):
		return reduce(int.__mul__, self.length)

	def __iter__(self):
		return self

	def next(self):
		if self.stop:
			raise StopIteration
		l = map(lambda l, i: l[i], self.lists, self.index)
		for i in xrange(len(self.index)):
			self.index[i] += 1
			if self.index[i] < self.length[i]:
				break
			self.index[i] = 0
		else:
			self.stop = True
		return l


def choices(vals, n):
	"""
	Return a list of all choices of n elements from the list vals.
	Order is preserved.

	Example:

	>>> choices(["a", "b", "c"], 2)
	[['a', 'b'], ['a', 'c'], ['b', 'c']]
	"""
	if n == len(vals):
		return [vals]
	if n == 1:
		return [[v] for v in vals]
	if n < 1:
		raise ValueError, n
	l = []
	for i in xrange(len(vals) - n + 1):
		for c in choices(vals[i+1:], n - 1):
			c.insert(0, vals[i])
			l.append(c)
	return l

