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


def MultiIter(*lists):
	"""
	An generator for iterating over the elements of multiple lists
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
		yield lists


def choices(vals, n):
	"""
	Iterate over all choices of n elements from the list vals.  In each
	result returned, the original order of the values is preserved.

	Example:

	>>> list(choices(["a", "b", "c"], 2))
	[['a', 'b'], ['a', 'c'], ['b', 'c']]
	"""
	if n == len(vals):
		yield vals
	elif n > 1:
		for i in xrange(len(vals) - n + 1):
			v = [vals[i]]
			for c in choices(vals[i+1:], n - 1):
				yield v + c
	elif n == 1:
		for v in vals:
			yield [v]
	else:
		# n < 1
		raise ValueError, n

