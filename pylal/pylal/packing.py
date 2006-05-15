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
This module provides bin packing utilities.
"""

__author__ = "Kipp Cannon <kipp@gravity.phys.uwm.edu>"
__version__ = "$Revision$"[11:-2]
__date__ = "$Date$"[7:-2]


#
# =============================================================================
#
#                                 Bin Packing
#
# =============================================================================
#

class Bin(object):
	"""
	Bin object for use in packing algorithm implementations.
	"""
	def __init__(self):
		self.size = 0
		self.objects = []

	def add_object(self, object, size):
		self.size += size
		self.objects.append(object)

	def get_objects(self):
		return self.objects

	def __cmp__(self, other):
		return cmp(self.size, other.size)

	def __repr__(self):
		return "(%d, %s)" % (self.size, self.objects)

	__str__ = __repr__


class Bins(list):
	"""
	List of Bins for use in packing algorithm implementations.
	"""
	def __new__(cls, n):
		l = list.__new__(cls, None)
		for i in xrange(n):
			l.append(Bin())
		return l

	def __init__(cls, n):
		pass


#
# =============================================================================
#
#                              Packing Algorithms
#
# =============================================================================
#

class PakingAlgorithm(object):
	"""
	Generic parent class for packing algorithms.
	"""
	def __init__(self, bins):
		"""
		Set the list of bins on which we shall operate
		"""
		self.bins = bins

	def pack(self, object, size):
		"""
		Pack an object of given size into the bins.
		"""
		raise NotImplemented

	def packlist(self, size_object_pairs):
		"""
		Pack the list of (size, object) tuples into the bins.
		"""
		raise NotImplemented


class BiggestIntoEmptiest(PakingAlgorithm):
	"""
	Packs the biggest object into the emptiest bin.
	"""
	def pack(self, object, size):
		min(self.bins).add_object(object, size)

	def packlist(self, size_object_pairs):
		l.sort()
		l.reverse()
		map(lambda (s, o): self.pack(o, s), l)
