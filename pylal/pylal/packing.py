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

import os

from glue.lal import CacheEntry
from glue import segments


#
# =============================================================================
#
#                                     Bins
#
# =============================================================================
#

class Bin(object):
	"""
	Bin object for use in packing algorithm implementations.
	"""
	def __init__(self):
		"""
		Initialize a new Bin instance.
		"""
		self.size = 0
		self.objects = []

	def add(self, object, size):
		"""
		Add the object, whose size is as given, to the bin.
		"""
		self.objects.append(object)
		self.size += size
		return self

	def __iadd__(self, other):
		"""
		Add the contents of another Bin object to this one.
		"""
		self.objects.extend(other.objects)
		self.size += other.size
		return self

	def __cmp__(self, other):
		"""
		Compare two Bin objects by their sizes.
		"""
		return cmp(self.size, other.size)

	def __repr__(self):
		"""
		A representation of the Bin object.
		"""
		return "(%s, %s)" % (str(self.size), str(self.objects))

	__str__ = __repr__


#
# =============================================================================
#
#                              Packing Algorithms
#
# =============================================================================
#

def bin_list(n, bintype = Bin):
	"""
	Convenience function for constructing a list of bins.
	"""
	# [bintype()] * n doesn't work because it produces a list of
	# references to the same object, rather than a list of references
	# to distinct objects.
	l = []
	for i in xrange(n):
		l.append(bintype())
	return l


class Packer(object):
	"""
	Generic parent class for packing algorithms.
	"""
	def __init__(self, bins):
		"""
		Set the list of bins on which we shall operate
		"""
		self.bins = bins

	def pack(self, size, object):
		"""
		Pack an object of given size into the bins.
		"""
		raise NotImplementedError

	def packlist(self, size_object_pairs):
		"""
		Pack the list of (size, object) tuples into the bins.
		"""
		raise NotImplementedError


class BiggestIntoEmptiest(Packer):
	"""
	Packs the biggest object into the emptiest bin.
	"""
	def pack(self, size, object):
		min(self.bins).add(object, size)

	def packlist(self, size_object_pairs):
		size_object_pairs.sort()
		size_object_pairs.reverse()
		for size, object in size_object_pairs:
			self.pack(size, object)
