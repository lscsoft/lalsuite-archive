# $Id$
#
# Copyright (C) 2006  Kipp C. Cannon
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
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
#                                    ILWDs
#
# =============================================================================
#


"""
The ilwd:char and ilwd:char_u types are used as IDs for objects within LIGO
Light-Weight XML files.  This module provides the ILWD class used to
represent a sequence of IDs.  An instance of the class can be used to
generate unique IDs.
"""


import re

import table


__author__ = "Kipp Cannon <kipp@gravity.phys.uwm.edu>"
__date__ = "$Date$"[7:-2]
__version__ = "$Revision$"[11:-2]


#
# =============================================================================
#
#                              ILWD Manipulation
#
# =============================================================================
#


#
# Regular expression to extract the parts of a row ID according to the LIGO
# LW naming conventions.
#


ILWDPattern = re.compile(r"(?P<Table>\w+):(?P<Column>\w+):(?P<ID>\d+)")


#
# Utility functions to extra ID parts.
#


def ILWDTableName(ilwdchar):
	"""
	Return the table name part of the row ID ilwdchar.  ValueError is
	raised if the ID cannot be parsed.
	"""
	try:
		return ILWDPattern.search(ilwdchar).group("Table")
	except AttributeError:
		raise ValueError, "unrecognized ID '%s'" % repr(ilwdchar)


def ILWDColumnName(ilwdchar):
	"""
	Return the column name part of the row ID ilwdchar.  ValueError is
	raised if the ID cannot be parsed.
	"""
	try:
		return ILWDPattern.search(ilwdchar).group("Column")
	except AttributeError:
		raise ValueError, "unrecognized ID '%s'" % repr(ilwdchar)


def ILWDID(ilwdchar):
	"""
	Return the ID part of the row ID ilwdchar.  ValueError is raised if
	the ID cannot be parsed.
	"""
	try:
		return int(ILWDPattern.search(ilwdchar).group("ID"))
	except AttributeError:
		raise ValueError, "unrecognized ID '%s'" % repr(ilwdchar)


#
# =============================================================================
#
#                               ID Parent Class
#
# =============================================================================
#


class ILWD(int):
	"""
	Row ID parent class.  This is only useful when subclassed in order
	to provide specific values of the class variables "table_name",
	"column_name", and "index_offset".
	"""
	__slots__ = ()
	table_name = None
	column_name = None
	index_offset = None

	# Wow, all these over rides suck

	def __add__(self, other):
		return self.__class__(int.__add__(self, other))

	def __and__(self, other):
		return self.__class__(int.__and__(self, other))

	def __cmp__(self, other):
		try:
			return cmp((self.table_name, self.column_name, int(self)), (other.table_name, other.column_name, int(other)))
		except:
			return NotImplemented

	def __coerce__(self, other):
		try:
			return self, self.__class__(other)
		except:
			return None

	def __div__(self, other):
		return self.__class__(int.__div__(self, other))

	def __divmod__(self, other):
		return self.__class__(int.__divmod__(self, other))

	def __floordiv__(self, other):
		return self.__class__(int.__floordiv__(self, other))

	def __invert__(self, other):
		return self.__class__(int.__invert__(self, other))

	def __lshift__(self, other):
		return self.__class__(int.__lshift__(self, other))

	def __mod__(self, other):
		return self.__class__(int.__mod__(self, other))

	def __mul__(self, other):
		return self.__class__(int.__mul__(self, other))

	def __neg__(self):
		return self.__class__(int.__neg__(self))

	def __or__(self, other):
		return self.__class__(int.__or__(self, other))

	def __pos__(self, other):
		return self.__class__(int.__pos__(self))

	def __pow__(self, *args):
		return self.__class__(int.__pow__(self, *args))

	def __radd__(self, other):
		return self.__class__(int.__radd__(self, other))

	def __rand__(self, other):
		return self.__class__(int.__rand__(self, other))

	def __rdiv__(self, other):
		return self.__class__(int.__rdiv__(self, other))

	def __rdivmod__(self, other):
		return self.__class__(int.__rdivmod__(self, other))

	def __rfloordiv__(self, other):
		return self.__class__(int.__rfloordiv__(self, other))

	def __rlshift__(self, other):
		return self.__class__(int.__rlshift__(self, other))

	def __rmod__(self, other):
		return self.__class__(int.__rmod__(self, other))

	def __rmul__(self, other):
		return self.__class__(int.__rmul__(self, other))

	def __ror__(self, other):
		return self.__class__(int.__ror__(self, other))

	def __rpow__(self, *args):
		return self.__class__(int.__rpow__(self, *args))

	def __rrshift__(self, other):
		return self.__class__(int.__rrshift__(self, other))

	def __rshift__(self, other):
		return self.__class__(int.__rshift__(self, other))

	def __rsub__(self, other):
		return self.__class__(int.__rsub__(self, other))

	def __rtruediv__(self, other):
		return self.__class__(int.__rtruediv__(self, other))

	def __rxor__(self, other):
		return self.__class__(int.__rxor__(self, other))

	def __str__(self):
		return "%s:%s:%d" % (self.table_name, self.column_name, self)

	def __sub__(self, other):
		return self.__class__(int.__sub__(self, other))

	def __truediv__(self, other):
		return self.__class__(int.__truediv__(self, other))

	def __xor__(self, other):
		return self.__class__(int.__xor__(self, other))


#
# =============================================================================
#
#                                Cached Classes
#
# =============================================================================
#


IDClassCache = {}


def get_id_class(tbl_name, col_name):
	tbl_name = str(tbl_name)
	col_name = str(col_name)

	#
	# if the class already exists, retrieve it
	#

	key = (tbl_name, col_name)
	if key in IDClassCache:
		return IDClassCache[key]

	#
	# define a new class
	#

	class cached_id_class(ILWD):
		__slots__ = ()
		table_name = tbl_name
		column_name = col_name
		index_offset = len("%s:%s:" % key)

	#
	# cache the new class and return it
	#

	IDClassCache[key] = cached_id_class

	return cached_id_class
