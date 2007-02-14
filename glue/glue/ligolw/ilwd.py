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

# Regular expression to extract the parts of a row ID according to the LIGO
# LW naming conventions.

ILWDPattern = re.compile(r"(?P<Table>\w+):(?P<Column>\w+):(?P<ID>\d+)")


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


class ILWD(object):
	"""
	Unique ILWD generator.
	"""
	def __init__(self, table_name, column_name, n = 0):
		"""
		Initialize an ILWD object.  table_name and column_name are
		the names of the table and column within the table for
		which these will be used as IDs, eg., "process" and
		"process_id".  The optional n parameter sets the starting
		value for the numeric suffix in the ilwd:char
		representation of ID (default is 0).  An initialized ILWD
		generator instance has the attribute index_offset
		containing the location of the numeric index portion of the
		ID within the string.
		"""
		self.table_name = table.StripTableName(table_name)
		self.column_name = table.StripColumnName(column_name)
		self.n = n
		self.index_offset = len(self.table_name) + len(self.column_name) + 2

	def __str__(self):
		"""
		Return an ilwd:char string representation of this ID.
		"""
		return "%s:%s:%d" % (self.table_name, self.column_name, self.n)

	def __cmp__(self, other):
		"""
		Compare IDs first by the base string, then by n.
		"""
		if isinstance(other, ILWD):
			return cmp((self.table_name, self.column_name, self.n), (other.table_name, other.column_name, other.n))
		return NotImplemented

	def __getitem__(self, n):
		return "%s:%s:%d" % (self.table_name, self.column_name, n)

	def __iter__(self):
		return self

	def set_next(self, n):
		"""
		Set the current value of the numeric suffix.
		"""
		self.n = n

	def next(self):
		"""
		Return the current value of the generator as a string, and
		increment the numeric suffix.
		"""
		s = str(self)
		self.n += 1
		return s

