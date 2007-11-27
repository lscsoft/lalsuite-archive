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
# Utility functions to extract ID parts.
#


def ILWDTableName(ilwdchar):
	"""
	Return the table name part of the row ID ilwdchar.  ValueError is
	raised if the ID cannot be parsed.
	"""
	try:
		return ILWDPattern.match(ilwdchar).group("Table")
	except AttributeError:
		raise ValueError, "unrecognized ID '%s'" % repr(ilwdchar)


def ILWDColumnName(ilwdchar):
	"""
	Return the column name part of the row ID ilwdchar.  ValueError is
	raised if the ID cannot be parsed.
	"""
	try:
		return ILWDPattern.match(ilwdchar).group("Column")
	except AttributeError:
		raise ValueError, "unrecognized ID '%s'" % repr(ilwdchar)


def ILWDID(ilwdchar):
	"""
	Return the ID part of the row ID ilwdchar.  ValueError is raised if
	the ID cannot be parsed.
	"""
	try:
		return int(ILWDPattern.match(ilwdchar).group("ID"))
	except AttributeError:
		raise ValueError, "unrecognized ID '%s'" % repr(ilwdchar)


#
# =============================================================================
#
#                               ID Parent Class
#
# =============================================================================
#


from __ilwd import *


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


def get_id(string):
	"""
	Convert an ilwd:char string into an instance of the matching
	subclass of ILWD.
	"""
	# try parsing the string as an ilwd:char formated string
	m = ILWDPattern.match(string)
	if m is None:
		# nope, didn't work
		raise ValueError, string
	# retrieve the matching class from the ID class cache, and return
	# an instance initialized to the desired value
	tbl, col, i = m.groups()
	return get_id_class(tbl, col)(int(i))
