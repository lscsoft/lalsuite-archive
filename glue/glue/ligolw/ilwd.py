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
Light-Weight XML files.  This module provides the ilwdchar class used as a
parent class for memory-efficient storage of ilwd:char strings.
"""


import re


__author__ = "Kipp Cannon <kipp@gravity.phys.uwm.edu>"
__date__ = "$Date$"[7:-2]
__version__ = "$Revision$"[11:-2]


#
# =============================================================================
#
#                            ilwd:char Manipulation
#
# =============================================================================
#


#
# Regular expression to extract the parts of an ilwd:char string according
# to the LIGO LW naming conventions.
#


ilwdchar_pattern = re.compile(r"(?P<Table>\w+):(?P<Column>\w+):(?P<ID>\d+)")


#
# =============================================================================
#
#                               ID Parent Class
#
# =============================================================================
#


#
# Load the C extension module that provides the ilwdchar parent class.
#


from __ilwd import *


#
# =============================================================================
#
#                                Cached Classes
#
# =============================================================================
#


#
# Cache of pre-defined ilwd:char subclasses, indexed by the string prefix
# used with each.
#


ilwdchar_class_cache = {}


#
# Functions for retrieving ilwdchar subclasses, and constructing subclass
# instances.
#


def get_ilwdchar_class(tbl_name, col_name):
	"""
	Searches the cache of pre-defined ilwdchar subclasses for a class
	whose table_name and column_name attributes match those provided.
	If a matching subclass is found it is returned.  Otherwise a new
	class is defined, added to the cache, and returned.
	"""
	#
	# if the class already exists, retrieve it
	#

	key = (str(tbl_name), str(col_name))
	try:
		return ilwdchar_class_cache[key]
	except KeyError:
		#
		# define a new class, and add it to the cache
		#

		class cached_ilwdchar_class(ilwdchar):
			__slots__ = ()
			table_name, column_name = key
			index_offset = len("%s:%s:" % key)

		ilwdchar_class_cache[key] = cached_ilwdchar_class

		return cached_ilwdchar_class


def get_ilwdchar(s):
	"""
	Convert an ilwd:char string into an instance of the matching
	subclass of ilwdchar.
	"""
	#
	# try parsing the string as an ilwd:char formated string
	#

	m = ilwdchar_pattern.match(s)
	if m is None:
		raise ValueError, repr(s)

	#
	# retrieve the matching class from the ID class cache, and return
	# an instance initialized to the desired value
	#

	table_name, column_name, i = m.groups()

	return get_ilwdchar_class(table_name, column_name)(int(i))
