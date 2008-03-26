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

LIGO Light Weight XML "ilwd:char" IDs are strings of the form
"table:column:integer", for example "process:process_id:10".  Large complex
documents can have many millions of these strings, and their storage
represents a significant RAM burden.  At the same time, however, while
there can be millions of ID strings in use there may be only a small number
(e.g. 10 or fewer) unique ID prefixes in use (the table name and column
name part).  The amount of RAM required to load a document can be
significantly reduced if the small number of unique string prefixes is
stored separately.  This module (and the __ilwd C extension module it uses)
implement the machinery used to do this.

The technique makes use of Python's ability to define new classes at
runtime.  Each unique string prefix, for example "process:process_id", is
mapped to a class.  The class definitions are stored in a look-up table
indexed by the (table name, column name) string pairs.  The table and
column strings are stored as class attributes (shared by all instances of
the class), while only the integer suffix unique to each ID is stored as an
instance attribute.  For those who have no idea what this means, an example
to illustrate:

	class Example(object):
		# a class attribute.  the value of this variable is shared
		# by all instances of the class, only one copy is stored in
		# memory
		prefix = "hello"

		def __init__(self):
			# an instance attribute.  each instance of this
			# class gets its own variable by this name, each
			# can have a different value
			suffix = 10

In detail, the implementation in this module begins with the ilwdchar
class, coded in C for speed, that acts as the parent class for all
prefix-specific ID classes.  The parent class implements all the methods
and features required by an ID object.  The dictionary ilwdchar_class_cache
is used to store subclasses of ilwd.  The helper function get_ilwdchar() is
used to convert an ID string of the form "table:column:integer" into an
instance of the appropriate subclass of ilwdchar.  This is done by parsing
the string and retrieving the class matching the table and column name
parts from the class dictionary, or creating a new class on the fly if a
matching class is not found.  Either way, a new instance of the class is
created initialized to the integer part of the ID.
"""


import re


__author__ = "Kipp Cannon <kipp@gravity.phys.uwm.edu>"
__date__ = "$Date$"[7:-2]
__version__ = "$Revision$"[11:-2]


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
	If a matching subclass is found it is returned;  otherwise a new
	class is defined, added to the cache, and returned.

	Example:

	>>> cls = get_ilwdchar_class("process", "process_id")
	>>> id = cls(10)
	>>> id
	<glue.ligolw.ilwd.cached_ilwdchar_class object at 0x2b8de0a186a8>
	>>> str(id)
	'process:process_id:10'
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

	Example:

	>>> id = get_ilwdchar("process:process_id:10")
	>>> str(id)
	'process:process_id:10'
	>>> id.table_name
	'process'
	>>> id.column_name
	'process_id'
	>>> int(id)
	10
	>>> id.index_offset
	19
	"""
	#
	# try parsing the string as an ilwd:char formated string
	#

	try:
		table_name, column_name, i = s.split(":")
	except ValueError, AttributeError:
		raise ValueError, "invalid ilwd:char %s" % repr(s)

	#
	# retrieve the matching class from the ID class cache, and return
	# an instance initialized to the desired value
	#

	return get_ilwdchar_class(table_name, column_name)(int(i))
