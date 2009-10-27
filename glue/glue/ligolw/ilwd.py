# Copyright (C) 2006  Kipp Cannon
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
The ilwd:char type is used to store ID strings for objects within LIGO
Light-Weight XML files.  This module provides the ilwdchar class used as a
parent class for memory-efficient storage of ilwd:char strings.

LIGO Light Weight XML "ilwd:char" IDs are strings of the form
"table:column:integer", for example "process:process_id:10".  Large complex
documents can have many millions of these strings, and their storage
represents a significant RAM burden.  However, while there can be millions
of ID strings in use there may be only a small number (e.g. 10 or fewer)
unique ID prefixes in use (the table name and column name part).  The
amount of RAM required to load a document can be significantly reduced if
the small number of unique string prefixes is stored separately.  This
module (and the __ilwd C extension module it uses) implement the machinery
used to do this.

To get started, the get_ilwdchar() function provided by this module is the
means by which an ID string is converted into an instance of the
appropriate subclass of ilwdchar.

Example:

>>> x = get_ilwdchar("process:process_id:10")

Like strings, the object resulting from this function call is immutable.
It provides two read-only attributes, "table_name" and "column_name", that
can be used to access the table and column parts of the original ID string.
The integer suffix can be retrieved by converting the object to an integer.

Example:

>>> x.table_name
'process'
>>> int(x)
10

The object also provides the read-only attribute "index_offset", giving the
length of the string preceding the interger suffix.

Example:

>>> x.index_offset
19

The objects support some arithmetic operations.

Example:

>>> y = x + 5
>>> str(y)
'process:process_id:15'
>>> int(y - x)
5
"""


import re


from glue import git_version


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date


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
# Special version for inspirals
#
# FIXME:  remove this when sngl_inspiral tables no longer contain 20
# bazillion digit long ids
#


class inspiral_ilwdchar(long):
	__slots__ = ()

	table_name = "sngl_inspiral"
	column_name = "event_id"
	index_offset = len("%s:%s:" % (table_name, column_name))

	def __hash__(self):
		return hash(self.table_name) ^ hash(self.column_name) ^ long.__hash__(self)

	def __add__(self, other):
		return self.__class__(long.__add__(self, other))

	def __cmp__(self, other):
		return cmp(self.table_name, other.table_name) or cmp(self.column_name, other.column_name) or long.__cmp__(self, other)

	def __str__(self):
		return "%s:%s:%d" % (self.table_name, self.column_name, long(self))

	def __conform__(self, protocol):
		return unicode(self)

	def __sub__(self, other):
		return long.__sub__(self, other)


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


# FIXME:  remove the special initializer when the inspiral ids no longer
# contain a bazillion digits.
ilwdchar_class_cache = {
	("sngl_inspiral", "event_id"): inspiral_ilwdchar
}


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

	>>> process_id = get_ilwdchar_class("process", "process_id")
	>>> x = process_id(10)
	>>> x
	<glue.ligolw.ilwd.cached_ilwdchar_class object at 0x2b8de0a186a8>
	>>> str(x)
	'process:process_id:10'

	Retrieving and storing the class provides a convenient mechanism
	for quickly constructing new ID objects.

	Example:

	>>> for i in range(10):
	...	print str(process_id(i))
	...
	process:process_id:0
	process:process_id:1
	process:process_id:2
	process:process_id:3
	process:process_id:4
	process:process_id:5
	process:process_id:6
	process:process_id:7
	process:process_id:8
	process:process_id:9
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

			def __conform__(self, protocol):
				# The presence of this method allows
				# ilwdchar sub-classes to be inserted
				# directly into SQLite databases as
				# strings. See
				#
				# http://www.python.org/dev/peps/pep-0246
				#
				# for more information.
				#
				# NOTE:  GvR has rejected that PEP, so this
				# mechanism is obsolete.  Be prepared to
				# fix this, replacing it with whatever
				# replaces it.
				#
				# NOTE:  The return should be inside an "if
				# protocol is sqlite3.PrepareProtocol:"
				# conditional, but that would require
				# importing sqlite3 which would break this
				# module on FC4 boxes, and I'm not going to
				# spend time fixing something that's
				# obsolete anyway.
				return unicode(self)

		ilwdchar_class_cache[key] = cached_ilwdchar_class

		return cached_ilwdchar_class


def get_ilwdchar(s):
	"""
	Convert an ilwd:char string into an instance of the matching
	subclass of ilwdchar.

	Example:

	>>> x = get_ilwdchar("process:process_id:10")
	>>> str(x)
	'process:process_id:10'
	>>> x.table_name
	'process'
	>>> x.column_name
	'process_id'
	>>> int(x)
	10
	>>> x.index_offset
	19
	>>> str(x)[x.index_offset:]
	'10'
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
