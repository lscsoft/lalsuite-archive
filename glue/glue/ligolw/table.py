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
While the ligolw module provides classes and parser support for reading and
writing LIGO Light Weight XML documents, this module supplements that code
with classes and parsers that add intelligence to the in-RAM document
representation.

In particular, the document tree associated with a Table element is
enhanced.  During parsing, the Stream element in this module converts the
character data contained within it into a list of objects.  The list
contains one object for each row of the table, and the objects' attributes
are the names of the table's columns.  When the document is written out
again, the Stream element serializes the row objects back into character
data.

The list of row objects is stored as an attribute of the Table element,
which itself exports a list-like interface to the rows.  The Column
elements also provide list-like access to the values in the corresponding
columns of the table.
"""

__author__ = "Kipp Cannon <kipp@gravity.phys.uwm.edu>"
__date__ = "$Date$"[7:-2]
__version__ = "$Revision$"[11:-2]

import copy
import re
import sys
from xml.sax.xmlreader import AttributesImpl

import ligolw
import tokenizer
import types


#
# =============================================================================
#
#                           Column Name Manipulation
#
# =============================================================================
#

# Regular expression to extract the significant part of a column name
# according to the LIGO LW naming conventions.

# FIXME: the pattern should be
#
# r"(?:\A[a-z0-9_]+:|\A)(?P<FullName>(?:[a-z0-9_]+:|\A)(?P<Name>[a-z0-9_]+))\Z"
#
# but people are putting upper case letters in names!!!!!  Someone is going
# to get the beats.

ColumnPattern = re.compile(r"(?:\A\w+:|\A)(?P<FullName>(?:(?P<Table>\w+):|\A)(?P<Name>\w+))\Z")


def StripColumnName(name):
	"""
	Return the significant portion of a column name according to LIGO
	LW naming conventions.
	"""
	try:
		return ColumnPattern.search(name).group("Name")
	except AttributeError:
		return name


def CompareColumnNames(name1, name2):
	"""
	Convenience function to compare two column names according to LIGO
	LW naming conventions.
	"""
	return cmp(StripColumnName(name1), StripColumnName(name2))


def getColumnsByName(elem, name):
	"""
	Return a list of columns with name name under elem.
	"""
	return elem.getElements(lambda e: (e.tagName == ligolw.Column.tagName) and (CompareColumnNames(e.getAttribute("Name"), name) == 0))


#
# =============================================================================
#
#                           Table Name Manipulation
#
# =============================================================================
#

# Regular expression used to extract the signifcant portion of a table or
# stream name, according to LIGO LW naming conventions.

TablePattern = re.compile(r"(?:\A[a-z0-9_]+:|\A)(?P<Name>[a-z0-9_]+):table\Z")


def StripTableName(name):
	"""
	Return the significant portion of a table name according to LIGO LW
	naming conventions.
	"""
	try:
		return TablePattern.search(name).group("Name")
	except AttributeError:
		return name


def CompareTableNames(name1, name2):
	"""
	Convenience function to compare two table names according to LIGO
	LW naming conventions.
	"""
	return cmp(StripTableName(name1), StripTableName(name2))


def getTablesByName(elem, name):
	"""
	Return a list of tables with name name under elem.
	"""
	return elem.getElements(lambda e: (e.tagName == ligolw.Table.tagName) and (CompareTableNames(e.getAttribute("Name"), name) == 0))


#
# =============================================================================
#
#                                  Utilities
#
# =============================================================================
#

def new_from_template(template):
	"""
	Construct a new Table document subtree whose structure is the same
	as the template table, that is it has the same columns etc..  The
	rows are not copied.  Note that a fair amount of metadata is shared
	between the original and new tables.  In particular, a copy of the
	Table object itself is created, copies of the child nodes are
	created, and a new row list is created.  All other object
	references are shared between the two instances, such as the
	RowType attribute on the Table object.
	"""
	new = copy.copy(template)
	new.childNodes = map(copy.copy, template.childNodes)
	for child in new.childNodes:
		child.parentNode = new
	new.rows = []
	return new


#
# =============================================================================
#
#                               Element Classes
#
# =============================================================================
#

class _ColumnIter(object):
	"""
	An iterator class for looping through the values in a column.
	"""
	def __init__(self, column):
		self.rowiter = iter(column.parentNode.rows)
		self.attr = column.asattribute

	def next(self):
		return getattr(self.rowiter.next(), self.attr)

class Column(ligolw.Column):
	"""
	High-level column element that provides list-like access to the
	values in a column.
	"""
	def __init__(self, attrs):
		ligolw.Column.__init__(self, attrs)
		self.asattribute = StripColumnName(self.getAttribute("Name"))

	def __len__(self):
		"""
		Return the number of values in this column.
		"""
		return len(self.parentNode.rows)

	def __getitem__(self, i):
		"""
		Retrieve the value in this column in row i.
		"""
		return getattr(self.parentNode.rows[i], self.asattribute)

	def __setitem__(self, i, value):
		"""
		Set the value in this column in row i.
		"""
		setattr(self.parentNode.rows[i], self.asattribute, value)

	def __iter__(self):
		"""
		Return an iterator object for iterating over values in this
		column.
		"""
		return _ColumnIter(self)

	def count(self, value):
		"""
		Return the number of rows with this column equal to value.
		"""
		n = 0
		for r in self.parentNode.rows:
			if getattr(r, self.asattribute) == value:
				n += 1
		return n

	def index(self, val):
		"""
		Return the smallest index of the row(s) with this column equal
		to value.
		"""
		for i in range(len(self.parentNode.rows)):
			if getattr(self.parentNode.rows[i], self.asattribute) == value:
				return i
		raise ValueError, "%s not found" % repr(val)

	def asarray(self):
		"""
		Construct a numarray array from this column.
		"""
		if self.getAttribute("Type") in types.StringTypes:
			raise TypeError, "Column does not have numeric type"
		# hack to work around bug in numarray:  numarray tests that
		# an object can be turned into an array, that is it is
		# "list like", by trying to retrieve element 0.  This fails
		# if the list like object has 0 length, causing numarray to
		# barf.  If the object is, in fact, a real Python list then
		# numarray is made happy.
		import numarray
		if not len(self):
			return numarray.array([], type = types.ToNumArrayType[self.getAttribute("Type")], shape = (len(self),))
		return numarray.array(self, type = types.ToNumArrayType[self.getAttribute("Type")], shape = (len(self),))


class TableStream(ligolw.Stream):
	"""
	High-level Stream element for use inside Tables.  This element
	knows how to parse the delimited character stream into rows in the
	parent element, and knows how to turn the parent's rows back into a
	character stream.
	"""
	def __init__(self, attrs):
		ligolw.Stream.__init__(self, attrs)
		self.tokenizer = tokenizer.Tokenizer(self.getAttribute("Delimiter"))
		self.__colnames = None
		self.__numcols = None
		self.__row = None
		self.__colindex = 0

	def appendData(self, content):
		# some initialization that can only be done once parentNode
		# has been set.
		if self.__row == None:
			self.__colnames = tuple(self.parentNode.columnnames)
			self.tokenizer.set_types(self.parentNode.columntypes)
			self.__numcols = len(self.__colnames)
			self.__loadcolumns = self.parentNode.loadcolumns
			self.__row = self.parentNode.RowType()

		# tokenize buffer, and construct row objects
		for token in self.tokenizer.add(content):
			colname = self.__colnames[self.__colindex]
			if self.__loadcolumns == None or colname in self.__loadcolumns:
				try:
					setattr(self.__row, colname, token)
				except AttributeError, e:
					# row object does not have an
					# attribute matching this column.
					# by allowing this, code can be
					# written that saves memory by
					# setting RowType to an object with
					# slots for only the desired
					# columns.
					pass
			self.__colindex += 1
			if self.__colindex >= self.__numcols:
				self.parentNode.append(self.__row)
				self.__row = self.parentNode.RowType()
				self.__colindex = 0

	def unlink(self):
		"""
		Break internal references within the document tree rooted
		on this element to promote garbage collection.
		"""
		self.__colnames = None
		self.__numcols = None
		self.__loadcolumns = None
		self.__row = None
		self.__colindex = 0
		ligolw.Stream.unlink(self)

	def _rowstr(self, row, columninfo):
		# FIXME: after calling getattr(), should probably check that
		# the result has the expected type.
		strs = []
		for isstring, name in columninfo:
			if isstring:
				strs.append("\"" + getattr(row, name) + "\"")
			else:
				strs.append(str(getattr(row, name)))
		return self.getAttribute("Delimiter").join(strs)

	def write(self, file = sys.stdout, indent = ""):
		columninfo = [(c.getAttribute("Type") in types.StringTypes, StripColumnName(c.getAttribute("Name"))) for c in self.parentNode.getElementsByTagName(ligolw.Column.tagName)]

		# loop over parent's rows.  This is complicated because we
		# need to not put a delimiter at the end of the last row.
		print >>file, self.start_tag(indent)
		rowiter = iter(self.parentNode)
		try:
			file.write(indent + ligolw.Indent + self._rowstr(rowiter.next(), columninfo))
			while True:
				file.write(self.getAttribute("Delimiter") + "\n" + indent + ligolw.Indent + self._rowstr(rowiter.next(), columninfo))
		except StopIteration:
			if len(self.parentNode) > 0:
				file.write("\n")
		print >>file, self.end_tag(indent)

	# FIXME: This function is for the metaio library:  metaio cares
	# what order the attributes of XML tags come in.  This function
	# will be removed when the metaio library is fixed.
	def start_tag(self, indent):
		"""
		Generate the element start tag.
		"""
		return indent + "<%s Name=\"%s\" Type=\"%s\" Delimiter=\"%s\">" % (self.tagName, self.getAttribute("Name"), self.getAttribute("Type"), self.getAttribute("Delimiter"))


class TableRow(object):
	"""
	Helpful parent class for row objects.
	"""
	pass


class Table(ligolw.Table):
	"""
	High-level Table element that knows about its columns and rows.
	"""
	validcolumns = None
	loadcolumns = None
	RowType = TableRow

	def __init__(self, *attrs):
		"""
		Initialize
		"""
		ligolw.Table.__init__(self, *attrs)
		self.columnnames = []
		self.columntypes = []
		self.rows = []

	#
	# Sequence methods
	#

	def __len__(self):
		"""
		Return the number of rows in this table.
		"""
		return len(self.rows)

	def __getitem__(self, key):
		"""
		Retrieve row(s).
		"""
		return self.rows[key]

	def __setitem__(self, key, value):
		"""
		Set row(s).
		"""
		self.rows[key] = value

	def __delitem__(self, key):
		"""
		Remove row(s).
		"""
		del self.rows[key]

	def __iter__(self):
		"""
		Return an iterator object for iterating over rows in this
		table.
		"""
		return iter(self.rows)

	def append(self, row):
		"""
		Append row to the list of rows for this table.
		"""
		self.rows.append(row)

	def extend(self, rows):
		"""
		Add a list of rows to the end of the table.
		"""
		self.rows.extend(rows)

	def pop(self, key):
		return self.rows.pop(key)

	def sort(self, *args):
		self.rows.sort(*args)

	def filterRows(self, func):
		"""
		Delete all rows for which func(row) evaluates to False.
		"""
		i = 0
		while i < len(self.rows):
			if not func(self.rows[i]):
				del self.rows[i]
			else:
				i += 1
		return self


	#
	# Column access
	#

	def getColumnByName(self, name):
		try:
			return getColumnsByName(self, name)[0]
		except IndexError:
			raise KeyError, name


	def appendColumn(self, name):
		"""
		Append a column named "name" to the table.  Returns the new
		child.
		"""
		if getColumnsByName(self, name):
			raise ValueError, "duplicate Column \"%s\"" % name
		column = Column(AttributesImpl({u"Name": "%s:%s" % (StripTableName(self.tableName), name), u"Type": self.validcolumns[name]}))
		streams = self.getElementsByTagName(ligolw.Stream.tagName)
		if streams:
			self.insertBefore(column, streams[0])
		else:
			self.appendChild(column)
		return column


	#
	# Element methods
	#

	def _updateColumninfo(self):
		self.columnnames = []
		self.columntypes = []
		for child in self.childNodes:
			if child.tagName != ligolw.Column.tagName:
				continue
			colname = StripColumnName(child.getAttribute("Name"))
			llwtype = child.getAttribute("Type")
			if self.validcolumns != None:
				if colname not in self.validcolumns.keys():
					raise ligolw.ElementError, "invalid Column name \"%s\" for Table \"%s\"" % (child.getAttribute("Name"), self.getAttribute("Name"))
				if self.validcolumns[colname] != llwtype:
					raise ligolw.ElementError, "invalid type \"%s\" for Column \"%s\"" % (llwtype, child.getAttribute("Name"))
			if colname in self.columnnames:
				raise ligolw.ElementError, "duplicate Column \"%s\"" % child.getAttribute("Name")
			self.columnnames.append(colname)
			if llwtype in types.StringTypes:
				self.columntypes.append(str)
			elif llwtype in types.IntTypes:
				self.columntypes.append(int)
			elif llwtype in types.FloatTypes:
				self.columntypes.append(float)
			else:
				raise ligolw.ElementError, "unrecognized Type attribute \"%s\" for Column \"%s\" in Table \"%s\"" % (llwtype, child.getAttribute("Name"), self.getAttribute("Name"))

	def _verifyChildren(self, i):
		ligolw.Table._verifyChildren(self, i)
		child = self.childNodes[i]
		if child.tagName == ligolw.Column.tagName:
			self._updateColumninfo()
		elif child.tagName == ligolw.Stream.tagName:
			if child.getAttribute("Name") != self.getAttribute("Name"):
				raise ligolw.ElementError, "Stream name \"%s\" does not match Table name \"%s\"" % (child.getAttribute("Name"), self.getAttribute("Name"))

	def removeChild(self, child):
		"""
		Remove a child from this element.  The child element is
		returned, and it's parentNode element is reset.
		"""
		ligolw.Table.removeChild(self, child)
		if child.tagName == ligolw.Column.tagName:
			self._updateColumninfo()
		return child

	def unlink(self):
		"""
		Break internal references within the document tree rooted
		on this element to promote garbage collection.
		"""
		ligolw.Table.unlink(self)
		self.rows = []


#
# =============================================================================
#
#                               Content Handler
#
# =============================================================================
#

#
# Override portions of the ligolw.LIGOLWContentHandler class
#

__parent_startStream = ligolw.LIGOLWContentHandler.startStream
__parent_endStream = ligolw.LIGOLWContentHandler.endStream

def startColumn(self, attrs):
	return Column(attrs)

def startStream(self, attrs):
	if self.current.tagName == ligolw.Table.tagName:
		return TableStream(attrs)
	return __parent_startStream(self, attrs)

def endStream(self):
	# stream tokenizer uses delimiter to identify end of each token, so
	# add a final delimiter to induce the last token to get parsed.
	if self.current.parentNode.tagName == ligolw.Table.tagName:
		self.current.appendData(self.current.getAttribute("Delimiter"))
	else:
		__parent_endStream(self)

def startTable(self, attrs):
	return Table(attrs)

ligolw.LIGOLWContentHandler.startColumn = startColumn
ligolw.LIGOLWContentHandler.startStream = startStream
ligolw.LIGOLWContentHandler.endStream = endStream
ligolw.LIGOLWContentHandler.startTable = startTable


#
# This class defined for backwards compatibility;  remove when nobody uses
# it.
#

import warnings
class LIGOLWContentHandler(ligolw.LIGOLWContentHandler):
	def __init__(*args):
		warnings.warn("table.LIGOLWContentHandler() class is deprecated:  use ligolw.LIGOLWContentHandler() instead", DeprecationWarning)
		ligolw.LIGOLWContentHandler.__init__(*args)
