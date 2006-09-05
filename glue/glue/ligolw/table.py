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
	del new[:]
	return new


def get_table(xmldoc, name):
	"""
	Scan xmldoc for a table named name.  Raises ValueError if not
	exactly 1 such table is found.
	"""
	tables = getTablesByName(xmldoc, name)
	if len(tables) != 1:
		raise ValueError, "document must contain exactly one %s table" % StripTableName(name)
	return tables[0]


#
# =============================================================================
#
#                               Element Classes
#
# =============================================================================
#

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
		return len(self.parentNode)

	def __getitem__(self, i):
		"""
		Retrieve the value in this column in row i.
		"""
		if type(i) == slice:
			return map(lambda r: getattr(r, self.asattribute), self.parentNode[i])
		else:
			return getattr(self.parentNode[i], self.asattribute)

	def __setitem__(self, i, value):
		"""
		Set the value in this column in row i.
		"""
		if type(i) == slice:
			map(lambda r: setattr(r, self.asattribute, value), self.parentNode[i])
		else:
			setattr(self.parentNode[i], self.asattribute, value)

	def __iter__(self):
		"""
		Return an iterator object for iterating over values in this
		column.
		"""
		for row in self.parentNode:
			yield getattr(row, self.asattribute)

	def count(self, value):
		"""
		Return the number of rows with this column equal to value.
		"""
		n = 0
		for row in self.parentNode:
			if getattr(row, self.asattribute) == value:
				n += 1
		return n

	def index(self, val):
		"""
		Return the smallest index of the row(s) with this column equal
		to value.
		"""
		for i in xrange(len(self.parentNode)):
			if getattr(self.parentNode[i], self.asattribute) == value:
				return i
		raise ValueError, "%s not found" % repr(val)

	def asarray(self):
		"""
		Construct a numarray array from this column.  Note that
		this creates a copy of the data, so modifications made to
		the array will not be recorded in the original document.
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
			self.tokenizer.set_types(self.parentNode.columnpytypes)
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
					import warnings
					warnings.warn("glue.ligolw.table.TableStream.appendData():  invalid attribute %s for %s;  in the future this will be a fatal error, use Table class' loadcolumns attribute to restrict parsed columns" % (colname, repr(self.__row)), DeprecationWarning)
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

	def write(self, file = sys.stdout, indent = ""):
		rowfmt = indent + ligolw.Indent + self.getAttribute("Delimiter").join([types.ToFormat[c.getAttribute("Type")] for c in self.parentNode.getElementsByTagName(ligolw.Column.tagName)])
		colnames = self.parentNode.columnnames

		# loop over parent's rows.  This is complicated because we
		# need to not put a delimiter at the end of the last row.
		print >>file, self.start_tag(indent)
		rowiter = iter(self.parentNode)
		try:
			row = rowiter.next()
			file.write(rowfmt % tuple([getattr(row, name) for name in colnames]))
			rowfmt = self.getAttribute("Delimiter") + "\n" + rowfmt
			while True:
				row = rowiter.next()
				file.write(rowfmt % tuple([getattr(row, name) for name in colnames]))
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
	Helpful parent class for row objects.  Also used as the default row
	class by Table instances.
	"""
	pass


class Table(ligolw.Table, list):
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
		self.columnpytypes = []

	#
	# Sequence methods
	#

	def filterRows(self, func):
		"""
		Delete all rows for which func(row) evaluates to False.
		"""
		for i in xrange(len(self) - 1, -1, -1):
			if not func(self[i]):
				del self[i]
		return self


	#
	# Column access
	#

	def getColumnByName(self, name):
		"""
		Retrieve and return the Column child element whose name is
		as given.  Raises KeyError if this table has no column by
		that name.
		"""
		try:
			return getColumnsByName(self, name)[0]
		except IndexError:
			raise KeyError, name


	def appendColumn(self, name):
		"""
		Append a column named "name" to the table.  Returns the new
		child.  Raises ValueError if the table already has a column
		by that name, and KeyError if the validcolumns attribute of
		this table does not contain an entry for a column by that
		name.
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

	def _update_column_info(self):
		"""
		Used for validation during parsing, and additional
		book-keeping.  For internal use only.
		"""
		self.columnnames = []
		self.columntypes = []
		self.columnpytypes = []
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
			self.columntypes.append(llwtype)
			if llwtype in types.StringTypes:
				self.columnpytypes.append(str)
			elif llwtype in types.IntTypes:
				self.columnpytypes.append(int)
			elif llwtype in types.FloatTypes:
				self.columnpytypes.append(float)
			else:
				raise ligolw.ElementError, "unrecognized Type attribute \"%s\" for Column \"%s\" in Table \"%s\"" % (llwtype, child.getAttribute("Name"), self.getAttribute("Name"))

	def _verifyChildren(self, i):
		"""
		Used for validation during parsing, and additional
		book-keeping.  For internal use only.
		"""
		ligolw.Table._verifyChildren(self, i)
		child = self.childNodes[i]
		if child.tagName == ligolw.Column.tagName:
			self._update_column_info()
		elif child.tagName == ligolw.Stream.tagName:
			if child.getAttribute("Name") != self.getAttribute("Name"):
				raise ligolw.ElementError, "Stream name \"%s\" does not match Table name \"%s\"" % (child.getAttribute("Name"), self.getAttribute("Name"))

	def _end_of_columns(self):
		"""
		Called during parsing to indicate that the last Column
		child element has been added.
		"""
		pass

	def removeChild(self, child):
		"""
		Remove a child from this element.  The child element is
		returned, and it's parentNode element is reset.
		"""
		ligolw.Table.removeChild(self, child)
		if child.tagName == ligolw.Column.tagName:
			self._update_column_info()
		return child

	def unlink(self):
		"""
		Break internal references within the document tree rooted
		on this element to promote garbage collection.
		"""
		ligolw.Table.unlink(self)
		del self[:]


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
		self.current._end_of_columns()
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

def endTable(self):
	# Table elements are allowed to contain 0 Stream children, but
	# _end_of_columns() needs to be called regardless.
	if self.current.childNodes[-1].tagName != ligolw.Stream.tagName:
		self.current._end_of_columns()

ligolw.LIGOLWContentHandler.startColumn = startColumn
ligolw.LIGOLWContentHandler.startStream = startStream
ligolw.LIGOLWContentHandler.endStream = endStream
ligolw.LIGOLWContentHandler.startTable = startTable
ligolw.LIGOLWContentHandler.endTable = endTable


#
# This class defined for backwards compatibility;  remove when nobody uses
# it.
#

import warnings
class LIGOLWContentHandler(ligolw.LIGOLWContentHandler):
	def __init__(*args):
		warnings.warn("table.LIGOLWContentHandler() class is deprecated:  use ligolw.LIGOLWContentHandler() instead", DeprecationWarning)
		ligolw.LIGOLWContentHandler.__init__(*args)
