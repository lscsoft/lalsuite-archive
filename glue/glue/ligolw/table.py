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

The Table element exports a list-like interface to the rows.  The Column
elements also provide list-like access to the values in the corresponding
columns of the table.
"""


import copy
import re
import sys
from xml.sax.saxutils import escape as xmlescape
from xml.sax.xmlreader import AttributesImpl


from glue import git_version
from glue.ligolw import ligolw
from glue.ligolw import tokenizer
from glue.ligolw import types as ligolwtypes


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date


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
	Table object itself is created (but with no rows), and copies of
	the child nodes are created.  All other object references are
	shared between the two instances, such as the RowType attribute on
	the Table object.
	"""
	new = copy.copy(template)
	new.childNodes = map(copy.copy, template.childNodes)
	for child in new.childNodes:
		child.parentNode = new
	del new[:]
	new._end_of_columns()
	new._end_of_rows()
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


def reassign_ids(elem):
	"""
	Recurses over all tables below elem whose next_id attributes are
	not none, and uses the next_id attribute to assign new IDs to the
	rows in each table.  The modifications are recorded in a mapping of
	old row keys to new row keys.  Finally, applies the mapping to all
	rows of all tables to update cross references.
	"""
	mapping = {}
	for tbl in elem.getElementsByTagName(ligolw.Table.tagName):
		if tbl.next_id is not None:
			tbl.updateKeyMapping(mapping)
	for tbl in elem.getElementsByTagName(ligolw.Table.tagName):
		tbl.applyKeyMapping(mapping)


def reset_next_ids(classes):
	"""
	For each class in the list, if the .next_id attribute is not None
	(meaning the table has an ID generator associated with it), set
	.next_id to 0.  This has the effect of reset the ID generators, and
	is useful in applications that process multiple documents and wish
	to reset all the ID generators between documents so that the
	assigned IDs don't grow without bound as each document is
	processed.

	Example:

	>>> reset_next_ids(lsctables.TableByName.values())
	"""
	for cls in classes:
		if cls.next_id is not None:
			cls.set_next_id(type(cls.next_id)(0))


#
# =============================================================================
#
#                                Column Element
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
		if isinstance(i, slice):
			return map(lambda r: getattr(r, self.asattribute), self.parentNode[i])
		else:
			return getattr(self.parentNode[i], self.asattribute)

	def __setitem__(self, i, value):
		"""
		Set the value in this column in row i.
		"""
		if isinstance(i, slice):
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

	def index(self, value):
		"""
		Return the smallest index of the row(s) with this column
		equal to value.
		"""
		for i in xrange(len(self.parentNode)):
			if getattr(self.parentNode[i], self.asattribute) == value:
				return i
		raise ValueError, value

	def __contains__(self, value):
		"""
		Returns True or False if there is or is not, respectively,
		a row containing val in this column.
		"""
		for i in xrange(len(self.parentNode)):
			if getattr(self.parentNode[i], self.asattribute) == value:
				return True
		return False

	def asarray(self):
		"""
		Construct a numpy array from this column.  Note that this
		creates a copy of the data, so modifications made to the
		array will not be recorded in the original document.
		"""
		if self.getAttribute("Type") not in ligolwtypes.NumericTypes:
			raise TypeError, "Column does not have numeric type"
		import numpy
		return numpy.fromiter(self, dtype = ligolwtypes.ToNumPyType[self.getAttribute("Type")])


#
# =============================================================================
#
#                                Stream Element
#
# =============================================================================
#


#
# A subclass of tokenizer.RowBuilder that interns strings.
#


class InterningRowBuilder(tokenizer.RowBuilder):
	"""
	This subclass of the tokenizer.RowBuilder class respects the
	"interning" hints provided by table definitions, and attempts to
	replace the values of row attributes associated with interned
	columns with references to shared instances of those values.  This
	results in a reduction in memory use which is small for most
	documents, but can be subtantial when dealing with poorly-designed
	tables containing large volumes of repeated information.
	
	The values are stored in a dictionary that is shared between all
	instances of this class, and which survives forever.  Nothing is
	ever naturally "uninterned", so the string dictionary grows without
	bound as more documents are processed.  This can be a problem in
	some use cases, and the work-around is to run

	>>> InterningRowBuilder.strings.clear()

	to reset the dictionary and appropriate points in the application.
	"""
	strings = {}
	def append(self, tokens):
		for row in tokenizer.RowBuilder.append(self, tokens):
			for col in self.interns:
				val = getattr(row, col)
				setattr(row, col, self.strings.setdefault(val, val))
			yield row


#
# Select the RowBuilder class to use when parsing tables.
#


RowBuilder = tokenizer.RowBuilder


#
# Stream class
#


class TableStream(ligolw.Stream):
	"""
	High-level Stream element for use inside Tables.  This element
	knows how to parse the delimited character stream into rows in the
	parent element, and knows how to turn the parent's rows back into a
	character stream.
	"""
	def __init__(self, attrs):
		ligolw.Stream.__init__(self, attrs)
		self._tokenizer = tokenizer.Tokenizer(self.getAttribute("Delimiter"))
		self._rowbuilder = None

	def config(self, parentNode):
		# some initialization that requires access to the
		# parentNode, and so cannot be done inside the __init__()
		# function.
		loadcolumns = set(parentNode.columnnames)
		if parentNode.loadcolumns is not None:
			# FIXME:  convert loadcolumns attributes to sets to
			# avoid the conversion.
			loadcolumns &= set(parentNode.loadcolumns)
		self._tokenizer.set_types([(colname in loadcolumns) and pytype or None for pytype, colname in zip(parentNode.columnpytypes, parentNode.columnnames)])
		columnnames = [name for name in parentNode.columnnames if name in loadcolumns]
		# FIXME:  convert interncolumns attributes to sets to
		# simplify computing the intersection
		interncolumns = [name for name in (parentNode.interncolumns or set()) if name in columnnames]
		self._rowbuilder = RowBuilder(parentNode.RowType, columnnames, interncolumns)
		return self

	def appendData(self, content):
		# tokenize buffer, pack into row objects, and append to
		# table
		for row in self._rowbuilder.append(self._tokenizer.append(content)):
			self.parentNode.append(row)

	def unlink(self):
		"""
		Break internal references within the document tree rooted
		on this element to promote garbage collection.
		"""
		self._tokenizer = None
		self._rowbuilder = None
		ligolw.Stream.unlink(self)

	def write(self, file = sys.stdout, indent = u""):
		# loop over parent's rows.  This is complicated because we
		# need to not put a delimiter at the end of the last row
		# unless it ends with a null token
		file.write(self.start_tag(indent))
		rowdumper = tokenizer.RowDumper(self.parentNode.columnnames, [ligolwtypes.FormatFunc[coltype] for coltype in self.parentNode.columntypes], self.getAttribute("Delimiter"))
		rowdumper.dump(self.parentNode)
		try:
			line = rowdumper.next()
		except StopIteration:
			# table is empty
			pass
		else:
			# write first row
			newline = u"\n" + indent + ligolw.Indent
			file.write(newline)
			# the xmlescape() call replaces things like "<"
			# with "&lt;" so that the string will not confuse
			# an XML parser when the file is read.  turning
			# "&lt;" back into "<" during file reading is
			# handled by the XML parser, so there is no code
			# in Glue related to that.
			file.write(xmlescape(line))
			# now add delimiter and write the remaining rows
			newline = rowdumper.delimiter + newline
			for line in rowdumper:
				file.write(newline)
				file.write(xmlescape(line))
			if rowdumper.tokens and rowdumper.tokens[-1] == u"":
				# the last token of the last row was null:
				# add a final delimiter to indicate that a
				# token is present
				file.write(rowdumper.delimiter)
		file.write(u"\n" + self.end_tag(indent) + u"\n")

	# FIXME: This function is for the metaio library:  metaio cares
	# what order the attributes of XML tags come in.  This function
	# will be removed when the metaio library is fixed.
	def start_tag(self, indent):
		"""
		Generate the element start tag.
		"""
		return indent + u"<%s Name=\"%s\" Type=\"%s\" Delimiter=\"%s\">" % (self.tagName, self.getAttribute("Name"), self.getAttribute("Type"), self.getAttribute("Delimiter"))


#
# =============================================================================
#
#                                Table Element
#
# =============================================================================
#


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
	interncolumns = None
	constraints = None
	how_to_index = None
	RowType = TableRow
	next_id = None

	def __init__(self, *attrs):
		"""
		Initialize
		"""
		ligolw.Table.__init__(self, *attrs)
		self.columnnames = []
		self.columntypes = []
		self.columnpytypes = []

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
			col, = getColumnsByName(self, name)
		except ValueError:
			# did not find exactly 1 matching child
			raise KeyError, name
		return col


	def appendColumn(self, name):
		"""
		Append a column named "name" to the table.  Returns the new
		child.  Raises ValueError if the table already has a column
		by that name, and KeyError if the validcolumns attribute of
		this table does not contain an entry for a column by that
		name.
		"""
		try:
			self.getColumnByName(name)
			# if we get here the table already has that column
			raise ValueError, "duplicate Column '%s'" % name
		except KeyError:
			pass
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
		del self.columnnames[:]
		del self.columntypes[:]
		del self.columnpytypes[:]
		for child in self.getElementsByTagName(ligolw.Column.tagName):
			colname = StripColumnName(child.getAttribute("Name"))
			llwtype = child.getAttribute("Type")
			if self.validcolumns is not None:
				try:
					if self.validcolumns[colname] != llwtype:
						raise ligolw.ElementError, "invalid type '%s' for Column '%s' in Table '%s', expected type '%s'" % (llwtype, child.getAttribute("Name"), self.getAttribute("Name"), self.validcolumns[colname])
				except KeyError:
					raise ligolw.ElementError, "invalid Column '%s' for Table '%s'" % (child.getAttribute("Name"), self.getAttribute("Name"))
			if colname in self.columnnames:
				raise ligolw.ElementError, "duplicate Column '%s' in Table '%s'" % (child.getAttribute("Name"), self.getAttribute("Name"))
			self.columnnames.append(colname)
			self.columntypes.append(llwtype)
			try:
				self.columnpytypes.append(ligolwtypes.ToPyType[llwtype])
			except KeyError:
				raise ligolw.ElementError, "unrecognized Type '%s' for Column '%s' in Table '%s'" % (llwtype, child.getAttribute("Name"), self.getAttribute("Name"))

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
				raise ligolw.ElementError, "Stream name '%s' does not match Table name '%s'" % (child.getAttribute("Name"), self.getAttribute("Name"))

	def _end_of_columns(self):
		"""
		Called during parsing to indicate that the last Column
		child element has been added.
		"""
		pass

	def _end_of_rows(self):
		"""
		Called during parsing to indicate that the last row has
		been added.
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
	# Row ID manipulation
	#

	def get_next_id(cls):
		"""
		Returns the current value of the next_id class attribute,
		and increments the next_id class attribute by 1.
		"""
		id = cls.next_id
		cls.next_id += 1
		return id
	get_next_id = classmethod(get_next_id)

	def set_next_id(cls, id):
		"""
		Sets the value of the next_id class attribute.  This is a
		convenience function to help prevent accidentally assigning
		a value to an instance attribute instead of the class
		attribute.
		"""
		cls.next_id = id
	set_next_id = classmethod(set_next_id)

	def sync_next_id(self):
		"""
		Determines the highest-numbered ID in this table, and sets
		the table's next_id attribute to the next highest ID in
		sequence.  If the next_id attribute is already set to a
		value greater than the max found, then it is left
		unmodified.  The return value is the ID identified by this
		method.  If the table's next_id attribute is None, then
		this function is a no-op.

		Note that tables of the same name typically share a common
		next_id attribute (it is a class attribute, not an
		attribute of each instance).  This is enforced so that IDs
		can be generated that are unique across all tables in the
		document.  Running sync_next_id() on all the tables in a
		document that are of the same type will have the effect of
		setting the ID to the next ID higher than any ID in any of
		those tables.
		"""
		if self.next_id is not None:
			if len(self):
				n = max(self.getColumnByName(self.next_id.column_name)) + 1
			else:
				n = type(self.next_id)(0)
			if n > self.next_id:
				self.set_next_id(n)
		return self.next_id

	def updateKeyMapping(self, mapping):
		"""
		Used as the first half of the row key reassignment
		algorithm.  Accepts a dictionary mapping old key --> new
		key.  Iterates over the rows in this table, using the
		table's next_id attribute to assign a new ID to each row,
		recording the changes in the mapping.  Returns the mapping.
		Raises ValueError if the table's next_id attribute is None.
		"""
		if self.next_id is None:
			raise ValueError, self
		try:
			column = self.getColumnByName(self.next_id.column_name)
		except KeyError:
			# table is missing its ID column, this is a no-op
			return mapping
		for i, old in enumerate(column):
			if old in mapping:
				column[i] = mapping[old]
			else:
				column[i] = mapping[old] = self.get_next_id()
		return mapping

	def applyKeyMapping(self, mapping):
		"""
		Used as the second half of the key reassignment algorithm.
		Loops over each row in the table, replacing references to
		old row keys with the new values from the mapping.
		"""
		for coltype, colname in zip(self.columntypes, self.columnnames):
			if coltype in ligolwtypes.IDTypes and (self.next_id is None or colname != self.next_id.column_name):
				column = self.getColumnByName(colname)
				for i, old in enumerate(column):
					if old in mapping:
						column[i] = mapping[old]


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
		return TableStream(attrs).config(self.current)
	return __parent_startStream(self, attrs)


def endStream(self):
	# stream tokenizer uses delimiter to identify end of each token, so
	# add a final delimiter to induce the last token to get parsed but
	# only if there's something other than whitespace left in the
	# tokenizer's buffer.  Also call _end_of_rows() hook.
	if self.current.parentNode.tagName == ligolw.Table.tagName:
		if not self.current._tokenizer.data.isspace():
			self.current.appendData(self.current.getAttribute("Delimiter"))
		self.current.parentNode._end_of_rows()
	else:
		__parent_endStream(self)


def startTable(self, attrs):
	return Table(attrs)


def endTable(self):
	# Table elements are allowed to contain 0 Stream children, but
	# _end_of_columns() and _end_of_rows() hooks must be called
	# regardless, so we do that here if needed.
	if self.current.childNodes[-1].tagName != ligolw.Stream.tagName:
		self.current._end_of_columns()
		self.current._end_of_rows()


ligolw.LIGOLWContentHandler.startColumn = startColumn
ligolw.LIGOLWContentHandler.startStream = startStream
ligolw.LIGOLWContentHandler.endStream = endStream
ligolw.LIGOLWContentHandler.startTable = startTable
ligolw.LIGOLWContentHandler.endTable = endTable
