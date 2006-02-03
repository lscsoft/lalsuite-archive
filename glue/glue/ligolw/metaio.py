try:
	import numarray
except:
	pass
import re
import sys
from xml import sax

import ligolw


#
# =============================================================================
#
#                               Type Information
#
# =============================================================================
#

StringTypes = ["char_s", "ilwd:char", "ilwd:char_u", "lstring", "string"]
IntTypes = ["int_2s", "int_2u", "int_4s", "int_4u", "int_8s", "int_8u", "int"]
FloatTypes = ["real_4", "real_8", "float", "double"]

Types = StringTypes + IntTypes + FloatTypes

ToNumArrayType = {
	"int_2s": "Int16",
	"int_2u": "UInt16",
	"int_4s": "Int32",
	"int_4u": "UInt32",
	"int_8s": "Int64",
	"int_8u": "UInt64",
	"int": "Int32",
	"real_4": "Float32",
	"real_8": "Float64",
	"float": "Float64",
	"double": "Float64"
}


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

ColumnPattern = re.compile(r"(?:\A\w+:|\A)(?P<FullName>(?:\w+:|\A)(?P<Name>\w+))\Z")


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


#
# =============================================================================
#
# Table Name Manipulation
#
# =============================================================================
#

# Regular expression used to extract the signifcant portion of a table or
# stream name, according to LIGO LW naming conventions.

TablePattern = re.compile(r"(?:\A[a-z0-9_]+:|\A)(?P<Name>[a-z0-9_]+:table)\Z")


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


#
# =============================================================================
#
#                             Unique ID Generator
#
# =============================================================================
#

class ILWD(object):
	def __init__(self, base, n = 0):
		self.base = base
		self.n = n

	def __str__(self):
		return "%s:%d" % (self.base, self.n)

	def __iter__(self):
		return self

	def next(self):
		self.n += 1
		return self


#
# =============================================================================
#
#                               Element Classes
#
# =============================================================================
#

class Column(ligolw.Column):
	"""
	High-level column element that knows how to turn column data from
	the table into an array.
	"""
	def asarray(self):
		attr = StripColumnName(self.getAttribute("Name"))
		if self.getAttribute("Type") in StringTypes:
			return [getattr(row, attr) for row in self.parentNode]
		else:
			return numarray.asarray([getattr(row, attr) for row in self.parentNode], type = ToNumArrayType[self.getAttribute("Type")])


class Stream(ligolw.Stream):
	"""
	High-level Stream element for use inside Tables.  This element
	knows how to parse the delimited character stream into rows in the
	parent element, and knows how to turn the parent's rows back into a
	character stream.
	"""
	def __init__(self, attrs):
		ligolw.Stream.__init__(self, attrs)
		self.tokenizer = re.compile(r"""\s*(?:"([^"]*)")|(?:([^""" + self.getAttribute("Delimiter") + r"""\s]+))\s*""" + self.getAttribute("Delimiter"))
		self.tokens = []

	def appendData(self, content):
		# append new data to buffer
		ligolw.Stream.appendData(self, content)

		# make sure we are inside a Table
		if self.parentNode.tagName != ligolw.Table.tagName:
			return

		# move tokens from buffer to token list
		match = None
		for match in self.tokenizer.finditer(self.pcdata):
			self.tokens.append(match.group(match.lastindex))
		if match != None:
			self.pcdata = self.pcdata[match.end():]

		# construct row objects from tokens, and append to parent
		while len(self.tokens) >= len(self.parentNode.columninfo):
			row = self.parentNode.RowType()
			for i, (colname, pytype) in enumerate(self.parentNode.columninfo):
				try:
					setattr(row, colname, pytype(self.tokens[i]))
				except ValueError, e:
					raise ligolw.ElementError, "Stream parsing error near tokens %s: %s" % (str(self.tokens), str(e))
				except AttributeError, e:
					pass
			self.tokens = self.tokens[i+1:]
			self.parentNode.appendRow(row)

	def _rowstr(self, row, columns):
		strs = []
		for column in columns:
			if column.getAttribute("Type") in StringTypes:
				strs.append("\"" + getattr(row, StripColumnName(column.getAttribute("Name"))) + "\"")
			else:
				strs.append(str(getattr(row, StripColumnName(column.getAttribute("Name")))))
		return self.getAttribute("Delimiter").join(strs)

	def write(self, file = sys.stdout, indent = ""):
		columns = self.parentNode.getElementsByTagName(ligolw.Column.tagName)

		# loop over parent's rows.  This is complicated because we
		# need to not put a comma at the end of the last row.
		print >>file, indent + self.start_tag()
		rowiter = iter(self.parentNode)
		try:
			file.write(indent + ligolw.Indent + self._rowstr(rowiter.next(), columns))
			while True:
				file.write(self.getAttribute("Delimiter") + "\n" + indent + ligolw.Indent + self._rowstr(rowiter.next(), columns))
		except StopIteration:
			file.write("\n")
		print >>file, indent + self.end_tag()


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
	RowType = TableRow

	def __init__(self, *attrs):
		ligolw.Table.__init__(self, *attrs)
		self.columninfo = []
		self.rows = []

	def __len__(self):
		return len(self.rows)

	def __iter__(self):
		return iter(self.rows)

	def getColumnByName(self, name):
		try:
			return self.getElements(lambda e: (e.tagName == ligolw.Column.tagName) and (CompareColumnNames(e.getAttribute("Name"), name) == 0))[0]
		except IndexError:
			raise KeyError, "no Column matching name %s" % name

	def appendChild(self, child):
		if child.tagName == ligolw.Column.tagName:
			colname = StripColumnName(child.getAttribute("Name"))
			llwtype = child.getAttribute("Type")
			if self.validcolumns != None:
				if colname not in self.validcolumns.keys():
					raise ligolw.ElementError, "invalid Column name %s for Table" % child.getAttribute("Name")
				if self.validcolumns[colname] != llwtype:
					raise ligolw.ElementError, "invalid type %s for Column %s" % (llwtype, child.getAttribute("Name"))
			if colname in [c[0] for c in self.columninfo]:
				raise ligolw.ElementError, "duplicate Column %s" % child.getAttribute("Name")
			if llwtype in StringTypes:
				self.columninfo.append((colname, str))
			elif llwtype in IntTypes:
				self.columninfo.append((colname, int))
			elif llwtype in FloatTypes:
				self.columninfo.append((colname, float))
			else:
				raise ligolw.ElementError, "unrecognized Type attribute %s for Column element" % llwtype
		elif child.tagName == ligolw.Stream.tagName:
			if child.getAttribute("Name") != self.getAttribute("Name"):
				raise ligolw.ElementError, "Stream name %s does not match Table name %s" % (child.getAttribute("Name"), self.getAttribute("Name"))
		ligolw.Table.appendChild(self, child)

	def removeChild(self, child):
		"""
		Remove a child from this element.  The child element is
		returned, and it's parentNode element is reset.
		"""
		ligolw.Table.removeChild(self, child)
		if child.tagName == ligolw.Column.tagName:
			for n in [n for n, item in enumerate(self.columninfo) if item[0] == StripColumnName(child.getAttribute("Name"))]:
				del self.columinfo[n]
		return child

	def appendRow(self, row):
		"""
		Append row to the list of rows for this table.
		"""
		self.rows.append(row)

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
# =============================================================================
#
#                               Content Handler
#
# =============================================================================
#

class LIGOLWContentHandler(ligolw.LIGOLWContentHandler):
	"""
	ContentHandler that redirects Column, Stream and Table elements to
	those defined in this module.
	"""
	def startColumn(self, attrs):
		return Column(attrs)

	def startStream(self, attrs):
		return Stream(attrs)

	def endStream(self):
		# stream tokenizer uses comma to identify end of each
		# token, so add a final comma to induce the last token to
		# get parsed.  FIXME: this is a hack, and hacks are the
		# source of bugs.
		self.current.appendData(self.current.getAttribute("Delimiter"))

	def startTable(self, attrs):
		return Table(attrs)
