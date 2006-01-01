try:
	import numarray
except:
	pass
import re
import sys
from xml import sax

import ligolw


StringTypes = ["ilwd:char", "ilwd:char_u", "lstring", "string"]
IntTypes = ["int_2s", "int_4s", "int_8s", "int"]
FloatTypes = ["real_4", "real_8", "float"]
Types = StringTypes + IntTypes + FloatTypes
ToNumArrayType = {
	"int_2s": "Int16",
	"int_4s": "Int32",
	"int_8s": "Int64",
	"int": "Int32",
	"real_4": "Float32",
	"real_8": "Float64",
	"float": "Float64"
}


class Column(ligolw.Column):
	"""
	High-level column element that knows how to turn column data from
	the table into an array.
	"""
	def asarray(self):
		name = self.getAttribute("Name").split(":")[-1]
		return numarray.asarray([getattr(row, name) for row in self.parentNode.rows], type = ToNumArrayType[self.getAttribute("Type")])


class Stream(ligolw.Stream):
	"""
	High-level Stream element for use inside Tables.  This element
	knows how to parse the delimited character stream into rows in the
	parent element, and knows how to turn the parent's rows back into a
	character stream.
	"""
	def __init__(self, attrs):
		ligolw.Stream.__init__(self, attrs)

		# initialize the tokenizer.
		self.token_pattern = re.compile(r"""\s*(?:"([^"]*)")|(?:([^""" + self.getAttribute("Delimiter") + r"""\s]+))\s*""" + self.getAttribute("Delimiter"))
		self.tokens = []

	def appendData(self, content):
		# append new data to buffer
		ligolw.Stream.appendData(self, content)

		# make sure we are inside a Table
		if self.parentNode.tagName != "Table":
			return

		# move tokens from buffer to token list
		match = None
		for match in self.token_pattern.finditer(self.pcdata):
			self.tokens.append(match.group(match.lastindex))
		if match != None:
			self.pcdata = self.pcdata[match.end():]

		# construct row objects from tokens
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
				strs += ["\"" + getattr(row, column.getAttribute("Name").split(":")[-1]) + "\""]
			else:
				strs += [str(getattr(row, column.getAttribute("Name").split(":")[-1]))]
		return self.getAttribute("Delimiter").join(strs)

	def write(self, file = sys.stdout, indent = ""):
		columns = self.parentNode.getElementsByTagName("Column")

		# loop over parent's rows.  This is complicated because we
		# need to not put a comma at the end of the last row.
		print >>file, indent + self.start_tag()
		rowiter = iter(self.parentNode.rows)
		try:
			row = rowiter.next()
			file.write(indent + ligolw.Indent + self._rowstr(row, columns))
			while True:
				row = rowiter.next()
				file.write(self.getAttribute("Delimiter") + "\n" + indent + ligolw.Indent + self._rowstr(row, columns))
		except StopIteration:
			file.write("\n")
		print >>file, indent + self.end_tag()


class TableRow(object):
	"""
	Place-holder, default row type.
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

	def columnName(self, name):
		return ":".join(self.tableName.split(":")[:-1] + [name])

	def getColumnsByName(self, name):
		return self.getChildrenByAttributes({"Name": self.columnName(name)})

	def appendChild(self, child):
		if child.tagName == "Column":
			colname = child.getAttribute("Name").split(":")[-1]
			llwtype = child.getAttribute("Type")
			if self.validcolumns != None:
				if colname not in self.validcolumns.keys():
					raise ligolw.ElementError, "invalid Column name %s for Table" % child.getAttribute("Name")
				if self.validcolumns[colname] != llwtype:
					raise ligolw.ElementError, "invalid type %s for Column %s" % (llwtype, child.getAttribute("Name"))
			if llwtype in StringTypes:
				self.columninfo.append((colname, str))
			elif llwtype in IntTypes:
				self.columninfo.append((colname, int))
			elif llwtype in FloatTypes:
				self.columninfo.append((colname, float))
			else:
				raise ligolw.ElementError, "unrecognized Type attribute %s for Column element" % llwtype
		elif child.tagName == "Stream":
			if child.getAttribute("Name") != self.getAttribute("Name"):
				raise ligolw.ElementError, "Stream name %s does not match Table name %s" % (child.getAttribute("Name"), self.getAttribute("Name"))
		ligolw.Table.appendChild(self, child)

	def appendRow(self, row):
		self.rows.append(row)

	def filterRows(self, func):
		i = 0
		while i < len(self.rows):
			if not func(self.rows[i]):
				del self.rows[i]
			else:
				i += 1
		return self


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
