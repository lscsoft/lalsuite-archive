import re
import sys
from xml import sax

import ligolw


StringTypes = ["ilwd:char", "ilwd:char_u", "lstring"]
IntTypes = ["int_2s", "int_4s", "int_8s"]
FloatTypes = ["real_4", "real_8"]
Types = StringTypes + IntTypes + FloatTypes


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
		self.token_pattern = re.compile(r"""\s*(?:"([^"]*)")|(?:([^""" + self.attributes["Delimiter"] + r"""\s]+))\s*""" + self.attributes["Delimiter"])
		self.tokens = []

	def appendData(self, content):
		# append new data to buffer
		ligolw.Stream.appendData(self, content)

		# make sure we are inside a Table
		if self.parent.tagName != "Table":
			return

		# move tokens from buffer to token list
		match = None
		for match in self.token_pattern.finditer(self.pcdata):
			self.tokens.append(match.group(match.lastindex))
		if match != None:
			self.pcdata = self.pcdata[match.end():]

		# construct row objects from tokens
		while len(self.tokens) >= self.parent.ncolumns:
			row = self.parent.RowType()
			for i, column in enumerate(self.parent.getElementsByTagName("Column")):
				key = column.attributes["Name"].split(":")[-1]
				value = column.attributes["Type"]
				if value in StringTypes:
					setattr(row, key, self.tokens[i])
					#type(row).__dict__[key].__set__(row, self.tokens[i])
				elif value in IntTypes:
					try:
						setattr(row, key, int(self.tokens[i]))
					except ValueError, e:
						raise ligolw.ElementError, "Stream parsing error near tokens %s: %s" % (str(self.tokens), str(e))
				elif value in FloatTypes:
					try:
						setattr(row, key, float(self.tokens[i]))
					except ValueError, e:
						raise ligolw.ElementError, "Stream parsing error near tokens %s: %s" % (str(self.tokens), str(e))
				else:
					pass
			self.tokens = self.tokens[i+1:]
			self.parent.rows.append(row)

	def _rowstr(self, row, columns):
		strs = []
		for column in columns:
			if column.attributes["Type"] in StringTypes:
				strs += ["\"" + getattr(row, column.attributes["Name"].split(":")[-1]) + "\""]
			else:
				strs += [str(getattr(row, column.attributes["Name"].split(":")[-1]))]
		return self.attributes["Delimiter"].join(strs)

	def write(self, file = sys.stdout, indent = ""):
		columns = self.parent.getElementsByTagName("Column")

		# loop over parent's rows.  This is complicated because we
		# need to not put a comma at the end of the last row.
		print >>file, indent + self.start_tag()
		rowiter = iter(self.parent.rows)
		try:
			row = rowiter.next()
			file.write(indent + ligolw.Indent + self._rowstr(row, columns))
			while True:
				row = rowiter.next()
				file.write(self.attributes["Delimiter"] + "\n" + indent + ligolw.Indent + self._rowstr(row, columns))
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
		self.ncolumns = 0
		self.rows = []

	def appendChild(self, child):
		if child.tagName == "Column":
			self.ncolumns += 1
			if self.validcolumns != None:
				key = child.attributes["Name"].split(":")[-1]
				if key not in self.validcolumns.keys():
					raise ligolw.ElementError, "invalid Column name %s for Table" % child.attributes["Name"]
				if self.validcolumns[key] != child.attributes["Type"]:
					raise ligolw.ElementError, "invalid type %s for Column %s" % (child.attributes["Type"], child.attributes["Name"])
		elif child.tagName == "Stream":
			if child.attributes["Name"] != self.attributes["Name"]:
				raise ligolw.ElementError, "Stream name %s does not match Table name %s" % (child.attributes["Name"], self.attributes["Name"])
		ligolw.Table.appendChild(self, child)


class LIGOLWContentHandler(ligolw.LIGOLWContentHandler):
	"""
	ContentHandler that redirects Stream and Table elements to those
	defined in this module.
	"""
	def startStream(self, attrs):
		return Stream(attrs)

	def endStream(self):
		# stream tokenizer uses comma to identify end of each
		# token, so add a final comma to induce the last token to
		# get parsed.  FIXME: this is a hack, and hacks are the
		# source of bugs.
		self.current.appendData(self.current.attributes["Delimiter"])

	def startTable(self, attrs):
		return Table(attrs)
