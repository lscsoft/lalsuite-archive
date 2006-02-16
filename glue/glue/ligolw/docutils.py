"""
High-level document manipulation utilities.
"""

__author__ = "Kipp Cannon <kipp@gravity.phys.uwm.edu>"
__date__ = "$Date$"
__version__ = "$Revision$"

import ligolw
import metaio
import lsctables


#
# General utilities.
#

def MergeElements(elem1, elem2):
	"""
	Move the children of elem2 to elem1, and unlink elem2 from its
	parent.  The return value is elem1.
	
	If the two elements are tables, then more the rows of the second
	table into the first table, and unlink the second table from the
	document tree.  The table, column, and stream names of the first
	table are retained, as well as the (optional) comment child
	element.
	"""
	if elem1.tagName != elem2.tagName:
		raise ligolw.ElementError, "MergeElements(): elements must have same names"
	if elem1.tagName == ligolw.LIGO_LW.tagName:
		# copy children;  LIGO_LW elements have no attributes
		map(elem1.appendChild, elem2.childNodes)
	elif elem1.tagName == ligolw.Table.tagName:
		# copy rows
		elem1.rows.extend(elem2.rows)
	else:
		raise ligolw.ElementError, "MergeElements(): can't merge %s elements." % elem1.tagName
	if elem2.parentNode:
		elem2.parentNode.removeChild(elem2)
	return elem1


#
# Table manipulation utilities.
#

def HasNonLSCTables(elem):
	"""
	Return True if the document tree below elem contains non-LSC
	tables, otherwise return False.
	"""
	for table in elem.getElementsByTagName(ligolw.Table.tagName):
		if metaio.StripTableName(table.getAttribute("Name")) not in lsctables.TableByName.keys():
			return True
	return False


def TablesCanBeMerged(a, b):
	"""
	Return True if the two tables a and b can be merged.  This means
	they have equivalent names, and equivalent columns according to
	LIGO LW name conventions.
	"""
	if metaio.CompareTableNames(a.getAttribute("Name"), b.getAttribute("Name")) != 0:
		return False
	acols = [(metaio.StripColumnName(col.getAttribute("Name")), col.getAttribute("Type")) for col in a.getElementsByTagName(ligolw.Column.tagName)]
	bcols = [(metaio.StripColumnName(col.getAttribute("Name")), col.getAttribute("Type")) for col in b.getElementsByTagName(ligolw.Column.tagName)]
	for acol in acols:
		if acol not in bcols:
			return False
	return True


def MergeCompatibleTables(elem):
	"""
	Below the given element, find all Tables whose structure is
	described in lsctables, and merge compatible ones of like type.
	That is, merge all SnglBurstTables that have the same columns into
	a single table, etc..
	"""
	for tname in lsctables.TableByName.keys():
		tables = metaio.getTablesByName(elem, tname)
		for i in range(1, len(tables)):
			if TablesCanBeMerged(tables[0], tables[i]):
				MergeElements(tables[0], tables[i])
	return elem


#
# Table row cross-reference utilities.
#

def NewILWDs(table, column_name):
	"""
	From the LSC table, return a compatible ILWD iterator object,
	initialized to the next unique ID following those found in the
	table.
	"""
	try:
		n = max(map(lambda id: int(lsctables.ILWDID(id)), table.dict.keys()))
	except ValueError:
		n = -1
	return lsctables.ILWD(metaio.StripTableName(table.getAttribute("Name")), column_name, n + 1)


def makeReference(elem):
	"""
	Run the makeReference() method on all LSC tables below elem,
	constructing references to other tables under elem.
	"""
	for table in lsctables.getLSCTables(elem):
		try:
			table.makeReference(elem)
		except AttributeError:
			# table is missing a cross-reference column.
			pass


def deReference(elem):
	"""
	Run the deReference() method on all LSC tables below elem.
	"""
	for table in lsctables.getLSCTables(elem):
		try:
			table.deReference()
		except AttributeError:
			# table is missing a cross-reference column.
			pass


def NewIDs(elem, ilwditers):
	"""
	Using the dictionary of table name and ILWD iterator object pairs,
	recurse over all tables below elem whose names are in the
	dictionary, and use the corresponding ILWD iterator object to
	construct a mapping of old row keys to new row keys.  Finally,
	apply the mapping to all rows.
	"""
	for tablename, ilwditer in ilwditers.iteritems():
		for table in metaio.getTablesByName(elem, tablename):
			keymap = {}
			try:
				for oldkey in table.dict.keys():
					keymap[oldkey] = ilwditer.next()
			except AttributeError:
				# table is missing its ID column
				continue
			if not len(keymap):
				continue
			for row in table.rows:
				try:
					row._set_key(keymap[row._get_key()])
				except KeyError:
					# row has a key not listed in the
					# table's dictionary.
					pass


#
# Utilities for partial document loading.
#

class PartialLIGOLWContentHandler(lsctables.LIGOLWContentHandler):
	"""
	LIGO LW content handler object with the ability to strip unwanted
	portions of the document from the input stream.  Useful, for
	example, when one wishes to read only a single table from the XML.
	"""
	def __init__(self, document, filter):
		"""
		Only those elements for which filter(name, attrs) evaluates
		to True, and the children of those elements, will be
		loaded.
		"""
		lsctables.LIGOLWContentHandler.__init__(self, document)
		self.filter = filter
		self.filtered_depth = 0

	def startElement(self, name, attrs):
		if self.filtered_depth > 0 or self.filter(name, attrs):
			lsctables.LIGOLWContentHandler.startElement(self, name, attrs)
			self.filtered_depth += 1

	def endElement(self, name):
		if self.filtered_depth > 0:
			self.filtered_depth -= 1
			lsctables.LIGOLWContentHandler.endElement(self, name)
