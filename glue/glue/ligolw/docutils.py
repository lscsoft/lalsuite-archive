"""
High-level document manipulation utilities.
"""

import ligolw
import lsctables


#
# General utilities.
#

def MergeElements(elem1, elem2):
	"""
	Move the children of elem2 to elem1, and unlink elem2 from its
	parent.  The return value is elem1.  If the two elements are
	tables, then the contents of the second table into the first table,
	and unlink the second table from the document tree.
	"""
	if elem1.tagName != elem2.tagName:
		raise ligolw.ElementError, "MergeElements(): elements must have same names"
	if elem1.tagName == ligolw.LIGO_LW.tagName:
		# copy children;  LIGO_LW elements have no attributes
		map(elem1.appendChild, elem2.childNodes)
	elif elem1.tagName == ligolw.Table.tagName:
		# FIXME: this doesn't work if one table has no stream element,
		# which should probably still be OK.
		if elem1.compare(elem2):
			raise ligolw.ElementError, "MergeElements(): Table elements must compare as equivalent."
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

def MergeCompatibleTables(elem):
	"""
	Below the given element, find all Tables whose structure is
	described in lsctables, and merge compatible ones of like type.
	That is, merge all SnglBurstTables that have the same columns into
	a single table, etc..
	"""
	for tname in lsctables.TableByName.keys():
		tables = [t for t in elem.getElementsByTagName(ligolw.Table.tagName) if t.getAttribute("Name") == tname]
		for i in range(1, len(tables)):
			if not tables[0].compare(tables[i]):
				MergeElements(tables[0], tables[i])


def RemapProcessIDs(elem, mapping):
	"""
	Recurse over all Table elements whose structure is described in the
	module lsctables, and remap the process IDs of all rows according
	to mapping.  Rows with process IDs not named in the mapping are
	ignored.
	"""
	for table in [t for t in elem.getElementsByTagName(ligolw.Table.tagName) if t.getAttribute("Name") in lsctables.TableByName.keys()]:
		for row in table.rows:
			try:
				row.process_id = mapping[row.process_id]
			except KeyError:
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
