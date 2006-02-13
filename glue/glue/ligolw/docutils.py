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
		tables = metaio.getTablesByName(elem, tname)
		for i in range(1, len(tables)):
			if not tables[0].compare(tables[i]):
				MergeElements(tables[0], tables[i])
	return elem


def RemapProcessIDs(elem, mapping):
	"""
	Recurse over all Table elements whose structure is described in the
	module lsctables, and remap the process IDs of all rows according
	to mapping.  Rows with process IDs not named in the mapping are
	ignored.
	"""
	for table in lsctables.getLSCTables(elem):
		for row in table:
			try:
				row.process_id = mapping[row.process_id]
			except KeyError:
				pass


#
# Process manipulation utilities.
#

class Process(object):
	def __init__(self, proc, params, summary):
		self.process = proc
		self.processparams = params
		self.searchsummary = summary

	def __cmp__(self, other):
		def listcmp(a, b):
			for i in range(min(len(a), len(b))):
				result = a[i].cmp(b[i])
				if result:
					return result
			return len(a) - len(b)
		result = self.process.cmp(other.process)
		if not result:
			result = listcmp(self.processparams, other.processparams)
			if not result:
				result = listcmp(self.searchsummary, other.searchsummary)
		return result


class ProcessList(object):
	def __init__(self, doc):
		def gettable(doc, Type):
			tables = lsctables.getTablesByType(doc, Type)
			if len(tables) != 1:
				raise Exception, "ProcessList(doc): doc must contain exactly 1 %s Table." % Type.tableName
			return tables[0]

		self.doc = doc
		self.processtable = gettable(doc, lsctables.ProcessTable)
		self.processparamstable = gettable(doc, lsctables.ProcessParamsTable)
		self.searchsummarytable = gettable(doc, lsctables.SearchSummaryTable)

	def __len__(self):
		return len(self.processtable)

	def __getitem__(self, key):
		proc = self.processtable.dict[key]
		try:
			params = self.processparamstable.dict[key]
		except KeyError:
			params = []
		try:
			summary = self.searchsummarytable.dict[key]
		except KeyError:
			summary = []
		return Process(proc, params, summary)

	def __setitem__(self, key, value):
		self.processtable.dict[key] = value.process
		self.processparamstable.dict[key] = value.processparams
		self.searchsummarytable.dict[key] = value.searchsummary

	def __delitem__(self, key):
		del self.processtable.dict[key]
		del self.processparamstable.dict[key]
		del self.searchsummarytable.dict[key]

	def __contains__(self, key):
		return key in self.processtable

	def keys(self):
		return self.processtable.dict.keys()


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
