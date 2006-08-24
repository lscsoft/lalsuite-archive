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
LSC Table definitions.  These must be kept synchronized with the official
definitions in the LDAS CVS repository at
http://www.ldas-sw.ligo.caltech.edu/cgi-bin/cvsweb.cgi/ldas/dbms/db2/sql.
Maintainership of the table definitions is left as an excercise to
interested users.
"""

__author__ = "Kipp Cannon <kipp@gravity.phys.uwm.edu>"
__date__ = "$Date$"[7:-2]
__version__ = "$Revision$"[11:-2]

import re
from xml import sax

from glue import lal
from glue import segments
import ligolw
import table
import types


#
# =============================================================================
#
#                            Convenience Functions
#
# =============================================================================
#

def New(Type, columns = None):
	"""
	Convenience function for constructing pre-defined LSC tables.  The
	optional columns argument is a list of the names of the columns the
	table should be constructed with.  If columns = None, then the
	table is constructed with all valid columns included (pass columns
	= [] to create a table with no columns).

	Example:

	>>> import lsctables
	>>> new = lsctables.New(lsctables.ProcessTable)
	"""
	new = Type(sax.xmlreader.AttributesImpl({u"Name": Type.tableName}))
	if columns != None:
		for key in columns:
			if key not in new.validcolumns.keys():
				raise ligolw.ElementError, "New(): invalid column \"%s\" for table \"%s\"." % (key, new.tableName)
	for key, value in new.validcolumns.items():
		if (columns == None) or (key in columns):
			new.appendChild(table.Column(sax.xmlreader.AttributesImpl({u"Name": ":".join(Type.tableName.split(":")[:-1]) + ":" + key, u"Type": value})))
	new.appendChild(table.TableStream(sax.xmlreader.AttributesImpl({u"Name": Type.tableName})))
	return new


def IsTableElement(Type, elem):
	"""
	Convenience function to check that an element is a Table of type
	Type.
	"""
	if elem.tagName != ligolw.Table.tagName:
		return False
	return table.CompareTableNames(elem.getAttribute("Name"), Type.tableName) == 0


def IsTableProperties(Type, tagname, attrs):
	"""
	Convenience function to check that the given tag name and
	attributes match those of a Table of type Type.
	"""
	if tagname != ligolw.Table.tagName:
		return False
	return table.CompareTableNames(attrs["Name"], Type.tableName) == 0


def getTablesByType(elem, Type):
	"""
	Return a list of tables of type Type under elem.
	"""
	return table.getTablesByName(elem, Type.tableName)


def getLSCTables(elem):
	"""
	Return a list of all LSC tables under elem.
	"""
	return elem.getElements(lambda e: (e.tagName == ligolw.Table.tagName) and (table.StripTableName(e.getAttribute("Name")) in TableByName.keys()))


def HasNonLSCTables(elem):
	"""
	Return True if the document tree below elem contains non-LSC
	tables, otherwise return False.
	"""
	for t in elem.getElementsByTagName(ligolw.Table.tagName):
		if table.StripTableName(t.getAttribute("Name")) not in TableByName:
			return True
	return False


#
# =============================================================================
#
#                              ILWD Manipulation
#
# =============================================================================
#

# Regular expression to extract the parts of a row ID according to the LIGO
# LW naming conventions.

ILWDPattern = re.compile(r"(?P<Table>\w+):(?P<Column>\w+):(?P<ID>\d+)")


def ILWDTableName(ilwdchar):
	"""
	Return the table name part of the row ID ilwdchar.  ValueError is
	raised if the ID cannot be parsed.
	"""
	try:
		return ILWDPattern.search(ilwdchar).group("Table")
	except AttributeError:
		raise ValueError, "unrecognized ID \"%s\"" % repr(ilwdchar)


def ILWDColumnName(ilwdchar):
	"""
	Return the column name part of the row ID ilwdchar.  ValueError is
	raised if the ID cannot be parsed.
	"""
	try:
		return ILWDPattern.search(ilwdchar).group("Column")
	except AttributeError:
		raise ValueError, "unrecognized ID \"%s\"" % repr(ilwdchar)


def ILWDID(ilwdchar):
	"""
	Return the ID part of the row ID ilwdchar.  ValueError is raised if
	the ID cannot be parsed.
	"""
	try:
		return int(ILWDPattern.search(ilwdchar).group("ID"))
	except AttributeError:
		raise ValueError, "unrecognized ID \"%s\"" % repr(ilwdchar)


def getTablesFromILWD(elem, ilwdchar):
	"""
	Return a list of all tables below elem which can contain rows
	having IDs matching ilwdchar.
	"""
	return elem.getElements(lambda e: (e.tagName == ligolw.Table.tagName) and (table.CompareTableNames(e.getAttribute("Name"), IDTableName(ilwdchar)) == 0))


def FindILWD(tables, ilwdchar):
	"""
	Search all rows in the list of tables for one with the given ID,
	and return a reference to the row object.  If the ID is not unique,
	the first row found is returned.  KeyError is raised if no rows are
	found or the ID cannot be parsed.
	"""
	for t in tables:
		try:
			return t.dict[ilwdchar]
		except KeyError:
			pass
	raise KeyError, repr(ilwdchar)


class ILWD(object):
	"""
	Unique ILWD generator.
	"""
	def __init__(self, table_name, column_name, n = 0):
		"""
		Initialize an ILWD object.  table_name and column_name are
		the names of the table and column within the table for
		which these will be used as IDs, eg., "process" and
		"process_id".  The optional n parameter sets the starting
		value for the numeric suffix in the ilwd:char
		representation of ID.
		"""
		self.table_name = table.StripTableName(table_name)
		self.column_name = table.StripColumnName(column_name)
		self.n = n

	def __str__(self):
		"""
		Return an ilwd:char string representation of this ID.
		"""
		return "%s:%s:%d" % (self.table_name, self.column_name, self.n)

	def __cmp__(self, other):
		"""
		Compare IDs first by the base string, then by n.
		"""
		return cmp((self.table_name, self.column_name, self.n), (other.table_name, other.column_name, other.n))

	def __getitem__(self, n):
		return "%s:%s:%d" % (self.table_name, self.column_name, n)

	def __iter__(self):
		return self

	def next(self):
		s = str(self)
		self.n += 1
		return s


def NewILWDs(table_elem, column_name):
	"""
	From the LSC table, return a compatible ILWD iterator object,
	initialized to the next unique ID following those found in the
	table.
	"""
	try:
		n = max(map(ILWDID, table_elem.dict.keys()))
	except ValueError:
		n = -1
	return ILWD(table.StripTableName(table_elem.getAttribute("Name")), column_name, n + 1)


def makeReference(elem):
	"""
	Run the makeReference() method on all LSC tables below elem,
	constructing references to other tables under elem.
	"""
	for table_elem in getLSCTables(elem):
		try:
			table_elem.makeReference(elem)
		except AttributeError:
			# table_elem is missing a cross-reference column.
			pass


def deReference(elem):
	"""
	Run the deReference() method on all LSC tables below elem.
	"""
	for table_elem in getLSCTables(elem):
		try:
			table_elem.deReference()
		except AttributeError:
			# table_elem is missing a cross-reference column.
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
		for table_elem in table.getTablesByName(elem, tablename):
			keymap = {}
			try:
				for oldkey in table_elem.dict.keys():
					keymap[oldkey] = ilwditer.next()
			except AttributeError:
				# table_elem is missing its ID column
				continue
			if not len(keymap):
				continue
			for row in table_elem:
				try:
					row._set_key(keymap[row._get_key()])
				except KeyError:
					# row has a key not listed in the
					# table's dictionary.
					pass


#
# =============================================================================
#
#                 Table and Row Class With a Mapping Protocol
#
# =============================================================================
#

class LSCTableUniqueItemIter(object):
	def __init__(self, dict):
		self.iter = iter(dict.rows)

	def __iter__(self):
		return self

	def next(self):
		row = self.iter.next()
		return (row._get_key(), row)


class LSCTableUniqueKeyIter(object):
	def __init__(self, dict):
		self.iter = iter(dict.rows)

	def __iter__(self):
		return self

	def next(self):
		return self.iter.next()._get_key()


class LSCTableUniqueDict(object):
	"""
	Class for implementing the Python mapping protocol on a list of table
	rows when each row has a unique key.
	"""
	def __init__(self, table_elem):
		"""
		Initialize the mapping on the list of rows.
		"""
		self.table = table_elem

	def __len__(self):
		"""
		Return the number of rows.
		"""
		return len(self.table)

	def __getitem__(self, key):
		"""
		Retrieve a row by key.
		"""
		for row in self.table:
			if row._has_key(key):
				return row
		raise KeyError, repr(key)

	def __setitem__(self, key, value):
		"""
		Set a row by key.  Note: the row key carried by value need
		not equal key, but this might not be allowed in the future.
		"""
		for i in xrange(len(self.table)):
			if self.table[i]._has_key(key):
				self.table[i] = value
				return
		# FIXME: should we call _set_key() on value to force it to have
		# the key that was searched for?
		self.append(value)

	def __delitem__(self, key):
		"""
		Delete a row by key.
		"""
		for i in xrange(len(self.table)):
			if self.table[i]._has_get(key):
				del self.table[i]
				return
		raise KeyError, repr(key)

	def __iter__(self):
		"""
		Return an iterator over the keys.
		"""
		return LSCTableUniqueKeyIter(self)

	iterkeys = __iter__

	def __contains__(self, key):
		"""
		Return True if a row has key equal to key, otherwise return
		False.
		"""
		for row in self.table:
			if row._has_key(key):
				return True
		return False

	has_key = __contains__

	def keys(self):
		"""
		Return a list of the keys.
		"""
		return [row._get_key() for row in self.table]

	def iteritems(self):
		"""
		Return an iterator over (key, value) pairs.
		"""
		return LSCTableUniqueItemIter(self)

	def itervalues(self):
		"""
		Return an iterator over rows.
		"""
		return iter(self.table)


class LSCTableMultiItemIter(object):
	def __init__(self, dict):
		self.dict = dict
		self.iter = iter(dict.keys())

	def __iter__(self):
		return self

	def next(self):
		key = self.iter.next()
		return (key, self.dict[key])


class LSCTableMultiDict(LSCTableUniqueDict):
	"""
	Class for implementing the Python mapping protocol on a list of table
	rows when multiple rows share the same key.
	"""
	def __init__(self, table_elem):
		"""
		Initialize the mapping on the list of rows.
		"""
		self.table = table_elem

	def __getitem__(self, key):
		"""
		Retrieve rows by key.
		"""
		l = [row for row in self.table if row._has_key(key)]
		if not len(l):
			raise KeyError, repr(key)
		return l

	def __setitem__(self, key, values):
		"""
		Replace the rows having key with the rows in the list values,
		appending to the table if there are no rows with that ID.
		"""
		# FIXME: should we call _set_key() on value to force it to have
		# the key that was searched for?
		del self[key]
		map(self.table.append, params)

	def __delitem__(self, key):
		"""
		Delete rows by key.
		"""
		self.table.filterRows(lambda row: not row._has_key(key))

	def __iter__(self):
		"""
		Return an iterator over the keys.
		"""
		return iter(self.keys())

	iterkeys = __iter__

	def keys(self):
		"""
		Return a list of the keys.
		"""
		keys = []
		for row in self.table:
			key = row._get_key()
			if key not in keys:
				keys.append(key)
		return keys

	def iteritems(self):
		"""
		Return an iterator over (key, value) pairs.
		"""
		return LSCTableMultiItemIter(self)


class LSCTableUnique(table.Table):
	def __init__(self, attrs):
		table.Table.__init__(self, attrs)
		self.dict = LSCTableUniqueDict(self)

	def makeReference(self, elem):
		"""
		Convert ilwd:char strings into object references.
		"""
		pass

	def deReference(self):
		"""
		Convert object references into ilwd:char strings.
		"""
		pass


class LSCTableMulti(table.Table):
	def __init__(self, attrs):
		table.Table.__init__(self, attrs)
		self.dict = LSCTableMultiDict(self)

	def makeReference(self, elem):
		"""
		Convert ilwd:char strings into object references.
		"""
		pass

	def deReference(self):
		"""
		Convert object references into ilwd:char strings.
		"""
		pass


# We don't subclass table.TableRow because that defeats the __slots__
# feature.
class LSCTableRow(object):
	__slots__ = []

	# Prefix with underscores to avoid collision with column names
	def _get_key(self):
		"""
		Get the unique ID for this row.
		"""
		raise KeyError, "row object does not define a key column"

	def _set_key(self, key):
		"""
		Set the unique ID for this row.
		"""
		raise KeyError, "row object does not define a key column"

	def _has_key(self, key):
		"""
		Check if this row's unique ID is equal to key.
		"""
		raise KeyError, "row object does not define a key column"


#
# =============================================================================
#
#                                process:table
#
# =============================================================================
#

class ProcessTable(LSCTableUnique):
	tableName = "process:table"
	validcolumns = {
		"program": "lstring",
		"version": "lstring",
		"cvs_repository": "lstring",
		"cvs_entry_time": "int_4s",
		"comment": "lstring",
		"is_online": "int_4s",
		"node": "lstring",
		"username": "lstring",
		"unix_procid": "int_4s",
		"start_time": "int_4s",
		"end_time": "int_4s",
		"jobid": "int_4s",
		"domain": "lstring",
		"ifos": "lstring",
		"process_id": "ilwd:char"
	}

class Process(LSCTableRow):
	__slots__ = ProcessTable.validcolumns.keys()

	def _get_key(self):
		return self.process_id

	def _set_key(self, key):
		self.process_id = key

	def _has_key(self, key):
		return key == self.process_id

	def cmp(self, other):
		# FIXME: this is a hack, but I need something so I can move
		# forward.
		for key in ProcessTable.validcolumns.keys():
			if key == "process_id":
				continue
			result = cmp(getattr(self, key), getattr(other, key))
			if result:
				return result
		return 0

ProcessTable.RowType = Process

class ProcessIDs(ILWD):
	def __init__(self, n = 0):
		ILWD.__init__(self, "process", "process_id", n)

#
# =============================================================================
#
#                                lfn:table
#
# =============================================================================
#

class LfnTable(LSCTableUnique):
	tableName = "lfn:table"
	validcolumns = {
		"process_id": "ilwd:char",
		"lfn_id": "ilwd:char",
		"name": "lstring",
		"comment": "lstring",
		"start_time": "int_4s",
		"end_time": "int_4s",
	}

	def makeReference(self, elem):
		"""
		Convert ilwd:char strings into object references.
		"""
		tables = table.getTablesByName(elem, ProcessTable.tableName)
		for row in self:
			row.process_id = FindILWD(tables, row.process_id)

	def deReference(self):
		"""
		Resolve object references back to ilwd:char strings.
		"""
		for row in self:
			row.process_id = row.process_id._get_key()

class Lfn(LSCTableRow):
	__slots__ = LfnTable.validcolumns.keys()

	def _get_key(self):
		return self.lfn_id

	def _set_key(self, key):
		self.lfn_id = key

	def _has_key(self, key):
		return key == self.lfn_id

	def cmp(self, other):
		# FIXME: this is a hack, but I need something so I can move
		# forward.
		for key in LfnTable.validcolumns.keys():
			if key == "lfn_id":
				continue
			result = cmp(getattr(self, key), getattr(other, key))
			if result:
				return result
		return 0

LfnTable.RowType = Lfn

class LfnIDs(ILWD):
	def __init__(self, n = 0):
		ILWD.__init__(self, "lfn", "lfn_id", n)


#
# =============================================================================
#
#                             process_params:table
#
# =============================================================================
#

class ProcessParamsTable(LSCTableMulti):
	tableName = "process_params:table"
	validcolumns = {
		"program": "lstring",
		"process_id": "ilwd:char",
		"param": "lstring",
		"type": "lstring",
		"value": "lstring"
	}

	def append(self, row):
		if row.type not in types.Types:
			raise ligolw.ElementError, "ProcessParamsTable.append():  unrecognized type \"%s\"" % row.type
		LSCTableMulti.append(self, row)

	def get_program(self, key):
		"""
		Return the name of the program associated with process ID
		key.
		"""
		for row in self:
			if row.process_id == key:
				return row.program
		raise KeyError, repr(key)

	def set_program(self, key, value):
		"""
		Set the program for all entries with process ID key to
		value.
		"""
		for row in self:
			if row.process_id == key:
				row.program = value

	def makeReference(self, elem):
		"""
		Convert ilwd:char strings into object references.
		"""
		tables = table.getTablesByName(elem, ProcessTable.tableName)
		for row in self:
			row.process_id = FindILWD(tables, row.process_id)

	def deReference(self):
		"""
		Resolve object references back to ilwd:char strings.
		"""
		for row in self:
			row.process_id = row.process_id._get_key()

class ProcessParams(LSCTableRow):
	__slots__ = ProcessParamsTable.validcolumns.keys()

	def _get_key(self):
		return self.process_id

	def _set_key(self, key):
		self.process_id = key

	def _has_key(self, key):
		return self.process_id == key

	def cmp(self, other):
		# FIXME: this is a hack, but I need something so I can move
		# forward.
		for key in ProcessParamsTable.validcolumns.keys():
			if key == "process_id":
				continue
			result = cmp(getattr(self, key), getattr(other, key))
			if result:
				return result
		return 0

ProcessParamsTable.RowType = ProcessParams


#
# =============================================================================
#
#                             search_summary:table
#
# =============================================================================
#

class SearchSummaryTable(LSCTableMulti):
	tableName = "search_summary:table"
	validcolumns = {
		"process_id": "ilwd:char",
		"shared_object": "lstring",
		"lalwrapper_cvs_tag": "lstring",
		"lal_cvs_tag": "lstring",
		"comment": "lstring",
		"ifos": "lstring",
		"in_start_time": "int_4s",
		"in_start_time_ns": "int_4s",
		"in_end_time": "int_4s",
		"in_end_time_ns": "int_4s",
		"out_start_time": "int_4s",
		"out_start_time_ns": "int_4s",
		"out_end_time": "int_4s",
		"out_end_time_ns": "int_4s",
		"nevents": "int_4s",
		"nnodes": "int_4s"
	}

	def makeReference(self, elem):
		"""
		Convert ilwd:char strings into object references.
		"""
		tables = table.getTablesByName(elem, ProcessTable.tableName)
		for row in self:
			row.process_id = FindILWD(tables, row.process_id)

	def deReference(self):
		"""
		Resolve object references back to ilwd:char strings.
		"""
		for row in self:
			row.process_id = row.process_id._get_key()

	def get_inlist(self):
		"""
		Return a segmentlist object describing the times spanned by
		the input segments of all rows in the table.
		"""
		return segments.segmentlist([row.get_in() for row in self])

	def get_outlist(self):
		"""
		Return a segmentlist object describing the times spanned by
		the output segments of all rows in the table.
		"""
		return segments.segmentlist([row.get_out() for row in self])

	def get_inprocs(self, seglist):
		"""
		Return a list of the process IDs for the processes whose
		input segments intersect some part of the segmentlist
		seglist.
		"""
		return [row.process_id for row in self if segments.segmentlist([row.get_in()]) & seglist]

	def get_outprocs(self, seglist):
		"""
		Return a list of the process IDs for the processes whose
		output segments intersect some part of the segmentlist
		seglist.
		"""
		return [row.process_id for row in self if segments.segmentlist([row.get_out()]) & seglist]

class SearchSummary(LSCTableRow):
	__slots__ = SearchSummaryTable.validcolumns.keys()

	def _get_key(self):
		return self.process_id

	def _set_key(self, key):
		self.process_id = key

	def _has_key(self, key):
		return self.process_id == key

	def cmp(self, other):
		# FIXME: this is a hack, but I need something so I can move
		# forward.
		for key in SearchSummaryTable.validcolumns.keys():
			if key == "process_id":
				continue
			result = cmp(getattr(self, key), getattr(other, key))
			if result:
				return result
		return 0

	def get_in(self):
		"""
		Return the input segment.
		"""
		return segments.segment(lal.LIGOTimeGPS(self.in_start_time, self.in_start_time_ns), lal.LIGOTimeGPS(self.in_end_time, self.in_end_time_ns))

	def set_in(self, seg):
		"""
		Set the input segment.
		"""
		self.in_start_time, self.in_start_time_ns = seg[0].seconds, seg[0].nanoseconds
		self.in_end_time, self.in_end_time_ns = seg[1].seconds, seg[1].nanoseconds

	def get_out(self):
		"""
		Get the output segment.
		"""
		return segments.segment(lal.LIGOTimeGPS(self.out_start_time, self.out_start_time_ns), lal.LIGOTimeGPS(self.out_end_time, self.out_end_time_ns))

	def set_out(self, seg):
		"""
		Set the output segment.
		"""
		self.out_start_time, self.out_start_time_ns = seg[0].seconds, seg[0].nanoseconds
		self.out_end_time, self.out_end_time_ns = seg[1].seconds, seg[1].nanoseconds

SearchSummaryTable.RowType = SearchSummary


#
# =============================================================================
#
#                            search_summvars:table
#
# =============================================================================
#

class SearchSummVarsTable(LSCTableMulti):
	tableName = "search_summvars:table"
	validcolumns = {
		"process_id": "ilwd:char",
		"name": "lstring",
		"string": "lstring",
		"value": "real_8"
	}

	def makeReference(self, elem):
		"""
		Convert ilwd:char strings into object references.
		"""
		tables = table.getTablesByName(elem, ProcessTable.tableName)
		for row in self:
			row.process_id = FindILWD(tables, row.process_id)

	def deReference(self):
		"""
		Resolve object references back to ilwd:char strings.
		"""
		for row in self:
			row.process_id = row.process_id._get_key()

class SearchSummVars(LSCTableRow):
	__slots__ = SearchSummVarsTable.validcolumns.keys()

	def _get_key(self):
		return self.process_id

	def _set_key(self, key):
		self.process_id = key

	def _has_key(self, key):
		return self.process_id == key

SearchSummVarsTable.RowType = SearchSummVars


#
# =============================================================================
#
#                               sngl_burst:table
#
# =============================================================================
#

class SnglBurstTable(LSCTableUnique):
	tableName = "sngl_burst:table"
	validcolumns = {
		"creator_db": "int_4s",
		"process_id": "ilwd:char",
		"filter_id": "ilwd:char",
		"ifo": "lstring",
		"search": "lstring",
		"channel": "lstring",
		"start_time": "int_4s",
		"start_time_ns": "int_4s",
		"stop_time": "int_4s",
		"stop_time_ns": "int_4s",
		"duration": "real_4",
		"flow": "real_4",
		"fhigh": "real_4",
		"central_freq": "real_4",
		"bandwidth": "real_4",
		"amplitude": "real_4",
		"snr": "real_4",
		"confidence": "real_4",
		"tfvolume": "real_4",
		"hrss": "real_4",
		"time_lag": "real_4",
		"peak_time": "int_4s",
		"peak_time_ns": "int_4s",
		"peak_frequency": "real_4",
		"peak_strain": "real_4",
		"peak_time_error": "real_4",
		"peak_frequency_error": "real_4",
		"peak_strain_error": "real_4",
		"ms_start_time": "int_4s",
		"ms_start_time_ns": "int_4s",
		"ms_stop_time": "int_4s",
		"ms_stop_time_ns": "int_4s",
		"ms_duration": "real_4",
		"ms_flow": "real_4",
		"ms_fhigh": "real_4",
		"ms_bandwidth": "real_4",
		"ms_hrss": "real_4",
		"ms_snr": "real_4",
		"ms_confidence": "real_4",
		"param_one_name": "lstring",
		"param_one_value": "real_8",
		"param_two_name": "lstring",
		"param_two_value": "real_8",
		"param_three_name": "lstring",
		"param_three_value": "real_8",
		"event_id": "ilwd:char"
	}

	def makeReference(self, elem):
		"""
		Convert ilwd:char strings into object references.
		"""
		tables = table.getTablesByName(elem, ProcessTable.tableName)
		for row in self:
			row.process_id = FindILWD(tables, row.process_id)

	def deReference(self):
		"""
		Resolve object references back to ilwd:char strings.
		"""
		for row in self:
			row.process_id = row.process_id._get_key()

class SnglBurst(LSCTableRow):
	__slots__ = SnglBurstTable.validcolumns.keys()

	def _get_key(self):
		return self.event_id

	def _set_key(self, key):
		self.event_id = key

	def _has_key(self, key):
		return self.event_id == key

	def get_start(self):
		return lal.LIGOTimeGPS(self.start_time, self.start_time_ns)

	def set_start(self, gps):
		self.start_time, self.start_time_ns = gps.seconds, gps.nanoseconds

	def get_stop(self):
		return lal.LIGOTimeGPS(self.stop_time, self.stop_time_ns)

	def set_stop(self, gps):
		self.stop_time, self.stop_time_ns = gps.seconds, gps.nanoseconds

	def get_peak(self):
		return lal.LIGOTimeGPS(self.peak_time, self.peak_time_ns)

	def set_peak(self, gps):
		self.peak_time, self.peak_time_ns = gps.seconds, gps.nanoseconds

	def get_period(self):
		start = lal.LIGOTimeGPS(self.start_time, self.start_time_ns)
		return segments.segment(start, start + self.duration)

	def set_period(self, period):
		# FIXME: should duration be forced to type float?
		self.start_time, self.start_time_ns = period[0].seconds, period[0].nanoseconds
		self.duration = float(period.duration())

	def get_band(self):
		return segments.segment(self.central_freq - self.bandwidth/2.0, self.central_freq + self.bandwidth/2.0)

	def set_band(self, band):
		self.central_freq = (band[0] + band[1])/2.0
		self.bandwidth = band.duration()

SnglBurstTable.RowType = SnglBurst

class SnglBurstIDs(ILWD):
	def __init__(self, n = 0):
		ILWD.__init__(self, "sngl_burst", "event_id", n)


#
# =============================================================================
#
#                             sngl_inspiral:table
#
# =============================================================================
#

class SnglInspiralTable(LSCTableUnique):
	tableName = "sngl_inspiral:table"
	validcolumns = {
		"process_id": "ilwd:char",
		"ifo": "lstring",
		"search": "lstring",
		"channel": "lstring",
		"end_time": "int_4s",
		"end_time_ns": "int_4s",
		"end_time_gmst": "real_8",
		"impulse_time": "int_4s",
		"impulse_time_ns": "int_4s",
		"template_duration": "real_8",
		"event_duration": "real_8",
		"amplitude": "real_4",
		"eff_distance": "real_4",
		"coa_phase": "real_4",
		"mass1": "real_4",
		"mass2": "real_4",
		"mchirp": "real_4",
		"mtotal": "real_4",
		"eta": "real_4",
		"tau0": "real_4",
		"tau2": "real_4",
		"tau3": "real_4",
		"tau4": "real_4",
		"tau5": "real_4",
		"ttotal": "real_4",
		"psi0": "real_4",
		"psi3": "real_4",
		"alpha": "real_4",
		"alpha1": "real_4",
		"alpha2": "real_4",
		"alpha3": "real_4",
		"alpha4": "real_4",
		"alpha5": "real_4",
		"alpha6": "real_4",
		"beta": "real_4",
		"f_final": "real_4",
		"snr": "real_4",
		"chisq": "real_4",
		"chisq_dof": "int_4s",
		"sigmasq": "real_8",
		"rsqveto_duration": "real_4",
		"Gamma0": "real_4",
		"Gamma1": "real_4",
		"Gamma2": "real_4",
		"Gamma3": "real_4",
		"Gamma4": "real_4",
		"Gamma5": "real_4",
		"Gamma6": "real_4",
		"Gamma7": "real_4",
		"Gamma8": "real_4",
		"Gamma9": "real_4",
		"event_id": "int_8s"
	}

	def makeReference(self, elem):
		"""
		Convert ilwd:char strings into object references.
		"""
		tables = table.getTablesByName(elem, ProcessTable.tableName)
		for row in self:
			row.process_id = FindILWD(tables, row.process_id)

	def deReference(self):
		"""
		Resolve object references back to ilwd:char strings.
		"""
		for row in self:
			row.process_id = row.process_id._get_key()

	def get_column(self,column):
		if column == 'effective_snr':
			return self.get_effective_snr()
		if column == 'snr_over_chi':
			return self.get_snr_over_chi()
		elif column == 'chirp_distance':
			return self.get_chirp_dist()
		else:
			return self.getColumnByName(column).asarray()

	def get_effective_snr(self):    
		snr = self.get_column('snr')
		chisq = self.get_column('chisq')
		chisq_dof = self.get_column('chisq_dof')
		return snr/ (1 + snr**2/250)**(0.25) / (chisq/(2*chisq_dof - 2) )**(0.25)
    
	def get_chirp_distance(self,ref_mass = 1.40):
		mchirp = self.get_column('mchirp')
		eff_dist = self.get_column('eff_distance')
		return eff_dist * (2.**(-1./5) * ref_mass / mchirp)**(5./6)

	def get_snr_over_chi(self):
		return self.get_column('snr')/self.get_column('chisq')**(1./2)
		
	def ifocut(self,ifo):
		ifoTrigs = table.new_from_template(self)
		for row in self:
			if row.ifo == ifo:
				ifoTrigs.append(row)
		return ifoTrigs

	def veto(self,seglist):
		vetoed = table.new_from_template(self)
		keep = table.new_from_template(self)
		for row in self:
			time = row.get_end()
			if time in seglist:
				vetoed.append(event)
			else:
				keep.append(event)

		return keep
	
	def getslide(self,slide_num):
		"""
		Return the triggers with a specific slide number.
		@param slide_num: the slide number to recover (contained in the event_id)
		"""
		slideTrigs = table.new_from_template(self)
		for row in self:
			if ( (row.event_id % 1000000000) / 100000 ) == slide_num:
				slideTrigs.append(row)
     
		return slideTrigs
	
		

class SnglInspiral(LSCTableRow):
	__slots__ = SnglInspiralTable.validcolumns.keys()

	def _get_key(self):
		return self.event_id

	def _set_key(self, key):
		self.event_id = key

	def _has_key(self, key):
		return self.event_id == key

	def get_end(self):
		return lal.LIGOTimeGPS(self.end_time, self.end_time_ns)

	def set_end(self, gps):
		self.end_time, self.end_time_ns = gps.seconds, gps.nanoseconds

	def get_effective_snr(self):
		return self.snr/ (1 + self.snr**2/250)**(0.25)/(self.chisq/(2*self.chisq_dof - 2) )**(0.25) 


SnglInspiralTable.RowType = SnglInspiral

# FIXME: class definition removed until LAL inspiral code generates ilwd:char
# event_ids.  Re-enable when LAL is fixed.
#
#class SnglInspiralIDs(ILWD):
#	def __init__(self, n = 0):
#		ILWD.__init__(self, "sngl_inspiral", "event_id", n)


#
# =============================================================================
#
#                             sngl_ringdown:table
#
# =============================================================================
#

class SnglRingDownTable(LSCTableUnique):
	tableName = "sngl_ringdown:table"
	validcolumns = {
		"process_id": "ilwd:char",
		"ifo": "lstring",
		"search": "lstring",
		"channel": "lstring",
		"start_time": "int_4s",
		"start_time_ns": "int_4s",
		"start_time_gmst": "real_8",
		"frequency": "real_4",
		"quality": "real_4",
		"mass": "real_4",
		"spin": "real_4",
		"snr": "real_4",
		"eff_distance": "real_4",
		"sigma_sq": "real_8",
		"event_id": "int_8s"
	}

	def makeReference(self, elem):
		"""
		Convert ilwd:char strings into object references.
		"""
		tables = table.getTablesByName(elem, ProcessTable.tableName)
		for row in self:
			row.process_id = FindILWD(tables, row.process_id)

	def deReference(self):
		"""
		Resolve object references back to ilwd:char strings.
		"""
		for row in self:
			row.process_id = row.process_id._get_key()

class SnglRingDown(LSCTableRow):
	__slots__ = SnglRingDownTable.validcolumns.keys()

	def _get_key(self):
		return self.event_id

	def _set_key(self, key):
		self.event_id = key

	def _has_key(self, key):
		return self.event_id == key

SnglRingDownTable.RowType = SnglRingDown

class SnglRingDownIDs(ILWD):
	def __init__(self, n = 0):
		ILWD.__init__(self, "sngl_ringdown", "event_id", n)


#
# =============================================================================
#
#                             multi_inspiral:table
#
# =============================================================================
#

class MultiInspiralTable(LSCTableMulti):
	tableName = "multi_inspiral:table"
	validcolumns = {
		"process_id": "ilwd:char",
		"ifos": "lstring",
		"search": "lstring",
		"end_time": "int_4s",
		"end_time_ns": "int_4s",
		"end_time_gmst": "real_8",
		"impulse_time": "int_4s",
		"impulse_time_ns": "int_4s",
		"amplitude": "real_4",
		"ifo1_eff_distance": "real_4",
		"ifo2_eff_distance": "real_4",
		"eff_distance": "real_4",
		"coa_phase": "real_4",
		"mass1": "real_4",
		"mass2": "real_4",
		"mchirp": "real_4",
		"eta": "real_4",
		"tau0": "real_4",
		"tau2": "real_4",
		"tau3": "real_4",
		"tau4": "real_4",
		"tau5": "real_4",
		"ttotal": "real_4",
		"ifo1_snr": "real_4",
		"ifo2_snr": "real_4",
		"snr": "real_4",
		"chisq": "real_4",
		"chisq_dof": "real_4",
		"sigmasq": "real_4",
		"ligo_axis_ra": "real_4",
		"ligo_axis_dec": "real_4",
		"ligo_angle": "real_4",
		"ligo_angle_sig": "real_4",
		"inclination": "real_4",
		"polarization": "real_4",
		"event_id": "ilwd:char"
	}

	def makeReference(self, elem):
		"""
		Convert ilwd:char strings into object references.
		"""
		tables = table.getTablesByName(elem, ProcessTable.tableName)
		for row in self:
			row.process_id = FindILWD(tables, row.process_id)

	def deReference(self):
		"""
		Resolve object references back to ilwd:char strings.
		"""
		for row in self:
			row.process_id = row.process_id._get_key()

class MultiInspiral(LSCTableRow):
	__slots__ = MultiInspiralTable.validcolumns.keys()

	def _get_key(self):
		return self.process_id

	def _set_key(self, key):
		self.process_id = key

	def _has_key(self, key):
		return self.process_id == key

MultiInspiralTable.RowType = MultiInspiral


#
# =============================================================================
#
#                              sim_inspiral:table
#
# =============================================================================
#

class SimInspiralTable(LSCTableUnique):
	tableName = "sim_inspiral:table"
	validcolumns = {
		"process_id": "ilwd:char",
		"waveform": "lstring",
		"geocent_end_time": "int_4s",
		"geocent_end_time_ns": "int_4s",
		"h_end_time": "int_4s",
		"h_end_time_ns": "int_4s",
		"l_end_time": "int_4s",
		"l_end_time_ns": "int_4s",
		"g_end_time": "int_4s",
		"g_end_time_ns": "int_4s",
		"t_end_time": "int_4s",
		"t_end_time_ns": "int_4s",
		"v_end_time": "int_4s",
		"v_end_time_ns": "int_4s",
		"end_time_gmst": "real_8",
		"source": "lstring",
		"mass1": "real_4",
		"mass2": "real_4",
		"mchirp": "real_4",
		"eta": "real_4",
		"distance": "real_4",
		"longitude": "real_4",
		"latitude": "real_4",
		"inclination": "real_4",
		"coa_phase": "real_4",
		"polarization": "real_4",
		"psi0": "real_4",
		"psi3": "real_4",
		"alpha": "real_4",
		"alpha1": "real_4",
		"alpha2": "real_4",
		"alpha3": "real_4",
		"alpha4": "real_4",
		"alpha5": "real_4",
		"alpha6": "real_4",
		"beta": "real_4",
		"spin1x": "real_4",
		"spin1y": "real_4",
		"spin1z": "real_4",
		"spin2x": "real_4",
		"spin2y": "real_4",
		"spin2z": "real_4",
		"theta0": "real_4",
		"phi0": "real_4",
		"f_lower": "real_4",
		"f_final": "real_4",
		"eff_dist_h": "real_4",
		"eff_dist_l": "real_4",
		"eff_dist_g": "real_4",
		"eff_dist_t": "real_4",
		"eff_dist_v": "real_4",
		"simulation_id": "ilwd:char"
	}

	def makeReference(self, elem):
		"""
		Convert ilwd:char strings into object references.
		"""
		tables = table.getTablesByName(elem, ProcessTable.tableName)
		for row in self:
			row.process_id = FindILWD(tables, row.process_id)

	def deReference(self):
		"""
		Resolve object references back to ilwd:char strings.
		"""
		for row in self:
			row.process_id = row.process_id._get_key()

	def get_column(self,column):
		if 'chirp_dist' in column:
			site = column[-1]
			return self.get_chirp_dist(site)
		elif column == 'spin1':
			return self.get_spin_mag(1)
		elif column == 'spin2':
			return self.get_spin_mag(2)
		else:
			return self.getColumnByName(column).asarray()

	def get_chirp_dist(self,site,ref_mass = 1.40):
		mchirp = self.get_column('mchirp')
		eff_dist = self.get_column('eff_dist_' + site)
		return eff_dist * (2.**(-1./5) * ref_mass / mchirp)**(5./6)

	def get_spin_mag(self,objectnumber):
		sx = self.get_column('spin' + str(objectnumber) + 'x')
		sy = self.get_column('spin' + str(objectnumber) + 'y')
		sz = self.get_column('spin' + str(objectnumber) + 'z')
		return (sx**2 + sy**2 + sz**2)**(0.5)

	def veto(self,seglist,site=None):
		keep = table.new_from_template(self)
		for row in self:
			time = row.get_end(site)
			if time not in seglist:
				keep.append(row)

		return keep

class SimInspiral(LSCTableRow):
	__slots__ = SimInspiralTable.validcolumns.keys()

	def _get_key(self):
		return self.simulation_id

	def _set_key(self, key):
		self.simulation_id = key

	def _has_key(self, key):
		return self.simulation_id == key

	def get_end(self,site = None):
		if not site:
			return lal.LIGOTimeGPS(self.geocent_end_time, self.geocent_end_time_ns)
		else:
			return lal.LIGOTimeGPS(getattr(self,site + 'end_time'), getattr(self,site + 'end_time_ns'))

SimInspiralTable.RowType = SimInspiral

class SimInspiralIDs(ILWD):
	def __init__(self, n = 0):
		ILWD.__init__(self, "sim_inspiral", "simulation_id", n)


#
# =============================================================================
#
#                               sim_burst:table
#
# =============================================================================
#

class SimBurstTable(LSCTableUnique):
	tableName = "sim_burst:table"
	validcolumns = {
		"process_id": "ilwd:char",
		"waveform": "lstring",
		"geocent_peak_time": "int_4s",
		"geocent_peak_time_ns": "int_4s",
		"h_peak_time": "int_4s",
		"h_peak_time_ns": "int_4s",
		"l_peak_time": "int_4s",
		"l_peak_time_ns": "int_4s",
		"peak_time_gmst": "real_8",
		"dtminus": "real_4",
		"dtplus": "real_4",
		"longitude": "real_4",
		"latitude": "real_4",
		"coordinates": "lstring",
		"polarization": "real_4",
		"hrss": "real_4",
		"hpeak": "real_4",
		"distance": "real_4",
		"freq": "real_4",
		"tau": "real_4",
		"zm_number": "int_4s",
		"simulation_id": "ilwd:char"
	}

	def makeReference(self, elem):
		"""
		Convert ilwd:char strings into object references.
		"""
		tables = table.getTablesByName(elem, ProcessTable.tableName)
		for row in self:
			row.process_id = FindILWD(tables, row.process_id)

	def deReference(self):
		"""
		Resolve object references back to ilwd:char strings.
		"""
		for row in self:
			row.process_id = row.process_id._get_key()

class SimBurst(LSCTableRow):
	__slots__ = SimBurstTable.validcolumns.keys()

	def _get_key(self):
		return self.simulation_id

	def _set_key(self, key):
		self.simulation_id = key

	def _has_key(self, key):
		return self.simulation_id == key

	def cmp(self, other):
		"""
		Return 0 if self and other describe the same injection,
		non-0 otherwise.
		"""
		a = (
			self.waveform,
			self.geocent_peak_time,
			self.geocent_peak_time_ns,
			self.h_peak_time,
			self.h_peak_time_ns,
			self.l_peak_time,
			self.l_peak_time_ns,
			self.peak_time_gmst,
			self.dtminus,
			self.dtplus,
			self.longitude,
			self.latitude,
			self.coordinates,
			self.polarization,
			self.hrss,
			self.hpeak,
			self.distance,
			self.freq,
			self.tau,
			self.zm_number
		)
		b = (
			other.waveform,
			other.geocent_peak_time,
			other.geocent_peak_time_ns,
			other.h_peak_time,
			other.h_peak_time_ns,
			other.l_peak_time,
			other.l_peak_time_ns,
			other.peak_time_gmst,
			other.dtminus,
			other.dtplus,
			other.longitude,
			other.latitude,
			other.coordinates,
			other.polarization,
			other.hrss,
			other.hpeak,
			other.distance,
			other.freq,
			other.tau,
			other.zm_number
		)
		return cmp(a, b)

	def get_geocent_peak(self):
		return lal.LIGOTimeGPS(self.geocent_peak_time, self.geocent_peak_time_ns)

	def set_geocent_peak(self, gps):
		self.geocent_peak_time, self.geocent_peak_time_ns = gps.seconds, gps.nanoseconds

	def get_peak(self, instrument):
		observatory = instrument[0]
		if observatory == "H":
			return lal.LIGOTimeGPS(self.h_peak_time, self.h_peak_time_ns)
		if observatory == "L":
			return lal.LIGOTimeGPS(self.l_peak_time, self.l_peak_time_ns)
		raise ValueError, instrument

	def set_peak(self, instrument, gps):
		observatory = instrument[0]
		if observatory == "H":
			self.h_peak_time, self.h_peak_time_ns = gps.seconds, gps.nanoseconds
		if observatory == "L":
			self.l_peak_time, self.l_peak_time_ns = gps.seconds, gps.nanoseconds
		raise ValueError, instrument

SimBurstTable.RowType = SimBurst

class SimBurstIDs(ILWD):
	def __init__(self, n = 0):
		ILWD.__init__(self, "sim_burst", "simulation_id", n)


#
# =============================================================================
#
#                              sim_ringdown:table
#
# =============================================================================
#

class SimRingDownTable(LSCTableUnique):
	tableName = "sim_ringdown:table"
	validcolumns = {
		"process_id": "ilwd:char",
		"waveform": "lstring",
		"coordinates": "lstring",
		"geocent_start_time": "int_4s",
		"geocent_start_time_ns": "int_4s",
		"h_start_time": "int_4s",
		"h_start_time_ns": "int_4s",
		"l_start_time": "int_4s",
		"l_start_time_ns": "int_4s",
		"start_time_gmst": "real_8",
		"mass": "real_4",
		"longitude": "real_4",
		"latitude": "real_4",
		"distance": "real_4",
		"inclination": "real_4",
		"polarization": "real_4",
		"epsilon": "real_4",
		"spin": "real_4",
		"frequency": "real_4",
		"quality": "real_4",
		"eff_dist_h": "real_4",
		"eff_dist_l": "real_4",
		"h0": "real_4",
		"hrss": "real_4",
		"hrss_h": "real_4",
		"hrss_l": "real_4",
		"simulation_id": "ilwd:char"
	}

	def makeReference(self, elem):
		"""
		Convert ilwd:char strings into object references.
		"""
		tables = table.getTablesByName(elem, ProcessTable.tableName)
		for row in self:
			row.process_id = FindILWD(tables, row.process_id)

	def deReference(self):
		"""
		Resolve object references back to ilwd:char strings.
		"""
		for row in self:
			row.process_id = row.process_id._get_key()

class SimRingDown(LSCTableRow):
	__slots__ = SimRingDownTable.validcolumns.keys()

	def _get_key(self):
		return self.simulation_id

	def _set_key(self, key):
		self.simulation_id = key

	def _has_key(self, key):
		return self.simulation_id == key

SimRingDownTable.RowType = SimRingDown

class SimRingDownIDs(ILWD):
	def __init__(self, n = 0):
		ILWD.__init__(self, "sim_ringdown", "simulation_id", n)


#
# =============================================================================
#
#                               summ_value:table
#
# =============================================================================
#

class SummValueTable(LSCTableMulti):
	tableName = "summ_value:table"
	validcolumns = {
		"program": "lstring",
		"process_id": "ilwd:char",
		"start_time": "int_4s",
		"start_time_ns": "int_4s",
		"end_time": "int_4s",
		"end_time_ns": "int_4s",
		"ifo": "lstring",
		"name": "lstring",
		"value": "real_4",
		"comment": "lstring"
	}

	def makeReference(self, elem):
		"""
		Convert ilwd:char strings into object references.
		"""
		tables = table.getTablesByName(elem, ProcessTable.tableName)
		for row in self:
			row.process_id = FindILWD(tables, row.process_id)

	def deReference(self):
		"""
		Resolve object references back to ilwd:char strings.
		"""
		for row in self:
			row.process_id = row.process_id._get_key()

class SummValue(LSCTableRow):
	__slots__ = SummValueTable.validcolumns.keys()

	def _get_key(self):
		return self.process_id

	def _set_key(self, key):
		self.process_id = key

	def _has_key(self, key):
		return self.process_id == key

SummValueTable.RowType = SummValue


#
# =============================================================================
#
#                            sim_inst_params:table
#
# =============================================================================
#

class SimInstParamsTable(LSCTableMulti):
	tableName = "sim_inst_params:table"
	validcolumns = {
		"simulation_id": "ilwd:char",
		"name": "lstring",
		"comment": "lstring",
		"value": "real_8"
	}

class SimInstParams(LSCTableRow):
	__slots__ = SimInstParamsTable.validcolumns.keys()

	def _get_key(self):
		return self.simulation_id

	def _set_key(self, key):
		self.simulation_id = key

	def _has_key(self, key):
		return self.simulation_id == key

SimInstParamsTable.RowType = SimInstParams

class SimInstParamsIDs(ILWD):
	def __init__(self, n = 0):
		ILWD.__init__(self, "sim_inst_params", "simulation_id", n)


#
# =============================================================================
#
#                               stochastic:table
#
# =============================================================================
#

class StochasticTable(LSCTableMulti):
	tableName = "stochastic:table"
	validcolumns = {
		"process_id": "ilwd:char",
		"ifo_one": "lstring",
		"ifo_two": "lstring",
		"channel_one": "lstring",
		"channel_two": "lstring",
		"start_time": "int_4s",
		"start_time_ns": "int_4s",
		"duration": "int_4s",
		"duration_ns": "int_4s",
		"f_min": "real_8",
		"f_max": "real_8",
		"cc_stat": "real_8",
		"cc_sigma": "real_8"
	}

	def makeReference(self, elem):
		"""
		Convert ilwd:char strings into object references.
		"""
		tables = table.getTablesByName(elem, ProcessTable.tableName)
		for row in self:
			row.process_id = FindILWD(tables, row.process_id)

	def deReference(self):
		"""
		Resolve object references back to ilwd:char strings.
		"""
		for row in self:
			row.process_id = row.process_id._get_key()

class Stochastic(LSCTableRow):
	__slots__ = StochasticTable.validcolumns.keys()

	def _get_key(self):
		return self.process_id

	def _set_key(self, key):
		self.process_id = key

	def _has_key(self, key):
		return self.process_id == key

StochasticTable.RowType = Stochastic


#
# =============================================================================
#
#                               stochsumm:table
#
# =============================================================================
#

class StochSummTable(LSCTableMulti):
	tableName = "stochsumm:table"
	validcolumns = {
		"process_id": "ilwd:char",
		"ifo_one": "lstring",
		"ifo_two": "lstring",
		"channel_one": "lstring",
		"channel_two": "lstring",
		"start_time": "int_4s",
		"start_time_ns": "int_4s",
		"end_time": "int_4s",
		"end_time_ns": "int_4s",
		"f_min": "real_8",
		"f_max": "real_8",
		"y_opt": "real_8",
		"error": "real_8"
	}

	def makeReference(self, elem):
		"""
		Convert ilwd:char strings into object references.
		"""
		tables = table.getTablesByName(elem, ProcessTable.tableName)
		for row in self:
			row.process_id = FindILWD(tables, row.process_id)

	def deReference(self):
		"""
		Resolve object references back to ilwd:char strings.
		"""
		for row in self:
			row.process_id = row.process_id._get_key()

class StochSumm(LSCTableRow):
	__slots__ = StochSummTable.validcolumns.keys()

	def _get_key(self):
		return self.process_id

	def _set_key(self, key):
		self.process_id = key

	def _has_key(self, key):
		return self.process_id == key

StochSummTable.RowType = StochSumm


#
# =============================================================================
#
#                            external_trigger:table
#
# =============================================================================
#

# FIXME: this table looks broken to me.  There is no unique ID column, thus
# it is not possible to refer to entries in this table from other tables.
# There is a "notice_id" column, but that is for recording the native
# identifier as used by the source of the trigger.  It cannot be relied
# upon to be unique within this table (two different sources might *happen*
# to use the same identifier format, like "event001").

class ExtTriggersTable(LSCTableMulti):
	tableName = "external_trigger:table"
	validcolumns = {
		"process_id": "ilwd:char",
		"det_alts": "lstring",
		"det_band": "lstring",
		"det_fluence": "lstring",
		"det_fluence_int": "lstring",
		"det_name": "lstring",
		"det_peak": "lstring",
		"det_peak_int": "lstring",
		"det_snr": "lstring",
		"email_time": "int_4s",
		"event_dec": "real_4",
		"event_dec_err": "real_4",
		"event_epoch": "lstring",
		"event_err_type": "lstring",
		"event_ra": "real_4",
		"event_ra_err": "real_4",
		"start_time": "int_4s",
		"start_time_ns": "int_4s",
		"event_type": "lstring",
		"event_z": "real_4",
		"event_z_err": "real_4",
		"notice_comments": "lstring",
		"notice_id": "lstring",
		"notice_sequence": "lstring",
		"notice_time": "int_4s",
		"notice_type": "lstring",
		"notice_url": "lstring",
		"obs_fov_dec": "real_4",
		"obs_fov_dec_width": "real_4",
		"obs_fov_ra": "real_4",
		"obs_fov_ra_width": "real_4",
		"obs_loc_ele": "real_4",
		"obs_loc_lat": "real_4",
		"obs_loc_long": "real_4",
		"ligo_fave_lho": "real_4",
		"ligo_fave_llo": "real_4",
		"ligo_delay": "real_4",
		"event_number_gcn": "int_4s",
		"event_number_grb": "lstring",
		"event_status": "int_4s"
	}

class ExtTriggers(LSCTableRow):
	__slots__ = ExtTriggersTable.validcolumns.keys()

	def _get_key(self):
		return self.process_id

	def _set_key(self, key):
		self.process_id = key

	def _has_key(self, key):
		return self.process_id == key

ExtTriggersTable.RowType = ExtTriggers


#
# =============================================================================
#
#                                 filter:table
#
# =============================================================================
#

class FilterTable(LSCTableMulti):
	tableName = "filter:table"
	validcolumns = {
		"process_id": "ilwd:char",
		"program": "lstring",
		"start_time": "int_4s",
		"filter_name": "lstring",
		"comment": "lstring"
	}

	def makeReference(self, elem):
		"""
		Convert ilwd:char strings into object references.
		"""
		tables = table.getTablesByName(elem, ProcessTable.tableName)
		for row in self:
			row.process_id = FindILWD(tables, row.process_id)

	def deReference(self):
		"""
		Resolve object references back to ilwd:char strings.
		"""
		for row in self:
			row.process_id = row.process_id._get_key()

class Filter(LSCTableRow):
	__slots__ = FilterTable.validcolumns.keys()

	def _get_key(self):
		return self.process_id

	def _set_key(self, key):
		self.process_id = key

	def _has_key(self, key):
		return self.process_id == key

FilterTable.RowType = Filter


#
# =============================================================================
#
#                                segment:table
#
# =============================================================================
#

class SegmentTable(LSCTableUnique):
	tableName = "segment:table"
	validcolumns = {
		"creator_db": "int_4s",
		"process_id": "ilwd:char",
		"segment_id": "ilwd:char",
		"start_time": "int_4s",
		"start_time_ns": "int_4s",
		"end_time": "int_4s",
		"end_time_ns": "int_4s",
		"active": "int_4s",
		"segnum": "int_4s",
		"insertion_time": "int_4s"
	}

	def makeReference(self, elem):
		"""
		Convert ilwd:char strings into object references.
		"""
		tables = table.getTablesByName(elem, ProcessTable.tableName)
		for row in self:
			row.process_id = FindILWD(tables, row.process_id)

	def deReference(self):
		"""
		Resolve object references back to ilwd:char strings.
		"""
		for row in self:
			row.process_id = row.process_id._get_key()

class Segment(LSCTableRow):
	__slots__ = SegmentTable.validcolumns.keys()

	def _get_key(self):
		return self.segment_id

	def _set_key(self, key):
		self.segment_id = key

	def _has_key(self, key):
		return self.segment_id == key

	def get(self):
		"""
		Return the segment described by this row.
		"""
		return segments.segment(lal.LIGOTimeGPS(self.start_time, self.start_time_ns), lal.LIGOTimeGPS(self.end_time, self.end_time_ns))

	def set(self, segment):
		"""
		Set the segment described by this row.
		"""
		self.start_time, self.start_time_ns = segment[0].seconds, segment[0].nanoseconds
		self.end_time, self.end_time_ns = segment[1].seconds, segment[1].nanoseconds

	def get_active(self):
		"""
		Return True if the segment is active, False if the segment
		is inactive and None if neither is the case.
		"""
		if self.active > 0:
			return True
		if self.active < 0:
			return False
		return None

	def set_active(self, active):
		"""
		Sets the segment to active if active is True, to inactive
		if active if False, and undefined if active is None.
		"""
		if active == None:
			self.active = 0
		elif active:
			self.active = 1
		else:
			self.active = -1
		return self

SegmentTable.RowType = Segment

class SegmentIDs(ILWD):
	def __init__(self, n = 0):
		ILWD.__init__(self, "segment", "segment_id", n)


#
# =============================================================================
#
#                            segment_def_map:table
#
# =============================================================================
#

class SegmentDefMapTable(LSCTableUnique):
	tableName = "segment_def_map:table"
	validcolumns = {
		"creator_db": "int_4s",
		"process_id": "ilwd:char",
		"seg_def_map_id": "ilwd:char",
		"segment_cdb": "int_4s",
		"segment_id": "ilwd:char",
		"segment_def_cdb": "int_4s",
		"segment_def_id": "ilwd:char",
		"state_vec_map": "int_4s",
		"insertion_time": "int_4s"
	}

	def makeReference(self, elem):
		"""
		Convert ilwd:char strings into object references.
		"""
		proctabs = table.getTablesByName(elem, ProcessTable.tableName)
		segtabs = table.getTablesByName(elem, SegmentTable.tableName)
		deftabs = table.getTablesByName(elem, SegmentDefTable.tableName)
		for row in self:
			row.process_id = FindILWD(proctabs, row.process_id)
			row.segment_id = FindILWD(segtabs, row.segment_id)
			row.segment_def_id = FindILWD(deftabs, row.segment_def_id)

	def deReference(self):
		"""
		Resolve object references back to ilwd:char strings.
		"""
		for row in self:
			row.process_id = row.process_id._get_key()
			row.segment_id = row.segment_id._get_key()
			row.segment_def_id = row.segment_def_id._get_key()

class SegmentDefMap(LSCTableRow):
	__slots__ = SegmentDefMapTable.validcolumns.keys()

	def _get_key(self):
		return self.seg_def_map_id

	def _set_key(self, key):
		self.seg_def_map_id = key

	def _has_key(self, key):
		return self.seg_def_map_id == key

SegmentDefMapTable.RowType = SegmentDefMap

class SegmentDefMapIDs(ILWD):
	def __init__(self, n = 0):
		ILWD.__init__(self, "segment_def_map", "segment_def_id", n)


#
# =============================================================================
#
#                            segment_definer:table
#
# =============================================================================
#

class SegmentDefTable(LSCTableUnique):
	tableName = "segment_definer:table"
	validcolumns = {
		"creator_db": "int_4s",
		"process_id": "ilwd:char",
		"segment_def_id": "ilwd:char",
		"run": "lstring",
		"ifos": "lstring",
		"name": "lstring",
		"version": "int_4s",
		"comment": "lstring",
		"state_vec_major": "int_4s",
		"state_vec_minor": "int_4s",
		"insertion_time": "int_4s"
	}

	def makeReference(self, elem):
		"""
		Convert ilwd:char strings into object references.
		"""
		tables = table.getTablesByName(elem, ProcessTable.tableName)
		for row in self:
			row.process_id = FindILWD(tables, row.process_id)

	def deReference(self):
		"""
		Resolve object references back to ilwd:char strings.
		"""
		for row in self:
			row.process_id = row.process_id._get_key()

class SegmentDef(LSCTableRow):
	__slots__ = SegmentDefTable.validcolumns.keys()

	def _get_key(self):
		return self.segment_def_id

	def _set_key(self, key):
		self.segment_def_id = key

	def _has_key(self, key):
		return self.segment_def_id == key

SegmentDefTable.RowType = SegmentDef


#
# =============================================================================
#
#                               time_slide:table
#
# =============================================================================
#

class TimeSlideTable(LSCTableMulti):
	tableName = "time_slide:table"
	validcolumns = {
		"process_id": "ilwd:char",
		"time_slide_id": "ilwd:char",
		"instrument": "lstring",
		"offset": "real_8"
	}

	def makeReference(self, elem):
		"""
		Convert ilwd:char strings into object references.
		"""
		tables = table.getTablesByName(elem, ProcessTable.tableName)
		for row in self:
			row.process_id = FindILWD(tables, row.process_id)

	def deReference(self):
		"""
		Resolve object references back to ilwd:char strings.
		"""
		for row in self:
			row.process_id = row.process_id._get_key()

	def get_offset_dict(self, id):
		"""
		Return a dictionary of offsets indexed by instrument for
		the given time slide ID.
		"""
		d = {}
		for row in self.dict[id]:
			if d.has_key(row.instrument):
				raise KeyError, "%s: duplicate instrument %s" % (repr(id), row.instrument)
			d[row.instrument] = row.offset
		return d


class TimeSlide(LSCTableRow):
	__slots__ = TimeSlideTable.validcolumns.keys()

	def _get_key(self):
		return self.time_slide_id

	def _set_key(self, key):
		self.time_slide_id = key

	def _has_key(self, key):
		return self.time_slide_id == key

TimeSlideTable.RowType = TimeSlide

class TimeSlideIDs(ILWD):
	def __init__(self, n = 0):
		ILWD.__init__(self, "time_slide", "time_slide_id", n)


#
# =============================================================================
#
#                             coinc_definer:table
#
# =============================================================================
#

class CoincDefTable(LSCTableMulti):
	tableName = "coinc_definer:table"
	validcolumns = {
		"coinc_def_id": "ilwd:char",
		"table_name": "lstring"
	}

	def get_contributors(self, id):
		"""
		Return a list of contributing table names for the given ID.
		"""
		l = [row.table_name for row in self.dict[id]]
		l.sort()
		return l

class CoincDef(LSCTableRow):
	__slots__ = CoincDefTable.validcolumns.keys()

	def _get_key(self):
		return self.coinc_def_id

	def _set_key(self, key):
		self.coinc_def_id = key

	def _has_key(self, key):
		return self.coinc_def_id == key

CoincDefTable.RowType = CoincDef

class CoincDefIDs(ILWD):
	def __init__(self, n = 0):
		ILWD.__init__(self, "coinc_definer", "coinc_def_id", n)


#
# =============================================================================
#
#                              coinc_event:table
#
# =============================================================================
#

class CoincTable(LSCTableUnique):
	tableName = "coinc_event:table"
	validcolumns = {
		"process_id": "ilwd:char",
		"coinc_def_id": "ilwd:char",
		"coinc_event_id": "ilwd:char",
		"time_slide_id": "ilwd:char",
		"nevents": "int_4u"
	}

	def makeReference(self, elem):
		"""
		Convert ilwd:char strings into object references.
		"""
		proctab = table.getTablesByName(elem, ProcessTable.tableName)
		slidetab = table.getTablesByName(elem, TimeSlideTable.tableName)
		deftab = table.getTablesByName(elem, CoincDefTable.tableName)
		for row in self:
			row.process_id = FindILWD(proctab, row.process_id)
			row.time_slide_id = FindILWD(slidetab, row.time_slide_id)
			row.coinc_def_id = FindILWD(deftab, row.coinc_def_id)

	def deReference(self):
		"""
		Resolve object references back to ilwd:char strings.
		"""
		for row in self:
			row.process_id = row.process_id._get_key()
			row.time_slide_id = row.time_slide_id[0]._get_key()
			row.coinc_def_id = row.coinc_def_id[0]._get_key()

class Coinc(LSCTableRow):
	__slots__ = CoincTable.validcolumns.keys()

	def _get_key(self):
		return self.coinc_event_id

	def _set_key(self, key):
		self.coinc_event_id = key

	def _has_key(self, key):
		return self.coinc_event_id == key

CoincTable.RowType = Coinc

class CoincIDs(ILWD):
	def __init__(self, n = 0):
		ILWD.__init__(self, "coinc_event", "coinc_event_id", n)


#
# =============================================================================
#
#                            coinc_event_map:table
#
# =============================================================================
#


# Tables that can provide "events" for the coinc_event_map table
# Wow:  it's annoying this has to be done by hand.
#
# FIXME: sngl_inspiral table cannot participate in the coinc_event table
# infrastructure until the LAL code generates event_ids of type ilwd:char.
# Re-list SnglInspiralTable in here when LAL is fixed.
CoincEventMapSourceNames = [
	table.StripTableName(SnglBurstTable.tableName),
#	table.StripTableName(SnglInspiralTable.tableName),
	table.StripTableName(SnglRingDownTable.tableName),
	table.StripTableName(MultiInspiralTable.tableName),
	table.StripTableName(SimInspiralTable.tableName),
	table.StripTableName(SimBurstTable.tableName),
	table.StripTableName(SimRingDownTable.tableName),
	table.StripTableName(CoincTable.tableName)
]

class CoincMapTable(LSCTableUnique):
	tableName = "coinc_event_map:table"
	validcolumns = {
		"coinc_event_id": "ilwd:char",
		"event_id": "ilwd:char"
	}

	def makeReference(self, elem):
		"""
		Convert ilwd:char strings into object references.
		"""
		coinctab = table.getTablesByName(elem, CoincTable.tableName)
		eventtab = []
		for tablename in CoincEventMapSourceNames:
			eventtab.extend(table.getTablesByName(elem, tablename))
		for row in self:
			row.coinc_event_id = FindILWD(coinctab, row.coinc_event_id)
			row.event_id = FindILWD(eventtab, row.event_id)

	def deReference(self):
		"""
		Resolve object references back to ilwd:char strings.
		"""
		for row in self:
			row.coinc_event_id = row.coinc_event_id._get_key()
			row.event_id = row.event_id._get_key()

class CoincMap(LSCTableRow):
	__slots__ = CoincMapTable.validcolumns.keys()

CoincMapTable.RowType = CoincMap


#
# =============================================================================
#
#                               ligolw_mon:table
#
# =============================================================================
#

class LIGOLWMonTable(LSCTableUnique):
	tableName = "ligolw_mon:table"
	validcolumns = {
		"creator_db": "int_4s",
		"process_id": "ilwd:char",
		"time": "int_4s",
		"time_ns": "int_4s",
		"amplitude": "real_8",
		"confidence": "real_8",
		"frequency": "real_8",
		"event_id": "ilwd:char",
		"insertion_time": "int_4s"
	}

	def makeReference(self, elem):
		"""
		Convert ilwd:char strings into object references.
		"""
		tables = table.getTablesByName(elem, ProcessTable.tableName)
		for row in self:
			row.process_id = FindILWD(tables, row.process_id)

	def deReference(self):
		"""
		Resolve object references back to ilwd:char strings.
		"""
		for row in self:
			row.process_id = row.process_id._get_key()

class LIGOLWMon(LSCTableRow):
	__slots__ = LIGOLWMonTable.validcolumns.keys()

	def _get_key(self):
		return self.event_id

	def _set_key(self, key):
		self.event_id = key

	def _has_key(self, key):
		return self.event_id == key

	def get_time(self):
		return lal.LIGOTimeGPS(self.time, self.time_ns)

	def set_time(self, gps):
		self.time, self.time_ns = gps.seconds, gps.nanoseconds

LIGOLWMonTable.RowType = LIGOLWMon

class LIGOLWMonIDs(ILWD):
	def __init__(self, n = 0):
		ILWD.__init__(self, "ligolw_mon", "event_id", n)


#
# =============================================================================
#
#                               Content Handler
#
# =============================================================================
#

# Table name ---> table type mapping.
TableByName = {
	table.StripTableName(ProcessTable.tableName): ProcessTable,
	table.StripTableName(LfnTable.tableName): LfnTable,
	table.StripTableName(ProcessParamsTable.tableName): ProcessParamsTable,
	table.StripTableName(SearchSummaryTable.tableName): SearchSummaryTable,
	table.StripTableName(SearchSummVarsTable.tableName): SearchSummVarsTable,
	table.StripTableName(SnglBurstTable.tableName): SnglBurstTable,
	table.StripTableName(SnglInspiralTable.tableName): SnglInspiralTable,
	table.StripTableName(SnglRingDownTable.tableName): SnglRingDownTable,
	table.StripTableName(MultiInspiralTable.tableName): MultiInspiralTable,
	table.StripTableName(SimInspiralTable.tableName): SimInspiralTable,
	table.StripTableName(SimBurstTable.tableName): SimBurstTable,
	table.StripTableName(SimRingDownTable.tableName): SimRingDownTable,
	table.StripTableName(SummValueTable.tableName): SummValueTable,
	table.StripTableName(SimInstParamsTable.tableName): SimInstParamsTable,
	table.StripTableName(StochasticTable.tableName): StochasticTable,
	table.StripTableName(StochSummTable.tableName): StochSummTable,
	table.StripTableName(ExtTriggersTable.tableName): ExtTriggersTable,
	table.StripTableName(FilterTable.tableName): FilterTable,
	table.StripTableName(SegmentTable.tableName): SegmentTable,
	table.StripTableName(SegmentDefMapTable.tableName): SegmentDefMapTable,
	table.StripTableName(SegmentDefTable.tableName): SegmentDefTable,
	table.StripTableName(TimeSlideTable.tableName): TimeSlideTable,
	table.StripTableName(CoincDefTable.tableName): CoincDefTable,
	table.StripTableName(CoincTable.tableName): CoincTable,
	table.StripTableName(CoincMapTable.tableName): CoincMapTable,
	table.StripTableName(LIGOLWMonTable.tableName): LIGOLWMonTable
}


# Table name ---> ILWD generator mapping.
#
# FIXME: sngl_inspiral table cannot participate in generic ilwd infrastructure
# until LAL code generates event_id columns with the correct type.  Re-enable
# when LAL is fixed.
ILWDGeneratorByTableName = {
	table.StripTableName(ProcessTable.tableName): ProcessIDs,
	table.StripTableName(LfnTable.tableName): LfnIDs,
	table.StripTableName(SnglBurstTable.tableName): SnglBurstIDs,
#	table.StripTableName(SnglInspiralTable.tableName): SnglInspiralIDs,
	table.StripTableName(SimRingDownTable.tableName): SnglRingDownIDs,
	table.StripTableName(SimInspiralTable.tableName): SimInspiralIDs,
	table.StripTableName(SimBurstTable.tableName): SimBurstIDs,
	table.StripTableName(SimRingDownTable.tableName): SimRingDownIDs,
	table.StripTableName(SimInstParamsTable.tableName): SimInstParamsIDs,
	table.StripTableName(SegmentTable.tableName): SegmentIDs,
	table.StripTableName(SegmentDefMapTable.tableName): SegmentDefMapIDs,
	table.StripTableName(TimeSlideTable.tableName): TimeSlideIDs,
	table.StripTableName(CoincDefTable.tableName): CoincDefIDs,
	table.StripTableName(CoincTable.tableName): CoincIDs,
	table.StripTableName(LIGOLWMonTable.tableName): LIGOLWMonIDs
}


#
# Override portions of the ligolw.LIGOLWContentHandler class
#

__parent_startTable = ligolw.LIGOLWContentHandler.startTable

def startTable(self, attrs):
	try:
		return TableByName[table.StripTableName(attrs["Name"])](attrs)
	except KeyError:
		return __parent_startTable(self, attrs)

ligolw.LIGOLWContentHandler.startTable = startTable
