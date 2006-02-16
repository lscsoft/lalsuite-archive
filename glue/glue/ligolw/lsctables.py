"""
LSC Table definitions.  These must be kept synchronized with the official
definitions in the LDAS CVS repository at
http://www.ldas-sw.ligo.caltech.edu/cgi-bin/cvsweb.cgi/ldas/dbms/db2/sql.
Maintainership of the table definitions is left as an excercise to
interested users.
"""

__author__ = "Kipp Cannon <kipp@gravity.phys.uwm.edu>"
__date__ = "$Date$"
__version__ = "$Revision$"

import re
from xml import sax

from glue import lal
from glue import segments
import ligolw
import metaio


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
		import lsctables

		table = lsctables.New(lsctables.ProcessTable)
	"""
	table = Type(sax.xmlreader.AttributesImpl({u"Name": Type.tableName}))
	if columns != None:
		for key in columns:
			if key not in table.validcolumns.keys():
				raise ligolw.ElementError, "New(): invalid column \"%s\" for table \"%s\"." % (key, table.tableName)
	for key, value in table.validcolumns.items():
		if (columns == None) or (key in columns):
			table.appendChild(metaio.Column(sax.xmlreader.AttributesImpl({u"Name": ":".join(Type.tableName.split(":")[:-1]) + ":" + key, u"Type": value})))
	table.appendChild(metaio.TableStream(sax.xmlreader.AttributesImpl({u"Name": Type.tableName})))
	return table


def IsTableElement(Type, elem):
	"""
	Convenience function to check that an element is a Table of type
	Type.
	"""
	if elem.tagName != ligolw.Table.tagName:
		return False
	return metaio.CompareTableNames(elem.getAttribute("Name"), Type.tableName) == 0


def IsTableProperties(Type, tagname, attrs):
	"""
	Convenience function to check that the given tag name and
	attributes match those of a Table of type Type.
	"""
	if tagname != ligolw.Table.tagName:
		return False
	return metaio.CompareTableNames(attrs["Name"], Type.tableName) == 0


def getTablesByType(elem, Type):
	"""
	Return a list of tables of type Type under elem.
	"""
	return metaio.getTablesByName(elem, Type.tableName)


def getLSCTables(elem):
	"""
	Return a list of all LSC tables under elem.
	"""
	return elem.getElements(lambda e: (e.tagName == ligolw.Table.tagName) and (metaio.StripTableName(e.getAttribute("Name")) in TableByName.keys()))


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
		return ILWDPattern.search(ilwdchar).group("ID")
	except AttributeError:
		raise ValueError, "unrecognized ID \"%s\"" % repr(ilwdchar)


def getTablesFromILWD(elem, ilwdchar):
	"""
	Return a list of all tables below elem which can contain rows
	having IDs matching ilwdchar.
	"""
	return elem.getElements(lambda e: (e.tagName == ligolw.Table.tagName) and (metaio.CompareTableNames(e.getAttribute("Name"), IDTableName(ilwdchar)) == 0))


def FindILWD(tables, ilwdchar):
	"""
	Search all rows in the list of tables for one with the given ID,
	and return a reference to the row object.  If the ID is not unique,
	the first row found is returned.  KeyError is raised if no rows are
	found or the ID cannot be parsed.
	"""
	for table in tables:
		try:
			return table.dict[ilwdchar]
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
		self.table_name = metaio.StripTableName(table_name)
		self.column_name = metaio.StripColumnName(column_name)
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

	def next(self):
		row = self.iter.next()
		return (row._get_key(), row)


class LSCTableUniqueKeyIter(object):
	def __init__(self, dict):
		self.iter = iter(dict.rows)

	def next(self):
		return self.iter.next()._get_key()


class LSCTableUniqueDict(object):
	"""
	Class for implementing the Python mapping protocol on a list of table
	rows when each row has a unique key.
	"""
	def __init__(self, table):
		"""
		Initialize the mapping on the list of rows.
		"""
		self.rows = table.rows

	def __len__(self):
		"""
		Return the number of rows.
		"""
		return len(self.rows)

	def __getitem__(self, key):
		"""
		Retrieve a row by key.
		"""
		for row in self.rows:
			if row._has_key(key):
				return row
		raise KeyError, repr(key)

	def __setitem__(self, key, value):
		"""
		Set a row by key.  Note: the row key carried by value need
		not equal key, but this might not be allowed in the future.
		"""
		for i in range(len(self.rows)):
			if self.rows[i]._has_key(key):
				self.rows[i] = value
				return
		# FIXME: should we call _set_key() on value to force it to have
		# the key that was searched for?
		self.append(value)

	def __delitem__(self, key):
		"""
		Delete a row by key.
		"""
		for i in range(len(self.rows)):
			if self.rows[i]._has_get(key):
				del self.rows[i]
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
		for row in self.rows:
			if row._has_key(key):
				return True
		return False

	has_key = __contains__

	def keys(self):
		"""
		Return a list of the keys.
		"""
		return [row._get_key() for row in self.rows]

	def iteritems(self):
		"""
		Return an iterator over (key, value) pairs.
		"""
		return LSCTableUniqueItemIter(self)

	def itervalues(self):
		"""
		Return an iterator over rows.
		"""
		return iter(self.rows)


class LSCTableMultiItemIter(object):
	def __init__(self, dict):
		self.dict = dict
		self.iter = iter(dict.keys())

	def next(self):
		key = self.iter.next()
		return (key, self.dict[key])


class LSCTableMultiDict(LSCTableUniqueDict):
	"""
	Class for implementing the Python mapping protocol on a list of table
	rows when multiple rows share the same key.
	"""
	def __init__(self, table):
		"""
		Initialize the mapping on the list of rows.
		"""
		self.table = table
		self.rows = table.rows

	def __getitem__(self, key):
		"""
		Retrieve rows by key.
		"""
		l = [row for row in self.rows if row._has_key(key)]
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
		for row in self.rows:
			key = row._get_key()
			if key not in keys:
				keys.append(key)
		return keys

	def iteritems(self):
		"""
		Return an iterator over (key, value) pairs.
		"""
		return LSCTableMultiItemIter(self)


class LSCTableUnique(metaio.Table):
	def __init__(self, attrs):
		metaio.Table.__init__(self, attrs)
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


class LSCTableMulti(metaio.Table):
	def __init__(self, attrs):
		metaio.Table.__init__(self, attrs)
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


# We don't subclass metaio.TableRow because that defeats the __slots__
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
		tables = metaio.getTablesByName(elem, ProcessTable.tableName)
		for row in self.rows:
			row.process_id = FindILWD(tables, row.process_id)

	def deReference(self):
		"""
		Resolve object references back to ilwd:char strings.
		"""
		for row in self.rows:
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
		if row.type not in metaio.Types:
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

	def __getitem__(self, key):
		"""
		Return a sorted list of rows matching the process ID key.
		"""
		params = LSCTableMulti.__getitem__(self, key)
		# sort by process ID, then parameter name (all rows should
		# be unique by this measure).
		params.sort(lambda a, b: cmp((a.process_id, a.param), (b.process_id, b.param)))
		return params

	def makeReference(self, elem):
		"""
		Convert ilwd:char strings into object references.
		"""
		tables = metaio.getTablesByName(elem, ProcessTable.tableName)
		for row in self.rows:
			row.process_id = FindILWD(tables, row.process_id)

	def deReference(self):
		"""
		Resolve object references back to ilwd:char strings.
		"""
		for row in self.rows:
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

	def __getitem__(self, key):
		"""
		Return a sorted list of rows matching the process ID key.
		"""
		summaries = LSCTableMulti.__getitem__(self, key)
		# sort by process ID, then output segment (all rows should
		# be unique by this measure).
		summaries.sort(lambda a, b: cmp((a.process_id, a.get_out()), (b.process_id, b.get_out())))
		return summaries

	def makeReference(self, elem):
		"""
		Convert ilwd:char strings into object references.
		"""
		tables = metaio.getTablesByName(elem, ProcessTable.tableName)
		for row in self.rows:
			row.process_id = FindILWD(tables, row.process_id)

	def deReference(self):
		"""
		Resolve object references back to ilwd:char strings.
		"""
		for row in self.rows:
			row.process_id = row.process_id._get_key()

	def get_inlist(self):
		return segments.segmentlist([row.get_in() for row in self])

	def get_outlist(self):
		return segments.segmentlist([row.get_out() for row in self])

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
		tables = metaio.getTablesByName(elem, ProcessTable.tableName)
		for row in self.rows:
			row.process_id = FindILWD(tables, row.process_id)

	def deReference(self):
		"""
		Resolve object references back to ilwd:char strings.
		"""
		for row in self.rows:
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
		tables = metaio.getTablesByName(elem, ProcessTable.tableName)
		for row in self.rows:
			row.process_id = FindILWD(tables, row.process_id)

	def deReference(self):
		"""
		Resolve object references back to ilwd:char strings.
		"""
		for row in self.rows:
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
		"event_id": "int_8s"
	}

	def makeReference(self, elem):
		"""
		Convert ilwd:char strings into object references.
		"""
		tables = metaio.getTablesByName(elem, ProcessTable.tableName)
		for row in self.rows:
			row.process_id = FindILWD(tables, row.process_id)

	def deReference(self):
		"""
		Resolve object references back to ilwd:char strings.
		"""
		for row in self.rows:
			row.process_id = row.process_id._get_key()

class SnglInspiral(LSCTableRow):
	__slots__ = SnglInspiralTable.validcolumns.keys()

	def _get_key(self):
		return self.event_id

	def _set_key(self, key):
		self.event_id = key

	def _has_key(self, key):
		return self.event_id == key

SnglInspiralTable.RowType = SnglInspiral

class SnglInspiralIDs(ILWD):
	def __init__(self, n = 0):
		ILWD.__init__(self, "sngl_inspiral", "event_id", n)


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
		tables = metaio.getTablesByName(elem, ProcessTable.tableName)
		for row in self.rows:
			row.process_id = FindILWD(tables, row.process_id)

	def deReference(self):
		"""
		Resolve object references back to ilwd:char strings.
		"""
		for row in self.rows:
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
		tables = metaio.getTablesByName(elem, ProcessTable.tableName)
		for row in self.rows:
			row.process_id = FindILWD(tables, row.process_id)

	def deReference(self):
		"""
		Resolve object references back to ilwd:char strings.
		"""
		for row in self.rows:
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
		tables = metaio.getTablesByName(elem, ProcessTable.tableName)
		for row in self.rows:
			row.process_id = FindILWD(tables, row.process_id)

	def deReference(self):
		"""
		Resolve object references back to ilwd:char strings.
		"""
		for row in self.rows:
			row.process_id = row.process_id._get_key()

class SimInspiral(LSCTableRow):
	__slots__ = SimInspiralTable.validcolumns.keys()

	def _get_key(self):
		return self.simulation_id

	def _set_key(self, key):
		self.simulation_id = key

	def _has_key(self, key):
		return self.simulation_id == key

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
		tables = metaio.getTablesByName(elem, ProcessTable.tableName)
		for row in self.rows:
			row.process_id = FindILWD(tables, row.process_id)

	def deReference(self):
		"""
		Resolve object references back to ilwd:char strings.
		"""
		for row in self.rows:
			row.process_id = row.process_id._get_key()

class SimBurst(LSCTableRow):
	__slots__ = SimBurstTable.validcolumns.keys()

	def _get_key(self):
		return self.simulation_id

	def _set_key(self, key):
		self.simulation_id = key

	def _has_key(self, key):
		return self.simulation_id == key

	def get_geocent_peak(self):
		return lal.LIGOTimeGPS(self.geocent_peak_time, self.geocent_peak_time_ns)

	def set_geocent_peak(self, gps):
		self.geocent_peak_time, self.geocent_peak_time_ns = gps.seconds, gps.nanoseconds

	def get_h_peak(self):
		return lal.LIGOTimeGPS(self.h_peak_time, self.h_peak_time_ns)

	def set_h_peak(self, gps):
		self.h_peak_time, self.h_peak_time_ns = gps.seconds, gps.nanoseconds

	def get_l_peak(self):
		return lal.LIGOTimeGPS(self.l_peak_time, self.l_peak_time_ns)

	def set_l_peak(self, gps):
		self.l_peak_time, self.l_peak_time_ns = gps.seconds, gps.nanoseconds

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
		tables = metaio.getTablesByName(elem, ProcessTable.tableName)
		for row in self.rows:
			row.process_id = FindILWD(tables, row.process_id)

	def deReference(self):
		"""
		Resolve object references back to ilwd:char strings.
		"""
		for row in self.rows:
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
		tables = metaio.getTablesByName(elem, ProcessTable.tableName)
		for row in self.rows:
			row.process_id = FindILWD(tables, row.process_id)

	def deReference(self):
		"""
		Resolve object references back to ilwd:char strings.
		"""
		for row in self.rows:
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
		tables = metaio.getTablesByName(elem, ProcessTable.tableName)
		for row in self.rows:
			row.process_id = FindILWD(tables, row.process_id)

	def deReference(self):
		"""
		Resolve object references back to ilwd:char strings.
		"""
		for row in self.rows:
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
		tables = metaio.getTablesByName(elem, ProcessTable.tableName)
		for row in self.rows:
			row.process_id = FindILWD(tables, row.process_id)

	def deReference(self):
		"""
		Resolve object references back to ilwd:char strings.
		"""
		for row in self.rows:
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
		tables = metaio.getTablesByName(elem, ProcessTable.tableName)
		for row in self.rows:
			row.process_id = FindILWD(tables, row.process_id)

	def deReference(self):
		"""
		Resolve object references back to ilwd:char strings.
		"""
		for row in self.rows:
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
		tables = metaio.getTablesByName(elem, ProcessTable.tableName)
		for row in self.rows:
			row.process_id = FindILWD(tables, row.process_id)

	def deReference(self):
		"""
		Resolve object references back to ilwd:char strings.
		"""
		for row in self.rows:
			row.process_id = row.process_id._get_key()

	def segmentlist(self, key, active = 1):
		"""
		Return the segment list having process_id equal to key, and
		with the activity flag having the same sign as active.  The
		list represents exactly the rows in the table, in order;
		in other words it has not been coalesced.
		"""
		if not active:
			raise ValueError, "SegmentTable.segmentlist(): activity flag must != 0."
		return segments.segmentlist([row.get_segment() for row in self if (row.process_id == key) and (row.active * active > 0)])

class Segment(LSCTableRow):
	__slots__ = SegmentTable.validcolumns.keys()

	def _get_key(self):
		return self.segment_id

	def _set_key(self, key):
		self.segment_id = key

	def _has_key(self, key):
		return self.segment_id == key

	def get_segment(self):
		"""
		Return the segment described by this row.
		"""
		return segments.segment(lal.LIGOTimeGPS(self.start_time, self.start_time_ns), lal.LIGOTimeGPS(self.end_time, self.end_time_ns))

	def set_segment(self, segment):
		"""
		Set the segment described by this row.
		"""
		self.start_time, self.start_time_ns = segment[0].seconds, segment[0].nanoseconds
		self.end_time, self.end_time_ns = segment[1].seconds, segment[1].nanoseconds

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
		proctabs = metaio.getTablesByName(elem, ProcessTable.tableName)
		segtabs = metaio.getTablesByName(elem, SegmentTable.tableName)
		deftabs = metaio.getTablesByName(elem, SegmentDefTable.tableName)
		for row in self.rows:
			row.process_id = FindILWD(proctabs, row.process_id)
			row.segment_id = FindILWD(segtabs, row.segment_id)
			row.segment_def_id = FindILWD(deftabs, row.segment_def_id)

	def deReference(self):
		"""
		Resolve object references back to ilwd:char strings.
		"""
		for row in self.rows:
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
		tables = metaio.getTablesByName(elem, ProcessTable.tableName)
		for row in self.rows:
			row.process_id = FindILWD(tables, row.process_id)

	def deReference(self):
		"""
		Resolve object references back to ilwd:char strings.
		"""
		for row in self.rows:
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
#                                 slide:table
#
# =============================================================================
#

class TimeSlideTable(LSCTableMulti):
	tableName = "time_slide:table"
	validcolumns = {
		"process_id": "ilwd:char",
		"time_slide_id": "ilwd:char",
		"ifo": "lstring",
		"offset": "int_4s",
		"offset_ns": "int_4s"
	}

	def makeReference(self, elem):
		"""
		Convert ilwd:char strings into object references.
		"""
		tables = metaio.getTablesByName(elem, ProcessTable.tableName)
		for row in self.rows:
			row.process_id = FindILWD(tables, row.process_id)

	def deReference(self):
		"""
		Resolve object references back to ilwd:char strings.
		"""
		for row in self.rows:
			row.process_id = row.process_id._get_key()

class TimeSlide(LSCTableRow):
	__slots__ = TimeSlideTable.validcolumns.keys()

	def _get_key(self):
		return self.slide_id

	def _set_key(self, key):
		self.slide_id = key

	def _has_key(self, key):
		return self.slide_id == key

TimeSlideTable.RowType = TimeSlide

class TimeSlideIDs(ILWD):
	def __init__(self, n = 0):
		ILWD.__init__(self, "time_slide", "time_slide_id", n)


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
		"coinc_event_id": "ilwd:char",
		"time_slide_id": "ilwd:char"
	}

	def makeReference(self, elem):
		"""
		Convert ilwd:char strings into object references.
		"""
		proctab = metaio.getTablesByName(elem, ProcessTable.tableName)
		slidetab = metaio.getTablesByName(elem, TimeSlideTable.tableName)
		for row in self.rows:
			row.process_id = FindILWD(proctab, row.process_id)
			row.time_slide_id = FindILWD(slidetab, row.time_slide_id)

	def deReference(self):
		"""
		Resolve object references back to ilwd:char strings.
		"""
		for row in self.rows:
			row.process_id = row.process_id._get_key()
			row.time_slide_id = row.time_slide_id[0]._get_key()

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
CoincEventMapSourceNames = [
	metaio.StripTableName(SnglBurstTable.tableName),
	metaio.StripTableName(SnglInspiralTable.tableName),
	metaio.StripTableName(SnglRingDownTable.tableName),
	metaio.StripTableName(MultiInspiralTable.tableName),
	metaio.StripTableName(SimInspiralTable.tableName),
	metaio.StripTableName(SimBurstTable.tableName),
	metaio.StripTableName(SimRingDownTable.tableName),
	metaio.StripTableName(CoincTable.tableName)
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
		coinctab = metaio.getTablesByName(elem, CoincTable.tableName)
		eventtab = []
		for tablename in CoincEventMapSourceNames:
			eventtab.extend(metaio.getTablesByName(elem, tablename))
		for row in self.rows:
			row.coinc_event_id = FindILWD(coinctab, row.coinc_event_id)
			row.event_id = FindILWD(eventtab, row.event_id)

	def deReference(self):
		"""
		Resolve object references back to ilwd:char strings.
		"""
		for row in self.rows:
			row.coinc_event_id = row.coinc_event_id._get_key()
			row.event_id = row.event_id._get_key()

class CoincMap(LSCTableRow):
	__slots__ = CoincMapTable.validcolumns.keys()

CoincMapTable.RowType = CoincMap


#
# =============================================================================
#
#                               Content Handler
#
# =============================================================================
#

# Table name ---> table type mapping.
TableByName = {
	metaio.StripTableName(ProcessTable.tableName): ProcessTable,
	metaio.StripTableName(LfnTable.tableName): LfnTable,
	metaio.StripTableName(ProcessParamsTable.tableName): ProcessParamsTable,
	metaio.StripTableName(SearchSummaryTable.tableName): SearchSummaryTable,
	metaio.StripTableName(SearchSummVarsTable.tableName): SearchSummVarsTable,
	metaio.StripTableName(SnglBurstTable.tableName): SnglBurstTable,
	metaio.StripTableName(SnglInspiralTable.tableName): SnglInspiralTable,
	metaio.StripTableName(SnglRingDownTable.tableName): SnglRingDownTable,
	metaio.StripTableName(MultiInspiralTable.tableName): MultiInspiralTable,
	metaio.StripTableName(SimInspiralTable.tableName): SimInspiralTable,
	metaio.StripTableName(SimBurstTable.tableName): SimBurstTable,
	metaio.StripTableName(SimRingDownTable.tableName): SimRingDownTable,
	metaio.StripTableName(SummValueTable.tableName): SummValueTable,
	metaio.StripTableName(SimInstParamsTable.tableName): SimInstParamsTable,
	metaio.StripTableName(StochasticTable.tableName): StochasticTable,
	metaio.StripTableName(StochSummTable.tableName): StochSummTable,
	metaio.StripTableName(ExtTriggersTable.tableName): ExtTriggersTable,
	metaio.StripTableName(FilterTable.tableName): FilterTable,
	metaio.StripTableName(SegmentTable.tableName): SegmentTable,
	metaio.StripTableName(SegmentDefMapTable.tableName): SegmentDefMapTable,
	metaio.StripTableName(SegmentDefTable.tableName): SegmentDefTable,
	metaio.StripTableName(TimeSlideTable.tableName): TimeSlideTable,
	metaio.StripTableName(CoincTable.tableName): CoincTable,
	metaio.StripTableName(CoincMapTable.tableName): CoincMapTable
}


# Table name ---> ILWD generator mapping.
ILWDGeneratorByTableName = {
	metaio.StripTableName(ProcessTable.tableName): ProcessIDs,
	metaio.StripTableName(LfnTable.tableName): LfnIDs,
	metaio.StripTableName(SnglBurstTable.tableName): SnglBurstIDs,
	metaio.StripTableName(SnglInspiralTable.tableName): SnglInspiralIDs,
	metaio.StripTableName(SimRingDownTable.tableName): SnglRingDownIDs,
	metaio.StripTableName(SimInspiralTable.tableName): SimInspiralIDs,
	metaio.StripTableName(SimBurstTable.tableName): SimBurstIDs,
	metaio.StripTableName(SimRingDownTable.tableName): SimRingDownIDs,
	metaio.StripTableName(SimInstParamsTable.tableName): SimInstParamsIDs,
	metaio.StripTableName(SegmentTable.tableName): SegmentIDs,
	metaio.StripTableName(SegmentDefMapTable.tableName): SegmentDefMapIDs,
	metaio.StripTableName(TimeSlideTable.tableName): TimeSlideIDs,
	metaio.StripTableName(CoincTable.tableName): CoincIDs
}


class LIGOLWContentHandler(metaio.LIGOLWContentHandler):
	"""
	ContentHandler that redirects Table elements with known structure
	to the definitions in this module, using the Table element in
	metaio as a fall-back for unrecognized tables.
	"""
	def startTable(self, attrs):
		try:
			return TableByName[metaio.StripTableName(attrs["Name"])](attrs)
		except KeyError:
			return metaio.Table(attrs)
