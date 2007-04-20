# $Id$
#
# Copyright (C) 2007  Kipp C. Cannon
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
This module provides an implementation of the Table element that uses a
database engine for storage.  On top of that it then re-implements a number
of the tables from the lsctables module to provide versions of their
methods that work against the SQL database.

*** CAUTION *** the API exported by this module is NOT STABLE.  The code
works well, and it hugely simplifies my life and it can yours too, but if
you use it you will have to be willing to track changes as I make them.
I'm still figuring out how this should work.
"""


import re
import sys
from xml.sax.xmlreader import AttributesImpl
# Python 2.3 compatibility
try:
	from sets import Set as set
except:
	pass

import ligolw
import table
import lsctables
import types
from glue import segments


__author__ = "Kipp Cannon <kipp@gravity.phys.uwm.edu>"
__date__ = "$Date$"[7:-2]
__version__ = "$Revision$"[11:-2]


#
# =============================================================================
#
#                                  Connection
#
# =============================================================================
#


def DBTable_set_connection(connection):
	"""
	Set the Python DB-API 2.0 compatible connection the DBTable class
	will use.
	"""
	DBTable.connection = connection


def DBTable_get_connection():
	"""
	Return the current connection object.
	"""
	return DBTable.connection


def DBTable_commit():
	"""
	Run commit on the DBTable class' connection.
	"""
	DBTable.connection.commit()


#
# =============================================================================
#
#                                  ID Mapping
#
# =============================================================================
#


def DBTable_idmap_create():
	"""
	Create the _idmap_ table.  This table has columns "old" and "new"
	containing text strings mapping old IDs to new IDs.  The old column
	is a primary key (is indexed and must contain unique entries).  The
	table is created as a temporary table, so it will be automatically
	dropped when the database connection is closed.
	"""
	DBTable.connection.cursor().execute("CREATE TEMPORARY TABLE _idmap_ (old TEXT PRIMARY KEY, new TEXT)")


def DBTable_idmap_reset():
	"""
	Erase the contents of the _idmap_ table.
	"""
	DBTable.connection.cursor().execute("DELETE FROM _idmap_")


def DBTable_idmap_get_new(old, ids):
	"""
	From the old ID string, obtain a replacement ID string by either
	grabbing it from the _idmap_ table if one has already been assigned
	to the old ID, or by creating a new one with the given ILWD
	generator.  In the latter case, the newly-generated ID is recorded
	in the _idmap_ table.  For internal use only.
	"""
	cursor = DBTable.connection.cursor()
	new = cursor.execute("SELECT new FROM _idmap_ WHERE old == ?", (old,)).fetchone()
	if new is not None:
		return new[0]
	new = ids.next()
	cursor.execute("INSERT INTO _idmap_ VALUES (?, ?)", (old, new))
	return new


#
# =============================================================================
#
#                             Database Information
#
# =============================================================================
#


#
# SQL parsing
#

_sql_create_table_pattern = re.compile(r"CREATE\s+TABLE\s+(?P<name>\w+)\s*\((?P<coldefs>.*)\)")
_sql_coldef_pattern = re.compile(r"\s*(?P<name>\w+)\s+(?P<type>\w+)[^,]*")


#
# Database info extraction utils
#


def DBTable_table_names():
	"""
	Return a list of the table names in the database.
	"""
	return [name for (name,) in DBTable.connection.cursor().execute("SELECT name FROM sqlite_master WHERE type == 'table'")]


def DBTable_column_info(table_name):
	"""
	Return an in order list of (name, type) tuples describing the
	columns for the given table.
	"""
	statement, = DBTable.connection.cursor().execute("SELECT sql FROM sqlite_master WHERE type == 'table' AND name == ?", (table_name,)).fetchone()
	coldefs = re.match(_sql_create_table_pattern, statement.upper()).groupdict()["coldefs"]
	return [(coldef.groupdict()["name"].lower(), coldef.groupdict()["type"]) for coldef in re.finditer(_sql_coldef_pattern, coldefs) if coldef.groupdict()["name"] not in ("PRIMARY", "UNIQUE", "CHECK")]


def DBTable_get_xml():
	"""
	Construct an XML document tree wrapping around the contents of the
	currently connected database.  On success the return value is a
	ligolw.LIGO_LW element containing the tables as children.
	"""
	ligo_lw = ligolw.LIGO_LW()
	for table_name in DBTable_table_names():
		# build the table document tree.  copied from
		# lsctables.New()
		try:
			cls = TableByName[table_name]
		except KeyError:
			cls = DBTable
		table_elem = cls(AttributesImpl({u"Name": table_name + ":table"}))
		colnamefmt = table_name + ":%s"
		for column_name, column_type in DBTable_column_info(table_elem.dbtablename):
			if table_elem.validcolumns is not None:
				# use the pre-defined column type
				column_type = table_elem.validcolumns[column_name]
			else:
				# guess the column type
				column_type = types.FromSQLiteType[column_type]
			column_name = colnamefmt % column_name
			table_elem.appendChild(table.Column(AttributesImpl({u"Name": column_name, u"Type": column_type})))

		table_elem._end_of_columns()
		table_elem.appendChild(table.TableStream(AttributesImpl({u"Name": table_name + ":table"})))
		ligo_lw.appendChild(table_elem)
	return ligo_lw


#
# =============================================================================
#
#                            DBTable Element Class
#
# =============================================================================
#


class DBTable(table.Table):
	"""
	A special version of the Table class using an SQL database for
	storage.  Many of the features of the Table class are not available
	here, but instead the user can use SQL to query the table's
	contents.  Before use, the connection attribute must be set to a
	Python DB-API 2.0 "connection" object.  The constraints attribute
	can be set to a text string that will be added to the table's
	CREATE statement where constraints go, for example you might wish
	to set this to "PRIMARY KEY (event_id)" for a table with an
	event_id column.

	Note:  because the table is stored in an SQL database, the use of
	this class imposes the restriction that table names be unique
	within a document.

	Also note that at the present time there is really only proper
	support for the pre-defined tables in the lsctables module.  It is
	possible to load unrecognized tables into a database from LIGO
	Light Weight XML files, but without developer intervention there is
	no way to indicate the constraints that should be imposed on the
	columns, for example which columns should be used as primary keys
	and so on.  This can result in poor query performance.  It is also
	possible to write a database' contents to a LIGO Light Weight XML
	file even when the database contains unrecognized tables, but
	without developer intervention the column types will be guessed.
	"""
	#
	# Global Python DB-API 2.0 connection object shared by all code.
	#

	connection = None

	def __init__(self, *attrs):
		"""
		Initialize
		"""
		table.Table.__init__(self, *attrs)
		self.dbtablename = table.StripTableName(self.getAttribute(u"Name"))
		try:
			# copy metadata from lsctables
			cls = lsctables.TableByName[self.dbtablename]
			self.tableName = cls.tableName
			self.validcolumns = cls.validcolumns
			self.constraints = cls.constraints
			self.ids = cls.ids
			self.RowType = cls.RowType
		except:
			# unknown table
			pass
		if self.connection is None:
			raise ligolw.ElementError, "connection attribute not set"
		self.cursor = self.connection.cursor()

	def _end_of_columns(self):
		table.Table._end_of_columns(self)
		# dbcolumnnames and types have the "not loaded" columns removed
		if self.loadcolumns is not None:
			self.dbcolumnnames = [name for name in self.columnnames if name in self.loadcolumns]
			self.dbcolumntypes = [name for i, name in enumerate(self.columntypes) if self.columnnames[i] in self.loadcolumns]
		else:
			self.dbcolumnnames = self.columnnames
			self.dbcolumntypes = self.columntypes

		# create the table
		statement = "CREATE TABLE IF NOT EXISTS " + self.dbtablename + " (" + ", ".join(map(lambda n, t: "%s %s" % (n, types.ToSQLiteType[t]), self.dbcolumnnames, self.dbcolumntypes))
		if self.constraints is not None:
			statement += ", " + self.constraints
		statement += ")"
		self.cursor.execute(statement)

		# record the highest internal row ID
		self.last_maxrowid = self.maxrowid() or 0

		# construct the SQL to be used to insert new rows
		self.append_statement = "INSERT INTO " + self.dbtablename + " VALUES (" + ",".join("?" * len(self.dbcolumnnames)) + ")"

	def _end_of_rows(self):
		# FIXME:  is this needed?
		table.Table._end_of_rows(self)
		self.connection.commit()

	def sync_ids(self):
		if self.ids is not None:
			statement = "SELECT MAX(CAST(SUBSTR(%s, %d, 10) AS INTEGER)) FROM %s" % (self.ids.column_name, self.ids.index_offset + 1, self.dbtablename)
			last = self.cursor.execute(statement).fetchone()[0]
			if last is None:
				self.ids.set_next(0)
			else:
				self.ids.set_next(last + 1)
		return self.ids

	def maxrowid(self):
		return self.cursor.execute("SELECT MAX(ROWID) FROM %s" % self.dbtablename).fetchone()[0]

	def __len__(self):
		return self.cursor.execute("SELECT COUNT(*) FROM " + self.dbtablename).fetchone()[0]

	def __iter__(self):
		for values in self.connection.cursor().execute("SELECT * FROM " + self.dbtablename):
			yield self._row_from_cols(values)

	def _append(self, row):
		# FIXME: in Python 2.5 use attrgetter() for attribute
		# tuplization.
		self.cursor.execute(self.append_statement, map(lambda n: getattr(row, n), self.dbcolumnnames))

	def _remapping_append(self, row):
		"""
		Replacement for the standard append() method.  This version
		performs on the fly row ID reassignment, and so also
		performs the function of the updateKeyMapping() method.
		SQL does not permit the PRIMARY KEY of a row to be
		modified, so it needs to be done prior to insertion.  This
		method is intended for internal use only.
		"""
		if self.ids is not None:
			# assign (and record) a new ID before inserting the
			# row to avoid collisions with existing rows
			setattr(row, self.ids.column_name, DBTable_idmap_get_new(getattr(row, self.ids.column_name), self.ids))
		# FIXME: in Python 2.5 use attrgetter() for attribute
		# tuplization.
		self.cursor.execute(self.append_statement, map(lambda n: getattr(row, n), self.dbcolumnnames))

	append = _append

	def _row_from_cols(self, values):
		"""
		Given a tuple or list of values in the order of columns in
		the database, construct and return a row object.  This is a
		convenience function for turning the results of database
		queries into Python objects.
		"""
		row = self.RowType()
		for c, v in zip(self.dbcolumnnames, values):
			setattr(row, c, v)
		return row

	def unlink(self):
		table.Table.unlink(self)
		#self.cursor.execute("DROP TABLE " + self.dbtablename)
		self.cursor = None

	def applyKeyMapping(self):
		"""
		Used as the second half of the key reassignment algorithm.
		Loops over each row in the table, replacing references to
		old row keys with the new values from the _idmap_ table.
		"""
		assignments = []
		for colname in [colname for coltype, colname in zip(self.dbcolumntypes, self.dbcolumnnames) if coltype in types.IDTypes and (self.ids is None or colname != self.ids.column_name)]:
			assignments.append(" %s = (SELECT new FROM _idmap_ WHERE old == %s)" % (colname, colname))
		if not assignments:
			# table has no columns to update
			return
		# SQLite documentation says ROWID is monotonically
		# increasing starting at 1 for the first row unless it ever
		# wraps around, then it is randomly assigned.  ROWID is a
		# 64 bit integer, so the only way it will wrap is if
		# somebody sets it to a very high number manually.  This
		# library does not do that, so I don't bother checking.
		statement = "UPDATE " + self.dbtablename + " SET" + ",".join(assignments) + " WHERE ROWID > %d" % self.last_maxrowid
		self.cursor.execute(statement)
		self.last_maxrowid = self.maxrowid() or 0


#
# =============================================================================
#
#                                  LSC Tables
#
# =============================================================================
#


class ProcessTable(DBTable):
	tableName = lsctables.ProcessTable.tableName
	validcolumns = lsctables.ProcessTable.validcolumns
	constraints = lsctables.ProcessTable.constraints
	ids = lsctables.ProcessTable.ids
	RowType = lsctables.ProcessTable.RowType

	def get_ids_by_program(self, program):
		"""
		Return a set of the process IDs from rows whose program
		string equals the given program.
		"""
		return set(id for (id,) in self.cursor.execute("""
			SELECT process_id FROM
				process
			WHERE
				program == ?
		""", (program,)))


class ProcessParamsTable(DBTable):
	tableName = lsctables.ProcessParamsTable.tableName
	validcolumns = lsctables.ProcessParamsTable.validcolumns
	constraints = lsctables.ProcessParamsTable.constraints
	ids = lsctables.ProcessParamsTable.ids
	RowType = lsctables.ProcessParamsTable.RowType

	def append(self, row):
		if row.type not in types.Types:
			raise ligolw.ElementError, "unrecognized type '%s'" % row.type
		DBTable.append(self, row)


class SearchSummaryTable(DBTable):
	tableName = lsctables.SearchSummaryTable.tableName
	validcolumns = lsctables.SearchSummaryTable.validcolumns
	constraints = lsctables.SearchSummaryTable.constraints
	ids = lsctables.SearchSummaryTable.ids
	RowType = lsctables.SearchSummaryTable.RowType

	def get_out_segmentlistdict(self, process_ids = None):
		"""
		Return a segmentlistdict mapping instrument to out segment
		list.  If process_ids is a list of process IDs, then only
		rows with matching IDs are included otherwise all rows are
		included.
		"""
		seglistdict = segments.segmentlistdict()
		for row in self:
			if process_ids is None or row.process_id in process_ids:
				for ifo in row.ifos.split(","):
					if ifo in seglistdict:
						seglistdict[ifo].append(row.get_out())
					else:
						seglistdict[ifo] = segments.segmentlist([row.get_out()])
		return seglistdict.coalesce()


class SnglBurstTable(DBTable):
	tableName = lsctables.SnglBurstTable.tableName
	validcolumns = lsctables.SnglBurstTable.validcolumns
	constraints = lsctables.SnglBurstTable.constraints
	ids = lsctables.SnglBurstTable.ids
	RowType = lsctables.SnglBurstTable.RowType


class SimBurstTable(DBTable):
	tableName = lsctables.SimBurstTable.tableName
	validcolumns = lsctables.SimBurstTable.validcolumns
	constraints = lsctables.SimBurstTable.constraints
	ids = lsctables.SimBurstTable.ids
	RowType = lsctables.SimBurstTable.RowType


class TimeSlideTable(DBTable):
	tableName = lsctables.TimeSlideTable.tableName
	validcolumns = lsctables.TimeSlideTable.validcolumns
	constraints = lsctables.TimeSlideTable.constraints
	ids = lsctables.TimeSlideTable.ids
	RowType = lsctables.TimeSlideTable.RowType

	def __len__(self):
		return self.cursor.execute("SELECT COUNT(DISTINCT time_slide_id) FROM time_slide").fetchone()[0]

	def __getitem__(self, id):
		offsets = {}
		for instrument, offset in self.cursor.execute("SELECT instrument, offset FROM time_slide WHERE time_slide_id == ?", (id,)):
			offsets[instrument] = offset
		return offsets

	def iterkeys(self):
		for (id,) in self.connection.cursor().execute("SELECT DISTINCT time_slide_id FROM time_slide"):
			yield id

	def is_null(self, id):
		return not self.cursor.execute("SELECT EXISTS (SELECT * FROM time_slide WHERE time_slide_id == ? AND offset != 0.0)", (id,)).fetchone()[0]


class CoincDefTable(DBTable):
	tableName = lsctables.CoincDefTable.tableName
	validcolumns = lsctables.CoincDefTable.validcolumns
	constraints = lsctables.CoincDefTable.constraints
	ids = lsctables.CoincDefTable.ids
	RowType = lsctables.CoincDefTable.RowType


class CoincTable(DBTable):
	tableName = lsctables.CoincTable.tableName
	validcolumns = lsctables.CoincTable.validcolumns
	constraints = lsctables.CoincTable.constraints
	ids = lsctables.CoincTable.ids
	RowType = lsctables.CoincTable.RowType
	how_to_index = lsctables.CoincTable.how_to_index


class CoincMapTable(DBTable):
	tableName = lsctables.CoincMapTable.tableName
	validcolumns = lsctables.CoincMapTable.validcolumns
	constraints = lsctables.CoincMapTable.constraints
	ids = lsctables.CoincMapTable.ids
	RowType = lsctables.CoincMapTable.RowType
	how_to_index = lsctables.CoincMapTable.how_to_index


#
# =============================================================================
#
#                                Table Metadata
#
# =============================================================================
#


def build_indexes(verbose = False):
	"""
	Using the how_to_index annotations in the table class definitions,
	construct a set of indexes for the database at the current
	connection.
	"""
	cursor = DBTable_get_connection().cursor()
	for table_name in DBTable_table_names():
		indexes = TableByName[table_name].how_to_index
		if verbose and indexes:
			print >>sys.stderr, "indexing %s table ..." % table_name
		for index_name, cols in indexes.iteritems():
			cursor.execute("CREATE INDEX IF NOT EXISTS %s ON %s (%s)" % (index_name, table_name, ",".join(cols)))


#
# =============================================================================
#
#                                Table Metadata
#
# =============================================================================
#


#
# Table name ---> table type mapping.
#


TableByName = {
	table.StripTableName(ProcessTable.tableName): ProcessTable,
	table.StripTableName(ProcessParamsTable.tableName): ProcessParamsTable,
	table.StripTableName(SearchSummaryTable.tableName): SearchSummaryTable,
	table.StripTableName(SnglBurstTable.tableName): SnglBurstTable,
	table.StripTableName(SimBurstTable.tableName): SimBurstTable,
	table.StripTableName(TimeSlideTable.tableName): TimeSlideTable,
	table.StripTableName(CoincDefTable.tableName): CoincDefTable,
	table.StripTableName(CoincTable.tableName): CoincTable,
	table.StripTableName(CoincMapTable.tableName): CoincMapTable
}


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


def startTable(self, attrs):
	name = table.StripTableName(attrs[u"Name"])
	if name in TableByName:
		return TableByName[name](attrs)
	return DBTable(attrs)


ligolw.LIGOLWContentHandler.startTable = startTable
