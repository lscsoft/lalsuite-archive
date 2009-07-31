# $Id$
#
# Copyright (C) 2007  Kipp C. Cannon
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
This module provides an implementation of the Table element that uses a
database engine for storage.  On top of that it then re-implements a number
of the tables from the lsctables module to provide versions of their
methods that work against the SQL database.
"""


import os
import re
import shutil
import signal
import sys
import tempfile
from xml.sax.xmlreader import AttributesImpl
# Python 2.3 compatibility
try:
	set
except NameError:
	from sets import Set as set


import ligolw
import table
import lsctables
import types as ligolwtypes
from glue import segments


__author__ = "Kipp Cannon <kcannon@ligo.caltech.edu>"
__date__ = "$Date$"[7:-2]
__version__ = "$Revision$"[11:-2]


#
# =============================================================================
#
#                                  Connection
#
# =============================================================================
#


def connection_db_type(connection):
	"""
	A totally broken attempt to determine what type of database a
	connection object is attached to.  Don't use this.

	The input is a DB API 2.0 compliant connection object, the return
	value is one of the strings "sqlite3" or "mysql".  Raises TypeError
	when the database type cannot be determined.
	"""
	if "sqlite3" in repr(connection):
		return "sqlite3"
	if "mysql" in repr(connection):
		return "mysql"
	raise TypeError, connection


def DBTable_set_connection(connection):
	"""
	Set the Python DB-API 2.0 compatible connection the DBTable class
	will use.
	"""
	DBTable.connection = connection


def get_connection_filename(filename, tmp_path = None, replace_file = False, verbose = False):
	"""
	Utility code for moving database files to a (presumably local)
	working location for improved performance and reduced fileserver
	load.
	"""
	def mktmp(path, verbose = False):
		fd, filename = tempfile.mkstemp(suffix = ".sqlite", dir = path)
		os.close(fd)
		if verbose:
			print >>sys.stderr, "using '%s' as workspace" % filename
		return filename

	def truncate(filename, verbose = False):
		if verbose:
			print >>sys.stderr, "'%s' exists, truncating ..." % filename
		try:
			fd = os.open(filename, os.O_WRONLY | os.O_TRUNC)
		except:
			if verbose:
				print >>sys.stderr, "cannot truncate '%s': %s" % (filename, str(e))
			return
		os.close(fd)
		if verbose:
			print >>sys.stderr, "done."

	def cpy(srcname, dstname, verbose = False):
		if verbose:
			print >>sys.stderr, "copying '%s' to '%s' ..." % (srcname, dstname)
		shutil.copy(srcname, dstname)

	database_exists = os.access(filename, os.F_OK)

	if tmp_path is not None:
		target = mktmp(tmp_path, verbose)
		if database_exists:
			if replace_file:
				# truncate database so that if this job
				# fails the user won't think the database
				# file is valid
				truncate(filename, verbose = verbose)
			else:
				# need to copy existing database to work
				# space for modifications
				i = 1
				while True:
					try:
						cpy(filename, target, verbose)
					except IOError, e:
						import errno
						import time
						if e.errno == errno.ENOSPC:
							if i < 5:
								if verbose:
									print >>sys.stderr, "warning: attempt %d: no space left on device, sleeping and trying again ..." % i
								time.sleep(10)
								i += 1
								continue
							else:
								if verbose:
									print >>sys.stderr, "warning: attempt %d: no space left on device: working with original file" % i
								os.remove(target)
								target = filename
						else:
							raise e
					break
	else:
		target = filename
		if database_exists and replace_file:
			truncate(target, verbose = verbose)

	del mktmp
	del truncate
	del cpy

	return target


class IOTrappedSignal(Exception):
	"""
	Raised by put_connection_filename() upon completion if it trapped a
	signal during the operation
	"""
	def __init__(self, signum):
		self.signum = signum

	def __str__(self):
		return "trapped signal %d" % self.signum


def put_connection_filename(filename, working_filename, verbose = False):
	"""
	This function reverses the effect of a previous call to
	get_connection_filename(), restoring the working copy to its
	original location if the two are different.  This function should
	always be called after calling get_connection_filename() when the
	file is no longer in use.

	This function traps the signals used by Condor to evict jobs, which
	reduces the risk of corrupting a document by the job terminating
	part-way through the restoration of the file to its original
	location.
	"""
	if working_filename != filename:
		# initialize SIGTERM and SIGTSTP trap
		global __llwapp_write_filename_got_sig
		__llwapp_write_filename_got_sig = []
		def newsigterm(signum, frame):
			global __llwapp_write_filename_got_sig
			__llwapp_write_filename_got_sig.append(signum)
		oldhandlers = {}
		for sig in (signal.SIGTERM, signal.SIGTSTP):
			oldhandlers[sig] = signal.getsignal(sig)
			signal.signal(sig, newsigterm)

		# replace document
		if verbose:
			print >>sys.stderr, "moving '%s' to '%s' ..." % (working_filename, filename)
		shutil.move(working_filename, filename)

		# restore original handlers, and report the most recently
		# trapped signal if any were
		for sig, oldhandler in oldhandlers.iteritems():
			signal.signal(sig, oldhandler)
		if __llwapp_write_filename_got_sig:
			raise
			IOTrappedSignal(__llwapp_write_filename_got_sig.pop())


def discard_connection_filename(filename, working_filename, verbose = False):
	"""
	Like put_connection_filename(), but the working copy is simply
	deleted instead of being copied back to its original location.
	This is a useful performance boost if it is known that no
	modifications were made to the file, for example if queries were
	performed but no updates.

	Note that the file is not deleted if the working copy and original
	file are the same, so it is always safe to call this function after
	a call to get_connection_filename() even if a separate working copy
	is not created.
	"""
	if working_filename != filename:
		if verbose:
			print >>sys.stderr, "removing '%s' ..." % working_filename
		os.remove(working_filename)


#
# =============================================================================
#
#                                  ID Mapping
#
# =============================================================================
#


def idmap_create(connection):
	"""
	Create the _idmap_ table.  This table has columns "old" and "new"
	containing text strings mapping old IDs to new IDs.  The old column
	is a primary key (is indexed and must contain unique entries).  The
	table is created as a temporary table, so it will be automatically
	dropped when the database connection is closed.

	This function is for internal use, it forms part of the code used
	to re-map row IDs when merging multiple documents.
	"""
	connection.cursor().execute("CREATE TEMPORARY TABLE _idmap_ (old TEXT PRIMARY KEY, new TEXT)")


def idmap_reset(connection):
	"""
	Erase the contents of the _idmap_ table, but leave the table in
	place.

	This function is for internal use, it forms part of the code used
	to re-map row IDs when merging multiple documents.
	"""
	connection.cursor().execute("DELETE FROM _idmap_")


def idmap_get_new(connection, old, tbl):
	"""
	From the old ID string, obtain a replacement ID string by either
	grabbing it from the _idmap_ table if one has already been assigned
	to the old ID, or by using the current value of the Table
	instance's next_id class attribute.  In the latter case, the new ID
	is recorded in the _idmap_ table, and the class attribute
	incremented by 1.

	This function is for internal use, it forms part of the code used
	to re-map row IDs when merging multiple documents.
	"""
	cursor = connection.cursor()
	cursor.execute("SELECT new FROM _idmap_ WHERE old == ?", (old,))
	new = cursor.fetchone()
	if new is not None:
		# a new ID has already been created for this old ID
		return new[0]
	# this ID was not found in _idmap_ table, assign a new ID and
	# record it
	new = unicode(tbl.get_next_id())
	cursor.execute("INSERT INTO _idmap_ VALUES (?, ?)", (old, new))
	return new


def idmap_get_max_id(connection, id_class):
	"""
	Given an ilwd:char ID class, return the highest ID from the table
	for whose IDs that is the class.

	Example:

	>>> id = ilwd.get_ilwdchar("sngl_burst:event_id:0")
	>>> print id
	sngl_inspiral:event_id:0
	>>> max = get_max_id(connection, type(id))
	>>> print max
	sngl_inspiral:event_id:1054
	"""
	cursor = connection.cursor()
	cursor.execute("SELECT MAX(CAST(SUBSTR(%s, %d, 10) AS INTEGER)) FROM %s" % (id_class.column_name, id_class.index_offset + 1, id_class.table_name))
	max = cursor.fetchone()[0]
	if max is None:
		return None
	return id_class(max)


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


_sql_create_table_pattern = re.compile(r"CREATE\s+TABLE\s+(?P<name>\w+)\s*\((?P<coldefs>.*)\)", re.IGNORECASE)
_sql_coldef_pattern = re.compile(r"\s*(?P<name>\w+)\s+(?P<type>\w+)[^,]*")


#
# Database info extraction utils
#


def get_table_names(connection):
	"""
	Return a list of the table names in the database.
	"""
	cursor = connection.cursor()
	cursor.execute("SELECT name FROM sqlite_master WHERE type == 'table'")
	return [name for (name,) in cursor]


def get_column_info(connection, table_name):
	"""
	Return an in order list of (name, type) tuples describing the
	columns in the given table.
	"""
	cursor = connection.cursor()
	cursor.execute("SELECT sql FROM sqlite_master WHERE type == 'table' AND name == ?", (table_name,))
	statement, = cursor.fetchone()
	coldefs = re.match(_sql_create_table_pattern, statement).groupdict()["coldefs"]
	return [(coldef.groupdict()["name"], coldef.groupdict()["type"]) for coldef in re.finditer(_sql_coldef_pattern, coldefs) if coldef.groupdict()["name"].upper() not in ("PRIMARY", "UNIQUE", "CHECK")]


def get_xml(connection, table_names = None):
	"""
	Construct an XML document tree wrapping around the contents of the
	database.  On success the return value is a ligolw.LIGO_LW element
	containing the tables as children.  Arguments are a connection to
	to a database, and an optional list of table names to dump.  If
	table_names is not provided the set is obtained from get_table_names()
	"""
	ligo_lw = ligolw.LIGO_LW()

	if table_names is None:
		table_names = get_table_names(connection)

	for table_name in table_names:
		# build the table document tree.  copied from
		# lsctables.New()
		try:
			cls = TableByName[table_name]
		except KeyError:
			cls = DBTable
		table_elem = cls(AttributesImpl({u"Name": u"%s:table" % table_name}), connection = connection)
		for column_name, column_type in get_column_info(connection, table_elem.dbtablename):
			if table_elem.validcolumns is not None:
				# use the pre-defined column type
				column_type = table_elem.validcolumns[column_name]
			else:
				# guess the column type
				column_type = ligolwtypes.FromSQLiteType[column_type]
			table_elem.appendChild(table.Column(AttributesImpl({u"Name": u"%s:%s" % (table_name, column_name), u"Type": column_type})))
		table_elem._end_of_columns()
		table_elem.appendChild(table.TableStream(AttributesImpl({u"Name": u"%s:table" % table_name})))
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
	contents.

	The constraints attribute can be set to a text string that will be
	added to the table's CREATE statement where constraints go, for
	example you might wish to set this to "PRIMARY KEY (event_id)" for
	a table with an event_id column.

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
	without developer intervention the column types will be guessed
	using a generic mapping of SQL types to LIGO Light Weight types.

	Each instance of this class must be connected to a database, but
	because the __init__() method must have the same signature as the
	__init__() methods of other Element subclasses (so that this class
	can be used as a drop-in replacement during document parsing), it
	is not possible to pass a connection object into the __init__()
	method as a separate argument.  Instead, the connection object is
	passed into the __init__() method via a class attribute.  The
	procedure is:  set the class attribute, create instances of the
	class connected to that database, set the class attribute to a new
	value, create instances of the class connected to a different
	database, and so on.  A utility function is available for setting
	the class attribute, and should be used.

	Example:

	import dbtables

	# set class attribute
	dbtables.DBTable_set_connection(connection)

	# create a process table instance connected to that database
	table_elem = dbtables.DBTable(AttributesImpl({u"Name": u"process:table"}))
	"""
	#
	# When instances of this class are created, they initialize their
	# connection attributes from the value of this class attribute at
	# the time of their creation.  The value must be passed into the
	# __init__() method this way because it cannot be passed in as a
	# separate argument;  this class' __init__() method must have the
	# same signature as normal Element subclasses.
	#

	connection = None

	def __new__(cls, *args, **kwargs):
		# does this class already have table-specific metadata?
		if not hasattr(cls, "tableName"):
			# no, try to retrieve it from lsctables
			attrs, = args
			name = table.StripTableName(attrs[u"Name"])
			try:
				lsccls = lsctables.TableByName[name]
			except KeyError:
				# unknown table, give up
				pass
			else:
				# found metadata, construct custom
				# subclass.  NOTE:  this works because when
				# using SQL-backed tables there can only be
				# ONE of any table in a document, which
				# solves the problem of trying to share the
				# next_id attribute across multiple
				# instances of a table.
				class CustomDBTable(cls):
					tableName = lsccls.tableName
					validcolumns = lsccls.validcolumns
					loadcolumns = lsccls.loadcolumns
					constraints = lsccls.constraints
					next_id = lsccls.next_id
					RowType = lsccls.RowType
					how_to_index = lsccls.how_to_index

				# save for re-use (required for ID
				# remapping across multiple documents in
				# ligolw_sqlite)
				TableByName[name] = CustomDBTable

				# replace input argument with new class
				cls = CustomDBTable
		return table.Table.__new__(cls, *args)

	def __init__(self, *args, **kwargs):
		# chain to parent class
		table.Table.__init__(self, *args)

		# save the stripped name
		self.dbtablename = table.StripTableName(self.getAttribute(u"Name"))

		# replace connection class attribute with an instance
		# attribute
		if "connection" in kwargs:
			self.connection = kwargs.pop("connection")
		else:
			self.connection = self.connection

		# pre-allocate a cursor for internal queries
		self.cursor = self.connection.cursor()

	def _end_of_columns(self):
		table.Table._end_of_columns(self)
		# dbcolumnnames and types have the "not loaded" columns
		# removed
		if self.loadcolumns is not None:
			self.dbcolumnnames = [name for name in self.columnnames if name in self.loadcolumns]
			self.dbcolumntypes = [name for i, name in enumerate(self.columntypes) if self.columnnames[i] in self.loadcolumns]
		else:
			self.dbcolumnnames = self.columnnames
			self.dbcolumntypes = self.columntypes

		# create the table
		ToSQLType = {
			"sqlite3": ligolwtypes.ToSQLiteType,
			"mysql": ligolwtypes.ToMySQLType
		}[connection_db_type(self.connection)]
		statement = "CREATE TABLE IF NOT EXISTS " + self.dbtablename + " (" + ", ".join(map(lambda n, t: "%s %s" % (n, ToSQLType[t]), self.dbcolumnnames, self.dbcolumntypes))
		if self.constraints is not None:
			statement += ", " + self.constraints
		statement += ")"
		self.cursor.execute(statement)

		# record the highest internal row ID
		self.last_maxrowid = self.maxrowid() or 0

		# construct the SQL to be used to insert new rows
		params = {
			"sqlite3": ",".join("?" * len(self.dbcolumnnames)),
			"mysql": ",".join(["%s"] * len(self.dbcolumnnames))
		}[connection_db_type(self.connection)]
		self.append_statement = "INSERT INTO %s (%s) VALUES (%s)" % (self.dbtablename, ",".join(self.dbcolumnnames), params)

	def _end_of_rows(self):
		# FIXME:  is this needed?
		table.Table._end_of_rows(self)
		self.connection.commit()

	def sync_next_id(self):
		if self.next_id is not None:
			max_id = idmap_get_max_id(self.connection, type(self.next_id))
			if max_id is None:
				self.set_next_id(type(self.next_id)(0))
			else:
				self.set_next_id(max_id + 1)
		return self.next_id

	def maxrowid(self):
		self.cursor.execute("SELECT MAX(ROWID) FROM %s" % self.dbtablename)
		return self.cursor.fetchone()[0]

	def __len__(self):
		self.cursor.execute("SELECT COUNT(*) FROM %s" % self.dbtablename)
		return self.cursor.fetchone()[0]

	def __iter__(self):
		cursor = self.connection.cursor()
		cursor.execute("SELECT * FROM %s" % self.dbtablename)
		for values in cursor:
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
		SQLite does not permit the PRIMARY KEY of a row to be
		modified, so it needs to be done prior to insertion.  This
		method is intended for internal use only.
		"""
		if self.next_id is not None:
			# assign (and record) a new ID before inserting the
			# row to avoid collisions with existing rows
			setattr(row, self.next_id.column_name, idmap_get_new(self.connection, getattr(row, self.next_id.column_name), self))
		# FIXME: in Python 2.5 use attrgetter() for attribute
		# tuplization.
		self.cursor.execute(self.append_statement, map(lambda n: getattr(row, n), self.dbcolumnnames))

	append = _append

	def _row_from_cols(self, values):
		"""
		Given an iterable of values in the order of columns in the
		database, construct and return a row object.  This is a
		convenience function for turning the results of database
		queries into Python objects.
		"""
		row = self.RowType()
		for c, v in zip(self.dbcolumnnames, values):
			setattr(row, c, v)
		return row

	def unlink(self):
		table.Table.unlink(self)
		self.connection = None
		self.cursor = None

	def applyKeyMapping(self):
		"""
		Used as the second half of the key reassignment algorithm.
		Loops over each row in the table, replacing references to
		old row keys with the new values from the _idmap_ table.
		"""
		assignments = ", ".join("%s = (SELECT new FROM _idmap_ WHERE old == %s)" % (colname, colname) for coltype, colname in zip(self.dbcolumntypes, self.dbcolumnnames) if coltype in ligolwtypes.IDTypes and (self.next_id is None or colname != self.next_id.column_name))
		if assignments:
			# SQLite documentation says ROWID is monotonically
			# increasing starting at 1 for the first row unless
			# it ever wraps around, then it is randomly
			# assigned.  ROWID is a 64 bit integer, so the only
			# way it will wrap is if somebody sets it to a very
			# high number manually.  This library does not do
			# that, so I don't bother checking.
			self.cursor.execute("UPDATE %s SET %s WHERE ROWID > %d" % (self.dbtablename, assignments, self.last_maxrowid))
			self.last_maxrowid = self.maxrowid() or 0


#
# =============================================================================
#
#                                  LSC Tables
#
# =============================================================================
#


class ProcessTable(DBTable):
	# FIXME:  remove this class
	tableName = lsctables.ProcessTable.tableName
	validcolumns = lsctables.ProcessTable.validcolumns
	constraints = lsctables.ProcessTable.constraints
	next_id = lsctables.ProcessTable.next_id
	RowType = lsctables.ProcessTable.RowType
	how_to_index = lsctables.ProcessTable.how_to_index

	def get_ids_by_program(self, program):
		"""
		Return a set of the process IDs from rows whose program
		string equals the given program.
		"""
		return set(id for (id,) in self.cursor.execute("SELECT process_id FROM process WHERE program == ?", (program,)))


class ProcessParamsTable(DBTable):
	tableName = lsctables.ProcessParamsTable.tableName
	validcolumns = lsctables.ProcessParamsTable.validcolumns
	constraints = lsctables.ProcessParamsTable.constraints
	next_id = lsctables.ProcessParamsTable.next_id
	RowType = lsctables.ProcessParamsTable.RowType
	how_to_index = lsctables.ProcessParamsTable.how_to_index

	def append(self, row):
		if row.type is not None and row.type not in ligolwtypes.Types:
			raise ligolw.ElementError, "unrecognized type '%s'" % row.type
		DBTable.append(self, row)


class TimeSlideTable(DBTable):
	tableName = lsctables.TimeSlideTable.tableName
	validcolumns = lsctables.TimeSlideTable.validcolumns
	constraints = lsctables.TimeSlideTable.constraints
	next_id = lsctables.TimeSlideTable.next_id
	RowType = lsctables.TimeSlideTable.RowType
	how_to_index = lsctables.TimeSlideTable.how_to_index

	def __len__(self):
		raise NotImplementedError

	def __getitem__(*args):
		raise NotImplementedError

	def get_offset_dict(self, id):
		offsets = dict(self.cursor.execute("SELECT instrument, offset FROM time_slide WHERE time_slide_id == ?", (id,)))
		if not offsets:
			raise KeyError, id
		return offsets

	def as_dict(self):
		"""
		Return a ditionary mapping time slide IDs to offset
		dictionaries.
		"""
		d = {}
		for id, instrument, offset in self.cursor.execute("SELECT time_slide_id, instrument, offset FROM time_slide"):
			if id not in d:
				d[id] = {}
			d[id][instrument] = offset
		return d

	def get_time_slide_id(self, offsetdict, create_new = None):
		"""
		Return the time_slide_id corresponding to the time slide
		described by offsetdict, a dictionary of instrument/offset
		pairs.  If no matching time_slide_id is found, then
		KeyError is raised.  If, however, the optional create_new
		argument is set to an lsctables.Process object (or any
		other object with a process_id attribute), then new rows
		are added to the table to describe the desired time slide,
		and the ID of the new rows is returned.
		"""
		# look for the ID
		for id, slide in self.as_dict().iteritems():
			if offsetdict == slide:
				# found it
				return id

		# time slide not found in table
		if create_new is None:
			raise KeyError, offsetdict
		self.sync_next_id()
		id = self.get_next_id()
		for instrument, offset in offsetdict.iteritems():
			row = self.RowType()
			row.process_id = create_new.process_id
			row.time_slide_id = id
			row.instrument = instrument
			row.offset = offset
			self.append(row)

		# return new ID
		return id

	def iterkeys(self):
		raise NotImplementedError


#
# =============================================================================
#
#                                Table Metadata
#
# =============================================================================
#


def build_indexes(connection, verbose = False):
	"""
	Using the how_to_index annotations in the table class definitions,
	construct a set of indexes for the database at the given
	connection.
	"""
	cursor = connection.cursor()
	for table_name in get_table_names(connection):
		# FIXME:  figure out how to do this extensibly
		if table_name in TableByName:
			how_to_index = TableByName[table_name].how_to_index
		elif table_name in lsctables.TableByName:
			how_to_index = lsctables.TableByName[table_name].how_to_index
		else:
			continue
		if how_to_index is not None:
			if verbose:
				print >>sys.stderr, "indexing %s table ..." % table_name
			for index_name, cols in how_to_index.iteritems():
				cursor.execute("CREATE INDEX IF NOT EXISTS %s ON %s (%s)" % (index_name, table_name, ",".join(cols)))
	connection.commit()


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
	table.StripTableName(TimeSlideTable.tableName): TimeSlideTable
}


#
# The database-backed table implementation requires there to be no more
# than one table of each name in the document.  Some documents require
# multiple tables with the same name, and those tables cannot be stored in
# the database.  Use this list to set which tables are not to be stored in
# the database.
#


NonDBTableNames = []


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


__parent_startTable = ligolw.LIGOLWContentHandler.startTable


def startTable(self, attrs):
	name = table.StripTableName(attrs[u"Name"])
	if name in map(table.StripTableName, NonDBTableNames):
		return __parent_startTable(self, attrs)
	if name in TableByName:
		return TableByName[name](attrs)
	return DBTable(attrs)


ligolw.LIGOLWContentHandler.startTable = startTable
