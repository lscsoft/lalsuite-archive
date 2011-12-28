# Copyright (C) 2006  Kipp Cannon
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


try:
	import sqlite3
except ImportError:
	# pre 2.5.x
	from pysqlite2 import dbapi2 as sqlite3
import sys


from glue import git_version
from glue.ligolw import ligolw
from glue.ligolw import dbtables
from glue.ligolw import utils


# FIXME: remove this hack when the SnglInspiralTable class uses the
# standard ID generator by default.
dbtables.lsctables.SnglInspiralTable.next_id = dbtables.lsctables.SnglInspiralID(0)


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date


#
# =============================================================================
#
#                                 Library Code
#
# =============================================================================
#


#
# Open database
#


def setup(target, check_same_thread=True):
	connection = sqlite3.connect(target, check_same_thread=check_same_thread)
	dbtables.DBTable_set_connection(connection)

	dbtables.idmap_sync(connection)

	return connection


#
# How to insert
#


def update_ids(connection, xmldoc=None, verbose = False):
	"""
	For internal use only.
	"""
	if xmldoc: table_elems = xmldoc.getElementsByTagName(ligolw.Table.tagName)
	else: table_elems = dbtables.get_xml(connection).getElementsByTagName(ligolw.Table.tagName)
	for i, tbl in enumerate(table_elems):
		if verbose:
			print >>sys.stderr, "updating IDs: %d%%\r" % (100.0 * i / len(table_elems)),
		tbl.applyKeyMapping()
		tbl.unlink()
	if verbose:
		print >>sys.stderr, "updating IDs: 100%"

	# reset ID mapping for next document
	dbtables.idmap_reset(connection)


def insert_from_url(connection, url, preserve_ids = False, verbose = False):
	"""
	Parse and insert the LIGO Light Weight document at the URL into the
	database the at the given connection.
	"""
	#
	# load document.  this process inserts the document's contents into
	# the database.  the document is unlinked to delete database cursor
	# objects it retains
	#

	xmldoc = utils.load_url(url, verbose = verbose)

	#
	# update references to row IDs
	#

	if not preserve_ids:
		update_ids(connection, xmldoc, verbose)

	xmldoc.unlink()


def insert_from_xmldoc(connection, xmldoc, preserve_ids = False, verbose = False):
	"""
	Insert the tables from an in-ram XML document into the database at
	the given connection.
	"""
	#
	# iterate over tables in the XML tree, reconstructing each inside
	# the database
	#

	for tbl in xmldoc.getElementsByTagName(ligolw.Table.tagName):
		#
		# instantiate the correct table class
		#

		name = dbtables.table.StripTableName(tbl.getAttribute("Name"))
		if name in dbtables.TableByName:
			dbtab = dbtables.TableByName[name](tbl.attributes, connection = connection)
		else:
			dbtab = dbtables.DBTable(tbl.attributes, connection = connection)

		#
		# copy table element child nodes from source XML tree
		#

		for elem in tbl.childNodes:
			if elem.tagName == dbtables.table.TableStream.tagName:
				dbtab._end_of_columns()
			dbtab.appendChild(type(elem)(elem.attributes))

		#
		# copy table rows from source XML tree
		#

		for row in tbl:
			dbtab.append(row)
		dbtab._end_of_rows()

		#
		# unlink to delete cursor objects
		#

		dbtab.unlink()

	#
	# update references to row IDs
	#
	connection.commit()
	if not preserve_ids:
		update_ids(connection, None, verbose)


def insert_from_urls(connection, urls, preserve_ids = False, verbose = False):
	"""
	Iterate over a sequence of URLs, calling insert_from_url() on each,
	then build the indexes indicated by the metadata in lsctables.py.
	"""
	#
	# enable/disable ID remapping
	#

	orig_DBTable_append = dbtables.DBTable.append
	if not preserve_ids:
		dbtables.idmap_create(connection)
		dbtables.DBTable.append = dbtables.DBTable._remapping_append
	else:
		dbtables.DBTable.append = dbtables.DBTable._append

	#
	# load documents
	#

	for n, url in enumerate(urls):
		if verbose:
			print >>sys.stderr, "%d/%d:" % (n + 1, len(urls)),
		insert_from_url(connection, url, preserve_ids = preserve_ids, verbose = verbose)
	connection.commit()

	#
	# done.  build indexes, restore original .append() method
	#

	dbtables.build_indexes(connection, verbose)
	dbtables.DBTable.append = orig_DBTable_append


#
# How to extract
#


def extract(connection, filename, table_names = None, verbose = False, xsl_file = None):
	xmldoc = ligolw.Document()
	xmldoc.appendChild(dbtables.get_xml(connection, table_names))
	utils.write_filename(xmldoc, filename, gz = (filename or "stdout").endswith(".gz"), verbose = verbose, xsl_file = xsl_file)

	# delete cursors
	xmldoc.unlink()
