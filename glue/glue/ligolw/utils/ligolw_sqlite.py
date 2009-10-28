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


# so they can be inserted into a database
dbtables.ligolwtypes.ToPyType["ilwd:char"] = unicode


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


def setup(target):
	connection = sqlite3.connect(target)
	dbtables.DBTable_set_connection(connection)

	for tbl in dbtables.get_xml(connection).getElementsByTagName(ligolw.Table.tagName):
		tbl.sync_next_id()

	return connection


#
# How to insert
#


def insert(connection, urls, preserve_ids = False, verbose = False):
	"""
	Iterate over a sequence of URLs and parse and insert each one into
	the database the dbtables.DBTable class is currently connected to.
	"""
	if not preserve_ids:
		# enable ID remapping
		dbtables.idmap_create(connection)
		dbtables.DBTable.append = dbtables.DBTable._remapping_append
	for n, url in enumerate(urls):
		# load document (if enabled, row IDs are reassigned on
		# input)
		if verbose:
			print >>sys.stderr, "%d/%d:" % (n + 1, len(urls)),
		xmldoc = utils.load_url(url, verbose = verbose, gz = (url or "stdin").endswith(".gz"))

		# update references to row IDs
		if not preserve_ids:
			table_elems = xmldoc.getElementsByTagName(ligolw.Table.tagName)
			for i, tbl in enumerate(table_elems):
				if verbose:
					print >>sys.stderr, "updating IDs: %d%%\r" % (100 * i / len(table_elems)),
				tbl.applyKeyMapping()
			if verbose:
				print >>sys.stderr, "updating IDs: 100%"

			# reset ID mapping for next document
			dbtables.idmap_reset(connection)

		# delete cursors
		xmldoc.unlink()
	connection.commit()

	dbtables.build_indexes(connection, verbose)


#
# How to extract
#


def extract(connection, filename, table_names = None, verbose = False, xsl_file = None):
	xmldoc = ligolw.Document()
	xmldoc.appendChild(dbtables.get_xml(connection, table_names))
	utils.write_filename(xmldoc, filename, gz = (filename or "stdout").endswith(".gz"), verbose = verbose, xsl_file = xsl_file)

	# delete cursors
	xmldoc.unlink()
