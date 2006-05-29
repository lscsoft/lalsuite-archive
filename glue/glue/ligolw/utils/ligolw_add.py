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
Add (merge) LIGO LW XML files containing LSC tables.
"""

import os
import sys
import urllib
from urlparse import urlparse

from glue.ligolw import ligolw
from glue.ligolw import metaio
from glue.ligolw import lsctables

__author__ = "Kipp Cannon <kipp@gravity.phys.uwm.edu>"
__date__ = "$Date$"[7:-2]
__version__ = "$Revision$"[11:-2]


#
# =============================================================================
#
#                                    Input
#
# =============================================================================
#

ContentHandler = ligolw.LIGOLWContentHandler

def append_document(doc, file):
	"""
	Parse the contents of the file object file, appending to the
	document tree doc.
	"""
	ligolw.make_parser(ContentHandler(doc)).parse(file)
	return doc


def url2path(url):
	"""
	If url identifies a file on the local host, return the path to the
	file otherwise raise ValueError.
	"""
	(scheme, location, path, params, query, frag) = urlparse(url)
	if scheme.lower() in ("", "file") and location.lower() in ("", "localhost"):
		return path
	raise ValueError, url


def remove_input(urls, preserves, verbose = False):
	"""
	Attempt to delete all files identified by the URLs in urls except
	any that are the same as the files in the preserves list.
	"""
	for path in map(url2path, urls):
		if True in map(os.path.samefile, [path] * len(preserves), preserves):
			continue
		if verbose:
			print >>sys.stderr, "removing %s ..." % path
		try:
			os.remove(path)
		except:
			pass


#
# =============================================================================
#
#                                Document Merge
#
# =============================================================================
#

def reassign_ids(doc):
	"""
	Reassign IDs to all rows in all LSC tables in doc so that there are
	no collisions when the LIGO_LW elements are merged.
	"""
	ilwditers = {}
	for tablename, ilwdclass in lsctables.ILWDGeneratorByTableName.iteritems():
		ilwditers[tablename] = ilwdclass()
	for elem in doc.getElementsByTagName(ligolw.LIGO_LW.tagName):
		lsctables.makeReference(elem)
		lsctables.NewIDs(elem, ilwditers)
		lsctables.deReference(elem)
	return doc


def merge_elements(elem1, elem2):
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
		raise ligolw.ElementError, "merge_elements(): elements must have same names"
	if elem1.tagName == ligolw.LIGO_LW.tagName:
		# copy children;  LIGO_LW elements have no attributes
		map(elem1.appendChild, elem2.childNodes)
	elif elem1.tagName == ligolw.Table.tagName:
		# copy rows
		elem1.rows.extend(elem2.rows)
	else:
		raise ligolw.ElementError, "merge_elements(): can't merge %s elements." % elem1.tagName
	if elem2.parentNode:
		elem2.parentNode.removeChild(elem2)
	return elem1


def tables_can_be_merged(a, b):
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


def merge_compatible_tables(elem):
	"""
	Below the given element, find all Tables whose structure is
	described in lsctables, and merge compatible ones of like type.
	That is, merge all SnglBurstTables that have the same columns into
	a single table, etc..
	"""
	for tname in lsctables.TableByName.keys():
		tables = metaio.getTablesByName(elem, tname)
		for i in range(1, len(tables)):
			if tables_can_be_merged(tables[0], tables[i]):
				merge_elements(tables[0], tables[i])
	return elem


def merge_all_elements(doc):
	"""
	Combine the child elements of all top-level LIGO_LW elements under
	a single LIGO_LW element, then combine equivalent tables into
	single tables.
	"""
	# LIGO_LW elements
	reduce(merge_elements, doc.getElementsByTagName(ligolw.LIGO_LW.tagName))

	# Table elements
	merge_compatible_tables(doc)

	return doc


#
# =============================================================================
#
#                                 Library API
#
# =============================================================================
#

def ligolw_add(doc, urls, **kwargs):
	"""
	An implementation of the LIGO LW add algorithm.  urls is a list of
	URLs to load, doc is the XML document tree to which they should be
	added.
	"""
	# Input
	for n, url in enumerate(urls):
		if kwargs["verbose"]:
			print >>sys.stderr, "loading %d/%d: %s" % (n + 1, len(urls), url)
		append_document(doc, urllib.urlopen(url))

	# ID reassignment
	if not kwargs["non_lsc_tables_ok"] and lsctables.HasNonLSCTables(doc):
		raise ValueError, "non-LSC tables found.  Use --non-lsc-tables-ok to force"
	if kwargs["verbose"]:
		print >>sys.stderr, "reasigning row IDs ..."
	reassign_ids(doc)

	# Document merge
	if kwargs["verbose"]:
		print >>sys.stderr, "merging elements ..."
	merge_all_elements(doc)

	return doc
