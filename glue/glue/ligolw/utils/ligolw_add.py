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
from glue.ligolw import lsctables
from glue.ligolw import docutils

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

def append_document(doc, file):
	"""
	Parse the contents of the file object file, appending to the
	document tree doc.
	"""
	ligolw.make_parser(ligolw.LIGOLWContentHandler(doc)).parse(file)
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


def remove_input(urls, output, verbose = False):
	"""
	Attempt to delete all files identified by the URLs in urls except
	any that are the same as the file whose name is output.
	"""
	for path in map(url2path, urls):
		if output and os.path.samefile(path, output):
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
		docutils.makeReference(elem)
		docutils.NewIDs(elem, ilwditers)
		docutils.deReference(elem)
	return doc


def element_merge(doc):
	"""
	Combine the child elements of all top-level LIGO_LW elements under
	a single LIGO_LW element, then combine equivalent tables into
	single tables.
	"""
	# LIGO_LW elements
	reduce(docutils.MergeElements, doc.getElementsByTagName(ligolw.LIGO_LW.tagName))

	# Table elements
	docutils.MergeCompatibleTables(doc)

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
	if not kwargs["non_lsc_tables_ok"] and docutils.HasNonLSCTables(doc):
		raise ValueError, "non-LSC tables found.  Use --non-lsc-tables-ok to force"
	if kwargs["verbose"]:
		print >>sys.stderr, "reasigning row IDs ..."
	reassign_ids(doc)

	# Document merge
	if kwargs["verbose"]:
		print >>sys.stderr, "merging elements ..."
	element_merge(doc)

	return doc
