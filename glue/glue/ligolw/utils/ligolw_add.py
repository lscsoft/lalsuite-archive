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
	ligolw.make_parser(lsctables.LIGOLWContentHandler(doc)).parse(file)
	return doc


def url2path(url):
	(scheme, location, path, params, query, frag) = urlparse(url)
	if scheme.lower() in ("", "file") and location.lower() in ("", "localhost"):
		return path
	raise ValueError, url


def remove_input(urls, output, verbose = False):
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
	URLs to load.
	"""
	# Input
	for n, url in enumerate(urls):
		if kwargs["verbose"]:
			print >>sys.stderr, "loading %d/%d: %s" % (n + 1, len(urls), url)
		append_document(doc, urllib.urlopen(url))

	# ID reassignment
	if not kwargs["non_lsc_tables_ok"] and docutils.HasNonLSCTables(doc):
		print >>sys.stderr, "error:  non-LSC tables found.  Use --non-lsc-tables-ok to force"
		sys.exit(1)
	if kwargs["verbose"]:
		print >>sys.stderr, "reasigning row IDs ..."
	reassign_ids(doc)

	# Document merge
	if kwargs["verbose"]:
		print >>sys.stderr, "merging elements ..."
	element_merge(doc)

	return doc
