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
High-level document manipulation utilities.
"""

__author__ = "Kipp Cannon <kipp@gravity.phys.uwm.edu>"
__date__ = "$Date$"[7:-2]
__version__ = "$Revision$"[11:-2]

import copy

import ligolw
import metaio
import lsctables

#
# =============================================================================
#
#                                    Tables
#
# =============================================================================
#

def HasNonLSCTables(elem):
	"""
	Return True if the document tree below elem contains non-LSC
	tables, otherwise return False.
	"""
	for table in elem.getElementsByTagName(ligolw.Table.tagName):
		if metaio.StripTableName(table.getAttribute("Name")) not in lsctables.TableByName.keys():
			return True
	return False


#
# Row cross-reference utilities.
#

def NewILWDs(table, column_name):
	"""
	From the LSC table, return a compatible ILWD iterator object,
	initialized to the next unique ID following those found in the
	table.
	"""
	try:
		n = max(map(lsctables.ILWDID, table.dict.keys()))
	except ValueError:
		n = -1
	return lsctables.ILWD(metaio.StripTableName(table.getAttribute("Name")), column_name, n + 1)


def makeReference(elem):
	"""
	Run the makeReference() method on all LSC tables below elem,
	constructing references to other tables under elem.
	"""
	for table in lsctables.getLSCTables(elem):
		try:
			table.makeReference(elem)
		except AttributeError:
			# table is missing a cross-reference column.
			pass


def deReference(elem):
	"""
	Run the deReference() method on all LSC tables below elem.
	"""
	for table in lsctables.getLSCTables(elem):
		try:
			table.deReference()
		except AttributeError:
			# table is missing a cross-reference column.
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
		for table in metaio.getTablesByName(elem, tablename):
			keymap = {}
			try:
				for oldkey in table.dict.keys():
					keymap[oldkey] = ilwditer.next()
			except AttributeError:
				# table is missing its ID column
				continue
			if not len(keymap):
				continue
			for row in table.rows:
				try:
					row._set_key(keymap[row._get_key()])
				except KeyError:
					# row has a key not listed in the
					# table's dictionary.
					pass


#
# =============================================================================
#
#                                     I/O
#
# =============================================================================
#

class PartialLIGOLWContentHandler(ligolw.LIGOLWContentHandler):
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
		ligolw.LIGOLWContentHandler.__init__(self, document)
		self.filter = filter
		self.filtered_depth = 0

	def startElement(self, name, attrs):
		if self.filtered_depth > 0 or self.filter(name, attrs):
			ligolw.LIGOLWContentHandler.startElement(self, name, attrs)
			self.filtered_depth += 1

	def endElement(self, name):
		if self.filtered_depth > 0:
			self.filtered_depth -= 1
			ligolw.LIGOLWContentHandler.endElement(self, name)
