# $Id$
#
# Copyright (C) 2008  Kipp C. Cannon
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
Ask Kipp to document this!
"""


from glue import iterutils
from glue import segments
from glue import segmentsUtils
from glue.ligolw import table
from glue.ligolw import lsctables
from pylal.date import LIGOTimeGPS


__author__ = "Kipp Cannon <kipp@gravity.phys.uwm.edu>"
__version__ = "$Revision$"[11:-2]
__date__ = "$Date$"[7:-2]


#
# =============================================================================
#
#                                 Segment List
#
# =============================================================================
#


def segment_def(instruments, name, comment, id, process):
	row = lsctables.SegmentDef()
	row.process_id = process.process_id
	row.segment_def_id = id
	instruments = list(instruments)
	instruments.sort()
	row.ifos = ",".join(instruments)
	row.name = name
	row.comment = comment
	return row


def seg_def_map(seg, seg_def, tbl, process):
	row = lsctables.SegmentDefMap()
	row.process_id = process.process_id
	row.segment_id = seg.segment_id
	row.segment_def_id = seg_def.segment_def_id
	row.seg_def_map_id = tbl.get_next_id()
	return row


class LigolwSegmentList(object):
	"""
	A description of a class of segments.
	"""
	def __init__(self, active = (), inactive = (), unknown = (), instruments = set(), name = None, comment = None):
		self.active = segments.segmentlist(active)
		self.inactive = segments.segmentlist(inactive)
		self.unknown = segments.segmentlist(unknown)
		self.instruments = instruments
		self.name = name
		self.comment = comment

	def sort(self, *args):
		"""
		Sort the three internal segment lists in time order.
		"""
		self.active.sort(*args)
		self.inactive.sort(*args)
		self.unknown.sort(*args)

	def coalesce(self):
		"""
		Coalesce the three internal segment lists.
		"""
		self.active.coalesce()
		self.inactive.coalesce()
		self.unknown.coalesce()

	def time_order_rows(self, tbl, process):
		"""
		A generator, yielding in time order
		glue.ligolw.lsctables.Segment objects representing all the
		segments from the three internal segment lists.
		"""
		def seglist_iterator(seglist, activity):
			for seg in seglist:
				yield seg, activity

		self.sort()

		for seg, activity in iterutils.inorder(seglist_iterator(self.active, True), seglist_iterator(self.inactive, False), seglist_iterator(self.unknown, None)):
			row = lsctables.Segment()
			row.set(seg)
			row.set_active(activity)
			row.process_id = process.process_id
			row.segment_id = tbl.get_next_id()
			yield row


#
# =============================================================================
#
#                                 Library API
#
# =============================================================================
#


class LigolwSegments(object):
	"""
	A high-level interface to the segments tables in a LIGO Light
	Weight XML document.
	"""
	def __init__(self, xmldoc, process):
		#
		# The row in the process table on which we will blame our
		# work
		#

		self.process = process

		#
		# Find tables, and synchronize ID generators
		#

		try:
			self.segment_def_table = table.get_table(xmldoc, lsctables.SegmentDefTable.tableName)
		except ValueError:
			self.segment_def_table = xmldoc.childNodes[0].appendChild(lsctables.New(lsctables.SegmentDefTable, ("process_id", "segment_def_id", "ifos", "name", "comment")))
		self.segment_def_table.sync_next_id()

		try:
			self.segment_table = table.get_table(xmldoc, lsctables.SegmentTable.tableName)
		except ValueError:
			self.segment_table = xmldoc.childNodes[0].appendChild(lsctables.New(lsctables.SegmentTable, ("process_id", "segment_id", "start_time", "start_time_ns", "end_time", "end_time_ns", "active")))
		self.segment_table.sync_next_id()

		try:
			self.segment_def_map_table = table.get_table(xmldoc, lsctables.SegmentDefMapTable.tableName)
		except ValueError:
			self.segment_def_map_table = xmldoc.childNodes[0].appendChild(lsctables.New(lsctables.SegmentDefMapTable, ("process_id", "segment_id", "segment_def_id", "seg_def_map_id")))
		self.segment_def_map_table.sync_next_id()

		#
		# Transform segment tables into a collection of
		# LigolwSegmentList objects for more convenient
		# manipulation
		#

		# segment_def_id --> LigolwSegmentList object mapping
		self.segment_lists = {}
		for row in self.segment_def_table:
			if row.segment_def_id in self.segment_lists:
				raise ValueError, "segment_definer table contains duplicate segment_definer keys"
			self.segment_lists[row.segment_def_id] = LigolwSegmentList(instruments = row.get_ifos(), name = row.name, comment = row.comment)
		del self.segment_def_table[:]

		# segment_id --> LigolwSegmentList object mapping
		segment_def_map = {}
		for row in self.segment_def_map_table:
			if row.segment_id in segment_def_map:
				raise ValueError, "segment_def_map table contains duplicate segment keys"
			segment_def_map[row.segment_id] = self.segment_lists[row.segment_def_id]
		del self.segment_def_map_table[:]

		# populate LigolwSegmentList objects from segment table
		for row in self.segment_table:
			active = row.get_active()
			if active is True:
				segment_def_map[row.segment_id].active.append(row.get())
			elif row.active is False:
				segment_def_map[row.segment_id].inactive.append(row.get())
			elif row.active is None:
				segment_def_map[row.segment_id].unknown.append(row.get())
			else:
				raise ValueError, "execution should not get to this line of code"
		del self.segment_table[:]
		del segment_def_map

		#
		# Done
		#


	def coalesce(self):
		#
		# Coalesce the segment lists.
		#

		for ligolw_segment_list in self.segment_lists.values():
			ligolw_segment_list.coalesce()


	def sort(self):
		for ligolw_segment_list in self.segment_lists.values():
			ligolw_segment_list.sort()


	def optimize(self):
		"""
		Identifies segment lists that differ only in their
		instruments --- they have the same active, inactive, and
		unknown segments, the same name and the same comment ---
		and then deletes all but one of them, leaving just a single
		list having the union of the instruments.
		"""
		self.sort()
		mergelist = []
		for (def_id_a, seglist_a), (def_id_b, seglist_b) in iterutils.choices(self.segment_lists.items(), 2):
			if seglist_a.active == seglist_b.active and seglist_a.inactive == seglist_b.inactive and seglist_a.unknown == seglist_b.unknown and seglist_a.name == seglist_b.name and seglist_a.comment == seglist_b.comment:
				mergelist.append((def_id_a, def_id_b))
		for target, source in mergelist:
			if source in self.segment_lists:
				self.segment_lists[target].instruments |= self.segment_lists.pop(source).instruments


	def finalize(self):
		#
		# Restore the LigolwSegmentList objects to the XML tables
		# in preparation for output.  All segments from all segment
		# lists are inserted into the tables in time order, but
		# this is NOT behaviour external applications should rely
		# on.  This is done simply in the belief that it might
		# assist in constructing well balanced indexed databases
		# from the resulting files.  If that proves not to be the
		# case, or for some reason this behaviour proves
		# inconvenient to preserve, then it might be discontinued
		# without notice.  You've been warned.
		#

		def row_iterator(rows, segment_definition):
			for row in rows:
				yield row, segment_definition

		row_iterators = []
		for segment_def_id, ligolw_segment_list in self.segment_lists.items():
			segment_definition = segment_def(ligolw_segment_list.instruments, ligolw_segment_list.name, ligolw_segment_list.comment, segment_def_id, self.process)
			self.segment_def_table.append(segment_definition)
			row_iterators.append(row_iterator(ligolw_segment_list.time_order_rows(self.segment_table, self.process), segment_definition))

		for row, segment_definition in iterutils.inorder(*row_iterators):
			self.segment_table.append(row)
			self.segment_def_map_table.append(seg_def_map(row, segment_definition, self.segment_def_map_table, self.process))

		self.segment_lists.clear()


	def insert_from_segwizard(self, fileobj, instruments, name, comment):
		self.segment_lists[self.segment_def_table.get_next_id()] = LigolwSegmentList(active = segmentsUtils.fromsegwizard(fileobj, coltype = LIGOTimeGPS), instruments = instruments, name = name, comment = comment)
