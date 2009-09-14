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
from glue.lal import LIGOTimeGPS
from glue.ligolw import table
from glue.ligolw import lsctables


__author__ = "Kipp Cannon <kcannon@ligo.caltech.edu>"
__version__ = "$Revision$"[11:-2]
__date__ = "$Date$"[7:-2]


#
# =============================================================================
#
#                                 Segment List
#
# =============================================================================
#


class LigolwSegmentList(object):
	"""
	A description of a class of segments.
	"""
	def __init__(self, active = (), valid = (), instruments = set(), name = None, comment = None):
		self.valid = segments.segmentlist(valid)
		self.active = segments.segmentlist(active)
		self.instruments = instruments
		self.name = name
		self.comment = comment

	def sort(self, *args):
		"""
		Sort the internal segment lists.  The optional args are
		passed to the .sort() method of the segment lists.  This
		can be used to control the sort order by providing an
		alternate comparison function.  The default is to sort by
		start time with ties broken by end time.
		"""
		self.valid.sort(*args)
		self.active.sort(*args)

	def coalesce(self):
		"""
		Coalesce the internal segment lists.
		"""
		self.valid.coalesce()
		self.active.coalesce()


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
	def __init__(self, xmldoc):
		#
		# Find tables
		#

		try:
			self.segment_def_table = table.get_table(xmldoc, lsctables.SegmentDefTable.tableName)
		except ValueError:
			self.segment_def_table = xmldoc.childNodes[0].appendChild(lsctables.New(lsctables.SegmentDefTable, ("process_id", "segment_def_id", "ifos", "name", "comment")))

		try:
			self.segment_sum_table = table.get_table(xmldoc, lsctables.SegmentSumTable.tableName)
		except ValueError:
			self.segment_sum_table = xmldoc.childNodes[0].appendChild(lsctables.New(lsctables.SegmentSumTable, ("process_id", "segment_sum_id", "start_time", "start_time_ns", "end_time", "end_time_ns", "segment_def_id", "comment")))

		try:
			self.segment_table = table.get_table(xmldoc, lsctables.SegmentTable.tableName)
		except ValueError:
			self.segment_table = xmldoc.childNodes[0].appendChild(lsctables.New(lsctables.SegmentTable, ("process_id", "segment_id", "start_time", "start_time_ns", "end_time", "end_time_ns", "segment_def_id")))

		#
		# Transform segment tables into a collection of
		# LigolwSegmentList objects for more convenient
		# manipulation
		#

		# construct empty LigolwSegmentList objects, one for each
		# entry in the segment_definer table, indexed by
		# segment_definer id
		self.segment_lists = dict((row.segment_def_id, LigolwSegmentList(instruments = row.get_ifos(), name = row.name, comment = row.comment)) for row in self.segment_def_table)
		if len(self.segment_lists) != len(self.segment_def_table):
			raise ValueError, "duplicate segment_definer IDs detected in segment_definer table"
		del self.segment_def_table[:]

		# populate LigolwSegmentList objects from segment table and
		# segment_summary table
		for row in self.segment_sum_table:
			self.segment_lists[row.segment_def_id].valid.append(row.get())
		del self.segment_sum_table[:]
		for row in self.segment_table:
			self.segment_lists[row.segment_def_id].active.append(row.get())
		del self.segment_table[:]

		# replace segment_lists dictionary with a list because the
		# segment_definer IDs no longer have any meaning
		self.segment_lists = self.segment_lists.values()

		#
		# Synchronize ID generators
		#

		# FIXME:  why am I doing this!?  I've just deleted all the
		# rows.
		self.segment_def_table.sync_next_id()
		self.segment_table.sync_next_id()
		self.segment_sum_table.sync_next_id()

		#
		# Done
		#


	def coalesce(self):
		"""
		Coalesce the segment lists.
		"""
		for ligolw_segment_list in self.segment_lists:
			ligolw_segment_list.coalesce()


	def sort(self, *args):
		"""
		Sort the segment lists.  The optional args are passed to
		the .sort() methods of the segment lists.  This can be used
		to control the sort order by providing an alternate
		comparison function (the default is to sort all lists by
		segment start time with ties broken by end time).
		"""
		for ligolw_segment_list in self.segment_lists:
			ligolw_segment_list.sort(*args)


	def optimize(self):
		"""
		Identifies segment lists that differ only in their
		instruments --- they have the same valid and active
		segments, the same name and the same comment --- and then
		deletes all but one of them, leaving just a single list
		having the union of the instruments.
		"""
		self.sort()
		segment_lists = dict(enumerate(self.segment_lists))
		for target, source in [(idx_a, idx_b) for (idx_a, seglist_a), (idx_b, seglist_b) in iterutils.choices(segment_lists.items(), 2) if seglist_a.valid == seglist_b.valid and seglist_a.active == seglist_b.active and seglist_a.name == seglist_b.name and seglist_a.comment == seglist_b.comment]:
			try:
				segment_lists[target].instruments |= segment_lists.pop(source).instruments
			except KeyError:
				pass
		self.segment_lists[:] = segment_lists.values()


	def finalize(self, process_row):
		"""
		Restore the LigolwSegmentList objects to the XML tables in
		preparation for output.  All segments from all segment
		lists are inserted into the tables in time order, but this
		is NOT behaviour external applications should rely on.
		This is done simply in the belief that it might assist in
		constructing well balanced indexed databases from the
		resulting files.  If that proves not to be the case, or for
		some reason this behaviour proves inconvenient to preserve,
		then it might be discontinued without notice.  You've been
		warned.
		"""
		#
		# put all segment lists in time order
		#

		self.sort()

		#
		# generator function to return active and inactive segment
		# table rows from a LigolwSegmentList object in time order,
		# each paired with a handle to the matching row in the
		# segment_definer table
		#

		def row_generator(segs, target_table, process_row, segment_def_row):
			id_column = target_table.next_id.column_name
			for seg in segs:
				row = target_table.RowType()
				row.set(seg)
				row.process_id = process_row.process_id
				setattr(row, id_column, target_table.get_next_id())
				row.segment_def_id = segment_def_row.segment_def_id
				if hasattr(row, "comment"):
					row.comment = None
				yield row, target_table

		#
		# populate the segment_definer table from the list of
		# LigolwSegmentList objects and construct a matching list
		# of table row generators
		#

		row_generators = []
		for ligolw_segment_list in self.segment_lists:
			segment_def_row = self.segment_def_table.RowType()
			segment_def_row.process_id = process_row.process_id
			segment_def_row.segment_def_id = self.segment_def_table.get_next_id()
			segment_def_row.set_ifos(ligolw_segment_list.instruments)
			segment_def_row.name = ligolw_segment_list.name
			segment_def_row.comment = ligolw_segment_list.comment
			self.segment_def_table.append(segment_def_row)

			row_generators.append(row_generator(ligolw_segment_list.valid, self.segment_sum_table, process_row, segment_def_row))
			row_generators.append(row_generator(ligolw_segment_list.active, self.segment_table, process_row, segment_def_row))

		#
		# populate segment and segment_summary tables by pulling
		# rows from the generators in time order
		#

		for row, target_table in iterutils.inorder(*row_generators):
			target_table.append(row)

		#
		# empty ourselves to prevent this process from being repeated
		#

		del self.segment_lists[:]


#
# =============================================================================
#
#                                  Utilities
#
# =============================================================================
#


def insert_from_segwizard(ligolw_segments, fileobj, instruments, name, comment):
	"""
	Parse the contents of the file object fileobj as a segwizard-format
	segment list, and insert the result as a new list of "active"
	segments into the LigolwSegments object ligolw_segments.  A new
	entry will be created in the segment_definer table for the segment
	list, and instruments, name and comment are used to populate the
	entry's metadata.
	"""
	ligolw_segments.segment_lists.append(LigolwSegmentList(active = segmentsUtils.fromsegwizard(fileobj, coltype = LIGOTimeGPS), instruments = instruments, name = name, comment = comment))


def has_segment_tables(xmldoc, name = None):
	"""
	Return True if the document contains a complete set of segment
	tables.  Returns False otherwise.  If name is given and not None
	then the return value is True only if the document's segment
	tables, if present, contain a segment list by that name.
	"""
	try:
		def_table = table.get_table(xmldoc, lsctables.SegmentDefTable.tableName)
	except ValueError:
		return False
	try:
		table.get_table(xmldoc, lsctables.SegmentTable.tableName)
	except ValueError:
		return False
	if name is not None and name not in set(row.name for row in def_table):
		return False
	return True


def segmenttable_get_by_name(xmldoc, name):
	"""
	Retrieve the segments whose name matches those given.  The result
	is a segmentlistdict indexed by instrument.

	The output of this function is not coalesced, each segmentlist
	contains the segments as found in the segment table.
	"""
	#
	# find required tables
	#

	def_table = table.get_table(xmldoc, lsctables.SegmentDefTable.tableName)
	seg_table = table.get_table(xmldoc, lsctables.SegmentTable.tableName)

	#
	# segment_def_id --> instrument names mapping but only for
	# segment_definer entries bearing the requested name
	#

	instrument_index = dict((row.segment_def_id, row.get_ifos()) for row in def_table if row.name == name)

	#
	# populate result segmentlistdict object from segment_def_map table
	# and index
	#

	result = segments.segmentlistdict()

	for row in seg_table:
		try:
			instruments = instrument_index[row.segment_def_id]
		except KeyError:
			# not a segment we want
			continue
		seg = row.get()
		for instrument in instruments:
			try:
				result[instrument].append(seg)
			except KeyError:
				result[instrument] = segments.segmentlist([seg])

	#
	# done
	#

	return result
