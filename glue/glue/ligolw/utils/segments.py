# Copyright (C) 2008  Kipp Cannon
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


from glue import git_version
from glue import iterutils
from glue import segments
from glue import segmentsUtils
from glue.lal import LIGOTimeGPS
from .. import table
from .. import lsctables


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date


#
# =============================================================================
#
#                                 Segment List
#
# =============================================================================
#


class LigolwSegmentList(object):
	"""
	A description of a LIGO Light-Weight XML segment list.  Instances
	of this class carry all the metadata associated with a LIGO Light-
	Weight XML segment list including its name, version number, a
	comment, and so on.

	LIGO Light-Weight XML segment lists are three-state objects.  A
	segment list can be on, off, or undefined.  Two separate sequences
	of segments are used for this:  the "valid" list defines the
	intervals when the state of the segment list is known, and the
	"active" list defines the intervals when the segment list is on.
	It is not an error for the active list to be on during times when
	the segment lists state is unknown, this code does not impose any
	policy in that regard, but it should be expected that consumers of
	the segment list will treat all times when the segment list's state
	is unknown the same way.
	"""
	#
	# the columns the segment_definer, segment_summary and segment
	# tables need to have
	#

	segment_def_columns = (u"process_id", u"segment_def_id", u"ifos", u"name", u"version", u"insertion_time", u"comment")
	segment_sum_columns = (u"process_id", u"segment_sum_id", u"start_time", u"start_time_ns", u"end_time", u"end_time_ns", u"segment_def_id", u"comment")
	segment_columns = (u"process_id", u"segment_id", u"start_time", u"start_time_ns", u"end_time", u"end_time_ns", u"segment_def_id")

	def __init__(self, active = (), valid = (), instruments = set(), name = None, version = None, insertion_time = None, comment = None):
		"""
		Initialize a new LigolwSegmentList instance.  active and
		valid are sequences that will be cast to
		segments.segmentlist objects.  They can be generator
		expressions.  The "active" sequence is what is usually
		thought of as the segment list, the "valid" sequence
		identifies the intervals of time for which the segment
		list's state is defined.
		"""
		self.valid = segments.segmentlist(valid)
		self.active = segments.segmentlist(active)
		self.instruments = instruments
		self.name = name
		self.version = version
		self.insertion_time = insertion_time
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
		Coalesce the internal segment lists.  Returns self.
		"""
		self.valid.coalesce()
		self.active.coalesce()
		return self


#
# =============================================================================
#
#                                 Library API
#
# =============================================================================
#


class LigolwSegments(set):
	"""
	An interface shim between code that makes use of segments in
	glue.segments form, and LIGO Light-Weight XML I/O code.

	This class is "attached" to an XML document object, at which time
	it parses and extracts the segment lists from the document, and
	clears the document's segment tables (preventing a second
	LigolwSegments object from being meaningfully attached to the same
	document).  When the application is finished manipulating the
	segment lists, they can be inserted into the XML document at which
	time the contents of the LigolwSegments object are cleared
	(preventing any further manipulations).

	This class is a subclass of the Python set builtin.  Each element
	of the set is a LigolwSegmentList instance describing one of the
	segment lists in the original XML document.
	"""
	def __init__(self, xmldoc):
		#
		# Find tables
		#

		try:
			self.segment_def_table = table.get_table(xmldoc, lsctables.SegmentDefTable.tableName)
		except ValueError:
			self.segment_def_table = xmldoc.childNodes[0].appendChild(lsctables.New(lsctables.SegmentDefTable, LigolwSegmentList.segment_def_columns))

		try:
			self.segment_sum_table = table.get_table(xmldoc, lsctables.SegmentSumTable.tableName)
		except ValueError:
			self.segment_sum_table = xmldoc.childNodes[0].appendChild(lsctables.New(lsctables.SegmentSumTable, LigolwSegmentList.segment_sum_columns))

		try:
			self.segment_table = table.get_table(xmldoc, lsctables.SegmentTable.tableName)
		except ValueError:
			self.segment_table = xmldoc.childNodes[0].appendChild(lsctables.New(lsctables.SegmentTable, LigolwSegmentList.segment_columns))

		#
		# Transform segment tables into a collection of
		# LigolwSegmentList objects for more convenient
		# manipulation
		#

		# construct empty LigolwSegmentList objects, one for each
		# entry in the segment_definer table, indexed by
		# segment_definer id
		segment_lists = dict((row.segment_def_id, LigolwSegmentList(instruments = row.get_ifos(), name = row.name, version = row.version, insertion_time = row.insertion_time, comment = row.comment)) for row in self.segment_def_table)
		if len(segment_lists) != len(self.segment_def_table):
			raise ValueError("duplicate segment_definer IDs detected in segment_definer table")
		del self.segment_def_table[:]

		# populate LigolwSegmentList objects from segment table and
		# segment_summary table
		for row in self.segment_sum_table:
			segment_lists[row.segment_def_id].valid.append(row.get())
		del self.segment_sum_table[:]
		for row in self.segment_table:
			segment_lists[row.segment_def_id].active.append(row.get())
		del self.segment_table[:]

		#
		# transcribe LigolwSegmentList objects into self.  segment
		# definer IDs are now meaningless
		#

		self.update(segment_lists.values())

		#
		# Synchronize ID generators
		#

		# FIXME:  why am I doing this!?  I've just deleted all the
		# rows.  the id generators should probably be 0'ed here,
		# and these sync's moved to the .finalize() method to
		# prevent collisions
		self.segment_def_table.sync_next_id()
		self.segment_table.sync_next_id()
		self.segment_sum_table.sync_next_id()

		#
		# Done
		#


	def insert_from_segwizard(self, fileobj, instruments, name, version = None, insertion_time = None, comment = None):
		"""
		Parse the contents of the file object fileobj as a
		segwizard-format segment list, and insert the result as a
		new list of "active" segments into this LigolwSegments
		object.  A new entry will be created in the segment_definer
		table for the segment list, and instruments, name and
		comment are used to populate the entry's metadata.  Note
		that the "valid" segments are left empty, nominally
		indicating that there are no periods of validity.
		"""
		self.add(LigolwSegmentList(active = segmentsUtils.fromsegwizard(fileobj, coltype = LIGOTimeGPS), instruments = instruments, name = name, version = version, insertion_time = insertion_time, comment = comment))


	def insert_from_segmentlistdict(self, seglists, name, version = None, insertion_time = None, comment = None):
		"""
		Insert the segments from the segmentlistdict object
		seglists as a new list of "active" segments into this
		LigolwSegments object.  The dictionary's keys are assumed
		to provide the instrument name for each segment list.  A
		new entry will be created in the segment_definer table for
		the segment lists, and the dictionary's keys, the name, and
		comment will be used to populate the entry's metadata.
		"""
		for instrument, segments in seglists.items():
			self.add(LigolwSegmentList(active = segments, instruments = set([instrument]), name = name, version = version, insertion_time = insertion_time, comment = comment))


	def coalesce(self):
		"""
		Coalesce the segment lists.  Returns self.
		"""
		for ligolw_segment_list in self:
			ligolw_segment_list.coalesce()
		return self


	def sort(self, *args):
		"""
		Sort the segment lists.  The optional args are passed to
		the .sort() methods of the segment lists.  This can be used
		to control the sort order by providing an alternate
		comparison function (the default is to sort all lists by
		segment start time with ties broken by end time).
		"""
		for ligolw_segment_list in self:
			ligolw_segment_list.sort(*args)


	def optimize(self):
		"""
		Identifies segment lists that differ only in their
		instruments --- they have the same valid and active
		segments, the same name, version and the same comment ---
		and then deletes all but one of them, leaving just a single
		list having the union of the instruments.
		"""
		self.sort()
		segment_lists = dict(enumerate(self))
		for target, source in [(idx_a, idx_b) for (idx_a, seglist_a), (idx_b, seglist_b) in iterutils.choices(segment_lists.items(), 2) if seglist_a.valid == seglist_b.valid and seglist_a.active == seglist_b.active and seglist_a.name == seglist_b.name and seglist_a.version == seglist_b.version and seglist_a.comment == seglist_b.comment]:
			try:
				source = segment_lists.pop(source)
			except KeyError:
				continue
			segment_lists[target].insertion_time = max(segment_lists[target].insertion_time, source.insertion_time)
			segment_lists[target].instruments |= source.instruments
		self.clear()
		self.update(segment_lists.values())


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
		# generator function to convert segments into row objects,
		# each paired with the table to which the row is to be
		# appended
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
		# of table row generators.  empty ourselves to prevent this
		# process from being repeated
		#

		row_generators = []
		while self:
			ligolw_segment_list = self.pop()
			segment_def_row = self.segment_def_table.RowType()
			segment_def_row.process_id = process_row.process_id
			segment_def_row.segment_def_id = self.segment_def_table.get_next_id()
			segment_def_row.set_ifos(ligolw_segment_list.instruments)
			segment_def_row.name = ligolw_segment_list.name
			segment_def_row.version = ligolw_segment_list.version
			segment_def_row.insertion_time = ligolw_segment_list.insertion_time
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
# =============================================================================
#
#                                  Utilities
#
# =============================================================================
#


def has_segment_tables(xmldoc, name = None):
	"""
	Return True if the document contains a complete set of segment
	tables.  Returns False otherwise.  If name is given and not None
	then the return value is True only if the document's segment
	tables, if present, contain a segment list by that name.
	"""
	try:
		def_table = table.get_table(xmldoc, lsctables.SegmentDefTable.tableName)
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

	instruments = set(instrument for instruments in instrument_index.values() for instrument in instruments)
	result = segments.segmentlistdict((instrument, segments.segmentlist()) for instrument in instruments)

	for row in seg_table:
		if row.segment_def_id in instrument_index:
			seg = row.get()
			for instrument in instrument_index[row.segment_def_id]:
				result[instrument].append(seg)

	#
	# done
	#

	return result
