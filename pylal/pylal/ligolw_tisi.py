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


from glue.ligolw import lsctables
from pylal import itertools
from pylal import llwapp


__author__ = "Kipp Cannon <kipp@gravity.phys.uwm.edu>"
__version__ = "$Revision$"[11:-2]
__date__ = "$Date$"[7:-2]


#
# =============================================================================
#
#                              Time Slide Parsing
#
# =============================================================================
#


def parse_slidespec(slidespec):
	"""
	Accepts a string in the format
	instrument=first:last:step[,first:last:step]...  and returns the
	tuple (instrument, [offset1, offset2, ....]) where the offsets are
	the sorted list of numbers described by the ranges.
	"""
	try:
		instrument, rangespec = slidespec.split("=")
	except ValueError:
		raise ValueError, "cannot parse time slide \"%s\"" % slidespec
	offsets = []
	for range in rangespec.strip().split(","):
		try:
			first, last, step = map(float, range.strip().split(":"))
		except ValueError:
			raise ValueError, "malformed range \"%s\" in \"%s\"" % (range, rangespec)
		if step == 0:
			if first != last:
				raise ValueError, "divide by zero in range \"%s\"" % range
			offsets.append(first)
			continue
		if (last - first) / step < 0.0:
			raise ValueError, "step has wrong sign in range \"%s\"" % range
		i = 0
		while True:
			offsets.append(first + i * step)
			if offsets[-1] >= last:
				if offsets[-1] > last:
					del offsets[-1]
				break
			i += 1
	offsets.sort()
	return instrument.strip(), offsets


def parse_slides(slides):
	"""
	Accepts a list of strings of the format understood by
	parse_slidespec() and returns a dictionary of instrument, offsets
	pairs
	"""
	d = {}
	for slidespec in slides:
		instrument, offsets = parse_slidespec(slidespec)
		if instrument in d:
			raise ValueError, "duplicate instrument \"%s\" in time slides" % instrument
		d[instrument] = offsets
	return d


#
# =============================================================================
#
#                           Add Process Information
#
# =============================================================================
#


def append_process(doc, **kwargs):
	process = llwapp.append_process(doc, program = "ligolw_tisi", version = __version__, cvs_repository = "lscsoft", cvs_entry_time = __date__, comment = kwargs["comment"])

	llwapp.append_process_params(doc, process, [("--instrument", "lstring", instrument) for instrument in kwargs["instrument"]])

	return process


#
# =============================================================================
#
#                              Build Time Slides
#
# =============================================================================
#


def SlidesIter(slides):
	"""
	Accepts a dictionary mapping instrument --> list-of-offsets, and
	iterates over the cartesian (outer) product of the offset lists,
	yielding all possible N-way instrument --> offset mappings.

	Example:

	>>> slides = {"H1": [-1, 0, 1], "H2": [-1, 0, 1], "L1": [0]}
	>>> list(SlidesIter(slides))
	[{'H2': -1, 'H1': -1, 'L1': 0}, {'H2': 0, 'H1': -1, 'L1': 0},
	{'H2': 1, 'H1': -1, 'L1': 0}, {'H2': -1, 'H1': 0, 'L1': 0}, {'H2':
	0, 'H1': 0, 'L1': 0}, {'H2': 1, 'H1': 0, 'L1': 0}, {'H2': -1, 'H1':
	1, 'L1': 0}, {'H2': 0, 'H1': 1, 'L1': 0}, {'H2': 1, 'H1': 1, 'L1':
	0}]
	"""
	instruments = slides.keys()
	for slide in itertools.MultiIter(slides.values()):
		offsetdict = {}
		for instrument, offset in zip(instruments, slide):
			offsetdict[instrument] = offset
		yield offsetdict


def RowsFromOffsetDict(offsetdict, time_slide_id, process):
	"""
	Accepts a dictionary mapping instrument --> offset, and a
	time_slide ID, and yields a sequence of rows to append to the
	time_slide table.  process must be the row in the process table on
	which the newly-constructed time_slide table rows are to be blamed.
	"""
	for instrument, offset in offsetdict.iteritems():
		row = lsctables.TimeSlide()
		row.process_id = process.process_id
		row.time_slide_id = time_slide_id
		row.instrument = instrument
		row.offset = offset
		yield row


#
# =============================================================================
#
#                            Time Slide Comparison
#
# =============================================================================
#


def time_slide_cmp(offsetdict1, offsetdict2):
	"""
	Compare two offset dictionaries mapping instrument --> offset.  The
	dictionaries compare as equal (return value is 0) if the relative
	offsets are all equal.

	Example:

	>>> d1 = {"H1": 0.0, "H2": 0.0, "L1": 0.0}
	>>> d2 = {"H1": 10.0, "H2": 10.0, "L1": 10.0}
	>>> time_slide_cmp(d1, d2)
	0

	because although the absolute offsets are not equal in the two
	dictionaries, the relative offsets are.
	"""
	if offsetdict2:
		# there is at least one entry in offsetdict2, pick an
		# instrument at random
		inst = offsetdict2.iterkeys().next()
		if inst in offsetdict1:
			# it is listed in offsetdict1, so make a working
			# copy of offsetdict2
			offsetdict2 = offsetdict2.copy()
			# compute the offset difference for the common
			# instrument
			delta = offsetdict1[inst] - offsetdict2[inst]
			# add it to the offsets in offsetdict2
			for inst in offsetdict2.keys():
				offsetdict2[inst] += delta
	# either the offsets have now been normalized to one another, or it
	# was discovered that the two offset dictionaries have different
	# instrument lists;  either way we can now use the built-in cmp
	# method
	return cmp(offsetdict1, offsetdict2)


def time_slides_vacuum(time_slides):
	"""
	Given a dictionary mapping time slide IDs to instrument-->offset
	mappings, for example as returned by the get_offsets() method of
	the TimeSlideTable class in glue.ligolw.lsctables, construct and
	return a mapping indicating time slide equivalences.  This can be
	used to delete redundant time slides from a time slide table, and
	then also used via the applyKeyMapping() method of
	glue.ligolw.table.Table instances to update cross references.

	Example:

	>>> slides = {"time_slide_id:0": {"H1": 0, "H2": 0},
	"time_slide_id:1": {"H1": 10, "H2": 10}, "time_slide_id:2": {"H1":
	0, "H2": 10}}
	>>> time_slides_vacuum(slides)
	{'time_slide_id:1': 'time_slide_id:0'}

	indicating that time_slide_id:1 describes a time slide that is
	equivalent to time_slide_id:0.  The calling code could use this
	information to delete time_slide_id:1 from the time_slide table,
	and replace occurances of that ID with time_slide_id:0 in other
	tables.
	"""
	# so we can modify it
	time_slides = time_slides.copy()
	# old --> new mapping
	mapping = {}
	# while there are time slide offset dictionaries remaining
	while time_slides:
		# pull out an ID/offset dictionary pair at random
		id1, offsetdict1 = time_slides.popitem()
		# for every other ID/offset dictionary pair in the time
		# slides
		for id2, offsetdict2 in time_slides.items():
			# if the offset dictionaries are equivalent
			if not time_slide_cmp(offsetdict1, offsetdict2):
				# remove it, and record in the old --> new
				# mapping
				time_slides.pop(id2)
				mapping[id2] = id1
	# done
	return mapping

