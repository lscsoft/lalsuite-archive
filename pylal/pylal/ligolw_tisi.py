# Copyright (C) 2006--2013  Kipp Cannon
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


import itertools
import sys


from glue import iterutils
from glue import offsetvector
from glue.ligolw import lsctables
from glue.ligolw import utils as ligolw_utils
from pylal import git_version


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date


#
# =============================================================================
#
#                             Command Line Parsing
#
# =============================================================================
#


def parse_slidespec(slidespec):
	"""
	Accepts a string in the format
	instrument=first:last:step[,first:last:step]...  and returns the
	tuple (instrument, [offset1, offset2, ....]) where the offsets are
	the sorted list of unique numbers described by the ranges.  A range
	with a positive step describes the offsets (first + n * step) where
	n is an integer such that first <= offset <= last.  A range with a
	negative step describes the offsets (first + n * step) where n is
	an integer such that last <= offset <= first.

	Example:

	>>> parse_slidespec("H1=-5:+5:0.5")
	('H1', [-5.0, -4.5, -4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5,
	0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0])
	"""
	try:
		instrument, rangespec = slidespec.split("=")
	except ValueError:
		raise ValueError("cannot parse time slide '%s'" % slidespec)
	offsets = set()
	for range in [s.strip() for s in rangespec.split(",")]:
		try:
			first, last, step = map(float, range.split(":"))
		except ValueError:
			raise ValueError("malformed range '%s' in '%s'" % (range, rangespec))
		if step == 0:
			if first != last:
				raise ValueError("divide by zero in range '%s'" % range)
			offsets.add(first)
			continue
		if (last - first) / step < 0.0:
			raise ValueError("step has wrong sign in range '%s'" % range)

		for i in itertools.count():
			x = first + i * step
			if step > 0:
				if not (first <= x <= last):
					break
			elif not (last <= x <= first):
				break
			offsets.add(x)
	return instrument.strip(), sorted(offsets)


def parse_slides(slides):
	"""
	Accepts a list of strings of the format understood by
	parse_slidespec(), and returns a dictionary mapping instrument
	names to sorted lists of unique offsets.

	Example:

	>>> parse_slides(["H1=-1:+1:+1", "H2=-1:+1:+1", "L1=0:0:0"])
	{'H2': [-1.0, 0.0, 1.0], 'H1': [-1.0, 0.0, 1.0], 'L1': [0.0]}
	"""
	d = {}
	# store the offsets for each instrument as sets to uniquify the
	# numbers
	for slidespec in slides:
		instrument, offsets = parse_slidespec(slidespec)
		try:
			d[instrument] |= set(offsets)
		except KeyError:
			d[instrument] = set(offsets)
	# convert offsets back to sorted lists
	return dict((instrument, sorted(offsets)) for instrument, offsets in d.items())


def parse_inspiral_num_slides_slidespec(slidespec):
	"""
	Accepts a string in the format
	count:instrument=offset[,instrument=offset...] and returns the
	tuple (count, {instrument: offset, ...})

	Example:

	>>> parse_inspiral_num_slides_slidespec("3:H1=0,H2=5,L1=10")
	(3, {'H2': 5.0, 'H1': 0.0, 'L1': 10.0})
	"""
	count, offsets = slidespec.strip().split(":")
	offsetvect = offsetvector.offsetvector((instrument.strip(), float(offset)) for instrument, offset in (token.strip().split("=") for token in offsets.strip().split(",")))
	return int(count), offsetvect


#
# =============================================================================
#
#                                     I/O
#
# =============================================================================
#


class DefaultContentHandler(lsctables.ligolw.LIGOLWContentHandler):
	pass
lsctables.use_in(DefaultContentHandler)


def load_time_slides(filename, verbose = False, gz = None, contenthandler = DefaultContentHandler):
	"""
	Load a time_slide table from the LIGO Light Weight XML file named
	filename, or stdin if filename is None.  See
	glue.ligolw.utils.load_filename() for a description of the verbose
	and gz parameters.  The return value is a dictionary mapping each
	time slide ID to a dictionary of instrument/offset pairs for that
	time slide.

	Note that a side effect of this function is that the ID generator
	associated with the TimeSlideTable class in glue.ligolw.lsctables
	is synchronized with the result, so that the next ID it generates
	will be immediately following the IDs listed in the dictionary
	returned by this function.

	Note also that this utility function should not be how applications
	that perform multiple manipulations with an XML file obtain the
	time slide table's contents since this function re-parses the file
	from scratch.  Instead, from the glue.ligolw package use
	table.get_table(...).as_dict().
	"""
	time_slide_table = lsctables.TimeSlideTable.get_table(ligolw_utils.load_filename(filename, verbose = verbose, gz = gz, contenthandler = contenthandler))
	time_slide_table.sync_next_id()
	return time_slide_table.as_dict()


def get_time_slide_id(xmldoc, time_slide, create_new = None, superset_ok = False, nonunique_ok = False):
	"""
	Return the time_slide_id corresponding to the offset vector
	described by time_slide, a dictionary of instrument/offset pairs.

	Example:

	>>> get_time_slide_id(xmldoc, {"H1": 0, "L1": 0})
	'time_slide:time_slide_id:10'

	This function is a wrapper around the .get_time_slide_id() method
	of the glue.ligolw.lsctables.TimeSlideTable class.  See the
	documentation for that class for the meaning of the create_new,
	superset_ok and nonunique_ok keyword arguments.

	This function requires the document to contain exactly one
	time_slide table.  If the document does not contain exactly one
	time_slide table then ValueError is raised, unless the optional
	create_new argument is not None.  In that case a new table is
	created.  This effect of the create_new argument is in addition to
	the affects described by the TimeSlideTable class.
	"""
	try:
		tisitable = lsctables.TimeSlideTable.get_table(xmldoc)
	except ValueError:
		# table not found
		if create_new is None:
			raise
		tisitable = lsctables.New(lsctables.TimeSlideTable)
		xmldoc.childNodes[0].appendChild(tisitable)
	# make sure the next_id attribute is correct
	tisitable.sync_next_id()
	# get the id
	return tisitable.get_time_slide_id(time_slide, create_new = create_new, superset_ok = superset_ok, nonunique_ok = nonunique_ok)


#
# =============================================================================
#
#                              Build Time Slides
#
# =============================================================================
#


def SlidesIter(slides):
	"""
	Accepts a dictionary mapping instrument --> list-of-offsets (for
	example, as returned by parse_slides()), and iterates over the
	cartesian (outer) product of the offset lists, yielding all
	possible N-way instrument --> offset mappings.

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
	for slide in iterutils.MultiIter(*slides.values()):
		yield offsetvector.offsetvector(zip(instruments, slide))


def Inspiral_Num_Slides_Iter(count, offsets):
	"""
	This generator yields a sequence of time slide dictionaries in the
	style of lalapps_thinca's time slides.  Each resulting dictionary
	maps instrument to offset.  The input is a count of time slides (an
	integer), and a dictionary mapping instrument to offset.  The
	output dictionaries describe time slides that are integer multiples
	of the input time shifts.

	Example:

	>>> list(Inspiral_Num_Slides_Iter(3, {"H1": 0.0, "H2": 5.0,"L1": 10.0}))
	[{'H2': -15.0, 'H1': -0.0, 'L1': -30.0}, {'H2': -10.0, 'H1': -0.0,
	'L1': -20.0}, {'H2': -5.0, 'H1': -0.0, 'L1': -10.0}, {'H2': 0.0,
	'H1': 0.0, 'L1': 0.0}, {'H2': 5.0, 'H1': 0.0, 'L1': 10.0}, {'H2':
	10.0, 'H1': 0.0, 'L1': 20.0}, {'H2': 15.0, 'H1': 0.0, 'L1': 30.0}]

	The output time slides are all integer multiples of the input time
	shift vector in the range [-count, +count], and are returned in
	increasing order of mupltiplier.
	"""
	offsets = offsets.items()
	for n in range(-count, +count + 1):
		yield offsetvector.offsetvector((instrument, offset * n) for instrument, offset in offsets)


def RowsFromOffsetDict(offsetvect, time_slide_id, process):
	"""
	Accepts a dictionary mapping instrument --> offset, and a
	time_slide ID, and yields a sequence of rows to append to the
	time_slide table.  process must be the row in the process table on
	which the newly-constructed time_slide table rows are to be blamed.
	"""
	# FIXME:  remove when all occurances are replaced with
	# glue.ligolw.lsctables.TimeSlideTable.append_offsetvector()
	for instrument, offset in offsetvect.items():
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


def time_slides_vacuum(time_slides, verbose = False):
	"""
	Given a dictionary mapping time slide IDs to instrument-->offset
	mappings, for example as returned by the as_dict() method of the
	TimeSlideTable class in glue.ligolw.lsctables or by the
	load_time_slides() function in this module, construct and return a
	mapping indicating time slide equivalences.  This can be used to
	delete redundant time slides from a time slide table, and then also
	used via the applyKeyMapping() method of glue.ligolw.table.Table
	instances to update cross references (for example in the
	coinc_event table).

	Example:

	>>> slides = {"time_slide_id:0": {"H1": 0, "H2": 0},
	"time_slide_id:1": {"H1": 10, "H2": 10}, "time_slide_id:2": {"H1":
	0, "H2": 10}}
	>>> time_slides_vacuum(slides)
	{'time_slide_id:1': 'time_slide_id:0'}

	indicating that time_slide_id:1 describes a time slide that is
	equivalent to time_slide_id:0.  The calling code could use this
	information to delete time_slide_id:1 from the time_slide table,
	and replace references to that ID in other tables with references
	to time_slide_id:0.
	"""
	# convert offsets to deltas
	time_slides = dict((time_slide_id, offsetvect.deltas) for time_slide_id, offsetvect in time_slides.items())
	N = len(time_slides)
	# old --> new mapping
	mapping = {}
	# while there are time slide offset dictionaries remaining
	while time_slides:
		n = N - len(time_slides)
		if verbose and not (n % 10):
			print >>sys.stderr, "\t%.1f%%\r" % (100.0 * n / N),
		# pick an ID/offset dictionary pair at random
		id1, deltas1 = time_slides.popitem()
		# for every other ID/offset dictionary pair in the time
		# slides
		ids_to_delete = []
		for id2, deltas2 in time_slides.items():
			# if the relative offset dictionaries are
			# equivalent record in the old --> new mapping
			if deltas2 == deltas1:
				mapping[id2] = id1
				ids_to_delete.append(id2)
		for id2 in ids_to_delete:
			time_slides.pop(id2)
	# done
	if verbose:
		print >>sys.stderr, "\t100.0%"
	return mapping


def time_slide_list_merge(slides1, slides2):
	"""
	Merges two lists of offset dictionaries into a single list with
	no duplicate (equivalent) time slides.
	"""
	deltas1 = set(frozenset(offsetvect1.deltas.items()) for offsetvect1 in slides1)
	return slides1 + [offsetvect2 for offsetvect2 in slides2 if frozenset(offsetvect2.deltas.items()) not in deltas1]
