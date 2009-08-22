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


import itertools
import sys
# Python 2.3 compatibility
try:
	set
except NameError:
	from sets import Set as set


from glue import iterutils
from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw import utils
from glue.ligolw.utils import process as ligolw_process
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
		raise ValueError, "cannot parse time slide '%s'" % slidespec
	offsets = set()
	for range in [s.strip() for s in rangespec.split(",")]:
		try:
			first, last, step = map(float, range.split(":"))
		except ValueError:
			raise ValueError, "malformed range '%s' in '%s'" % (range, rangespec)
		if step == 0:
			if first != last:
				raise ValueError, "divide by zero in range '%s'" % range
			offsets.add(first)
			continue
		if (last - first) / step < 0.0:
			raise ValueError, "step has wrong sign in range '%s'" % range
		i = 0
		while True:
			x = first + i * step
			if step > 0:
				if not (first <= x <= last):
					break
			elif not (last <= x <= first):
				break
			offsets.add(x)
			i += 1
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
	offsets = dict([(instrument.strip(), float(offset)) for instrument, offset in [token.strip().split("=") for token in offsets.strip().split(",")]])
	return int(count), offsets


def load_time_slides(filename, verbose = False, gz = False):
	"""
	Load a time_slide table from the LIGO Light Weight XML file named
	filename, or stdin if filename is None.  Extra verbosity is printed
	if verbose is True, and the file is gzip decompressed while reading
	if gz is Tue.  The output is returned as a dictionary, mapping each
	time slide ID to a dictionary providing a mapping of instrument to
	offset for that time slide.

	Note that a side effect of this function is that the ID generator
	associated with the TimeSlideTable class in glue.ligolw.lsctables
	is synchronized with the result, so that the next ID it generates
	will be immediately following the IDs listed in the dictionary
	returned by this function.
	"""
	time_slide_table = table.get_table(utils.load_filename(filename, verbose = verbose, gz = (filename or "stdin")[-3:] == ".gz"), lsctables.TimeSlideTable.tableName)
	time_slides = time_slide_table.as_dict()
	time_slide_table.sync_next_id()
	return time_slides


#
# =============================================================================
#
#                           Add Process Information
#
# =============================================================================
#


def append_process(doc, **kwargs):
	process = llwapp.append_process(doc, program = u"ligolw_tisi", version = __version__, cvs_repository = u"lscsoft", cvs_entry_time = __date__, comment = kwargs["comment"])

	ligolw_process.append_process_params(doc, process, [(u"--instrument", u"lstring", instrument) for instrument in kwargs["instrument"]])

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
		yield dict(zip(instruments, slide))


def Inspiral_Num_Slides_Iter(count, offsets):
	"""
	This generator yields a sequence of time slide dictionaries in the
	style of the inspiral pipeline's time slides.  Each resulting
	dictionary maps instrument to offset.  The input is a count of time
	slides (an integer), and a dictionary mapping instrument to offset.
	The output dictionaries describe time slides that are integer
	multiples of the input time shifts.

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
		yield dict([(instrument, offset * n) for instrument, offset in offsets])


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


def offset_vector_to_deltas(offset_vector):
	"""
	Construct a dictionary of relative offsets from a dictionary of
	absolute offsets.

	Example:

	>>> offset_vector_to_deltas({"H1": 0, "L1": 10, "V1": 20})
	{('H1', 'L1'): 10, ('L1', 'V1'): 10}

	The keys in the result are instrument pairs, (a, b), and the values
	are the relative time offsets, (offset[b] - offset[a]).
	"""
	instruments = sorted(offset_vector)
	return dict(((a, b), offset_vector[b] - offset_vector[a]) for a, b in zip(instruments[:-1], instruments[1:]))


def time_slide_cmp(offsetdict1, offsetdict2):
	"""
	Compare two offset dictionaries mapping instrument --> offset.  The
	dictionaries compare as equal (return value is 0) if the relative
	offsets (not absolute offsets) are all equal.

	Example:

	>>> offsets1 = {"H1": 0.0, "H2": 0.0, "L1": 0.0}
	>>> offsets2 = {"H1": 10.0, "H2": 10.0, "L1": 10.0}
	>>> time_slide_cmp(offsets1, offsets2)
	0

	because although the absolute offsets are not equal in the two
	dictionaries, all relative offsets are.
	"""
	return cmp(offset_vector_to_deltas(offsetdict1), offset_vector_to_deltas(offsetdict2))


def time_slide_contains(offset_vector1, offset_vector2):
	"""
	Returns True if offset vector 2 can be found in offset vector 1,
	False otherwise.  An offset vector is "found in" another offset
	vector if the latter contains all of the former's instruments and
	the relative offsets among those instruments agree (the absolute
	offsets need not).
	"""
	return offset_vector_to_deltas(dict((instrument, offset) for instrument, offset in offset_vector1.items() if instrument in offset_vector2)) == offset_vector_to_deltas(offset_vector2)


def time_slides_vacuum(time_slides, verbose = False):
	"""
	Given a dictionary mapping time slide IDs to instrument-->offset
	mappings, for example as returned by the as_dict() method of the
	TimeSlideTable class in glue.ligolw.lsctables, construct and return
	a mapping indicating time slide equivalences.  This can be used to
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
	time_slides = dict((id, offset_vector_to_deltas(offset_vector)) for id, offset_vector in time_slides.items())
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
	new = []
	for offsetdict2 in slides2:
		for offsetdict1 in slides1:
			if not time_slide_cmp(offsetdict1, offsetdict2):
				# these are the same slides, discard
				break
		else:
			# loop completed without finding a match
			new.append(offsetdict2)
	return slides1 + new


#
# =============================================================================
#
#                           Time Slide Manipulation
#
# =============================================================================
#


def time_slide_component_vectors(offset_vectors, n):
	"""
	Given an iterable of time slide vectors, return the shortest list
	of the unique n-instrument time slide vectors from which all the
	vectors in the input list can be consructed.  This can be used to
	determine the minimal set of n-instrument coincs required to
	construct all of the coincs for all of the requested instrument and
	offset combinations in the time slide list.

	It is assumed that the coincs for the vector {"H1": 0, "H2": 10,
	"L1": 20} can be constructed from the coincs for the vectors {"H1":
	0, "H2": 10} and {"H2": 0, "L1": 10}, that is only the relative
	offsets are significant in determining if two events are
	coincident, not the absolute offsets.  This assumption is not true
	for the standard inspiral pipeline, where the absolute offsets are
	significant.
	"""
	#
	# collect unique instrument set / deltas combinations
	#

	delta_sets = {}
	for offset_vector in offset_vectors:
		for instruments in iterutils.choices(sorted(offset_vector), n):
			delta_sets.setdefault(instruments, set()).add(tuple(offset_vector[instrument] - offset_vector[instruments[0]] for instrument in instruments[1:]))

	#
	# translate into a list of n-instrument offset vectors
	#

	return [dict(zip(instruments, (0.0,) + deltas)) for instruments, delta_set in delta_sets.items() for deltas in delta_set]


def time_slide_normalize(time_slide, **kwargs):
	"""
	The time slide, a mapping of instrument --> offset, is adjusted so
	that a particular instrument in the slide has the desired offset.
	All other instruments have their offsets adjusted so that the
	relative offsets are preserved.  The instrument to noramlize, and
	the offset one wishes it to have, are provided as a key-word
	argument.  More than one key-word argument can be given, in which
	case they are considered in order until one is found that is
	applicable, that is the instrument is in the time slide.  This
	function is a no-op if no key-word argument is found that applies.
	The return value is the time slide dictionary, which is modified in
	place.

	Example:

	>>> time_slide_normalize({"H1": -10, "H2": -10, "L1": -10}, L1 = 0)
	{'H2': 0, 'H1': 0, 'L1': 0}
	>>> time_slide_normalize({"H1": -10, "H2": -10}, L1 = 0, H2 = 5)
	{'H2': 5, 'H1': 5}
	"""
	for instrument, offset in kwargs.iteritems():
		if instrument in time_slide:
			delta = offset - time_slide[instrument]
			for instrument in time_slide.keys():
				time_slide[instrument] += delta
			break
	return time_slide


#
# =============================================================================
#
#                                    Other
#
# =============================================================================
#


def display_component_offsets(component_offset_vectors, fileobj = sys.stderr):
	"""
	Print a summary of the output of time_slide_component_offsets().
	"""
	#
	# organize the information
	#
	# groupby requires its input to be grouped (= sorted) by the
	# grouping key (the instruments), so we have to do this first.
	# after constructing the strings, we make sure the lists of offset
	# strings are all the same length by appending empty strings as
	# needed.  finally we transpose the whole mess so that it's stored
	# as rows instead of columns.
	#

	l = sorted(component_offset_vectors, lambda a, b: cmp(sorted(a), sorted(b)))
	l = [[", ".join("%s-%s" % (b, a) for a, b in zip(instruments[:-1], instruments[1:]))] + [", ".join("%.17g s" % (offset_vector[b] - offset_vector[a]) for a, b in zip(instruments[:-1], instruments[1:])) for offset_vector in offset_vectors] for instruments, offset_vectors in itertools.groupby(l, lambda v: sorted(v))]
	n = max(len(offsets) for offsets in l)
	for offsets in l:
		offsets += [""] * (n - len(offsets))
	l = zip(*l)

	#
	# find the width of the columns
	#

	width = max(max(len(s) for s in line) for line in l)
	format = "%%%ds" % width

	#
	# print the offsets
	#

	lines = iter(l)
	print >>fileobj, " | ".join(format % s for s in lines.next())
	print >>fileobj, "-+-".join(["-" * width] * len(l[0]))
	for line in lines:
		print >>fileobj, " | ".join(format % s for s in line)
