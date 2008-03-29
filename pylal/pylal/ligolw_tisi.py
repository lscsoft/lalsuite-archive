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
	offsets = list(offsets)
	offsets.sort()
	return instrument.strip(), offsets


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
		if instrument not in d:
			d[instrument] = set()
		d[instrument] |= set(offsets)
	# convert offsets back to lists, and sort
	for instrument, offsets in d.items():
		d[instrument] = list(offsets)
		d[instrument].sort()
	return d


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


def Inspiral_Num_Slides_Iter(slidespec):
	"""
	Accepts a string in the format understood by
	parse_inspiral_slidespec(), and generates a sequence of time slide
	dictionaries mapping instrument to offset.

	Example:

	>>> list(Inspiral_Num_Slides_Iter("3:H1=0,H2=5,L1=10"))
	[{'H2': -15.0, 'H1': -0.0, 'L1': -30.0}, {'H2': -10.0, 'H1': -0.0, 'L1': -20.0}, {'H2': -5.0, 'H1': -0.0, 'L1': -10.0}, {'H2': 5.0, 'H1': 0.0, 'L1': 10.0}, {'H2': 10.0, 'H1': 0.0, 'L1': 20.0}, {'H2': 15.0, 'H1': 0.0, 'L1': 30.0}]

	The instrument/offset pairs in the input string form a vector of
	offsets, and the output time slides are all integer multiples of
	that vector with multiples in the range [-count, +count] excluding
	0.
	"""
	count, offsets = parse_inspiral_num_slides_slidespec(slidespec)
	offsets = offsets.items()
	for n in range(-count, +count + 1):
		if n == 0:
			continue
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
	if offsetdict2:
		# offsetdict2 is not empty, pick an instrument at random
		instrument, offset = offsetdict2.iteritems().next()
		if instrument in offsetdict1:
			# the instrument is listed in offsetdict1, so make
			# a working copy of offsetdict2
			offsetdict2 = offsetdict2.copy()
			# compute the offset difference for the common
			# instrument
			delta = offsetdict1[instrument] - offset
			# add it to the offsets in the working copy of
			# offsetdict2
			for instrument in offsetdict2.keys():
				offsetdict2[instrument] += delta
	# either the offsets have now been normalized to one another, or it
	# was discovered that the two offset dictionaries have different
	# instrument lists;  either way we can now use the built-in cmp
	# method
	return cmp(offsetdict1, offsetdict2)


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
	# so we can modify it
	time_slides = time_slides.copy()
	N = len(time_slides)
	# old --> new mapping
	mapping = {}
	# while there are time slide offset dictionaries remaining
	while time_slides:
		n = N - len(time_slides)
		if verbose and not (n % 10):
			print >>sys.stderr, "\t%.1f%%\r" % (100.0 * n / N),
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
