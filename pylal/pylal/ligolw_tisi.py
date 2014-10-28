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
import warnings


from glue import iterutils
from glue import offsetvector
from glue.ligolw import lsctables
from glue.ligolw import utils as ligolw_utils
from glue.ligolw.utils import time_slide as ligolw_time_slide
from pylal import git_version


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date


def get_time_slide_id(*args, **kwargs):
	warnings.warn("pylal.ligolw_tisi.get_time_slide_id() is deprecated.  use glue.ligolw.utils.time_slide.get_time_slide_id() instead.", DeprecationWarning)
	return ligolw_time_slide.get_time_slide_id(*args, **kwargs)


def time_slides_vacuum(*args, **kwargs):
	warnings.warn("pylal.ligolw_tisi.time_slides_vacuum() is deprecated.  use glue.ligolw.utils.time_slide.time_slides_vacuum() instead.", DeprecationWarning)
	return ligolw_time_slide.time_slides_vacuum(*args, **kwargs)


def time_slide_list_merge(*args, **kwargs):
	warnings.warn("pylal.ligolw_tisi.time_slide_list_merge() is deprecated.  use glue.ligolw.utils.time_slide.time_slide_list_merge() instead.", DeprecationWarning)
	return ligolw_time_slide.time_slide_list_merge(*args, **kwargs)


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


@lsctables.use_in
class DefaultContentHandler(lsctables.ligolw.LIGOLWContentHandler):
	pass


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
	will not conflict with the IDs listed in the dictionary returned by
	this function.

	Note also that this utility function should not be how applications
	that perform multiple manipulations with an XML file obtain the
	time slide table's contents since this function re-parses the file
	from scratch.  Instead, from the glue.ligolw package use
	table.get_table(...).as_dict().
	"""
	warnings.warn("pylal.ligolw_tisi.load_time_slides() is deprecated.  use glue.ligolw library directly to load document instead.", DeprecationWarning)
	time_slide_table = lsctables.TimeSlideTable.get_table(ligolw_utils.load_filename(filename, verbose = verbose, gz = gz, contenthandler = contenthandler))
	time_slide_table.sync_next_id()
	return time_slide_table.as_dict()


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
