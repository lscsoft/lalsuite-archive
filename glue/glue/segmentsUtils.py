# Copyright (C) 2006  Kipp Cannon
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
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
This module provides additional utilities for use with segments.segmentlist
objects.
"""


import re


from glue import git_version
from glue import iterutils
from glue import lal
from glue import segments


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date


#
# =============================================================================
#
#                                     I/O
#
# =============================================================================
#


#
# A list of file names
#


def fromfilenames(filenames, coltype=int):
	"""
	Return a segmentlist describing the intervals spanned by the files
	whose names are given in the list filenames.  The segmentlist is
	constructed by parsing the file names, and the boundaries of each
	segment are coerced to type coltype.

	The file names are parsed using a generalization of the format
	described in Technical Note LIGO-T010150-00-E, which allows the
	start time and duration appearing in the file name to be
	non-integers.

	NOTE:  the output is a segmentlist as described by the file names;
	if the file names are not in time order, or describe overlaping
	segments, then thusly shall be the output of this function.  It is
	recommended that this function's output be coalesced before use.
	"""
	pattern = re.compile(r"-([\d.]+)-([\d.]+)\.[\w_+#]+\Z")
	l = segments.segmentlist()
	for name in filenames:
                if name.endswith('.gz'):
                  name = name.rstrip('.gz')
		[(s, d)] = pattern.findall(name.strip())

		s = coltype(s)
		d = coltype(d)
		l.append(segments.segment(s, s + d))
	return l


#
# LAL cache files
#


def fromlalcache(cachefile, coltype = int):
	"""
	Construct a segmentlist representing the times spanned by the files
	identified in the LAL cache contained in the file object file.  The
	segmentlist will be created with segments whose boundaries are of
	type coltype, which should raise ValueError if it cannot convert
	its string argument.
	"""
	return segments.segmentlist([lal.CacheEntry(l, coltype = coltype).segment for l in cachefile])


#
# Segwizard-formated segment list text files
#


def fromsegwizard(file, coltype=int, strict=True):
	"""
	Read a segmentlist from the file object file containing a segwizard
	compatible segment list.  Parsing stops on the first line that
	cannot be parsed (which is consumed).  The segmentlist will be
	created with segment whose boundaries are of type coltype, which
	should raise ValueError if it cannot convert its string argument.
	Both two-column and four-column segwizard files are recognized, but
	the entire file must be in the same format, which is decided by the
	first parsed line.  If strict is True and the file is in
	four-column format, then each segment's duration is checked against
	that column in the input file.

	NOTE:  the output is a segmentlist as described by the file;  if
	the segments in the input file are not coalesced or out of order,
	then thusly shall be the output of this function.  It is
	recommended that this function's output be coalesced before use.
	"""
	commentpat = re.compile(r"\s*([#;].*)?\Z", re.DOTALL)
	twocolsegpat = re.compile(r"\A\s*([\d.+-eE]+)\s+([\d.+-eE]+)\s*\Z")
	fourcolsegpat = re.compile(r"\A\s*([\d]+)\s+([\d.+-eE]+)\s+([\d.+-eE]+)\s+([\d.+-eE]+)\s*\Z")
	format = None
	l = segments.segmentlist()
	for line in file:
		line = commentpat.split(line)[0]
		if not len(line):
			continue
		try:
			[tokens] = fourcolsegpat.findall(line)
			num = int(tokens[0])
			seg = segments.segment(map(coltype, tokens[1:3]))
			duration = coltype(tokens[3])
			this_line_format = 4
		except ValueError:
			try:
				[tokens] = twocolsegpat.findall(line)
				seg = segments.segment(map(coltype, tokens[0:2]))
				duration = abs(seg)
				this_line_format = 2
			except ValueError:
				break
		if strict:
			if abs(seg) != duration:
				raise ValueError, "segment '%s' has incorrect duration" % line
			if format is None:
				format = this_line_format
			elif format != this_line_format:
				raise ValueError, "segment '%s' format mismatch" % line
		l.append(seg)
	return l


def tosegwizard(file, seglist, header=True, coltype=int):
	"""
	Write the segmentlist seglist to the file object file in a
	segwizard compatible format.  If header is True, then the output
	will begin with a comment line containing column names.  The
	segment boundaries will be coerced to type coltype before output.
	"""
	if header:
		print >>file, "# seg\tstart    \tstop     \tduration"
	for n, seg in enumerate(seglist):
		print >>file, "%d\t%s\t%s\t%s" % (n, str(coltype(seg[0])), str(coltype(seg[1])), str(coltype(abs(seg))))


#
# TAMA-formated segment list text files
#


def fromtama(file, coltype=lal.LIGOTimeGPS):
	"""
	Read a segmentlist from the file object file containing TAMA
	locked-segments data.  Parsing stops on the first line that cannot
	be parsed (which is consumed).  The segmentlist will be created
	with segments whose boundaries are of type coltype, which should
	raise ValueError if it cannot convert its string argument.

	NOTE:  TAMA locked-segments files contain non-integer start and end
	times, so the default column type is set to LIGOTimeGPS.

	NOTE:  the output is a segmentlist as described by the file;  if
	the segments in the input file are not coalesced or out of order,
	then thusly shall be the output of this function.  It is
	recommended that this function's output be coalesced before use.
	"""
	segmentpat = re.compile(r"\A\s*\S+\s+\S+\s+\S+\s+([\d.+-eE]+)\s+([\d.+-eE]+)")
	l = segments.segmentlist()
	for line in file:
		try:
			[tokens] = segmentpat.findall(line)
			l.append(segments.segment(map(coltype, tokens[0:2])))
		except ValueError:
			break
	return l


#
# Command line or config file strings
#


def from_range_strings(ranges, boundtype = int):
	"""
	Parse a list of ranges expressed as strings in the form "value" or
	"first:last" into an equivalent glue.segments.segmentlist.  In the
	latter case, an empty string for "first" and(or) "last" indicates a
	(semi)infinite range.  A typical use for this function is in
	parsing command line options or entries in configuration files.

	NOTE:  the output is a segmentlist as described by the strings;  if
	the segments in the input file are not coalesced or out of order,
	then thusly shall be the output of this function.  It is
	recommended that this function's output be coalesced before use.

	Example:

	>>> text = "0:10,35,100:"
	>>> from_range_strings(text.split(","))
	[segment(0, 10), segment(35, 35), segment(100, infinity)]
	"""
	# preallocate segmentlist
	segs = segments.segmentlist([None] * len(ranges))

	# iterate over strings
	for i, range in enumerate(ranges):
		parts = range.split(":")
		if len(parts) == 1:
			parts = boundtype(parts[0])
			segs[i] = segments.segment(parts, parts)
			continue
		if len(parts) != 2:
			raise ValueError, range
		if parts[0] == "":
			parts[0] = segments.NegInfinity
		else:
			parts[0] = boundtype(parts[0])
		if parts[1] == "":
			parts[1] = segments.PosInfinity
		else:
			parts[1] = boundtype(parts[1])
		segs[i] = segments.segment(parts[0], parts[1])

	# success
	return segs


def to_range_strings(seglist):
	"""
	Turn a segment list into a list of range strings as could be parsed
	by from_range_strings().  A typical use for this function is in
	machine-generating configuration files or command lines for other
	programs.

	Example:

	>>> from glue.segments import *
	>>> segs = segmentlist([segment(0, 10), segment(35, 35), segment(100, infinity())])
	>>> ",".join(to_range_strings(segs))
	'0:10,35,100:'
	"""
	# preallocate the string list
	ranges = [None] * len(seglist)

	# iterate over segments
	for i, seg in enumerate(seglist):
		if not seg:
			ranges[i] = str(seg[0])
		elif (seg[0] is segments.NegInfinity) and (seg[1] is segments.PosInfinity):
			ranges[i] = ":"
		elif (seg[0] is segments.NegInfinity) and (seg[1] is not segments.PosInfinity):
			ranges[i] = ":%s" % str(seg[1])
		elif (seg[0] is not segments.NegInfinity) and (seg[1] is segments.PosInfinity):
			ranges[i] = "%s:" % str(seg[0])
		elif (seg[0] is not segments.NegInfinity) and (seg[1] is not segments.PosInfinity):
			ranges[i] = "%s:%s" % (str(seg[0]), str(seg[1]))
		else:
			raise ValueError, seg

	# success
	return ranges


def segmentlistdict_to_short_string(seglists):
	return "/".join(["%s=%s" % (str(key), ",".join(to_range_strings(value))) for key, value in seglists.items()])


def segmentlistdict_from_short_string(s):
	d = segments.segmentlistdict()
	for token in s.strip().split("/"):
		key, ranges = token.strip().split("=")
		d[key.strip()] = from_range_strings(ranges.strip().split(","))
	return d


#
# =============================================================================
#
#                    Pre-defined Segments and Segment Lists
#
# =============================================================================
#


def S2playground(extent):
	"""
	Return a segmentlist identifying the S2 playground times within the
	interval defined by the segment extent.
	"""
	lo = int(extent[0])
	lo -= (lo - 729273613) % 6370
	hi = int(extent[1]) + 1
	return segments.segmentlist(segments.segment(t, t + 600) for t in range(lo, hi, 6370)) & segments.segmentlist([extent])


def segmentlist_range(start, stop, period):
	"""
	Analogous to Python's range() builtin, this generator yields a
	sequence of continuous adjacent segments each of length "period"
	with the first starting at "start" and the last ending not after
	"stop".  Note that the segments generated do not form a coalesced
	list (they are not disjoint).  start, stop, and period can be any
	objects which support basic arithmetic operations.

	Example:

	>>> from glue.segments import *
	>>> segmentlist(segmentlist_range(0, 15, 5))
	[segment(0, 5), segment(5, 10), segment(10, 15)]
	>>> segmentlist(segmentlist_range('', 'xxx', 'x'))
	[segment('', 'x'), segment('x', 'xx'), segment('xx', 'xxx')]
	"""
	n = 1
	b = start
	while True:
		a, b = b, start + n * period
		if b > stop:
			break
		yield segments.segment(a, b)
		n += 1


#
# =============================================================================
#
#                         Extra Manipulation Routines
#
# =============================================================================
#


def Fold(seglist1, seglist2):
	"""
	An iterator that generates the results of taking the intersection
	of seglist1 with each segment in seglist2 in turn.  In each result,
	the segment start and stop values are adjusted to be with respect
	to the start of the corresponding segment in seglist2.  See also
	the segmentlist_range() function.

	This has use in applications that wish to convert ranges of values
	to ranges relative to epoch boundaries.  Below, a list of time
	intervals in hours is converted to a sequence of daily interval
	lists with times relative to midnight.

	Example:

	>>> from glue.segments import *
	>>> x = segmentlist([segment(0, 13), segment(14, 20), segment(22, 36)])
	>>> for y in Fold(x, segmentlist_range(0, 48, 24)): print y
	...
	[segment(0, 13), segment(14, 20), segment(22, 24)]
	[segment(0, 12)]
	"""
	for seg in seglist2:
		yield (seglist1 & segments.segmentlist([seg])).shift(-seg[0])


def vote(seglists, n):
	"""
	Given a sequence of segmentlists, returns the intervals during
	which at least n of them intersect.

	Example:

	>>> from glue.segments import *
	>>> w = segmentlist([segment(0, 15)])
	>>> x = segmentlist([segment(5, 20)])
	>>> y = segmentlist([segment(10, 25)])
	>>> z = segmentlist([segment(15, 30)])
	>>> vote((w, x, y, z), 3)
	[segment(10, 20)]
	"""
	return reduce(lambda x, y: x | y, (reduce(lambda a, b: a & b, votes) for votes in iterutils.choices(seglists, n)))
