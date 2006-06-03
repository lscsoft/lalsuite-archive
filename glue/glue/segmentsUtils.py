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
This module provides additional utilities for use with segments.segmentlist
objects.
"""

__author__ = "Kipp Cannon <kipp@gravity.phys.uwm.edu>"
__date__ = "$Date$"[7:-2]
__version__ = "$Revision$"[11:-2]

import re
import segments
import lal


#
# =============================================================================
#
#                                     I/O
#
# =============================================================================
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
		[(s, d)] = pattern.findall(name.strip())
		s = coltype(s)
		d = coltype(d)
		l.append(segments.segment(s, s + d))
	return l


def fromlalcache(cachefile, coltype = int):
	"""
	Construct a segmentlist representing the times spanned by the files
	identified in the LAL cache contained in the file object file.  The
	segmentlist will be created with segments whose boundaries are of
	type coltype, which should raise ValueError if it cannot convert
	its string argument.
	"""
	return segments.segmentlist([c.segment for c in map(lambda l: lal.CacheEntry(l, coltype = coltype), cachefile)])


def fromsegwizard(file, coltype=int, strict=True):
	"""
	Read a segmentlist from the file object file containing a segwizard
	compatible segment list.  Parsing stops on the first line that
	cannot be parsed (which is consumed).  The segmentlist will be
	created with segments whose boundaries are of type coltype, which
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
	twocolsegpat = re.compile(r"\A\s*([\d.]+)\s+([\d.]+)\s*\Z")
	fourcolsegpat = re.compile(r"\A\s*([\d]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s*\Z")
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
				duration = seg.duration()
				this_line_format = 2
			except ValueError:
				break
		if strict:
			if seg.duration() != duration:
				raise ValueError, "segment \"" + line + "\" has incorrect duration"
			if not format:
				format = this_line_format
			elif format != this_line_format:
				raise ValueError, "segment \"" + line + "\" format mismatch"
		l.append(seg)
	return l


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
	segmentpat = re.compile(r"\A\s*\S+\s+\S+\s+\S+\s+([\d.]+)\s+([\d.]+)")
	l = segments.segmentlist()
	for line in file:
		try:
			[tokens] = segmentpat.findall(line)
			l.append(segments.segment(map(coltype, tokens[0:2])))
		except ValueError:
			break
	return l


def tosegwizard(file, seglist, header=True, coltype=int):
	"""
	Write the segmentlist seglist to the file object file in a
	segwizard compatible format.  If header is True, then the output
	will begin with a comment line containing column names.  The
	segment boundaries will be coerced to type coltype before output.
	"""
	if header:
		print >>file, "# segment\tstart\tstop\tduration"
	for n, seg in enumerate(seglist):
		print >>file, "%d\t%s\t%s\t%s" % (n, coltype(seg[0]), coltype(seg[1]), coltype(seg.duration()))


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
	return segments.segmentlist([segments.segment(t, t + 600) for t in range(extent[0] - ((extent[0] - 729273613) % 6370), extent[1], 6370)]) & segments.segmentlist([extent])


def segment_range(start, stop, period):
	"""
	Analogous to Python's range() builtin, returns a segmentlist of
	continuous adjacent segments each of length "period" with the first
	starting at "start" and the last ending not after "stop".  Note
	that the output of this function is not a coalesced list.  start,
	stop, and period can be any objects which support basic arithmetic
	operations.
	"""
	new = segments.segmentlist()
	for n in xrange((stop - start) / period):
		new.append(segments.segment(start + n * period, start + (n + 1) * period))
	return new


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
	of seglist1 with each segment in seglist2.  In each result, the
	segment start and stop values are adjusted to be with respect to
	the start of the corresponding segment in seglist2.  See also the
	range() function.
	"""
	for seg in seglist2:
		yield (seglist1 & segments.segmentlist([seg])).shift(-seg[0])
