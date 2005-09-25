"""
This module provides additional utilities for use with segments.segmentlist
objects.
"""

__author__ = "Kipp Cannon <kipp@gravity.phys.uwm.edu>"
__date__ = "$Date$"
__version__ = "$Revision$"

import re
import shlex
import segments


def fromfilenames(filenames, coltype=long):
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
		[(s, d)] = pattern.findall(name)
		s = coltype(s)
		d = coltype(d)
		l.append(segments.segment(s, s + d))
	return l


def fromsegwizard(file, coltype=long, strict=True):
	"""
	Read a segmentlist from the file object file containing a segwizard
	compatible segment list.  Parsing stops on the first line that
	cannot be parsed (which is consumed).  The segmentlist will be
	created with segments whose boundaries are of type coltype, which
	should raise ValueError if it cannot convert its string argument.
	If strict is True, then each segment's duration is checked against
	the input file.

	NOTE:  the output is a segmentlist as described by the file;  if
	the segments in the input file are not coalesced or out of order,
	then thusly shall be the output of this function.  It is
	recommended that this function's output be coalesced before use.
	"""
	l = segments.segmentlist()
	for line in file:
		tokens = shlex.split(line, comments=True)
		if len(tokens) == 0:
			continue
		if len(tokens) != 4:
			break
		try:
			seg = segments.segment(map(coltype, tokens[1:2]))
			duration = coltype(tokens[3])
		except ValueError:
			break
		if strict and seg.duration() != duration:
			raise ValueError, "segment \"" + line + "\" has incorrect duration"
		l.append(seg)
	return l


def tosegwizard(file, seglist, header=True, coltype=long):
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
