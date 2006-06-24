#!/usr/bin/python
#
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
	tuple (instruct, [offset1, offset2, ....]) where the offsets are
	the sorted list of numbers described by the ranges.
	"""
	try:
		[instrument, rangespec] = slidespec.split("=")
	except ValueError:
		raise ValueError, "cannot parse time slide \"%s\"" % slidespec
	offsets = []
	for range in rangespec.split(","):
		try:
			[first, last, step] = range.split(":")
			first, last, step = float(first), float(last), float(step)
		except ValueError:
			raise ValueError, "malformed range \"%s\" in \"%s\"" % (range, rangespec)
		if step == 0.0:
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
	return instrument, offsets


def parse_slides(slides):
	"""
	Accepts a list of strings of the format understood by
	parse_slidespec() and returns a dictionary of instrument, offsets
	pairs
	"""
	d = {}
	for slidespec in slides:
		instrument, offsets = parse_slidespec(slidespec)
		if d.has_key(instrument):
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

class SlidesIter(object):
	def __init__(self, slides):
		self.instruments = slides.keys()
		self.slides = itertools.MultiIter(slides.values())

	def __iter__(self):
		return self

	def next(self):
		offsetdict = {}
		for n, offset in enumerate(self.slides.next()):
			offsetdict[self.instruments[n]] = offset
		return offsetdict


#
# =============================================================================
#
#                                 Library API
#
# =============================================================================
#

def ligolw_tisi(doc, **kwargs):
	timeslidetable = llwapp.get_table(doc, lsctables.TimeSlideTable.tableName)
	process = append_process(doc, **kwargs)
	ids = lsctables.TimeSlideIDs()
	for offsetdict in SlidesIter(parse_slides(kwargs["instrument"])):
		id = ids.next()
		for instrument, offset in offsetdict.iteritems():
			row = lsctables.TimeSlide()
			row.process_id = process.process_id
			row.time_slide_id = id
			row.instrument = instrument
			row.offset = offset
			timeslidetable.append(row)
	llwapp.set_process_end_time(process)
	return doc
