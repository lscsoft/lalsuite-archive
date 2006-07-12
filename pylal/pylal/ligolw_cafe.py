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

"""
LIGO Light-Weight XML Coincidence Analysis Front End.
"""

from math import log10
import sys

from glue import segments
from glue.lal import CacheEntry
from glue.ligolw import lsctables
from pylal import llwapp
from pylal import packing
from pylal.date import LIGOTimeGPS

__author__ = "Kipp Cannon <kipp@gravity.phys.uwm.edu>"
__date__ = "$Date$"[7:-2]
__version__ = "$Revision$"[11:-2]


#
# =============================================================================
#
#                                    Input
#
# =============================================================================
#

def load_cache(filename, verbose = False):
	if verbose:
		print >>sys.stderr, "reading %s ..." % (filename or "stdin")
	if filename:
		fileobj = file(filename)
	else:
		fileobj = sys.stdin
	return map(lambda l: CacheEntry(l, coltype = LIGOTimeGPS), fileobj)


def cache_to_seglistdict(cache):
	s = segments.segmentlistdict()
	for c in cache:
		try:
			s[c.observatory].append(c.segment)
		except KeyError:
			s[c.observatory] = segments.segmentlist([c.segment])
	s.coalesce()
	return s


def get_time_slides(filename, verbose = False):
	doc = llwapp.load_filename(filename, verbose)
	tisitable = llwapp.get_table(doc, lsctables.TimeSlideTable.tableName)
	doc.unlink()
	for row in tisitable:
		row.offset = LIGOTimeGPS(row.offset)
	return map(tisitable.get_offset_dict, tisitable.dict.keys())


def null_time_slides(cache):
	timeslide = {}
	for entry in cache:
		timeslide[entry.observatory] = LIGOTimeGPS(0)
	return [timeslide]


#
# =============================================================================
#
#                             Output Cache Packing
#
# =============================================================================
#

def cacheentry_to_seglistdict(cacheentry):
	return segments.segmentlistdict({cacheentry.observatory: segments.segmentlist([cacheentry.segment])})


class CafePacker(packing.Packer):
	def set_time_slides(self, offsetdictlist):
		self.timeslides = offsetdictlist

	def pack(self, size, object):
		# find all bins in which this object belongs.  the test is,
		# for each bin, to iterate through time slide dictionaries
		# applying the offsets to both the bin and the object and
		# checking if the object's segments intersect the segments
		# already in the bin, comparing only those that are
		# involved in the time slide.  FIXME: what's done here is
		# only an approximation of the correct test, that finds
		# false positives but no false negatives.
		matching_bins = []
		for n, bin in enumerate(self.bins):
			for offsetdict in self.timeslides:
				size.offsets.update(offsetdict)
				bin.size.offsets.update(offsetdict)
				a = segments.segmentlist()
				b = segments.segmentlist()
				for key in offsetdict.iterkeys():
					if key in size:
						a |= size[key]
					if key in bin.size:
						b |= bin.size[key]
				if a.intersects(b):
					matching_bins.append((n, bin))
					break

		# reset all offsets
		size.offsets.clear()
		for bin in self.bins:
			bin.size.offsets.clear()

		# add object by either adding a new bin or putting it into
		# the first bin that was found
		if not matching_bins:
			self.bins.append(packing.LALCache())
			self.bins[-1].add(object, size)
			return
		matching_bins[0][1].add(object, size)

		# if object belongs in more than one bin, merge bins.
		# reverse the list of matching bins so that we delete the
		# highest-numbered bin first (otherwise the bins would
		# shift as we go and we would delete the wrong ones).
		matching_bins.reverse()
		for n, bin in matching_bins[:-1]:
			matching_bins[-1][1].objects.extend(bin.objects)
			matching_bins[-1][1].size += bin.size
			del self.bins[n]


#
# =============================================================================
#
#                               Pack Input Files
#
# =============================================================================
#

#
# Find cache entries that intersect the surviving segments and pack into
# output caches.
#

def build_output_caches(cache, seglists, time_slides, verbose):
	outputcaches = []
	packer = CafePacker(outputcaches)
	packer.set_time_slides(time_slides)

	if verbose:
		print >>sys.stderr, "packing files..."
	for n, cacheentry in enumerate(cache):
		if verbose and not n % max(5, (len(cache)/1000)):
			print >>sys.stderr, "	%.1f%%	(%d files, %d caches)\r" % (100.0 * n / len(cache), n, len(outputcaches)),
		cache_seglistdict = cacheentry_to_seglistdict(cacheentry)
		if seglists.intersects(cache_seglistdict):
			packer.pack(cache_seglistdict, cacheentry)
	if verbose:
		print >>sys.stderr, "	100.0%%	(%d files, %d caches)" % (n, len(outputcaches))

	return outputcaches


#
# =============================================================================
#
#                                    Output
#
# =============================================================================
#

def write_caches(base, bins, instruments, verbose = False):
	filenames = []
	if len(bins):
		pattern = "%%s%%0%dd.cache" % int(log10(len(bins)) + 1)
	for n, bin in enumerate(bins):
		filename = pattern % (base, n)
		filenames.append(filename)
		if verbose:
			print >>sys.stderr, "writing %s" % filename
		f = file(filename, "w")
		for cacheentry in bin.objects:
			if cacheentry.observatory in instruments:
				print >>f, str(cacheentry)
	return filenames


def write_single_instrument_caches(base, bins, instruments, verbose = False):
	for instrument in instruments:
		write_caches("%s%s_" % (base, instrument), bins, [instrument], verbose)


#
# =============================================================================
#
#                                 Library API
#
# =============================================================================
#

def ligolw_cafe(cache, time_slides, verbose = False):
	if verbose:
		print >>sys.stderr, "computing segment list..."
	seglists = cache_to_seglistdict(cache)
	# For each instrument compute the times for which it will
	# contribute triggers to a coincidence analysis.  Times not spanned
	# by these lists are those times for which no time slide can
	# possibly produce coincident triggers.
	seglists = llwapp.get_coincident_segmentlistdict(seglists, time_slides)

	return seglists.keys(), build_output_caches(cache, seglists, time_slides, verbose)
