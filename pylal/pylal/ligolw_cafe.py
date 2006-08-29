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
from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw import utils
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
	cache = map(lambda l: CacheEntry(l, coltype = LIGOTimeGPS), fileobj)
	if verbose:
		print >>sys.stderr, "sorting by segment ..."
	cache.sort(lambda a, b: cmp(a.segment, b.segment))
	return cache


def cache_to_seglistdict(cache):
	s = segments.segmentlistdict()
	for c in cache:
		for instrument in c.observatory.split(","):
			if instrument in s:
				s[instrument].append(c.segment)
			else:
				s[instrument] = segments.segmentlist([c.segment])
	s.coalesce()
	return s


def get_time_slides(filename, verbose = False):
	doc = utils.load_filename(filename, verbose)
	tisitable = table.get_table(doc, lsctables.TimeSlideTable.tableName)
	for row in tisitable:
		row.offset = LIGOTimeGPS(row.offset)
	return map(tisitable.get_offset_dict, tisitable.dict.keys())


#
# =============================================================================
#
#                             Output Cache Packing
#
# =============================================================================
#

def cacheentry_to_seglistdict(cacheentry):
	s = segments.segmentlistdict()
	for instrument in cacheentry.observatory.split(","):
		s[instrument] = segments.segmentlist([cacheentry.segment])
	return s


class CafePacker(packing.Packer):
	def set_time_slides(self, offsetdictlist):
		self.timeslides = offsetdictlist

	def pack(self, cache_entry):
		# find all bins in which this cache_entry belongs.  the
		# test is, for each bin, to iterate through time slide
		# dictionaries applying the offsets to both the bin and the
		# cache_entry and checking if the cache_entry's segment
		# intersect the segments already in the bin, comparing only
		# those that are involved in the time slide.  FIXME: what's
		# done here is only an approximation of the correct test;
		# an approximation that finds false positives but no false
		# negatives.
		matching_bins = []
		for n, bin in enumerate(self.bins):
			for offsetdict in self.timeslides:
				if True in map(offsetdict.__contains__, cache_entry.observatory.split(",")):
					bin.size.offsets.update(offsetdict)
					b = segments.segmentlist()
					for key in offsetdict.iterkeys():
						if key in bin.size:
							b.extend(bin.size[key])
					# FIXME: does not yet handle the case
					# of a cache entry with multiple
					# instruments
					if b.coalesce().intersects_segment(cache_entry.segment.shift(offsetdict[cache_entry.observatory])):
						matching_bins.append((n, bin))
						break
			bin.size.offsets.clear()

		# add cache_entry by either adding a new bin or putting it
		# into the first bin that was found
		size = cacheentry_to_seglistdict(cache_entry)
		if not matching_bins:
			self.bins.append(packing.LALCache())
			self.bins[-1].add(cache_entry, size)
			return
		dest = matching_bins.pop(0)[1]
		dest.add(cache_entry, size)

		# if cache_entry belongs in more than one bin, merge bins.
		# reverse the list of matching bins so that we delete the
		# highest-numbered bin first (otherwise the bins would
		# shift as we go and we would delete the wrong ones).
		matching_bins.reverse()
		for n, bin in matching_bins:
			dest.objects.extend(bin.objects)
			dest.size += bin.size
			del self.bins[n]


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
			if True in map(instruments.__contains__, cacheentry.observatory.split(",")):
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
		print >>sys.stderr, "computing segment list ..."
	seglists = cache_to_seglistdict(cache)
	# For each instrument compute the times for which it will
	# contribute triggers to a coincidence analysis.  Times not spanned
	# by these lists are those times for which no time slide can
	# possibly produce coincident triggers.
	seglists = llwapp.get_coincident_segmentlistdict(seglists, time_slides)

	# Pack cache entries into output caches.
	outputcaches = []
	packer = CafePacker(outputcaches)
	packer.set_time_slides(time_slides)
	if verbose:
		print >>sys.stderr, "packing files..."
	for n, cacheentry in enumerate(cache):
		if verbose and not n % max(5, (len(cache)/1000)):
			print >>sys.stderr, "	%.1f%%	(%d files, %d caches)\r" % (100.0 * n / len(cache), n, len(outputcaches)),
		packer.pack(cacheentry)
	if verbose:
		print >>sys.stderr, "	100.0%%	(%d files, %d caches)" % (n, len(outputcaches))

	return seglists.keys(), outputcaches
