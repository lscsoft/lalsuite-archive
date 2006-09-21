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
	return map(lambda l: CacheEntry(l, coltype = LIGOTimeGPS), fileobj)


def cache_to_seglistdict(cache):
	s = segments.segmentlistdict()
	for c in cache:
		s |= c.to_segmentlistdict()
	return s


def get_time_slides(filename, verbose = False, gz = False):
	"""
	Load the XML document contained in filename and convert the time
	slide table therein to a list of dictionaries of instrument/offset
	pairs.  The optional verbose and gz parameters are the same as for
	the glue.ligolw.utils.load_filename() function.  Raises ValueError
	if the document does not contain exactly 1 time slide table, or if
	one or more offsets in the table cannot be expressed as LIGOTimeGPS
	objects.
	"""
	tisitable = table.get_table(utils.load_filename(filename, verbose = verbose, gz = gz), lsctables.TimeSlideTable.tableName)
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

class LALCache(packing.Bin):
	"""
	Bin object representing a LAL file cache.  The objects attribute
	contains a list of glue.lal.CacheEntry objects, and the size
	attribute holds a glue.segments.segmentlistdict object summarizing
	the times spanned by the files in the cache.
	"""
	def __init__(self):
		packing.Bin.__init__(self)
		self.size = segments.segmentlistdict()

	def __str__(self):
		return "\n".join(map(str, self.objects))


class CafePacker(packing.Packer):
	def set_time_slides(self, offsetdictlist):
		self.timeslides = offsetdictlist

	def pack(self, cache_entry):
		# find all bins in which this cache_entry belongs.  the
		# test is, for each bin, to iterate through time slide
		# dictionaries applying the offsets to both the bin and the
		# cache_entry and checking if the cache_entry's segment
		# intersect the segments already in the bin, comparing only
		# those that are involved in the time slide.  FIXME: almost
		# correct:  still doesn't take care to compare only the
		# segments for the instruments identified in the offset
		# dictionary (so gives false positives but no false
		# negatives)
		size = cache_entry.to_segmentlistdict()
		matching_bins = []
		for n, bin in enumerate(self.bins):
			for offsetdict in self.timeslides:
				size.offsets.update(offsetdict)
				bin.size.offsets.update(offsetdict)
				if bin.size.is_coincident(size):
					matching_bins.append(n)
					break
			bin.size.offsets.clear()
		size.offsets.clear()

		# add cache_entry by either adding a new bin or putting it
		# into the first bin that was found
		if not matching_bins:
			self.bins.append(LALCache())
			self.bins[-1].add(cache_entry, size)
			return
		dest = self.bins[matching_bins.pop(0)]
		dest.add(cache_entry, size)

		# if cache_entry belongs in more than one bin, merge bins.
		# reverse the list of matching bins so that we delete the
		# highest-numbered bin first (otherwise the bins would
		# shift as we go and we would delete the wrong ones).
		matching_bins.reverse()
		for n in matching_bins:
			dest += self.bins.pop(n)


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

	# optimization: adding files to bins in time order keeps the number
	# of bins from growing larger than needed.
	if verbose:
		print >>sys.stderr, "sorting input cache ..."
	cache = list(cache)
	cache.sort(lambda a, b: cmp(a.segment, b.segment))

	# Pack cache entries into output caches.
	outputcaches = []
	packer = CafePacker(outputcaches)
	packer.set_time_slides(time_slides)
	if verbose:
		print >>sys.stderr, "packing files ..."
	for n, cacheentry in enumerate(cache):
		if verbose and not n % max(5, (len(cache)/1000)):
			print >>sys.stderr, "	%.1f%%	(%d files, %d caches)\r" % (100.0 * n / len(cache), n + 1, len(outputcaches)),
		packer.pack(cacheentry)
	if verbose:
		print >>sys.stderr, "	100.0%%	(%d files, %d caches)" % (n + 1, len(outputcaches))
		print >>sys.stderr, "sorting output caches ..."
	for cache in outputcaches:
		cache.objects.sort()

	return seglists.keys(), outputcaches
