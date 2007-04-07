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
LIGO Light-Weight XML coincidence analysis front end.
"""


import bisect
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
	"""
	Parses a LAL cache file named filenamed into a list of
	glue.lal.CacheEntry objects.
	"""
	if verbose:
		print >>sys.stderr, "reading %s ..." % (filename or "stdin")
	if filename:
		fileobj = file(filename)
	else:
		fileobj = sys.stdin
	return map(lambda l: CacheEntry(l, coltype = LIGOTimeGPS), fileobj)


def cache_to_seglistdict(cache):
	"""
	Constructs a coalesced segmentlistdict object from a list of
	glue.lal.CacheEntry objects.
	"""
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
	objects.  Raises KeyError if a time slide lists the same instrument
	more than once.
	"""
	tisitable = table.get_table(utils.load_filename(filename, verbose = verbose, gz = gz), lsctables.TimeSlideTable.tableName)
	slides = {}
	for row in tisitable:
		if row.time_slide_id not in slides:
			slides[row.time_slide_id] = {}
		if row.instrument in slides[row.time_slide_id]:
			raise KeyError, "repeated instruments in %s" % row.time_slide_id
		slides[row.time_slide_id][row.instrument] = LIGOTimeGPS(row.offset)
	return slides.values()


#
# =============================================================================
#
#                             Output Cache Packing
#
# =============================================================================
#


class LALCacheBin(packing.Bin):
	"""
	Bin object representing a LAL file cache.  The objects attribute
	contains a list of glue.lal.CacheEntry objects, and the size
	attribute holds a glue.segments.segmentlistdict object summarizing
	the times spanned by the files in the cache.
	"""
	def __init__(self):
		packing.Bin.__init__(self)
		self.size = segments.segmentlistdict()
		self.extent = None

	def add(self, *args):
		packing.Bin.add(self, *args)
		self.extent = self.size.extent_all()
		return self

	def __iadd__(self, *args):
		packing.Bin.__iadd__(self, *args)
		self.extent = self.size.extent_all()
		return self

	def __cmp__(self, other):
		return cmp(self.extent, other.extent)

	def __str__(self):
		return "\n".join(map(str, self.objects))


class CafePacker(packing.Packer):
	def set_time_slides(self, offsetdictlist):
		self.timeslides = offsetdictlist
		min_offset = min([min(timeslide.itervalues()) for timeslide in offsetdictlist])
		max_offset = max([max(timeslide.itervalues()) for timeslide in offsetdictlist])
		# largest gap that can conceivably be closed by the time
		# slides
		self.max_gap = abs(max_offset - min_offset)

	def pack(self, cache_entry):
		# find all bins in which this cache_entry belongs.  the
		# test is, for each bin, to iterate through time slide
		# dictionaries applying the offsets to both the bin and the
		# cache_entry and checking if the cache_entry's segment
		# intersects the segments already in the bin, comparing
		# only those that are involved in the time slide.  FIXME:
		# almost correct:  still doesn't take care to compare only
		# the segments for the instruments identified in the offset
		# dictionary (so gives false positives but no false
		# negatives)
		new = LALCacheBin()
		new.add(cache_entry, cache_entry.to_segmentlistdict())
		new.extent.protract(self.max_gap)
		matching_bins = []
		for n in xrange((bisect.bisect_left(self.bins, new) or 1) - 1, len(self.bins)):
			bin = self.bins[n]
			for offsetdict in self.timeslides:
				new.size.offsets.update(offsetdict)
				bin.size.offsets.update(offsetdict)
				if bin.size.is_coincident(new.size):
					matching_bins.append(n)
					break
			bin.size.offsets.clear()
		new.size.offsets.clear()
		new.extent = new.size.extent_all()

		if not matching_bins:
			# no matching bins, add a new one
			self.bins.append(new)
		else:
			# put file into first bin that was found to match
			dest = self.bins[matching_bins.pop(0)]
			dest += new

			# if cache_entry belongs in more than one bin,
			# merge bins.  reverse the list of matching bins so
			# that we delete the highest-numbered bin first
			# (otherwise the bins would shift as we go and we
			# would delete the wrong ones).
			matching_bins.reverse()
			for n in matching_bins:
				dest += self.bins.pop(n)
		self.bins.sort()


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
			print >>sys.stderr, "writing %s ..." % filename
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
	# Construct a segment list dictionary from the cache
	if verbose:
		print >>sys.stderr, "computing segment list ..."
	seglists = cache_to_seglistdict(cache)

	# For each instrument compute the times for which it will
	# contribute triggers to a coincidence analysis.  Times not spanned
	# by these lists are those times for which no time slide can
	# possibly produce coincident triggers.
	seglists = llwapp.get_coincident_segmentlistdict(seglists, time_slides)

	# Remove files that will not participate in a coincidence.  Rather
	# than modifying the list in place, we make a copy so as to not
	# modify the calling code's data.
	if verbose:
		print >>sys.stderr, "filtering input cache ..."
	cache = [c for c in cache if seglists.intersects(c.to_segmentlistdict())]

	# optimization: adding files to bins in time order keeps the number
	# of bins from growing larger than needed.
	if verbose:
		print >>sys.stderr, "sorting input cache ..."
	cache.sort(lambda a, b: cmp(a.segment, b.segment))

	# Pack cache entries into output caches.
	outputcaches = []
	packer = CafePacker(outputcaches)
	packer.set_time_slides(time_slides)
	if verbose:
		print >>sys.stderr, "packing files ..."
	for n, cacheentry in enumerate(cache):
		if verbose and not n % 7:
			print >>sys.stderr, "	%.1f%%	(%d files, %d caches)\r" % (100.0 * n / len(cache), n + 1, len(outputcaches)),
		packer.pack(cacheentry)
	if verbose:
		print >>sys.stderr, "	100.0%%	(%d files, %d caches)" % (n + 1, len(outputcaches))
		print >>sys.stderr, "sorting output caches ..."
	for cache in outputcaches:
		cache.objects.sort()

	return seglists, outputcaches
