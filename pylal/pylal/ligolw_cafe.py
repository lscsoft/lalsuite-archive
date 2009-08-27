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
# Python 2.3 compatibility
try:
	set
except NameError:
	from sets import Set as set


from glue import segments
from glue.lal import CacheEntry
from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw import utils
from pylal import llwapp
from pylal import packing


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
	Parses a LAL cache file named filename into a list of
	glue.lal.CacheEntry objects.
	"""
	if verbose:
		print >>sys.stderr, "reading %s ..." % (filename or "stdin")
	if filename is not None:
		f = file(filename)
	else:
		f = sys.stdin
	return [CacheEntry(line, coltype = lsctables.LIGOTimeGPS) for line in f]


def cache_to_seglistdict(cache):
	"""
	Constructs a coalesced segmentlistdict object from a list of
	glue.lal.CacheEntry objects.
	"""
	s = segments.segmentlistdict()
	for c in cache:
		s |= c.to_segmentlistdict()
	return s


#
# =============================================================================
#
#                             Performance Helpers
#
# =============================================================================
#


def segmentlistdict_normalize(seglistdict, origin):
	"""
	Convert the times in a segmentlist dictionary to floats relative to
	an origin.  The purpose is to allow segment lists stored as
	LIGOTimeGPS times to be manipulated more quickly without loss of
	precision.  The modification is done in place.
	"""
	for seglist in seglistdict.itervalues():
		for i, seg in enumerate(seglist):
			seglist[i] = segments.segment(float(seg[0] - origin), float(seg[1] - origin))


def segmentlistdict_unnormalize(seglistdict, origin):
	"""
	The opposite of segmentlistdict_normalize(), restores the times in
	a segmentlist dictionary to absolute times.  The modification is
	done in place.
	"""
	for seglist in seglistdict.itervalues():
		for i, seg in enumerate(seglist):
			seglist[i] = segments.segment(origin + seg[0], origin + seg[1])


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
	"""
	Packing algorithm implementing the ligolw_cafe event list file
	packing algorithm.
	"""
	def set_time_slides(self, offsetdictlist):
		"""
		Set the list of time slides to be considered when deciding
		the bins in which each file belongs.  Must be called before
		packing any files.  The input is a list of dictionaries,
		each mapping instruments to offsets.
		"""
		self.timeslides = offsetdictlist
		min_offset = min(min(timeslide.values()) for timeslide in offsetdictlist)
		max_offset = max(max(timeslide.values()) for timeslide in offsetdictlist)
		# largest gap that can conceivably be closed by the time
		# slides
		self.max_gap = max_offset - min_offset
		if self.max_gap < 0:
			raise Exception, "crash!"

	def pack(self, cache_entry):
		"""
		Find all bins in which this glue.lal.CacheEntry instance
		belongs.  The test is, for each bin, to iterate through the
		time slide dictionaries applying the offsets to both the
		bin and the CacheEntry, and checking if the (time shifted)
		CacheEntry's segment(s) intersects the (time shifted)
		segments already in the bin, comparing only those that are
		involved in the time slide.
		"""
		new = LALCacheBin()
		new.add(cache_entry, cache_entry.to_segmentlistdict())
		matching_bins = []
		for n in xrange(bisect.bisect_left([bin.extent[1] for bin in self.bins], new.extent[0] - self.max_gap), len(self.bins)):
			bin = self.bins[n]
			for offsetdict in self.timeslides:
				new.size.offsets.update(offsetdict)
				bin.size.offsets.update(offsetdict)
				if bin.size.is_coincident(new.size, keys = offsetdict.keys()):
					matching_bins.append(n)
					break
			bin.size.offsets.clear()
		new.size.offsets.clear()

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
			if instruments & set(cacheentry.to_segmentlistdict().keys()):
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
	"""
	cache is a sequence of glue.lal.CacheEntry objects.  time_slides is
	a sequence of dictionaries describing the time slides to consider.
	Set verbose to True for verbosity.  The output is a tuple, the
	first element of which is a glue.segments.segmentlistdict object
	describing the times for which coincident data is available.  The
	second element is a list of LALCacheBin objects, providing the file
	groups identified by the coincidence analysis.
	"""
	#
	# Construct a segment list dictionary from the cache
	#

	if verbose:
		print >>sys.stderr, "computing segment list ..."
	seglists = cache_to_seglistdict(cache)

	#
	# For each instrument compute the times for which it will (could)
	# contribute to a coincidence analysis.
	#

	segmentlistdict_normalize(seglists, lsctables.LIGOTimeGPS(800000000))
	seglists = llwapp.get_coincident_segmentlistdict(seglists, time_slides)
	segmentlistdict_unnormalize(seglists, lsctables.LIGOTimeGPS(800000000))

	#
	# Remove files that will not participate in a coincidence.  Take
	# care not to modify the calling code's data.  Note that because we
	# have established that this segment list describes exactly the
	# times spanned by the input files that are coincident under at
	# least one time slide, a file participates in a multi-instrument
	# coincidence if and only if it intersects these times.
	#

	if verbose:
		print >>sys.stderr, "filtering input cache ..."
	cache = [c for c in cache if seglists.intersects_all(c.to_segmentlistdict())]

	#
	# Optimization: adding files to bins in time order keeps the number
	# of bins from growing larger than needed.
	#

	if verbose:
		print >>sys.stderr, "sorting input cache ..."
	cache.sort(lambda a, b: cmp(a.segment, b.segment))

	#
	# Pack cache entries into output caches.  Having reduced the file
	# list to just those that participate in coincidences, it only
	# remains to determine which other files each must be grouped with.
	#

	outputcaches = []
	packer = CafePacker(outputcaches)
	packer.set_time_slides(time_slides)
	if verbose:
		print >>sys.stderr, "packing files ..."
	for n, cacheentry in enumerate(cache):
		if verbose and not n % 13:
			print >>sys.stderr, "\t%.1f%%\t(%d files, %d caches)\r" % (100.0 * n / len(cache), n + 1, len(outputcaches)),
		packer.pack(cacheentry)
	if verbose:
		print >>sys.stderr, "\t100.0%%\t(%d files, %d caches)" % (n + 1, len(outputcaches))
		print >>sys.stderr, "sorting output caches ..."
	for cache in outputcaches:
		cache.objects.sort()

	#
	# Done.
	#

	return seglists, outputcaches
