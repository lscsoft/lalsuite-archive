# Copyright (C) 2009  Nickolas Fotopoulos
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

from __future__ import division

__author__ = "Nickolas Fotopoulos <nvf@gravity.phys.uwm.edu>"

import os
import os.path
import shutil
import sys

from glue.segments import segment, segmentlist
from pylal.metaarray import TimeSeries, TimeSeriesList
from pylal.Fr import frgetvect1d

__all__ = ('__author__', 'FrameCache')

class FrameCache(object):
    """
FrameCache is a transparent interface to LSC data. The user provides a LAL-
formatted cache file and the returned FrameCache object allows repeated
queries for channels and time, even across frame files. It also supports
smart, lazy local caching. Limitations: It only works for one-dimensional
time-series data.

Constructor:
    FrameCache(cache_entries=None, scratchdir=None, verbose=False)

Inputs:
    cache is a list of glue.lal.CacheEntry objects or a glue.lal.Cache.
    Data will be retrieved from the frame files described within.

    Scratchdir determines where to locally cache frames.  If None, no
    caching is performed.

Example:
>>> from glue import lal
>>> from pylal import frutils
>>> c = lal.Cache.fromfile(open("test.cache"))
>>> d = frutils.FrameCache(c, scratchdir="/tmp", verbose=True)
>>> data = d.fetch("H1:LSC-STRAIN", 861417967, 861417969)
Copying /Users/nvf/temp/H-H1_RDS_C03_L2-861417967-128.gwf -->
          /tmp/H-H1_RDS_C03_L2-861417967-128.gwf.
>>> print data
[  1.68448009e-16   1.69713183e-16   1.71046196e-16 ...,   1.80974629e-16
   1.80911765e-16   1.80804879e-16] {'dt': 6.103515625e-05, 'segments': [segment(861417967, 861417969)], 'comments': [], 'name': 'H1:LSC-STRAIN'}
>>> exit()
Removing /tmp/H-H1_RDS_C03_L2-861417967-128.gwf.

"""

    def __init__(self, cache_entries=None, scratchdir=None, verbose=False):
        """ Initializes interface to frame data.  See .__class__.__doc__"""

        # Simple initializations
        # Use list of segments vs segmentlist to prevent merging.
        self._verbose = verbose
        self._scratchdir = scratchdir
        self._remotefiles = []               # filename list
        self._remotesegs = segmentlist()     # list of segments
        self._remotecoverage = segmentlist() # coalesced copy of remotesegs

        # if we have a scratchdir, maintain independent lists
        if scratchdir is not None:
            self._cachedfiles = []
            self._cachedsegs = segmentlist()
            self._cachecoverage = segmentlist()
        else:
            self._cachedfiles = self._remotefiles
            self._cachedsegs = self._remotesegs
            self._cachecoverage = self._remotecoverage

        if cache_entries is not None:
          self.add_cache(cache_entries)

    def add_cache(self, cache_entries):
        """
        Add information from some cache entries.
        """
        newfiles = [entry.path() for entry in cache_entries \
                    if entry.path() not in self._remotefiles]
        newsegs = segmentlist([entry.segment for entry in cache_entries])
        self._remotefiles.extend(newfiles)
        self._remotesegs.extend(newsegs)
        self._remotecoverage |= segmentlist(newsegs)
        self._remotecoverage.coalesce()

    def __del__(self):
        """
        Clear cache in local scratch.
        """
        if self._scratchdir is None:
            return
        for f,s in zip(self._cachedfiles, self._cachedsegs):
            self._unfetch(f, s)
        return

    def fetch(self, channel, start, end):
        """
        Retrieve data, caching file locations and the files themselves.
        """
        seg = segment(start, end)

        if seg not in self._remotecoverage:
            raise ValueError, "%s not found in cache" % repr(seg)

        # Need to cache files locally
        # Note: seg *will* be in self._cachecoverage if self.scratchdir is None.
        if seg not in self._cachecoverage:
            for f,s in zip(self._remotefiles, self._remotesegs):
                if seg.intersects(s) and s not in self._cachecoverage:
                    dest = os.path.join(self._scratchdir, os.path.split(f)[-1])
                    if self._verbose:
                        print "Copying %s -->\n          %s." % (f, dest)
                    shutil.copy(f, dest)
                    self._cachedfiles.append(dest)
                    self._cachedsegs.append(s)
                    self._cachecoverage |= segmentlist([s])
            assert seg in self._cachecoverage

        # Finally, return the cached data
        return self._fetch(channel, start, end)

    def _fetch(self, channel, start, end, comments=[]):
        """
        Internal method to actually retrieve and return data as TimeSeries,
        assuming that self._framefiles is all set.  Does not check boundaries.
        """
        toreturn = TimeSeriesList([])

        if start==end:
            return toreturn

        # Find first frame
        try:
            index = self._cachedsegs.find(start)
        except ValueError:
            print >>sys.stderr, "Couldn't find any frame files to cover",\
                str(start),"to",str(end),"among:"
            print >>sys.stderr, str(framefilelist)
            return toreturn

        # Get frames; an error probably means that the frames didn't cover
        # the whole period of time.  Cleanly handles frames of varying lengths.
        now = start
        while now < end:
            dur = min(end, self._cachedsegs[index][1]) - now
            data, GPS_start, t_low, dt, x_unit, y_unit = \
                frgetvect1d(self._cachedfiles[index], channel, now, dur, 0)
            meta = {"name": channel, "dt": dt,
                "segments": [segment(now, now+dur)], "comments": comments}
            toreturn.append(TimeSeries(data, meta))
            now += dur
            index += 1

        if len(toreturn) == 0:
            print >>sys.stderr, "This print statement should never execute."
            print >>sys.stderr,"Couldn't find all frame files needed to cover",\
                str(start), "to", str(end), "among:"
            print >>sys.stderr, str(self._cachedfiles)

        toreturn = toreturn.merge_list()
        toreturn.metadata.segments.coalesce()

        return toreturn

    def unfetch(self, start, end):
        """
        Removes files from local scratch space based on start, end
        pairs.  Silently ignores non-existent times.  Remove if file end
        is between start and end.  This is biased to prevent cache misses
        for future fetches being in the future.  (Processing frames in
        chronological order)
        """
        if self._scratchdir is None:
            return

        for f,s in zip(self._cachedfiles, self._cachedsegs):
            if start < s[1] <= end:
                self._unfetch(f,s)

    def _unfetch(self, filename, seg):
        """
        Internal method to actually remove a file from cache.
        """
        if self._scratchdir is None:
            return
        if filename not in self._cachedfiles:
            print >>sys.stderr, \
                "Cache inconsistency: Delete request for file not in cache."
            return
        if self._verbose: print "Removing %s." % filename
        os.remove(filename)
        self._cachedfiles.remove(filename)
        self._cachedsegs.remove(seg)
        self._cachecoverage -= segmentlist([seg])
        return
