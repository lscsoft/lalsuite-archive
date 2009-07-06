from __future__ import division

import os
import re
import sys
import urllib

### Speed hacks
import glue.segments
import glue.__segments
glue.segments.segment = glue.__segments.segment
glue.segments.segmentlist = glue.__segments.segmentlist
glue.segments.PosInfinity = glue.__segments.PosInfinity
glue.segments.NegInfinity = glue.__segments.NegInfinity
### End speed hacks

from glue.segments import segmentlistdict, segment, segmentlist
from glue import segmentsUtils

__author__ = "Nickolas Fotopoulos (nvf@gravity.phys.uwm.edu)"
__version__ = "$Revision$"[11:-2]
__date__ = "$Date$"

dq_url_pattern = "http://ldas-cit.ligo.caltech.edu/segments/S5/%s/dq_segments.txt"

##############################################################################
# I/O functions

def get_dq_file(ifo, dest=None, verbose=False):
    """
    Connect to a server and download the dq_segments.txt file for the given
    IFO.  Save it to dest.  If dest is None, a randomly named
    temporary file is created.  The path to the destination is returned.
    """
    if verbose:
        print "Downloading %s" % (dq_url_pattern % ifo)

    dest, info = urllib.urlretrieve(dq_url_pattern % ifo, filename=dest)

    if verbose >= 2:
        print info
    return dest

def fromDQsegments(fileobj, coltype=int, discard_disabled=True, keep_versions=None):
    """
    Given a data quality segment file, return a segmentlistdict containing a
    segmentlist for each DQ flag.  Column 0 is the flag name and channel 1
    is the version of that flag.  %(name)s_v%(version)s will be the key on the
    segmentlistdict and columns 2 and 3 will be cast to coltype and used as
    start and end times for segments.

    If discard_disabled==True, then check that the enabled flag (column) is 1.
    """
    d = segmentlistdict()

    for line in fileobj:
        if len(line.strip()) == 0 or line[0] == "#":
            continue
        try:
            tokens = line.split()
            name = tokens[0] + "_v" + tokens[1]
            seg = segment(coltype(tokens[2]), coltype(tokens[3]))
            if (discard_disabled and tokens[4] != "1") or \
               (keep_versions is not None and tokens[1] not in keep_versions):
                continue
        except IndexError:
            raise ValueError, "not enough values on this line:\n" + line
        target_seglist = d.setdefault(name, segmentlist())
        target_seglist.append(seg)
    return d

def fromDQsegments_fast(fileobj):
    """
    Same as fromDQsegments, but the following options are hardcoded for
    speed (17% improvement):
      coltype=int
      discard_disabled=True
      keep_versions=["99"]
    
    Also, any empty or malformed lines will throw an exception.
    """
    d = segmentlistdict()

    for line in fileobj:
        if line[0] == "#":
            continue
        tokens = line.split()
        name = tokens[0] + "_v" + tokens[1]
        seg = segment(int(tokens[2]), int(tokens[3]))
        if (tokens[1] != "99") or (tokens[4] != "1"):
            continue
        target_seglist = d.setdefault(name, segmentlist())
        target_seglist.append(seg)
    return d

def from_veto_file(fileobj,flags_have_ifo_colon=False):
    """
    Return a dictionary keyed by the flag name and mapping to windows to be
    added to the start and end times of DQ-flagged segments.
    """
    veto_window_lines = [line.split() for line in fileobj\
                         if len(line.strip()) > 0 \
                            and not line.startswith("#")]
    if flags_have_ifo_colon: 
      veto_windows = dict([(key.split(':')[1], (int(minus), int(plus))) \
                   for key, minus, plus in veto_window_lines])
    else:
      veto_windows = dict([(key, (int(minus), int(plus))) for key, minus, plus\
                         in veto_window_lines])
    return veto_windows

##############################################################################
# segmentlist/segmentdict utility functions

def apply_veto_windows(seg_dict, veto_window_dict):
    """
    For each DQ flag, widen each flagged segment by the amount specified
    in the veto windows.  veto_window_dict should be keyed by the
    flag name and map to windows to be added to the start and end times.
    """
    veto_dict = segmentlistdict()

    for flag, window in veto_window_dict.iteritems():
        # FIXME: For now we assume that we want the v99 flags
        v99_flag = flag + "_v99"
        if v99_flag not in seg_dict:
            print >>sys.stderr, "warning: %s not in DQ segments file" % v99_flag
            continue
        segs = segmentlist([segment(s[0]+window[0], s[1]+window[1]) 
                for s in seg_dict[v99_flag]])
        segs.coalesce()
        veto_dict[flag] = segs

    return veto_dict

def separate_contributing_flags(segdict):
    """
    Partition a segmentlistdict into zero-length and non-zero-length entries.
    Return a tuple of (keys to zero-length entries, keys to non-zero-length
    entries).
    """
    nocontrib_flags = []
    contrib_flags = []
    for key, val in segdict.iteritems():
        if len(val) > 0:
            contrib_flags.append(key)
        else:
            nocontrib_flags.append(key)
    return contrib_flags, nocontrib_flags

def subtract_segmentlist_from_segmentlistdict(seglist, segdict):
    """
    Return the a dictionary containing the subtraction of the segmentlist
    from every entry in the segmentlistdict.
    """
    def _subtract_from_sublist(sublist):
        return sublist - seglist
    return segmentlistdict(segdict.map(_subtract_from_sublist))

def segmentlistdict_at_point(segdict, point):
    """
    Return a segmentlistdict whose segmentlists contain the segments of
    segdict that contain point.  The output is a segmentlistdict, so entries are
    either an empty segmentlist or contain a length-one segmentlist.
    """
    def find(seglist):
        try:
            ind = seglist.find(point)
            return segmentlist([seglist[ind]])
        except ValueError:
            return segmentlist()
    return segdict.map(find)

##############################################################################
# higher-level functions

def download_and_parse_dq_segs(dq_segfile, ifo, verbose=False):
    """
    Download and parse v99 DQ segments into a segmentlistdict.  Rename DQ flags
    to have the prefix of the ifo and a colon.  Return said segmentlistdict.
    """
    dq_segdict = segmentlistdict()
    
    # download DQ segment dump
    if os.path.isfile(dq_segfile):
      if verbose:
        print "using existing DQ segments file", dq_segfile
    else:
      get_dq_file(ifo, dest=dq_segfile, verbose=verbose)
  
    # parse DQ segments; rename flags to keep track of IFO
    if verbose:
      print "parsing", ifo, "DQ segments..."
    new_dict = fromDQsegments_fast(open(dq_segfile))
    for flag, val in new_dict.iteritems():
        dq_segdict["%s:%s" % (ifo, flag)] = val
    
    # sort each list; important property for segment arithmetic to work
    # correctly
    for seglist in dq_segdict.itervalues():
        seglist.sort()
  
    return dq_segdict

def determine_vetoes(veto_windows, dq_segdict, analyzable):
    """
    Determine all vetoable segments, then reduce to actual segments vetoed.
    Return this list of segments.
    """
    dq_veto_segdict = apply_veto_windows(dq_segdict, veto_windows)
    
    mask = ~analyzable
    vetoed_segdict = subtract_segmentlist_from_segmentlistdict(mask,
      dq_veto_segdict)
    vetoed_segs = vetoed_segdict.union(vetoed_segdict.iterkeys())
    
    return vetoed_segs
