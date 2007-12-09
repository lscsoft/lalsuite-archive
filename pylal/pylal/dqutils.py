from __future__ import division

import re
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

    If discard_disabled==True, then check that the enabled flag (column) is >0.
    """
    fivecolsegpat = re.compile(r"\A(\w+)\s+([\d]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s*\Z")
    d = segmentlistdict({})

    for line in fileobj:
        if len(line) == 0 or line.startswith("#"):
            continue
        try:
            tokens = fivecolsegpat.match(line).groups()
            name = tokens[0] + "_v" + tokens[1]
            seg = segment(coltype(tokens[2]), coltype(tokens[3]))
            enabled = int(tokens[4])
        except ValueError:
            break
        if (discard_disabled and enabled <= 0) or \
           (keep_versions is not None and tokens[1] not in keep_versions):
            continue
        target_seglist = d.setdefault(name, segmentlist())
        target_seglist.append(seg)
    return d

def from_veto_file(fileobj):
    """
    Return a dictionary keyed by the flag name and mapping to windows to be
    added to the start and end times of DQ-flagged segments.
    """
    veto_window_lines = [line.split() for line in fileobj\
                         if len(line) > 0 and not line.startswith("#")]
    veto_windows = dict([(key, (minus, plus)) for key, minus, plus\
                         in veto_window_lines])
    return veto_windows

def apply_veto_windows(seg_dict, veto_window_dict):
    """
    For each DQ flag, widen each flagged segment by the amount specified
    in the veto windows.  veto_window_dict should be keyed by the
    flag name and map to windows to be added to the start and end times.
    """
    veto_dict = segmentlistdict()

    for flag, window in veto_window_dict.iteritems():
        # FIXME: For now we assume that we want the v99 flags
        segs = segmentlist([segment(s[0]+window[0], s[1]+window[1]) 
                for s in seg_dict[flag + "_v99"]])
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
    return segdict.map(_subtract_from_sublist)

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

def determine_vetoes(dq_segdict, analyzable_seglist):
    """
    Return a new segmentlistdict like dq_segdict, but whose segmentlists
    are within analyzable times.  These are generally the relevant vetoes.
    """
    mask = ~analyzable_seglist
    return subtract_segmentlist_from_segmentlistdict(mask, dq_segdict)