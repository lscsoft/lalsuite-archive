from __future__ import division

import os
import re
import socket
import sys
import tempfile
import time
import urlparse
itertools = __import__("itertools")  # absolute import of system-wide itertools

from lalapps import inspiral

from glue import lal
from glue import pipeline
from glue import segments, segmentsUtils
from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw import utils

# from pylal import webCondor

##############################################################################
# Custom classes
##############################################################################

class GRBSummaryDAG(pipeline.CondorDAG):
  def __init__(self, config_file, log_path):
    self.basename = config_file.replace(".ini", "") 
    logfile, logfilename = tempfile.mkstemp(prefix=self.basename, suffix=".dag.log", dir=log_path)
    os.close(logfile)
    pipeline.CondorDAG.__init__(self, logfilename)
    self.set_dag_file(self.basename)

##############################################################################
# Utility functions
##############################################################################

def compute_masked_segments(analyzable_seglist, on_source_segment,
    veto_seglist=None, quantization_time=None):
    """
    Return veto segmentlists for on-source and off-source regions,
    respectively.  Optionally, use vetos from veto_seglist.  Optionally,
    quantize the off-source with quantization_time (seconds).
    """
    analyzable_seglist = segments.segmentlist(analyzable_seglist[:]).coalesce()
    if veto_seglist is None:
        veto_seglist = segments.segmentlist()
    off_source_segs = analyzable_seglist - segments.segmentlist([on_source_segment])

    ## on-source mask
    on_source_mask = off_source_segs | veto_seglist

    ## off-source mask
    # first, assign its value without quantization
    off_source_mask = segments.segmentlist([on_source_segment]) | veto_seglist

    # then, quantize as necessary
    if quantization_time is not None:
        off_source_quantized = segments.segmentlist(
            [segments.segment(s[0], s[0] + abs(s)//quantization_time) \
             for s in (off_source_segs - off_source_mask)])
        off_source_mask = analyzable_seglist - off_source_quantized

    return on_source_mask, off_source_mask

def compute_offsource_segment(analyzable, on_source, padding_time=0,
    max_trials=None, symmetric=True):
    """
    Return the longest symmetric protraction of the on_source, constrained
    to fall inside analyzable times, with the proctraction time minus
    padding_time divisible by quantization_time (the length of on_source).
    Return a protraction of exactly (max_trials//2) * quantization_time if
    max_trials is specified. Return None if there are less than max_trials
    quantization_times available.
    """
    quantization_time = abs(on_source)
    
    try:
        super_seg = analyzable[analyzable.find(on_source)].contract(padding_time)
    except ValueError:
        return None
    
    # check again after taking padding into account
    if on_source not in super_seg:
        return None
    
    nplus = (super_seg[1] - on_source[1]) // quantization_time
    nminus = (on_source[0] - super_seg[0]) // quantization_time
    
    if (max_trials is not None) and (nplus + nminus > max_trials):
        # try to make this as centered as possible
        half_max = max_trials // 2
        if nplus < half_max:  # cut left
            remainder = max_trials - nplus
            nminus = min(remainder, nminus)
        elif nminus < half_max:  # cut right
            remainder = max_trials - nminus
            nplus = min(remainder, nplus)
        else:  # cut both
            nplus = min(half_max, nplus)
            nminus = min(half_max, nminus)
    
    if symmetric:
        nplus = nminus = min(nplus, nminus)
    
    return segments.segment((on_source[0] - nminus*quantization_time - padding_time,
                             on_source[1] + nplus*quantization_time + padding_time))

def ext_trigger_gpstimes_from_xml(doc):
    """
    Return a dictionary of GPS times of external triggers keyed by the GRB
    name found in the ExtTriggersTables present in doc.  If there are no
    ExtTriggersTables in doc, return None.
    """
    ext_triggers_tables = lsctables.getTablesByType(doc, lsctables.ExtTriggersTable)
    if ext_triggers_tables is None:
        return None
    ext_triggers = {}
    for tab in ext_triggers_tables:
        for name, time in itertools.izip(tab.getColumnByName("event_number_grb"),
                                         tab.getColumnByName("start_time")):
            if name in ext_triggers:
                print >>sys.stderr, "warning: GRB %s appears twice in document; taking second definition"
            ext_triggers[name] = time
    return ext_triggers
    