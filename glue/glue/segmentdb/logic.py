# $Id$
#
# Copyright (C) 2009  Larne Pekowsky
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
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
Utility methods for doing logical operations on sets of segments

TODO: The boundries between these methods and bin/ligolw_segment_union,
bin/ligolw_segment_intersect are not as clean as I would like.  Stuff
like registering the program to the prcess table should really go in the
respective mains, but because of the way the documents are handled it
is simpler to put them here.
"""


import sys
import os
import glue.segments

from glue.ligolw import ligolw
from glue.ligolw import table
from glue.ligolw import lsctables

from glue.ligolw.utils import ligolw_add
from glue.ligolw.utils import process

from glue.segmentdb.segmentdb_utils import add_to_segment_definer
from glue.segmentdb.segmentdb_utils import add_to_segment
from glue.segmentdb.segmentdb_utils import find_segments

PROGRAM_NAME = sys.argv[0].replace('./','')
PROGRAM_PID  = os.getpid()
try:
        USER_NAME = os.getlogin()
except:
        USER_NAME = pwd.getpwuid(os.getuid())[0]


__author__  = "Larne Pekowsky <lppekows@physics.syr.edu>"
__date__    = "$Date$"[7:-2]
__version__ = "$Revision$"[11:-2]


def run_file_operation(options, filenames, operation):
    """
    Performs an operation (intersect or union) across a set of files.
    That is, given a set of files each with segment definers DMT-FLAG1,
    DMT-FLAG2 etc the result is a file where 

    DMT-FLAG1 = (file 1's DMT-FLAG1 operation file 2's DMT-FLAG1 operation ...)
    DMT-FLAG2 = (file 1's DMT-FLAG2 operation file 2's DMT-FLAG2 operation ...)
    
    etc
    """
    # load it up with the input files
    if options.preserve:
        resultdoc = ligolw_add.ligolw_add(ligolw.Document(), filenames)
    else:
        resultdoc = ligolw.Document()
        resultdoc.appendChild(ligolw.LIGO_LW())

    # Register ourselves
    proc_id = process.register_to_xmldoc(resultdoc, PROGRAM_NAME, options.__dict__).process_id

    # load up the files into individual documents
    xmldocs = [ligolw_add.ligolw_add(ligolw.Document(), [fname]) for fname in filenames]

    # Get the list of dinstinct segment_definers across all docs
    segment_definers = {}

    def register_definer(seg_def):
        key = (seg_def.ifos, seg_def.name, seg_def.version)
        segment_definers[key] = True
        return key

    for xmldoc in xmldocs:
        seg_def_table = table.get_table(xmldoc, lsctables.SegmentDefTable.tableName)
        map (register_definer, seg_def_table)

    # For each unique segment definer, find the intersection
    for ifo, name, version in segment_definers:
        if operation == 'intersect':
            result = glue.segments.segmentlist([glue.segments.segment(-glue.segments.infinity(), glue.segments.infinity())])
        else:
            result = glue.segments.segmentlist([])

        for xmldoc in xmldocs:
            if operation == 'intersect':
                result &= find_segments(xmldoc, '%s:%s:%d' % (ifo, name, version))
            else:
                result |= find_segments(xmldoc, '%s:%s:%d' % (ifo, name, version))


        result.coalesce()

        # Add a segment definer for the result
        seg_def_id = add_to_segment_definer(resultdoc, proc_id, ifo, name, version)

        # Add the segments
        add_to_segment(resultdoc, proc_id, seg_def_id, result)

    return resultdoc





def run_segment_operation(options, filenames, operation):
    """
    Performs an operation (intersect or union) across a set of segments.
    That is, given a set of files each with segment definers DMT-FLAG1,
    DMT-FLAG2 etc and a list of segments DMT-FLAG1,DMT-FLAG1 this returns

    RESULT = (table 1's DMT-FLAG1 union table 2's DMT-FLAG1 union ...)
             operation
             (table 1's DMT-FLAG2 union table 2's DMT-FLAG2 union ...)
             operation
    etc
    """
    # Create a temporary doc.  This will also be the output if we're
    # preserving
    xmldoc = ligolw_add.ligolw_add(ligolw.Document(), filenames)

    # Start with a segment covering all of time, then
    # intersect with each of the fields of interest
    if operation == 'intersect':
        sgmntlist = glue.segments.segmentlist([glue.segments.segment(-glue.segments.infinity(), glue.segments.infinity())])
    else:
        sgmntlist = glue.segments.segmentlist([])

    for key in options.segments.split(','):
        if operation == 'intersect':
            sgmntlist &= find_segments(xmldoc, key)
        else:
            sgmntlist |= find_segments(xmldoc, key)

    # Now that we have the result, if we're not preserving start a fresh document
    if not options.preserve:
        xmldoc = ligolw.Document()
        xmldoc.appendChild(ligolw.LIGO_LW())
        # Could also reset the ids so we start counting at 0

    # Register ourselves
    proc_id = process.register_to_xmldoc(xmldoc, PROGRAM_NAME, options.__dict__).process_id

    # Add a segment definer for the result
    seg_def_id = add_to_segment_definer(xmldoc, proc_id, "", options.result_name or 'RESULT', 1)

    # Add the segments
    add_to_segment(xmldoc, proc_id, seg_def_id, sgmntlist)

    return xmldoc

