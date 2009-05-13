#!/usr/bin/env python
#
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

import os
import re

import glue.segments
from glue.ligolw import lsctables
from glue.ligolw import table
from glue.segmentdb import query_engine

def get_all_files_in_range(dirname, starttime, endtime):
    """Returns all files in dirname and all its subdirectories whose
    names indicate that they contain segments in the range starttime
    to endtime"""
    
    ret = []

    # Maybe the user just wants one file...
    if os.path.isfile(dirname):
        return [dirname]

    first_four_start = starttime / 100000
    first_four_end   = endtime   / 100000

    for filename in os.listdir(dirname):
        if re.match('.*-[0-9]{4}$', filename):
            dirtime = int(filename[-4:])
            if dirtime >= first_four_start and dirtime <= first_four_end:
                ret += get_all_files_in_range(os.path.join(dirname,filename), starttime, endtime)
        elif re.match('.*-[0-9]*-[0-9]*\.xml', filename):
            file_time = int(filename.split('-')[-2])
            if file_time >= (starttime-64) and file_time <= (endtime+64):
                ret.append(os.path.join(dirname,filename))

    return ret



def setup_database(host_and_port):
    """Opens a connection to a LDBD Server"""
    global PROGRAM_NAME

    port = 30020
    
    if host_and_port.find(':') < 0:
        host = host_and_port
    else:
        # server and port specified
        host, portString = host_and_port.split(':')
        port = int(portString)


    identity = "/DC=org/DC=doegrids/OU=Services/CN=ldbd/%s" % host

    # open connection to LDBD Server
    client = None

    try:
        client = LDBDClient.LDBDClient(host, port, identity)
    except Exception, e:
        print >>sys.stderr, \
              "Unable to connect to LDBD Server %s:%d" % (host, port)
        if gsiserverutils.checkCredentials():
            print >>sys.stderr, "Got the following error : " + str(e)
            print >>sys.stderr, "Enter '%s --help' for usage" % PROGRAM_NAME
        sys.exit(-1)

    return client


def build_segment_list(engine, gps_start_time, gps_end_time, ifo, segment_name, version = None, start_pad = 0, end_pad = 0):
    """Optains a list of segments for the given ifo, name and version between the
    specified times.  If a version is given the request is straightforward and is
    passed on to build_segment_list_one.  Otherwise more complex processing is
    performed (not yet implemented)"""
    if version is not None:
        return build_segment_list_one(engine, gps_start_time, gps_end_time, ifo, segment_name, version, start_pad, end_pad)

    # This needs more sophisticated logic, for the moment just return the latest
    # available version
    sql  = "SELECT max(version) FROM segment_definer "
    sql += "WHERE  segment_definer.ifos = '%s' " % ifo
    sql += "AND   segment_definer.name = '%s' " % segment_name

    rows = engine.query(sql)
    version = len(rows[0]) and rows[0][0] or 1

    return build_segment_list_one(engine, gps_start_time, gps_end_time, ifo, segment_name, version, start_pad, end_pad)


def build_segment_list_one(engine, gps_start_time, gps_end_time, ifo, segment_name, version = None, start_pad = 0, end_pad = 0):
    """Builds a list of segments satisfying the given criteria """
    seg_result = glue.segments.segmentlist([])
    sum_result = glue.segments.segmentlist([])

    # Is there any way to get segment and segement summary in one query?
    # Maybe some sort of outer join where we keep track of which segment
    # summaries we've already seen.
    sql = "SELECT segment_summary.start_time, segment_summary.end_time "
    sql += "FROM segment_definer, segment_summary "
    sql += "WHERE segment_summary.segment_def_id = segment_definer.segment_def_id "
    sql += "AND   segment_definer.ifos = '%s' " % ifo
    sql += "AND   segment_definer.name = '%s' " % segment_name
    sql += "AND   segment_definer.version = %s " % version
    sql += "AND NOT (%s > segment_summary.end_time OR segment_summary.start_time > %s)" % (gps_start_time, gps_end_time)

    rows = engine.query(sql)

    for sum_start_time, sum_end_time in rows:
        sum_start_time = (sum_start_time < gps_start_time) and gps_start_time or sum_start_time
        sum_end_time = (sum_end_time > gps_end_time) and gps_end_time or sum_end_time

        sum_result |= glue.segments.segmentlist([glue.segments.segment(sum_start_time, sum_end_time)])

    # We can't use queries paramaterized with ? since the ldbd protocol doesn't support it...
    sql = "SELECT segment.start_time + %d, segment.end_time + %d " % (start_pad, end_pad)
    sql += "FROM segment, segment_definer "
    sql += "WHERE segment.segment_def_id = segment_definer.segment_def_id "

    if engine.__class__ == query_engine.LdbdQueryEngine:
       sql += "AND segment.segment_def_cdb = segment_definer.creator_db "
    sql += "AND   segment_definer.ifos = '%s' " % ifo
    sql += "AND   segment_definer.name = '%s' " % segment_name
    sql += "AND   segment_definer.version = %s " % version
    sql += "AND NOT (%s > segment.end_time OR segment.start_time > %s)" % (gps_start_time, gps_end_time)

    rows = engine.query(sql)
    
    for seg_start_time, seg_end_time in rows:
        seg_start_time = (seg_start_time < gps_start_time) and gps_start_time or seg_start_time
        seg_end_time = (seg_end_time > gps_end_time) and gps_end_time or seg_end_time

        seg_result |= glue.segments.segmentlist([glue.segments.segment(seg_start_time, seg_end_time)])

    engine.close()

    return sum_result, seg_result



def run_query_segments(doc, proc_id, engine, gps_start_time, gps_end_time, included_segments_string, excluded_segments_string = None, write_segments = True, start_pad = 0, end_pad = 0):
    """Runs a segment query.  This was originally part of ligolw_query_segments, but now is also
    used by ligolw_segments_from_cats.

    The write_segments option is provided so callers can coalesce segments obtained over
    sever invocations (as segments_from_cats does).
    """
    
    # Find or create the tables we'll need
    try:
        seg_def_table  = table.get_table(doc, lsctables.SegmentDefTable.tableName)
        new_seg_def_id = seg_def_table[0].segment_def_id       
    except ValueError:
        seg_def_table = lsctables.New(lsctables.SegmentDefTable,
                                      columns = ["process_id", "segment_def_id", "ifos", "name", "version", "comment"])
        doc.childNodes[0].appendChild(seg_def_table)

        # Add ourselves as a segment definer
        new_seg_def_id                 = seg_def_table.get_next_id()
        segment_definer                = lsctables.SegmentDef()
        segment_definer.process_id     = proc_id
        segment_definer.segment_def_id = new_seg_def_id
        segment_definer.ifos           = split_segment_ids(included_segments_string.split(','))[0][0]
        segment_definer.name           = "result"
        segment_definer.version        = 0
        segment_definer.comment        = ''
        
        seg_def_table.append(segment_definer)

    try:
        seg_sum_table = table.get_table(doc, lsctables.SegmentSumTable.tableName)
    except:
        seg_sum_table = lsctables.New(lsctables.SegmentSumTable,
                                      columns = ['process_id','segment_def_id','segment_sum_id','start_time', 'end_time'])
        doc.childNodes[0].appendChild(seg_sum_table)

        # with a segment summary spanning the interval from the command line
        segment_summary                = lsctables.SegmentSumTable()
        segment_summary.process_id     = proc_id
        segment_summary.segment_def_id = new_seg_def_id
        segment_summary.segment_sum_id = seg_sum_table.get_next_id()
        segment_summary.start_time     = gps_start_time
        segment_summary.end_time       = gps_end_time
        
        seg_sum_table.append(segment_summary)


    try:
        segtable = table.get_table(doc, lsctables.SegmentTable.tableName)
    except:
        segtable = lsctables.New(lsctables.SegmentTable,
                                 columns = ["process_id", "segment_def_id", "segment_id", "start_time", "end_time"])
        doc.childNodes[0].appendChild(segtable)


    included_segments = glue.segments.segmentlist([])
    excluded_segments = glue.segments.segmentlist([])


    for ifo, segment_name, version in split_segment_ids(included_segments_string.split(',')):
        sum_segments, seg_segments = build_segment_list(engine, gps_start_time, gps_end_time, ifo, segment_name, version, start_pad, end_pad)

        seg_def_id                     = seg_def_table.get_next_id()
        segment_definer                = lsctables.SegmentDef()
        segment_definer.process_id     = proc_id
        segment_definer.segment_def_id = seg_def_id
        segment_definer.ifos           = ifo
        segment_definer.name           = segment_name
        segment_definer.version        = version
        segment_definer.comment        = ''

        seg_def_table.append(segment_definer)

        # Add segment summary entries
        for sum_segment in sum_segments:
            segment_summary                = lsctables.SegmentSumTable()
            segment_summary.process_id     = proc_id
            segment_summary.segment_def_id = seg_def_id
            segment_summary.segment_sum_id = seg_sum_table.get_next_id()
            segment_summary.start_time     = sum_segment[0]
            segment_summary.end_time       = sum_segment[1]

            seg_sum_table.append(segment_summary)

        # and accumulate segments 
        included_segments |= seg_segments

    # Excluded segments are not required
    if excluded_segments_string:
        for ifo, segment_name, version in split_segment_ids(excluded_segments_string.split(',')):
            sum_segments, seg_segments = build_segment_list(engine, gps_start_time, gps_end_time, ifo, segment_name, version)
            excluded_segments |= seg_segments


    result = included_segments - excluded_segments
    result.coalesce()
    
    # Add the segments
    if write_segments:
        write_segments_to_xml(doc, proc_id, result)

    return result


def write_segments_to_xml(doc, proc_id, result):
    """Write out the segments in result to the doc"""
    
    seg_def_table  = table.get_table(doc, lsctables.SegmentDefTable.tableName)
    new_seg_def_id = seg_def_table[0].segment_def_id       

    try:
        segtable = table.get_table(doc, lsctables.SegmentTable.tableName)
    except:
        segtable = lsctables.New(lsctables.SegmentTable,
                                 columns = ["process_id", "segment_def_id", "segment_id", "start_time", "end_time"])
        doc.childNodes[0].appendChild(segtable)

    for seg in result:
        segment                = lsctables.Segment()
        segment.process_id     = proc_id
        segment.segment_def_id = new_seg_def_id
        segment.segment_id     = segtable.get_next_id()
        segment.start_time     = seg[0]
        segment.end_time       = seg[1]

        segtable.append(segment)


def split_segment_ids(segment_ids):
    """Given an array of strings of the form ifo:name and
    ifo:name:version, returns an array of tuples of the form (ifo,
    name, version) where version may be None"""
    
    def split_segment_id(segment_id):
        temp = segment_id.split(':')
        if len(temp) == 2:
            temp.append(None)
        elif temp[2] == '*':
            temp[2] = None
        else:
            temp[2] = int(temp[2])
            
        return temp

    return map(split_segment_id, segment_ids)




def find_segments(doc, key):
    key_pieces = key.split(':')
    while len(key_pieces) < 3:
        key_pieces.append('*')

    filter_func = lambda x: str(x.ifos) == key_pieces[0] and (str(x.name) == key_pieces[1] or key_pieces[1] == '*') and (str(x.version) == key_pieces[2] or key_pieces[2] == '*') 

    # Find all segment definers matching the critieria
    seg_def_table = table.get_table(doc, lsctables.SegmentDefTable.tableName)
    seg_defs      = filter(filter_func, seg_def_table)
    seg_def_ids   = map(lambda x: str(x.segment_def_id), seg_defs)

    # Find all segments belonging to those definers
    seg_table     = table.get_table(doc, lsctables.SegmentTable.tableName)
    seg_entries   = filter(lambda x: str(x.segment_def_id) in seg_def_ids, seg_table)

    # Combine into a segmentlist
    ret = glue.segments.segmentlist(map(lambda x: glue.segments.segment(x.start_time, x.end_time), seg_entries))

    ret.coalesce()

    return ret


def add_to_segment_definer(xmldoc, proc_id, ifo, name, version):
    try:
        seg_def_table = table.get_table(xmldoc, lsctables.SegmentDefTable.tableName)
    except:
        seg_def_table = lsctables.New(lsctables.SegmentDefTable, columns = ["process_id", "segment_def_id", "ifos", "name", "version", "comment"])
        xmldoc.childNodes[0].appendChild(seg_def_table)

    seg_def_id                     = seg_def_table.get_next_id()
    segment_definer                = lsctables.SegmentDef()
    segment_definer.process_id     = proc_id
    segment_definer.segment_def_id = seg_def_id
    segment_definer.ifos           = ifo
    segment_definer.name           = name
    segment_definer.version        = version
    segment_definer.comment        = ''

    seg_def_table.append(segment_definer)

    return seg_def_id



def add_to_segment(xmldoc, proc_id, seg_def_id, sgmtlist):
    try:
        segtable = table.get_table(xmldoc, lsctables.SegmentTable.tableName)
    except:
        segtable = lsctables.New(lsctables.SegmentTable, columns = ["process_id", "segment_def_id", "segment_id", "start_time", "end_time"])
        xmldoc.childNodes[0].appendChild(segtable)

    for seg in sgmtlist:
        segment                = lsctables.Segment()
        segment.process_id     = proc_id
        segment.segment_def_id = seg_def_id
        segment.segment_id     = segtable.get_next_id()
        segment.start_time     = seg[0]
        segment.end_time       = seg[1]

        segtable.append(segment)


