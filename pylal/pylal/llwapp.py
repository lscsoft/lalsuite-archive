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
A collection of utilities to assist in writing applications that manipulate
data in LIGO Light-Weight format.
"""

import bisect
import os
import pickle
import socket
import time
import urllib

from glue import segments
from glue.ligolw import ligolw
from glue.ligolw import table
from glue.ligolw import param
from glue.ligolw import array
from glue.ligolw import lsctables
from glue.ligolw import utils
from pylal.date import XLALUTCToGPS

__author__ = "Kipp Cannon <kipp@gravity.phys.uwm.edu>"
__version__ = "$Revision$"[11:-2]
__date__ = "$Date$"[7:-2]


#
# =============================================================================
#
#                                    Tables
#
# =============================================================================
#

def segmentlistdict_fromsearchsummary(xmldoc, live_time_program = None):
	"""
	Convenience wrapper for a common case usage of the segmentlistdict
	class: searches the process table in xmldoc for occurances of a
	program named live_time_program, then scans the search summary
	table for matching process IDs and constructs a segmentlistdict
	object from those rows.
	"""
	stbl = table.get_table(xmldoc, lsctables.SearchSummaryTable.tableName)
	ptbl = table.get_table(xmldoc, lsctables.ProcessTable.tableName)
	return stbl.get_out_segmentlistdict(ptbl.get_ids_by_program(program))


def get_time_slide_id(xmldoc, time_slide, create_new = None):
	"""
	Return the time_slide_id corresponding to the time slide described
	by time_slide, a dictionary of instrument/offset pairs.  If the
	document does not contain exactly 1 time_slide table, or the table
	does not describe a matching time slide, then KeyError is raised.
	If, however, the optional create_new argument is set to an
	lsctables.Process object (or any other object with a process_id
	attribute), then a time slide table will be created if needed and
	or the missing time slide description added to the table, and
	indicated as having been created by the given process.
	"""
	try:
		tisitable = table.get_table(xmldoc, lsctables.TimeSlideTable.tableName)
	except ValueError:
		# table not found
		if create_new is None:
			raise KeyError, time_slide
		tisitable = lsctables.New(lsctables.TimeSlideTable)
		xmldoc.childNodes[0].appendChild(tisitable)
	for id in tisitable.dict.iterkeys():
		if tisitable.get_offset_dict(id) == time_slide:
			break
	else:
		# time slide not found in table
		if create_new is None:
			raise KeyError, time_slide
		id = tisitable.sync_ids().next()
		for instrument, offset in time_slide.iteritems():
			row = lsctables.TimeSlide()
			row.process_id = create_new.process_id
			row.time_slide_id = id
			row.instrument = instrument
			row.offset = offset
			tisitable.append(row)
	return id


def get_coinc_def_id(xmldoc, table_names, create_new = True):
	"""
	Return the coinc_def_id corresponding to coincidences consisting
	exclusively of events from the given table names.  If no matching
	coinc_def_id is found, then a new one is created and the ID
	returned.  If the document does not contain a coinc_definer table,
	then one is added, a new coind_def_id created, and the ID returned.
	If, however, create_new is False, and for any reason the ID isn't
	found then KeyError is raised.
	"""
	try:
		coincdeftable = table.get_table(xmldoc, lsctables.CoincDefTable.tableName)
	except ValueError:
		# table not found
		if not create_new:
			raise KeyError, table_names
		coincdeftable = lsctables.New(lsctables.CoincDefTable)
		xmldoc.childNodes[0].appendChild(coincdeftable)
	table_names.sort()
	for id in coincdeftable.dict.iterkeys():
		if coincdeftable.get_contributors(id) == table_names:
			break
	else:
		# contributor list not found in table
		if not create_new:
			raise KeyError, table_names
		id = coincdeftable.sync_ids().next()
		for name in table_names:
			row = lsctables.CoincDef()
			row.coinc_def_id = id
			row.table_name = name
			coincdeftable.append(row)
	return id


#
# =============================================================================
#
#                                    Params
#
# =============================================================================
#

def pickle_to_param(obj, name):
	"""
	Return the top-level element of a document sub-tree containing the
	pickled serialization of a Python object.
	"""
	# WARNING:  this is an abuse of most of the things involved;  but,
	# hey, it's a one-liner so it's hard to resist.  Do NOT use this in
	# files you care about.  If there proves to be a real desire for
	# this kind of thing, then we need to think of a proper way to do
	# it.  For example, there is at least one library that can
	# serialize a Python object into proper XML.
	return param.new_param("pickle:%s" % name, "lstring", urllib.quote(pickle.dumps(obj)))


def pickle_from_param(elem, name):
	"""
	Retrieve a pickled Python object from the document tree rooted at
	elem.
	"""
	return pickle.loads(urllib.unquote(param.get_param(elem, "pickle:%s:param" % name).pcdata))


#
# =============================================================================
#
#                               Process Metadata
#
# =============================================================================
#

def append_process(doc, program = "", version = "", cvs_repository = "", cvs_entry_time = "", comment = "", is_online = False, jobid = 0, domain = "", ifos = ""):
	"""
	Add an entry to the process table in doc.  program, version,
	cvs_repository, comment, domain, and ifos should all be strings.
	cvs_entry_time should be a string in the format "YYYY/MM/DD
	HH:MM:SS".  is_online should be a boolean, jobid an integer.
	"""
	proctable = table.get_table(doc, lsctables.ProcessTable.tableName)
	process = lsctables.Process()
	process.program = program
	process.version = version
	process.cvs_repository = cvs_repository
	process.cvs_entry_time = XLALUTCToGPS(time.strptime(cvs_entry_time, "%Y/%m/%d %H:%M:%S")).seconds
	process.comment = comment
	process.is_online = int(is_online)
	process.node = socket.gethostbyaddr(socket.gethostname())[0]
	process.username = os.environ["LOGNAME"]
	process.unix_procid = os.getpid()
	process.start_time = XLALUTCToGPS(time.gmtime()).seconds
	process.end_time = 0
	process.jobid = jobid
	process.domain = domain
	process.ifos = ifos
	process.process_id = proctable.sync_ids().next()
	proctable.append(process)
	return process


def set_process_end_time(process):
	"""
	Set the end time in a row in a process table to the current time.
	"""
	process.end_time = XLALUTCToGPS(time.gmtime()).seconds
	return process


def append_process_params(doc, process, params):
	"""
	doc is an XML document tree, process is the row in the process
	table for which these are the parameters, and params is a list of
	(name, type, value) tuples one for each parameter.
	"""
	paramtable = table.get_table(doc, lsctables.ProcessParamsTable.tableName)
	for name, type, value in params:
		param = lsctables.ProcessParams()
		param.program = process.program
		param.process_id = process.process_id
		param.param = str(name)
		param.type = str(type)
		param.value = str(value)
		paramtable.append(param)
	return process


def doc_includes_process(doc, program):
	"""
	Return True if the process table in doc includes entries for a
	program named program.
	"""
	return program in table.get_table(doc, lsctables.ProcessTable.tableName).getColumnByName("program")


#
# =============================================================================
#
#                               Search Metadata
#
# =============================================================================
#

def append_search_summary(doc, process, shared_object = "standalone", lalwrapper_cvs_tag = "", lal_cvs_tag = "", comment = None, ifos = None, inseg = None, outseg = None, nevents = 0, nnodes = 1):
	"""
	Append search summary information associated with the given process
	to the search summary table in doc.
	"""
	summary = lsctables.SearchSummary()
	summary.process_id = process.process_id
	summary.shared_object = shared_object
	summary.lalwrapper_cvs_tag = lalwrapper_cvs_tag
	summary.lal_cvs_tag = lal_cvs_tag
	summary.comment = comment or process.comment
	summary.ifos = ifos or process.ifos
	summary.set_in(inseg)
	summary.set_out(outseg)
	summary.nevents = nevents
	summary.nnodes = nnodes
	table.get_table(doc, lsctables.SearchSummaryTable.tableName).append(summary)
	return summary


#
# =============================================================================
#
#                                    Other
#
# =============================================================================
#

def bisect_contains(array, val):
	"""
	Uses a bisection search to determine if val is in array.  Returns
	True or False.
	"""
	try:
		return array[bisect.bisect_left(array, val)] == val
	except IndexError:
		return False


def get_coincident_segmentlistdict(seglistdict, offsetdictlist):
	"""
	This function answers the question "Given a set of segment lists,
	and a set of time slides to apply to those segment lists, what
	segments do I need to keep in the original lists so that I have all
	the segments that will participate in a coincidence analysis done
	over those time slides?

	This function constructs and returns a segmentlistdict object that,
	for each key in seglistdict, contains the segments from the
	corresponding list in seglistdict which are coincident under at
	least one of the time slides described by offsetdictlist.

	offsetlistdict is a list of dictionaries of instrument/offset
	pairs, with each dictionary describing a time slide and the
	instruments that participate in it.  Each element in the list is
	free to contain only subsets of the keys in seglistdict.  In those
	cases, the coincidence is computed only between the segment lists
	corresponding to the given keys.

	For example, let us say that "input" is a segmentlistdict object
	containing segment lists for three instruments, "H1", "H2" and
	"L1".  And let us say that "slides" is a list of dictionaries, and
	is equal to [{"H1":0, "H2":0, "L1":0}, {"H1":0, "H2":10}].  Then if

	output = get_coincident_segmentlistdict(input, slides)

	output will contain, for each instrument, the segments (or parts
	thereof) from the original lists that are required in order to
	perform a triple-coincident analysis at zero lag betwen the three
	instruments, *and* a double-coincident analysis between H1 and H2
	with H2 offset by 10 seconds.

	During the computations, the input segmentlistdict object will have
	offsets applied to it in place, but they will be restored to their
	original value upon exit.
	"""
	origoffsets = {}
	origoffsets.update(seglistdict.offsets)
	coincseglists = segments.segmentlistdict()
	for offsetdict in offsetdictlist:
		seglistdict.offsets.update(offsetdict)
		intersection = seglistdict.extract_common(offsetdict.keys())
		intersection.offsets.update(origoffsets)
		coincseglists |= intersection
	seglistdict.offsets.update(origoffsets)
	return coincseglists
