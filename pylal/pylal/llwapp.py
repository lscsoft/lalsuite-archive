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
data in LIGO Light-Weight XML format.
"""


import os
import pickle
import socket
import time
import warnings


from glue import segments
from glue.ligolw import ligolw
from glue.ligolw import table
from glue.ligolw import param
from glue.ligolw import lsctables
from glue.ligolw.utils import process as ligolwprocess
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


def get_time_slide_id(xmldoc, time_slide, create_new = None):
	"""
	Return the time_slide_id corresponding to the time slide described
	by time_slide, a dictionary of instrument/offset pairs.  If the
	document does not contain exactly 1 time_slide table then
	ValueError is raised.  If the table does not describe a matching
	time slide then KeyError is raised.  If, however, the optional
	create_new argument is set to an lsctables.Process object (or any
	other object with a process_id attribute), then a time slide table
	will be created if needed and or the missing time slide description
	added to the table, and indicated as having been created by the
	given process.
	"""
	try:
		tisitable = table.get_table(xmldoc, lsctables.TimeSlideTable.tableName)
	except ValueError, e:
		# table not found
		if create_new is None:
			raise e
		tisitable = lsctables.New(lsctables.TimeSlideTable)
		xmldoc.childNodes[0].appendChild(tisitable)
	# make sure the next_id attribute is correct
	tisitable.sync_next_id()
	# get the id
	return tisitable.get_time_slide_id(time_slide, create_new = create_new)


def get_coinc_def_id(xmldoc, search, coinc_type, create_new = True, description = u""):
	"""
	Wrapper for the get_coinc_def_id() method of the CoincDefiner table
	class in glue.ligolw.lsctables.  This wrapper will optionally
	create a new coinc_definer table in the document if one does not
	already exist.
	"""
	try:
		coincdeftable = table.get_table(xmldoc, lsctables.CoincDefTable.tableName)
	except ValueError, e:
		# table not found
		if not create_new:
			raise e
		# FIXME:  doesn't work if the document is stored in a
		# database.
		coincdeftable = lsctables.New(lsctables.CoincDefTable)
		xmldoc.childNodes[0].appendChild(coincdeftable)
	# make sure the next_id attribute is correct
	coincdeftable.sync_next_id()
	# get the id
	return coincdeftable.get_coinc_def_id(search, coinc_type, create_new = create_new, description = description)


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
	return param.from_pyvalue(u"pickle:%s" % name, unicode(pickle.dumps(obj)))


def pickle_from_param(elem, name):
	"""
	Retrieve a pickled Python object from the document tree rooted at
	elem.
	"""
	return pickle.loads(str(param.get_pyvalue(elem, u"pickle:%s" % name)))


def yaml_to_param(obj, name):
	"""
	Return the top-level element of a document sub-tree containing the
	YAML serialization of a Python object.
	"""
	import yaml
	return param.from_pyvalue(u"yaml:%s" % name, unicode(yaml.dump(obj)))


def yaml_from_param(elem, name):
	"""
	Retrieve a YAMLed Python object from the document tree rooted at
	elem.
	"""
	import yaml
	return yaml.load(param.get_pyvalue(elem, u"yaml:%s" % name))


#
# =============================================================================
#
#                               Process Metadata
#
# =============================================================================
#


def append_process(*args, **kwargs):
	"""
	Identical to the append_process() function in
	glue.ligolw.utils.process except uses LAL to convert UTC to GPS
	time to get the leap seconds correct.
	"""
	process = ligolwprocess.append_process(*args, **kwargs)
	# FIXME:  remove the "" case when the git metadata business is
	# sorted out
	if "cvs_entry_time" in kwargs and kwargs["cvs_entry_time"] is not None and kwargs["cvs_entry_time"] != "":
		process.cvs_entry_time = XLALUTCToGPS(time.strptime(kwargs["cvs_entry_time"], "%Y/%m/%d %H:%M:%S")).seconds
	process.start_time = XLALUTCToGPS(time.gmtime()).seconds
	return process


def set_process_end_time(process):
	"""
	Identical to the set_process_end_time() function in
	glue.ligolw.utils.process except uses LAL to convert UTC to GPS
	time to get the leap seconds correct.
	"""
	process.end_time = XLALUTCToGPS(time.gmtime()).seconds
	return process


def append_process_params(*args, **kwargs):
	"""
	Deprecated.  Use glue.ligolw.utils.process.append_process_params
	instead.
	"""
	warnings.warn("function pylal.llwapp.append_process_params is deprecated, use glue.ligolw.utils.process.append_process_params instead", DeprecationWarning, stacklevel = 2)
	return ligolwprocess.append_process_params(*args, **kwargs)


def get_process_params(*args, **kwargs):
	"""
	Deprecated.  Use glue.ligolw.utils.process.get_process_params()
	instead.
	"""
	warnings.warn("function pylal.llwapp.get_process_params is deprecated, use glue.ligolw.utils.process.get_process_params instead", DeprecationWarning, stacklevel = 2)
	return ligolwprocess.get_process_params(*args, **kwargs)


def doc_includes_process(*args, **kwargs):
	"""
	Deprecated.  Use glue.ligolw.utils.doc_includes_process() instead.
	"""
	warnings.warn("function pylal.llwapp.doc_includes_process is deprecated, use glue.ligolw.utils.process.doc_includes_process instead", DeprecationWarning, stacklevel = 2)
	return ligolwprocess.doc_includes_process(*args, **kwargs)


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
	if ifos is not None:
		summary.set_ifos(ifos)
	else:
		summary.set_ifos(process.get_ifos())
	summary.set_in(inseg)
	summary.set_out(outseg)
	summary.nevents = nevents
	summary.nnodes = nnodes
	table.get_table(doc, lsctables.SearchSummaryTable.tableName).append(summary)
	return summary


def segmentlistdict_fromsearchsummary(xmldoc, program = None):
	"""
	Convenience wrapper for a common case usage of the segmentlistdict
	class:  searches the process table in xmldoc for occurances of a
	program named program, then scans the search summary table for
	matching process IDs and constructs a segmentlistdict object from
	the out segments in those rows.

	Note:  the segmentlists in the segmentlistdict are not necessarily
	coalesced, they contain the segments as they appear in the
	search_summary table.
	"""
	stbl = table.get_table(xmldoc, lsctables.SearchSummaryTable.tableName)
	ptbl = table.get_table(xmldoc, lsctables.ProcessTable.tableName)
	return stbl.get_out_segmentlistdict(program and ptbl.get_ids_by_program(program))


#
# =============================================================================
#
#                                    Other
#
# =============================================================================
#


def get_coincident_segmentlistdict(seglistdict, offsetdictlist):
	"""
	Compute the segments for which data is required in order to perform
	a complete coincidence analysis given the segments for which data
	is available and the list of offset vectors to be applied to the
	data during the coincidence analysis.

	seglistdict is a segmentlistdict object defining the instruments
	and times for which data is available.  offsetdictlist is a list of
	offset vectors to be applied to the data --- dictionaries of
	instrument/offset pairs.

	The offset vectors in offsetdictlist are applied to the input
	segments one by one and the interesection of the shifted segments
	is computed.  The segments surviving the intersection are unshifted
	to their original positions and stored.  The return value is the
	union of the results of this operation.

	In all cases the full n-way intersection is computed, that is if an
	offset vector lists three instruments then this function returns
	the times when exactly all three of those isntruments are on.  If
	the calling code requires times when any two of the three are on
	the list of offset vectors should be pre-processed to indicate this
	by listing the allowed instrument combinations as separate offset
	vectors.  See ligolw_tisi.time_slide_component_vectors() for a
	function to assist in doing this.

	For example, let us say that "input" is a segmentlistdict object
	containing segment lists for three instruments, "H1", "H2" and
	"L1".  And let us say that "slides" is a list of dictionaries, and
	is equal to [{"H1":0, "H2":0, "L1":0}, {"H1":0, "H2":10}].  Then if

	output = get_coincident_segmentlistdict(input, slides)

	output will contain, for each of the three instruments, the
	segments (or parts thereof) from the original lists that are
	required in order to perform a triple-coincident analysis at zero
	lag betwen the three instruments, *and* a double-coincident
	analysis between H1 and H2 with H2 offset by 10 seconds.

	During the computations, the input segmentlistdict object will have
	offsets applied to it in place, but they will be restored to their
	original values upon exit.  The segmentlistdict object returned by
	this function has its offsets set to those of the input
	segmentlistdict.
	"""
	origoffsets = dict(seglistdict.offsets)
	coincseglists = segments.segmentlistdict()
	for offsetdict in offsetdictlist:
		seglistdict.offsets.update(offsetdict)
		intersection = seglistdict.extract_common(offsetdict.keys())
		intersection.offsets.clear()
		coincseglists |= intersection
	seglistdict.offsets.update(origoffsets)
	coincseglists.offsets.update(origoffsets)
	return coincseglists
