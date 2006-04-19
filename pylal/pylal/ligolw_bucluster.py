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

import sys

from glue.ligolw import lsctables
from pylal import llwapp
from pylal import SnglBurstUtils

__author__ = "Kipp Cannon <kipp@gravity.phys.uwm.edu>"
__version__ = "$Revision$"[11:-2]
__date__ = "$Date$"[7:-2]


#
# =============================================================================
#
#                                 Preparation
#
# =============================================================================
#

def get_tables(doc):
	searchsummtable = llwapp.get_table(doc, lsctables.SearchSummaryTable.tableName)
	snglbursttable = llwapp.get_table(doc, lsctables.SnglBurstTable.tableName)
	return searchsummtable.get_inlist().extent(), searchsummtable.get_outlist().extent(), snglbursttable


#
# =============================================================================
#
#                           Add Process Information
#
# =============================================================================
#

def append_process(doc, **kwargs):
	process = llwapp.append_process(doc, program = "ligolw_bucluster", version = __version__, cvs_repository = "lscsoft", cvs_entry_time = __date__, comment = kwargs["comment"])

	params = [("--cluster", "lstring", kwargs["cluster"])]
	if kwargs["input"] != None:
		params += [("--input", "lstring", kwargs["input"])]
	if kwargs["output"] != None:
		params += [("--output", "lstring", kwargs["output"])]
	llwapp.append_process_params(doc, process, params)

	return process


#
# =============================================================================
#
#                                 Library API
#
# =============================================================================
#

def ligolw_bucluster(doc, **kwargs):
	# Extract segments and tables
	inseg, outseg, snglbursttable = get_tables(doc)

	# Add process information
	process = append_process(doc, **kwargs)

	# Cluster
	if kwargs["verbose"]:
		print >>sys.stderr, "clustering..."
	SnglBurstUtils.ClusterSnglBurstTable(snglbursttable.rows, kwargs["testfunc"], kwargs["clusterfunc"], kwargs["bailoutfunc"])

	# Add search summary information
	llwapp.append_search_summary(doc, process, inseg = inseg, outseg = outseg, nevents = len(snglbursttable))
	llwapp.set_process_end_time(process)

	return doc
