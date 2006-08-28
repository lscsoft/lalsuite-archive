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

from glue import segments
from glue.ligolw import table
from glue.ligolw import lsctables
from pylal import llwapp

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
	searchsummtable = table.get_table(doc, lsctables.SearchSummaryTable.tableName)
	snglbursttable = table.get_table(doc, lsctables.SnglBurstTable.tableName)
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

	llwapp.append_process_params(doc, process, [("--cluster", "lstring", kwargs["cluster"])])

	return process


#
# =============================================================================
#
#                             Clustering Algorithm
#
# =============================================================================
#

def smallest_enclosing_seg(a, b):
	"""
	Return the smallest segment that contains both a and b.
	"""
	return segments.segment(min(a[0], b[0]), max(a[1], b[1]))


def SnglBurstCluster(a, b):
	"""
	Replace a with a cluster constructed from a and b.  The cluster's
	time-frequency tile is the smallest tile that contains the original
	two tiles, and the cluster's other properties are taken from which
	ever is the statistically more confident of the two triggers.
	"""
	# The cluster's frequency band is the smallest band containing the
	# bands of the two original events

	a.set_band(smallest_enclosing_seg(a.get_band(), b.get_band()))

	# The cluster's time interval is the smallest interval containing
	# the intervals of the two original events

	a.set_period(smallest_enclosing_seg(a.get_period(), b.get_period()))

	# The amplitude, SNR, confidence, and peak time of the cluster are
	# those of the most confident of the two events (more negative
	# confidence == more confident).

	if a.confidence > b.confidence:
		a.amplitude = b.amplitude
		a.snr = b.snr
		a.confidence = b.confidence
		a.set_peak(b.get_peak())
		if hasattr(a, "tfvolume"):
			a.tfvolume = b.tfvolume


def ClusterSnglBurstTable(triggers, testfunc, clusterfunc, bailoutfunc = None):
	"""
	Cluster the triggers in the list.  testfunc should accept a pair of
	triggers, and return 0 if they should be clustered.  clusterfunc
	should accept a pair of triggers, and replace the contents of the
	first with a cluster constructed from the two.  If bailoutfunc is
	provided, the triggers will be sorted using testfunc as a
	comparison operator, and then only pairs of triggers for which
	bailoutfunc returns 0 will be considered for clustering.
	"""
	while True:
		did_cluster = False

		if bailoutfunc:
			triggers.sort(testfunc)

		i = 0
		while i < len(triggers):
			j = i + 1
			while j < len(triggers):
				if not testfunc(triggers[i], triggers[j]):
					clusterfunc(triggers[i], triggers[j])
					del triggers[j]
					did_cluster = True
				else:
					if bailoutfunc:
						if bailoutfunc(triggers[i], triggers[j]):
							break
					j += 1
			i += 1

		if not did_cluster:
			return


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
	ClusterSnglBurstTable(snglbursttable, kwargs["testfunc"], kwargs["clusterfunc"], kwargs["bailoutfunc"])

	# Add search summary information
	llwapp.append_search_summary(doc, process, inseg = inseg, outseg = outseg, nevents = len(snglbursttable))
	llwapp.set_process_end_time(process)

	return doc
