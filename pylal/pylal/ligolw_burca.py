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
from pylal import snglcoinc

__author__ = "Kipp Cannon <kipp@gravity.phys.uwm.edu>"
__version__ = "$Revision$"[11:-2]
__date__ = "$Date$"[7:-2]


#
# =============================================================================
#
#                           Add Process Information
#
# =============================================================================
#

def append_process(xmldoc, **kwargs):
	process = llwapp.append_process(xmldoc, program = "ligolw_burca", version = __version__, cvs_repository = "lscsoft", cvs_entry_time = __date__, comment = kwargs["comment"])

	params = [("--program", "lstring", kwargs["program"])]
	for key, value in kwargs["window"].iteritems():
		if key[0] < key[1]:
			params += [("--window", "lstring", "%s,%s=%s" % (key[0], key[1], str(value)))]
	llwapp.append_process_params(xmldoc, process, params)

	return process


#
# =============================================================================
#
#                                 Library API
#
# =============================================================================
#

def ligolw_burca(xmldoc, **kwargs):
	# add an entry in the process table
	process = append_process(xmldoc, **kwargs)

	if kwargs["verbose"]:
		print >>sys.stderr, "indexing..."

	# prepare the coincidence table interface
	coinc_tables = snglcoinc.CoincTables(xmldoc, [lsctables.SnglBurstTable.tableName])

	# build the event list accessors, populated with events from those
	# processes that can participate in a coincidence
	eventlists = snglcoinc.make_eventlists(xmldoc, lsctables.SnglBurstTable.tableName, kwargs["window"], kwargs["program"])

	# iterate over time slides
	time_slide_ids = coinc_tables.time_slide_ids()
	for n, time_slide_id in enumerate(time_slide_ids):
		offsetdict = coinc_tables.get_time_slide(time_slide_id)
		if kwargs["verbose"]:
			print >>sys.stderr, "time slide %d/%d: %s" % (n + 1, len(time_slide_ids), str(offsetdict))
		if False in map(eventlists.__contains__, offsetdict.keys()):
			if kwargs["verbose"]:
				print >>sys.stderr, "\twarning: skipping due to insufficient data"
			continue
		if kwargs["verbose"]:
			print >>sys.stderr, "\tapplying time offsets ..."
		eventlists.set_offsetdict(offsetdict)
		if kwargs["verbose"]:
			print >>sys.stderr, "\tsearching ..."
		# search for and record coincidences
		for ntuple in snglcoinc.CoincidentNTuples(eventlists, offsetdict.keys(), kwargs["window"], kwargs["verbose"]):
			coinc_tables.append_coinc(process.process_id, time_slide_id, ntuple)

	# clean up and finish
	eventlists.remove_offsetdict()
	llwapp.set_process_end_time(process)
	return xmldoc
