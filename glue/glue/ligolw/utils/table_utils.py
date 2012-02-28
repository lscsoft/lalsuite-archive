# Copyright (C) 2012  Matthew West
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
#					Preamble
#
# =============================================================================
#

"""

"""


import sys
from glue import git_version
from glue.ligolw import lsctables

__author__ = "Matt West <matthew.west@ligo.org>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date



#
# =============================================================================
#
#		     Depopulate ligolw_xml tables 
#
# =============================================================================
#

def depopulate_sngl_inspiral(xmldoc, verbose = False):
	"""
	This function takes the lists of event_ids from the sngl_inspiral and coinc_event_map 
	tables and determine the difference, such that the newlist contains only  non-coinc 
	single-ifo triggers. Then it remove these non-coinc triggers from the sngl_inspiral table.
	"""
	sngls_tbl = lsctables.table.get_table(xmldoc, lsctables.SnglInspiralTable.tableName)
	sngls_tbl_eid = list(sngls_tbl.getColumnByName("event_id"))

	coinc_map_tbl = lsctables.table.get_table(xmldoc, lsctables.CoincMapTable.tableName)

	if len(coinc_map_tbl) == 0:
		for i in xrange(len(sngls_tbl_eid) - 1, -1, -1):
			 del sngls_tbl[i]
		if verbose:
			 print >> sys.stderr, "This file lacks any coincident events, so all triggers in the sngl_inspiral table have been removed.",
	else:
		coinc_map_tbl_eid = list(coinc_map_tbl.getColumnByName("event_id"))
		non_coincs = list(set(sngls_tbl_eid) - set(coinc_map_tbl_eid))
		non_coincs.sort(reverse = True)

		for event in non_coincs:
			del sngls_tbl[sngls_tbl_eid.index(event)]
		if verbose:
			print >> sys.stderr, "\n\tremoved %i single-ifo triggers that were not associated with a coincident found in the coinc_event_map table." %( len(non_coincs) )

	return xmldoc

