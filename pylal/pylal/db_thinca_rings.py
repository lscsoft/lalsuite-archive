# $Id$
#
# Copyright (C) 2009  Kipp Cannon, Chad Hanna
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


from glue import iterutils
from glue import segments
from glue.ligolw import lsctables
from glue.ligolw import dbtables
from glue.ligolw.utils import segments as ligolw_segments
from pylal import SnglInspiralUtils


# FIXME:  remove when Python >= 2.5 required
try:
	any
except NameError:
	any = iterutils.any


#
# =============================================================================
#
#                                  Functions
#
# =============================================================================
#


def get_thinca_rings(connection, program_name = "thinca"):
  """
  Return the thinca rings from the database at the given connection.  The
  rings are returned as a glue.segments.segmentlistdict indexed by the
  instruments that were analyzed in that ring.

  Example:

  >>> ring_sets = get_thinca_rings(connection)
  >>> print ring_sets.keys()
  [frozenset(['H1', 'L1'])]
  """
  ring_sets = segments.segmentlistdict()
  for row in map(lsctables.table.get_table(dbtables.get_xml(connection), lsctables.SearchSummaryTable.tableName)._row_from_cols, connection.cursor().execute("""
SELECT
  search_summary.*
FROM
  search_summary
  JOIN process ON (
    process.process_id == search_summary.process_id
  )
WHERE
  process.program == ?
  AND EXISTS (
    SELECT
      *
    FROM
      process_params
    WHERE
      process_params.process_id == search_summary.process_id
      AND process_params.param == '--num-slides'
  )
  """, (program_name,))):
    available_instruments = frozenset(row.get_ifos())
    try:
      ring_sets[available_instruments].append(row.get_out())
    except KeyError:
      ring_sets[available_instruments] = [row.get_out()]
    return ring_sets


def get_veto_segments(connection, name):
  """
  Return a coalesced glue.segments.segmentlistdict object containing the
  segments of the given name extracted from the database at the given
  connection.
  """
  return ligolw_segments.segmenttable_get_by_name(dbtables.get_xml(connection), name).coalesce()


def get_background_offset_vectors(connection):
  """
  Return a list of the non-zero offset vectors extracted from the database
  at the given connection.  Each offset vector is returned as a dictionary
  mapping instrument name to offset.
  """
  return [offsetvector for offsetvector in lsctables.table.get_table(dbtables.get_xml(connection), lsctables.TimeSlideTable.tableName).as_dict().values() if any(offsetvector.values())]


def get_thinca_livetimes(ring_sets, veto_segments, offset_vectors, verbose = False):
  # FIXME:  somebody should document this
  livetimes = {}
  for available_instruments, rings in ring_sets.items():
    for on_instruments in (combo for m in range(2, len(available_instruments) + 1) for combo in iterutils.choices(sorted(available_instruments), m)):
      if verbose:
        print >>sys.stderr, "%s/%s" % (",".join(on_instruments), ",".join(sorted(available_instruments))),
      on_instruments = frozenset(on_instruments)
      if on_instruments not in livetimes:
        livetimes[on_instruments] = 0
      livetimes[on_instruments] += SnglInspiralUtils.compute_thinca_livetime(on_instruments, available_instruments - on_instruments, rings, veto_segments, offset_vectors)
  return livetimes
