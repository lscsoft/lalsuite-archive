# $Id$
#
# Copyright (C) 2006  Duncan A. Brown
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

from pylal.date import LIGOTimeGPS
from glue.ligolw import ligolw
from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw import utils
from glue.ligolw.utils import ligolw_add

#
# =============================================================================
#
#                                   Input
#
# =============================================================================
#


def ReadSnglInspiralFromFiles(fileList, mangle_event_id=False, verbose=False):
  """
  Read the SnglInspiralTables from a list of files

  @param fileList: list of input files
  @param mangle_event_id: ID remapping is necessary in cases where multiple
    files might have event_id collisions (ex: exttrig injections)
  @param verbose: print ligolw_add progress
  """
  # turn on ID remapping if necessary
  if mangle_event_id:
    lsctables.SnglInspiralTable.ids = lsctables.SnglInspiralIDs_old()

  # ligolw_add will merge all tables, which is overkill, but merge time is
  # much less than I/O time.
  xmldoc = ligolw_add.ligolw_add(ligolw.Document(), fileList, verbose=verbose)

  # extract the SnglInspiral table
  try:
    snglInspiralTriggers = table.get_table(xmldoc, \
      lsctables.SnglInspiralTable.tableName)
  except:
    snglInspiralTriggers = None

  # return ID remapping to its normal state (off)
  lsctables.SnglInspiralTable.ids = None

  return snglInspiralTriggers


#
# =============================================================================
#
#                                  Clustering
#
# =============================================================================
#

def CompareSnglInspiralByEndTime(a, b):
  """
  Orders a and b by peak time.
  """
  return cmp(a.get_end(), b.get_end())


def CompareSnglInspiralBySnr(a, b):
  """
  Orders a and b by peak time.
  """
  return cmp(a.snr, b.snr)


def CompareSnglInspiral(a, b, twindow = LIGOTimeGPS(0)):
  """
  Returns 0 if a and b are less than twindow appart.
  """
  tdiff = abs(a.get_end() - b.get_end())
  if tdiff < twindow:
    return 0
  else:
    return cmp(a.get_end(), b.get_end())
