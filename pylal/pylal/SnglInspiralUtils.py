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
from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw import utils
from glue.ligolw import utils
#
# =============================================================================
#
#                                   Input
#
# =============================================================================
#


def ReadSnglInspiralFromFiles(fileList, mangle_event_id = False):
  """
  Read the snglInspiral tables from a list of files
  @param fileList: list of input files
  """
  snglInspiralTriggers = None
  ncoincs = 0
  for thisFile in fileList:
    doc = utils.load_filename(thisFile, gz=(thisFile or "stdin").endswith(".gz"))
    # extract the sngl inspiral table
    try: snglInspiralTable = \
      table.get_table(doc, lsctables.SnglInspiralTable.tableName)
    except: snglInspiralTable = None

    # if there is an inspiral table, update the event_id to make sure
    # they are unique across all input files
    if snglInspiralTable and mangle_event_id:
      row_dict = {}
      for trig in snglInspiralTable:  # N
        if trig.event_id not in row_dict: # log_2 N
          ncoincs += 1
          row_dict[trig.event_id] = ReassignEventId(trig, ncoincs)
        trig.event_id = row_dict[trig.event_id]

    if snglInspiralTriggers and snglInspiralTable: 
      snglInspiralTriggers.extend(snglInspiralTable)
    elif not snglInspiralTriggers:
      snglInspiralTriggers = snglInspiralTable

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

#
# =============================================================================
#
#                                  Those damn event ids
#
# =============================================================================
#

def DisectEventId(eventid):
  """
  return the three pieces of the event id
  """
  x = eventid // 1000000000
  slidenum = (eventid % 1000000000) // 100000
  y = eventid % 100000

  return x,slidenum,y

def ReassignEventId(a, newid):
  """
  Applies a new event id based on the new id, but preserving the time
  slide information.  Ugh.  This is a hack.  FIXME!  Please.
  """
  newx = 100000000 + (newid // 100000)
  newy = newid % 100000
  x,slidenum,y = DisectEventId(a.event_id)
  return (1000000000 * newx + slidenum * 100000 + newy)

