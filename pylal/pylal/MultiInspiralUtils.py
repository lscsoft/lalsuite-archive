# Copyright (C) 2006  Sukanta Bose
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

from pylal.xlal.datatypes.ligotimegps import LIGOTimeGPS
from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw import utils
#
# =============================================================================
#
#                                   Input
#
# =============================================================================
#

def ReadMultiInspiralFromFiles(fileList):
  """
  Read the multiInspiral tables from a list of files
  @param fileList: list of input files
  """
  if not fileList:
    return multiInspiralTable(), None

  multis = None

  for thisFile in fileList:
    doc = utils.load_filename(thisFile,
        gz=(thisFile or "stdin").endswith(".gz"))
    # extract the multi inspiral table
    try:
      multiInspiralTable = table.get_table(doc,
          lsctables.MultiInspiralTable.tableName)
      if multis: multis.extend(multiInspiralTable)
      else: multis = multiInspiralTable
    except: multiInspiralTable = None
  return multis

#
# =============================================================================
#
#                                  Clustering
#
# =============================================================================
#

def CompareMultiInspiralByEndTime(a, b):
  """
  Orders a and b by peak time.
  """
  return cmp(a.get_end(), b.get_end())


def CompareMultiInspiralBySnr(a, b):
  """
  Orders a and b by peak time.
  """
  return cmp(a.snr, b.snr)


def CompareMultiInspiral(a, b, twindow = LIGOTimeGPS(0)):
  """
  Returns 0 if a and b are less than twindow appart.
  """
  tdiff = abs(a.get_end() - b.get_end())
  if tdiff < twindow:
    return 0
  else:
    return cmp(a.get_end(), b.get_end())

