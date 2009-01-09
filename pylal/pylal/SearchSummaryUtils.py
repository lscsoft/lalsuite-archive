# $Id$
#
# Copyright (C) 2006  Craig Robinson and Anand Sengupta
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

from glue import segments

from glue.ligolw import utils
from pylal import llwapp

def GetSegListFromSearchSummaries(fileList, verbose=False):
  """
  Read segment lists from search summary tables
  @param fileList: list of input files.
  """
  segList = segments.segmentlistdict()

  for thisFile in fileList:
    doc = utils.load_filename(thisFile, gz = thisFile.endswith(".gz"),
        verbose = verbose)
    try: 
      segs = llwapp.segmentlistdict_fromsearchsummary(doc)
    except:
      raise ValueError, "Cannot extract segments from the SearchSummaryTable of %s" % thisFile

    #Now add these segments to the existing list
    segList.extend(segs)

  for value in segList.values():
    value.sort()

  return segList
