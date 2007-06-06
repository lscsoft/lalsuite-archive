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

from glue.ligolw import utils
from pylal import llwapp

def GetSegListFromSearchSummary(fileList):
  """
  Read segment lists from search summary tables
  @param fileList: list of input files.
  """
  segList = None

  for thisFile in fileList:
    doc = utils.load_filename(thisFile)
    try: segs = llwapp.segmentlistdict_fromsearchsummary(doc)
    except: segs = None

    #Now add these segments to the existing list
    if segs and segList:
      segList.union(segs)
    elif not segList:
      segList = segs

  return segList
