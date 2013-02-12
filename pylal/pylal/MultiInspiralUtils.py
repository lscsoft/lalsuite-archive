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
from glue.ligolw import ilwd

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

def ReadMultiInspiralTimeSlidesFromFiles(fileList):
  """
  Read time-slid multiInspiral tables from a list of files
  @param fileList: list of input files
  """
  if not fileList:
    return multiInspiralTable(), None

  multis = None
  timeSlides = []

  for thisFile in fileList:

    doc = utils.load_filename(thisFile,
        gz=(thisFile or "stdin").endswith(".gz"))
    # Extract the time slide table
    timeSlideTable = table.get_table(doc,
          lsctables.TimeSlideTable.tableName)
    slideMapping = {}
    currSlides = {}
    for slide in timeSlideTable:
      currID = int(slide.time_slide_id)
      if currID not in currSlides.keys():
        currSlides[currID] = {}
        currSlides[currID][slide.instrument] = slide.offset
      elif slide.instrument not in currSlides[currID].keys():
        currSlides[currID][slide.instrument] = slide.offset

    for slideID,offsetDict in currSlides.items():
      try:
        # Is the slide already in the list and where?
        offsetIndex = timeSlides.index(offsetDict)
        slideMapping[slideID] = offsetIndex
      except ValueError:
        # If not then add it
        timeSlides.append(offsetDict)
        slideMapping[slideID] = len(timeSlides) - 1
    
    # extract the multi inspiral table
    try:
      multiInspiralTable = table.get_table(doc,
          lsctables.MultiInspiralTable.tableName)
      # Remap the time slide IDs
      for multi in multiInspiralTable:
        newID = slideMapping[int(multi.time_slide_id)]
        multi.time_slide_id = ilwd.ilwdchar(\
                              "multi_inspiral:time_slide_id:%d" % (newID))
      if multis: multis.extend(multiInspiralTable)
      else: multis = multiInspiralTable
#    except: multiInspiralTable = None
    except: raise

  # Make a new time slide table
  timeSlideTab = lsctables.New(lsctables.TimeSlideTable)

  for slideID,offsetDict in enumerate(timeSlides):
    for instrument in offsetDict.keys():
      currTimeSlide = lsctables.TimeSlide()
      currTimeSlide.instrument = instrument
      currTimeSlide.offset = offsetDict[instrument]
      currTimeSlide.time_slide_id = ilwd.ilwdchar(\
                              "time_slide:time_slide_id:%d" % (slideID))
      currTimeSlide.process_id = ilwd.ilwdchar(\
                              "time_slide:process_id:%d" % (0))
      timeSlideTab.append(currTimeSlide)

  return multis,timeSlides,timeSlideTab


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

