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

"""
This module is a wrapper of the xlal.date module, supplementing the C code
in that module with additional features that are more easily implemented in
Python.  It is recommended that you import this module rather than
importing xlal.date directly.
"""

#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#

from xlal.date import *

__author__ = "Kipp Cannon <kipp@gravity.phys.uwm.edu>"
__date__ = "$Date$"[7:-2]
__version__ = "$Revision$"[11:-2]


#
# =============================================================================
#
#                              Function Wrappers
#
# =============================================================================
#

def XLALGreenwichMeanSiderealTime(gps):
	return XLALGreenwichSiderealTime(gps, 0.0)

def XLALTimeDelayFromEarthCenter(pos, ra, dec, gps):
	return XLALArrivalTimeDiff(pos, [0.0, 0.0, 0.0], ra, dec, gps)


#
# =============================================================================
#
#                                Plotting Tools
#
# =============================================================================
#

def utc_midnight(gps):
	"""
	Truncate a LIGOTimeGPS to UTC midnight.
	"""
	# convert to UTC (as list so we can edit it)
	tm = list(XLALGPSToUTC(gps))

	# truncate to midnight
	tm[3] = 0       # hours
	tm[4] = 0       # minutes
	tm[5] = 0       # seconds

	# convert back to LIGOTimeGPS
	return XLALUTCToGPS(tuple(tm))


class UTCMidnights(object):
	"""
	Iterator for generating LIGOTimeGPS objects for UTC midnights.
	"""
	def __init__(self, start, end):
		"""
		LIGOTimeGPS objects will be returned for each UTC midnight
		in the range [start, end).
		"""
		self.midnight = utc_midnight(start)
		if self.midnight < start:
			self.midnight = utc_midnight(self.midnight + 86402)
		self.end = end

	def __iter__(self):
		return self

	def next(self):
		if self.midnight >= self.end:
			raise StopIteration
		result = self.midnight
		self.midnight = utc_midnight(self.midnight + 86402)
		return result
