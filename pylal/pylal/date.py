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

__author__ = "Kipp Cannon <kipp@gravity.phys.uwm.edu>"
__date__ = "$Date$"[7:-2]
__version__ = "$Revision$"[11:-2]

import math

from glue import segmentsUtils
from xlal.date import *


#
# =============================================================================
#
#                              Function Wrappers
#
# =============================================================================
#

def XLALGreenwichMeanSiderealTime(gps):
	return XLALGreenwichSiderealTime(gps, 0.0)

def XLALGreenwichSiderealTimeToGPS(gmst, equation_of_equinoxes):
	return XLALGreenwichMeanSiderealTimeToGPS(gmst) - equation_of_equinoxes

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


def UTCMidnights(start, end):
	"""
	Iterator for generating LIGOTimeGPS objects for UTC midnights.
	"""
	midnight = utc_midnight(start)
	if midnight < start:
		midnight = utc_midnight(midnight + 86402)
	while midnight < end:
		yield midnight
		midnight = utc_midnight(midnight + 86402)


def gmst_0h(gps):
	"""
	Truncate a LIGOTimeGPS to Greenwich mean sidereal 0 rad.
	"""
	gmst = XLALGreenwichMeanSiderealTime(gps)
	residual = gmst % (2.0 * math.pi)
	if residual:
		gmst -= residual
	return XLALGreenwichMeanSiderealTimeToGPS(gmst)


def GMST_0hs(start, end):
	"""
	Iterator for generating LIGOTimeGPS objects for Greenwich Mean
	Sidereal 0h.
	"""
	gmst = XLALGreenwichMeanSiderealTime(start)
	residual = gmst % (2.0 * math.pi)
	if residual:
		gmst -= residual
	if gmst < start:
		gmst += 2.0 * math.pi
	end = XLALGreenwichMeanSiderealTime(end)
	while gmst < end:
		yield XLALGreenwichMeanSiderealTimeToGPS(gmst)
		gmst += 2.0 * math.pi


#
# =============================================================================
#
# Segment Lists
#
# =============================================================================
#

def gmst_days(gps_start, gps_stop):
	"""
	Generates a segmentlist whose segments are the Greenwich Mean
	Sidereal days spanning the given range of GPS times.  Input and
	output times are all GPS seconds.
	"""
	gmst_0hs = GMST_0hs(gmst_0h(gps_start), gps_stop + 86402)
	l = segments.segmentlist([segments.segment(gmst_0hs.next(), gmst_0hs.next())])
	while l[-1][1] < gps_stop:
		l.append(segments.segment(l[-1][1], gmst_0hs.next()))
	return l
