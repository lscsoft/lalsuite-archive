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

import xlal.date

__author__ = "Kipp Cannon <kipp@gravity.phys.uwm.edu>"
__date__ = "$Date$"[7:-2]
__version__ = "$Revision$"[11:-2]


#
# =============================================================================
#
#                               LIGOTimeGPS Type
#
# =============================================================================
#

class LIGOTimeGPS(xlal.date.LIGOTimeGPS):
	def __init__(self, seconds, nanoseconds = None):
		if type(seconds) in [int, long]:
			if nanoseconds:
				if type(nanoseconds) not in [int, long]:
					raise TypeError, repr(nanoseconds)
				xlal.date.LIGOTimeGPS.__init__(self, seconds, nanoseconds)
			else:
				xlal.date.LIGOTimeGPS.__init__(self, seconds, 0)
		elif nanoseconds != None:
			raise TypeError, "LIGOTimeGPS(x): function takes exactly 1 argument with non-integer x (2 given)"
		elif type(seconds) == float:
			t = XLALREAL8ToGPS(seconds)
			self.seconds, self.nanoseconds = t.seconds, t.nanoseconds
		elif type(seconds) == str:
			t = XLALStrToGPS(seconds)
			self.seconds, self.nanoseconds = t.seconds, t.nanoseconds
		else:
			try:
				self.seconds, self.nanoseconds = seconds.seconds, seconds.nanoseconds
			except:
				raise ValueError, repr(seconds)

	def __abs__(self):
		return LIGOTimeGPS(xlal.date.LIGOTimeGPS.__abs__(self))

	def __add__(self, other):
		return LIGOTimeGPS(xlal.date.LIGOTimeGPS.__add__(self, other))

	#def __coerce__(self, other):
	#	if type(other) != LIGOTimeGPS:
	#		try:
	#			return (self, LIGOTimeGPS(other))
	#		except:
	#			pass
	#	return (self, other)

	def __div__(self, other):
		return LIGOTimeGPS(0, XLALGPSToINT8NS(self) / other)

	def __mod__(self, other):
		return LIGOTimeGPS(0, XLALGPSToINT8NS(self) % (other * 1000000000L))

	def __mul__(self, other):
		return LIGOTimeGPS(0, XLALGPSToINT8NS(self) * other)

	def __neg__(self):
		return LIGOTimeGPS(xlal.date.LIGOTimeGPS.__neg__(self))

	def __pos__(self):
		return LIGOTImeGPS(self)

	def __reduce__(self):
		return (LIGOTimeGPS, (self.seconds, self.nanoseconds))

	def __repr__(self):
		return "LIGOTimeGPS(%d,%d)" % (self.seconds, self.nanoseconds)

	__rmul__ = __mul__

	def __str__(self):
		return ("%d.%09d" % (self.seconds, abs(self.nanoseconds))).rstrip("0").rstrip(".")

	def __sub__(self, other):
		return LIGOTimeGPS(xlal.date.LIGOTimeGPS.__sub__(self, other))


#
# =============================================================================
#
#                              Function Wrappers
#
# =============================================================================
#

def XLALREAL8ToGPS(*args):
	return LIGOTimeGPS(xlal.date.XLALREAL8ToGPS(*args))

XLALGPSToINT8NS = xlal.date.XLALGPSToINT8NS

def XLALINT8NSToGPS(*args):
	return LIGOTimeGPS(xlal.date.XLALINT8NSToGPS(*args))

def XLALStrToGPS(*args):
	return LIGOTimeGPS(xlal.date.XLALStrToGPS(*args))

XLALLeapSeconds = xlal.date.XLALLeapSeconds

XLALLeapSecondsUTC = xlal.date.XLALLeapSecondsUTC

XLALGPSToUTC = xlal.date.XLALGPSToUTC

def XLALUTCToGPS(*args):
	return LIGOTimeGPS(xlal.date.XLALUTCToGPS(*args))

XLALJulianDay = xlal.date.XLALJulianDay

XLALModifiedJulianDay = xlal.date.XLALModifiedJulianDay

XLALGreenwichSiderealTime = xlal.date.XLALGreenwichSiderealTime

XLALArrivalTimeDiff = xlal.date.XLALArrivalTimeDiff

def XLALGreenwichMeanSiderealTime(gps):
	return XLALGreenwichSiderealTime(gps, 0.0)

def XLALTimeDelayFromEarthCenter(pos, ra, dec, gps):
	return XLALArrivalTimeDiff(pos, [0.0, 0.0, 0.0], ra, dec, gps)
