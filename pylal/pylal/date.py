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
from xlal.date import LIGOTimeGPS as _LIGOTimeGPS

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

class LIGOTimeGPS(_LIGOTimeGPS):
	def __init__(self, seconds, nanoseconds = None):
		if type(seconds) in [int, long]:
			if nanoseconds:
				_LIGOTimeGPS.__init__(self, seconds, nanoseconds)
			else:
				_LIGOTimeGPS.__init__(self, seconds, 0)
		elif nanoseconds != None:
			raise TypeError, "cannot set nanoseconds in combination with non-integer seconds"
		elif type(seconds) == LIGOTimeGPS:
			self.seconds, self.nanoseconds = seconds.seconds, seconds.nanoseconds
		elif type(seconds) == float:
			XLALGPSSetREAL8(self, seconds)
		elif type(seconds) == str:
			t = XLALStrToGPS(seconds)
			self.seconds, self.nanoseconds = t.seconds, t.nanoseconds
		else:
			raise TypeError, seconds

	def __mul__(self, other):
		return LIGOTimeGPS(0, XLALGPSToINT8NS(self) * other)

	__rmul__ = __mul__

	def __div__(self, other):
		return LIGOTimeGPS(0, XLALGPSToINT8NS(self) / other)

	def __mod__(self, other):
		return LIGOTimeGPS(0, XLALGPSToINT8NS(self) % (other * 1000000000L))


#
# =============================================================================
#
#                               Time Conversion
#
# =============================================================================
#

def XLALGreenwichMeanSiderealTime(gps):
	return XLALGreenwichSiderealTime(gps, 0.0)

def XLALTimeDelayFromEarthCenter(pos, ra, dec, gps):
	return XLALArrivalTimeDiff(pos, [0.0, 0.0, 0.0], ra, dec, gps)
