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

#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#

from glue import segments
from pylal.date import LIGOTimeGPS


#
# =============================================================================
#
#                                  Clustering
#
# =============================================================================
#

def cmp_segs(a, b):
	"""
	Returns 1 if a covers an interval above b's interval, -1 if a
	covers an interval below b's, and 0 if the two intervals overlap
	(including if their edges touch).
	"""
	if a[0] > b[1]:
		return 1
	if a[1] < b[0]:
		return -1
	return 0


def CompareSnglBurstByPeakTime(a, b):
	"""
	Orders a and b by peak time.
	"""
	return cmp(a.get_peak(), b.get_peak())


def CompareSnglBurstByPeakTimeAndFreq(a, b):
	"""
	Orders a and b by peak time, then by frequency band.  Returns 0 if
	a and b have the same peak time, and their frequency bands
	intersect.
	"""
	return cmp(a.get_peak(), b.get_peak()) or cmp_segs(a.get_band(), b.get_band())


def CompareSnglBurst(a, b, twindow = LIGOTimeGPS(0)):
	"""
	Orders a and b by time interval, then by frequency band.  Returns 0
	if a and b's time-frequency tiles intersect.  A time window can be
	optionally applied, and the time-frequency tiles will continue to
	compare as equal if they do not overlap by as much as the window
	amount.
	"""
	return cmp_segs(a.get_period().protract(twindow), b.get_period()) or cmp_segs(a.get_band(), b.get_band())


#
# =============================================================================
#
#                              Injection Related
#
# =============================================================================
#

def CompareSimBurstAndSnglBurstByTime(sim, burst):
	"""
	Return True if the peak time of the injection sim lies within the
	time interval of burst.
	"""
	if sim.coordinates == "ZENITH":
		return sim.get_geocent_peak() in burst.get_period()
	else:
		return sim.get_peak(burst.ifo) in burst.get_period()

def CompareSimBurstAndSnglBurstByTimeandFreq(sim, burst):
	"""
	Return True if the peak time and centre frequency of sim lie within
	the time-frequency tile of burst.
	"""
	return CompareSimBurstAndSnglBurstByTime(sim, burst) and (sim.freq in burst.get_band())

