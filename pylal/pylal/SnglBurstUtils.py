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


def smallest_enclosing_seg(a, b):
	"""
	Return the smallest segment that contains both a and b.
	"""
	return segments.segment(min(a[0], b[0]), max(a[1], b[1]))


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


def SnglBurstCluster(a, b):
	"""
	Replace a with a cluster constructed from a and b.  The cluster's
	time-frequency tile is the smallest tile that contains the original
	two tiles, and the cluster's other properties are taken from which
	ever is the statistically more confident of the two triggers.
	"""
	# The cluster's frequency band is the smallest band containing the
	# bands of the two original events

	a.set_band(smallest_enclosing_seg(a.get_band(), b.get_band()))

	# The cluster's time interval is the smallest interval containing
	# the intervals of the two original events

	a.set_period(smallest_enclosing_seg(a.get_period(), b.get_period()))

	# The amplitude, SNR, confidence, and peak time of the cluster are
	# those of the most confident of the two events (more negative
	# confidence == more confident).

	if a.confidence > b.confidence:
		a.amplitude = b.amplitude
		a.snr = b.snr
		a.confidence = b.confidence
		a.set_peak(b.get_peak())
		if hasattr(a, "tfvolume"):
			a.tfvolume = b.tfvolume


def ClusterSnglBurstTable(triggers, testfunc, clusterfunc, bailoutfunc = None):
	"""
	Cluster the triggers in the list.  testfunc should accept a pair of
	triggers, and return 0 if they should be clustered.  clusterfunc
	should accept a pair of triggers, and replace the contents of the
	first with a cluster constructed from the two.  If bailoutfunc is
	provided, the triggers will be sorted using testfunc as a
	comparison operator, and then only pairs of triggers for which
	bailoutfunc returns 0 will be considered for clustering.
	"""
	while True:
		did_cluster = False

		if bailoutfunc:
			triggers.sort(testfunc)

		i = 0
		while i < len(triggers):
			j = i + 1
			while j < len(triggers):
				if not testfunc(triggers[i], triggers[j]):
					clusterfunc(triggers[i], triggers[j])
					del triggers[j]
					did_cluster = True
				else:
					if bailoutfunc:
						if bailoutfunc(triggers[i], triggers[j]):
							break
					j += 1
			i += 1

		if not did_cluster:
			return


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
		tsim = sim.get_geocent_peak()
	elif burst.ifo == "H1":
		tsim = sim.get_h_peak()
	elif burst.ifo == "H2":
		tsim = sim.get_h_peak()
	elif burst.ifo == "L1":
		tsim = sim.get_l_peak()
	else:
		raise Exception, "unrecognized sngl_burst IFO \"%s\"" % burst.ifo
	return tsim in burst.get_period()

def CompareSimBurstAndSnglBurstByTimeandFreq(sim, burst):
	"""
	Return True if the peak time and centre frequency of sim lie within
	the time-frequency tile of burst.
	"""
	return CompareSimBurstAndSnglBurstByTime(sim, burst) and (sim.freq in burst.get_band())

