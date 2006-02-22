from glue import segments

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
	covers an interval below b's, and 0 if the two intervals overlap.
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
	return cmp(a.get_peak(), b.get_peak())


def CompareSnglBurstByPeakTimeAndFreq(a, b):
	result = cmp(a.get_peak(), b.get_peak())
	if not result:
		result = cmp_segs(a.get_band(), b.get_band())
	return result


def CompareSnglBurst(a, b):
	result = cmp_segs(a.get_period(), b.get_period())
	if not result:
		result = cmp_segs(a.get_band(), b.get_band())
	return result


def SnglBurstCluster(a, b):
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
		a.tfvolume = b.tfvolume


def ClusterSnglBurstTable(triggers, testfunc, clusterfunc, bailoutfunc = None):
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
	return CompareSimBurstAndSnglBurstByTime(sim, burst) and (sim.freq in burst.get_band())
