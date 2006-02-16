from glue import segments

#
# =============================================================================
#
#                                  Clustering
#
# =============================================================================
#

def CompareSnglBurstByPeakTime(a, b):
	return cmp(a.get_peak(), b.get_peak())


def CompareSnglBurstByPeakTimeAndFreq(a, b):
	return cmp((a.get_peak(), a.get_band()), (b.get_peak(), b.get_band()))


def SnglBurstCluster(a, b):
	def smallest_enclosing_interval(seg1, seg2):
		return segments.segment(min(seg1[0], seg2[0]), max(seg1[1], seg2[1]))

	# The cluster's frequency band is the smallest band containing the
	# bands of the two original events

	a.set_band(smallest_enclosing_interval(a.get_band(), b.get_band()))

	# The cluster's time interval is the smallest interval containing
	# the intervals of the two original events

	a.set_period(smallest_enclosing_interval(a.get_period(), b.get_period()))

	# The amplitude, SNR, confidence, and peak time of the cluster are
	# those of the most confident of the two events (more negative
	# confidence == more confident).

	if a.confidence > b.confidence:
		a.amplitude = b.amplitude
		a.snr = b.snr
		a.confidence = b.confidence
		a.set_peak(b.get_peak())
		a.peak_dof = b.peak_dof


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
