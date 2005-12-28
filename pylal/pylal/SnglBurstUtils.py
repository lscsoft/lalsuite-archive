from glue import segments

def CompareSnglBurstByPeakTime(a, b):
	return cmp(a.get_peak(), b.get_peak())


def CompareSnglBurstByPeakTimeAndFreq(a, b):
	result = cmp(a.get_peak(), b.get_peak())
	if not result:
		result = cmp(a.get_band(), b.get_band())
	return result


def SnglBurstCluster(a, b):
	# The cluster's frequency band is the smallest band containing the
	# bands of the two original events

	band1 = a.get_band()
	band2 = b.get_band()
	a.set_band(segments.segment(min(band1[0], band2[0]), max(band1[1], band2[1])))

	# The cluster's time interval is the smallest interval containing
	# the intervals of the two original events

	period1 = a.get_period()
	period2 = b.get_period()
	a.set_period(segments.segment(min(period1[0], period2[0]), max(period1[1], period2[1])))

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
