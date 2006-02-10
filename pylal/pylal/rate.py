import math
import numarray
from numarray import convolve

#
# =============================================================================
#
#      Convolve impulse events with a Gaussian window to compute the rate
#
# =============================================================================
#

def GaussianWindow(halfwidth):
	"""
	Generate a normalized Gaussian window (integral = 1).
	"""
	bins_per_unit = 10.0 / halfwidth
	return numarray.exp(-numarray.arrayrange(-10.0 * halfwidth, +10.0 * halfwidth, 1.0/bins_per_unit)**2.0 / (2.0 * halfwidth**2.0)) / math.sqrt(2.0 * math.pi) / halfwidth


class Rate(object):
	"""
	An object for binning and smoothing data to compute a moving average
	rate.
	"""
	def __init__(self, segment, halfwidth):
		"""
		Initialize the bins for the given segment and width.
		"""
		self.halfwidth = halfwidth
		self.bins_per_unit = 10.0 / halfwidth
		self.start = segment[0]
		self.xvals = numarray.arrayrange(0.0, float(segment.duration()), 1.0/self.bins_per_unit) + float(segment[0])
		self.yvals = numarray.zeros(len(self.xvals), "Float32")

	def bin(self, x):
		"""
		Return the index for the bin corresponding to x.
		"""
		return int(float(x - self.start) * self.bins_per_unit)

	def __getitem__(self, x):
		"""
		Retrieve the weight in bin corresponding to x.
		"""
		return self.yvals[self.bin(x)]

	def __setitem__(self, x, weight):
		"""
		Add weight to the bin corresponding to x.
		"""
		self.yvals[self.bin(x)] += weight

	def convolve(self):
		"""
		Generate a window, and convolve it with the binned weights
		to generate rate data.
		"""
		self.yvals = convolve.convolve(self.yvals, GaussianWindow(self.halfwidth), mode=convolve.SAME)
		return self


def smooth(impulses, segment, width, weights = None):
	rate = Rate(segment, width / 2.0)
	if True:
		if weights == None:
			for x in impulses:
				if segment[0] <= x < segment[1]:
					rate[x] = 1.0
		else:
			for n, x in enumerate(impulses):
				if segment[0] <= x < segment[1]:
					rate[x] = weights[n]
	else:
		# inject pulses at regular intervals to test normalization
		for x in numarray.arrayrange(float(segment[0]), float(segment[1]), 0.5):
			rate[x] = 1.0
	rate.convolve()
	return rate.xvals, rate.yvals

