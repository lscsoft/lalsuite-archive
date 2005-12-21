#!/usr/bin/python

import math
from matplotlib.patches import Patch
import pylab

from glue import lal
from glue import segments

import webplot


#
# Rate plot description
#

class RatePlot(webplot.PlotDescription):
	def trig_segment(self):
		return self.segment.protract(5.0 * self.ratewidth)


#
# How to make a rate plot
#

def makeplot(desc, table):
	halfwidth = desc.ratewidth / 2.0
	bins_per_second = 10.0 / halfwidth
	fig = pylab.figure(1)
	fig.set_figsize_inches(16,8)
	axes = pylab.gca()

	def bin(x):
		return int(float(x - desc.trig_segment()[0]) * bins_per_second)

	window = pylab.exp(-pylab.arrayrange(-10.0 * halfwidth, +10.0 * halfwidth, 1.0/bins_per_second)**2.0 / (2.0 * halfwidth**2.0)) / math.sqrt(2.0 * math.pi) / halfwidth

	xvals = pylab.arrayrange(0.0, float(desc.trig_segment().duration()) + 1.0/bins_per_second, 1.0/bins_per_second)

	peaktimes = pylab.zeros(len(xvals), "Float32")
	#for t in pylab.arrayrange(float(desc.trig_segment()[0]), float(desc.trig_segment()[1]), 0.5):
	#	if desc.trig_segment()[0] <= t <= desc.trig_segment()[1]:
	#		peaktimes[bin(t)] += 1.0
	for t in [float(row.get_peak()) for row in table.rows]:
		if desc.trig_segment()[0] <= t <= desc.trig_segment()[1]:
			peaktimes[bin(t)] += 1.0

	pylab.plot(xvals + float(desc.trig_segment()[0]), pylab.convolve(peaktimes, window, mode=1))

	pylab.set(axes, xlim = list(desc.segment))

	for greyseg in ~desc.seglist & segments.segmentlist([desc.segment]):
		pylab.axvspan(greyseg[0], greyseg[1], facecolor = "k", alpha = 0.2)

	pylab.title(desc.instrument + " Excess Power Trigger Rate\n(%g s Average)" % (halfwidth * 2.0))
	pylab.xlabel("GPS Time (s)")
	pylab.ylabel("Trigger Rate (Hz)")

	pylab.savefig(desc.filename)
	pylab.clf()


#
# Make a plot and send to client
#

class Row(object):
	__slots__ = ["peak_time", "peak_time_ns"]
	def get_peak(self):
		return lal.LIGOTimeGPS(self.peak_time, self.peak_time_ns)

description = RatePlot().parse_form()

makeplot(description, webplot.gettriggers(description, Row))

webplot.SendOutput(description)
