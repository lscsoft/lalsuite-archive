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

class Plot(webplot.PlotDescription):
	def trig_segment(self):
		return self.segment.protract(5.0 * self.ratewidth)


#
# How to make a rate plot
#

def makeplot(desc, table):
	fig = pylab.figure(1)
	fig.set_figsize_inches(16,8)
	axes = pylab.gca()

	pylab.plot(*webplot.smooth([float(row.get_peak()) for row in table.rows], desc.trig_segment(), desc.ratewidth))

	pylab.set(axes, xlim = list(desc.segment))
	pylab.grid(True)

	for greyseg in ~desc.seglist & segments.segmentlist([desc.segment]):
		pylab.axvspan(greyseg[0], greyseg[1], facecolor = "k", alpha = 0.2)

	pylab.title(desc.instrument + " Excess Power Trigger Rate vs. Time\n(%g s Average)" % (desc.ratewidth))
	pylab.xlabel("GPS Time (s)")
	pylab.ylabel("Trigger Rate (Hz)")

	pylab.savefig(desc.filename)


#
# Make a plot and send to client
#

class Row(object):
	__slots__ = ["peak_time", "peak_time_ns"]
	def get_peak(self):
		return lal.LIGOTimeGPS(self.peak_time, self.peak_time_ns)

description = Plot().parse_form()

makeplot(description, webplot.gettriggers(description, Row))

webplot.SendPNG(description)
