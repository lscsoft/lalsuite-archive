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
	def __init__(self):
		webplot.PlotDescription.__init__(self)
		self.widthratio = 3.0
		self.glitchthreshold = 1.0 * math.sqrt(self.widthratio)

	def trig_segment(self):
		return self.segment.protract(5.0 * self.widthratio * self.ratewidth)


#
# How to make a rate plot
#

def makeplot(desc, table):
	fig = pylab.figure(1)
	fig.set_figsize_inches(16,8)
	axes = pylab.gca()

	peaktime = [float(row.get_peak()) for row in table]
	confidence = -table.getColumnByName("confidence").asarray()

	pylab.plot(*webplot.smooth(peaktime, desc.trig_segment(), desc.ratewidth, weights = confidence))
	xvals, yvals = webplot.smooth(peaktime, desc.trig_segment(), desc.widthratio * desc.ratewidth, weights = confidence)
	pylab.plot(xvals, yvals, xvals, desc.glitchthreshold * yvals)

	pylab.set(axes, xlim = list(desc.segment))
	pylab.grid(True)

	for greyseg in ~desc.seglist & segments.segmentlist([desc.segment]):
		pylab.axvspan(greyseg[0], greyseg[1], facecolor = "k", alpha = 0.2)

	pylab.title(desc.instrument + " Excess Power Confidence Rate vs. Time\n(%d Triggers, %g s and %g s Averages)" % (len(table), desc.ratewidth, desc.widthratio * desc.ratewidth))
	pylab.xlabel("GPS Time (s)")
	pylab.ylabel("Confidence Rate (trigger confidence / s)")

	pylab.savefig(desc.filename)


#
# Make a plot and send to client
#

description = Plot().parse_form()

makeplot(description, webplot.gettriggers(description))

webplot.SendImage(description)
