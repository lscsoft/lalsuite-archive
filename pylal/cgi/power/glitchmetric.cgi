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
		return self.segment.protract(10.0 * self.ratewidth)


#
# How to make a rate plot
#

def makeplot(desc, table):
	fig = pylab.figure(1)
	fig.set_figsize_inches(16,8)
	axes = pylab.gca()
	widthratio = 5.0

	pylab.semilogy(*webplot.smooth([float(row.get_peak()) for row in table], desc.trig_segment(), desc.ratewidth, weights = -table.getColumnByName("confidence").asarray()))
	pylab.semilogy(*webplot.smooth([float(row.get_peak()) for row in table], desc.trig_segment(), widthratio * desc.ratewidth, weights = -table.getColumnByName("confidence").asarray()))

	pylab.set(axes, xlim = list(desc.segment))
	pylab.set(axes, ylim = [1, 10000])
	pylab.grid(True)

	for greyseg in ~desc.seglist & segments.segmentlist([desc.segment]):
		pylab.axvspan(greyseg[0], greyseg[1], facecolor = "k", alpha = 0.2)

	pylab.title(desc.instrument + " Excess Power Glitch Metric vs. Time\n(%d Triggers, %g s Average)" % (len(table), desc.ratewidth))
	pylab.xlabel("GPS Time (s)")
	pylab.ylabel("Glitch Metric (trigger confidence / s)")

	pylab.savefig(desc.filename)


#
# Make a plot and send to client
#

description = Plot().parse_form()

makeplot(description, webplot.gettriggers(description))

webplot.SendImage(description)
