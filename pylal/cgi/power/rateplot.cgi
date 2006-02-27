#!/usr/bin/python

import math
from matplotlib.patches import Patch
import pylab

from glue import lal
from glue import segments
from pylal import rate

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

	pylab.plot(*rate.smooth([float(row.get_peak()) for row in table], desc.trig_segment(), desc.ratewidth))

	pylab.set(axes, xlim = list(desc.segment))
	pylab.grid(True)

	for seg in ~desc.seglist & segments.segmentlist([desc.segment]):
		pylab.axvspan(seg[0], seg[1], facecolor = "k", alpha = 0.2)

	pylab.title(desc.instrument + " Excess Power Trigger Rate vs. Time\n(%d Triggers, %g s Average)" % (len(table), desc.ratewidth))
	pylab.xlabel("GPS Time (s)")
	pylab.ylabel("Trigger Rate (Hz)")

	pylab.savefig(desc.filename)


#
# Make a plot and send to client
#

description = Plot().parse_form()

makeplot(description, webplot.gettriggers(description)[0])

webplot.SendImage(description)
