#!/usr/bin/python

from matplotlib.patches import Patch
import pylab

from glue import lal
from glue import segments

import webplot


#
# Confidence vs. Time plot description
#

class Plot(webplot.PlotDescription):
	pass


#
# How to make a confidence vs. time plot
#

def makeplot(desc, table):
	fig = pylab.figure(1)
	fig.set_figsize_inches(16,8)
	axes = pylab.gca()

	confidence = -table.getColumnByName("confidence").asarray()
	peak_time = pylab.asarray([float(row.get_peak()) for row in table])

	pylab.semilogy(peak_time, confidence, "b+")

	for seg in ~desc.seglist & segments.segmentlist([desc.segment]):
		pylab.axvspan(seg[0], seg[1], facecolor = "k", alpha = 0.2)

	pylab.set(axes, xlim = list(desc.segment))
	pylab.grid(True)

	pylab.title(desc.instrument + " Excess Power Trigger Confidence vs. Time\n(%d Triggers)" % (len(table)))
	pylab.xlabel("GPS Time (s)")
	pylab.ylabel("|Confidence|")

	pylab.savefig(desc.filename)


#
# Make a plot and send to client
#

description = Plot().parse_form()

makeplot(description, webplot.gettriggers(description))

webplot.SendImage(description)
