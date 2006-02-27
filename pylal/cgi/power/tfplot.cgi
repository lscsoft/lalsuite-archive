#!/usr/bin/python

from matplotlib.patches import Patch
import pylab

from glue import lal
from glue import segments
from pylal import viz

import webplot


#
# TF plot description
#

class Plot(webplot.PlotDescription):
	def trig_segment(self):
		return self.segment.protract(2)


#
# How to make a time-frequency plane plot
#


def makeplot(desc, table):
	fig = pylab.figure(1)
	fig.set_figsize_inches(16,8)
	axes = pylab.gca()

	bandwidth = table.getColumnByName("bandwidth").asarray()
	lo_freq = table.getColumnByName("central_freq").asarray() - 0.5 * bandwidth
	start_time = pylab.asarray([float(row.get_start()) for row in table])

	if len(table):
		viz.tfplot(start_time, table.getColumnByName("duration").asarray(), lo_freq, bandwidth, pylab.log(-table.getColumnByName("confidence").asarray()), alpha=0.3)

	pylab.set(axes, xlim = list(desc.segment), ylim = list(desc.band))

	for seg in ~desc.seglist & segments.segmentlist([desc.segment]):
		pylab.axvspan(seg[0], seg[1], facecolor = "k", alpha = 0.2)
	pylab.title(desc.instrument + " Excess Power Time-Frequency Plane\n(%d Triggers)" % (len(table)))
	pylab.xlabel("GPS Time (s)")
	pylab.ylabel("Frequency (Hz)")

	pylab.savefig(desc.filename)


#
# Make a plot and send to client
#

description = Plot().parse_form()

makeplot(description, webplot.gettriggers(description)[0])

webplot.SendImage(description)
