#!/usr/bin/python

import math
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

	bandwidth = table.getColumnsByName("bandwidth")[0].asarray()
	lo_freq = table.getColumnsByName("central_freq")[0].asarray() - 0.5 * bandwidth
	start_time = pylab.asarray([float(row.get_start()) for row in table.rows])

	if len(table.rows):
		viz.tfplot(start_time, table.getColumnsByName("duration")[0].asarray(), lo_freq, bandwidth, pylab.log(-table.getColumnsByName("confidence")[0].asarray()), alpha=0.3)

	pylab.set(axes, xlim = list(desc.segment), ylim = list(desc.band))

	for greyseg in ~desc.seglist & segments.segmentlist([desc.segment]):
		pylab.axvspan(greyseg[0], greyseg[1], facecolor = "k", alpha = 0.2)
	pylab.title(desc.instrument + " Excess Power Time-Frequency Plane")
	pylab.xlabel("GPS Time (s)")
	pylab.ylabel("Frequency (Hz)")

	pylab.savefig(desc.filename)


#
# Make a plot and send to client
#

class Row(object):
	__slots__ = ["start_time", "start_time_ns", "duration", "confidence", "central_freq", "bandwidth"]
	def get_start(self):
		return lal.LIGOTimeGPS(self.start_time, self.start_time_ns)

description = Plot().parse_form()

makeplot(description, webplot.gettriggers(description, Row))

webplot.SendOutput(description)
