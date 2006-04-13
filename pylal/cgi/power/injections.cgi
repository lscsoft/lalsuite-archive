#!/usr/bin/python

import pylab

from glue import segments

import webplot


#
# Injection plot description
#

class Plot(webplot.PlotDescription):
	def trig_segment(self):
		return self.segment.protract(2)


#
# How to make a time-frequency plane plot of injections
#


def makeplot(desc, table):
	fig = pylab.figure(1)
	fig.set_figsize_inches(16,8)
	axes = pylab.gca()

	if desc.instrument[:2] in ["H1", "H2"]:
		time = pylab.asarray([float(row.get_h_peak()) for row in table.rows])
	elif desc.instrument[:2] in ["L1"]:
		time = pylab.asarray([float(row.get_l_peak()) for row in table.rows])

	freq = table.getColumnByName("freq").asarray()

	pylab.plot(time, freq, "b+")

	for seg in ~desc.seglist & segments.segmentlist([desc.segment]):
		pylab.axvspan(seg[0], seg[1], facecolor = "k", alpha = 0.2)

	pylab.set(axes, xlim = list(desc.segment), ylim = list(desc.band))
	pylab.grid(True)

	pylab.title(desc.instrument + " Injection Locations\n(%d Injections)" % (len(table)))
	pylab.xlabel("GPS Time (s)")
	pylab.ylabel("Frequency (Hz)")

	pylab.savefig(desc.filename)


#
# Make a plot and send to client
#

description = Plot().parse_form()

makeplot(description, webplot.gettriggers(description)[1])

webplot.SendImage(description)
