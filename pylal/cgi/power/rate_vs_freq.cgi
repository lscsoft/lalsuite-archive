#!/usr/bin/python

from matplotlib.patches import Patch
import pylab

from glue import segments

import webplot


#
# Trigger rate vs. frequency plot description
#

class Plot(webplot.PlotDescription):
	pass


#
# How to make a trigger rate vs. frequency plot
#

def makeplot(desc, table):
	duration = float(segments.segmentlist.duration(desc.seglist & segments.segmentlist([desc.segment])))
	
	fig = pylab.figure(1)
	fig.set_figsize_inches(16,8)
	axes = pylab.gca()

	pylab.plot(*webplot.smooth(table.getColumnsByName("central_freq")[0].asarray(), desc.band, desc.freqwidth))

	pylab.set(axes, xlim = list(desc.band))
	pylab.grid(True)

	pylab.title(desc.instrument + " Excess Power Trigger Rate vs. Central Frequency\n(GPS Times %s ... %s, %d Triggers, %g Hz Average)" % (desc.segment[0], desc.segment[1], len(table.rows), desc.freqwidth))
	pylab.xlabel("Central Frequency (Hz)")
	pylab.ylabel("Rate (Triggers per Hz)")

	pylab.savefig(desc.filename)


#
# Make a plot and send to client
#

class Row(object):
	__slots__ = ["central_freq"]

description = Plot().parse_form()

makeplot(description, webplot.gettriggers(description, Row))

webplot.SendPNG(description)
