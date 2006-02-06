#!/usr/bin/python

from matplotlib.patches import Patch
import pylab

from glue import lal
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

	xvals, yvals = webplot.smooth(table.getColumnByName("central_freq").asarray(), desc.band, desc.freqwidth)
	pylab.plot(xvals, yvals / duration)

	pylab.set(axes, xlim = list(desc.band))
	pylab.xticks(pylab.arange(desc.band[0], desc.band[1], 100))
	pylab.grid(True)

	pylab.title(desc.instrument + " Excess Power Trigger Rate vs. Central Frequency\n(GPS Times %s ... %s, %d Triggers, %g Hz Average)" % (desc.segment[0], desc.segment[1], len(table), desc.freqwidth))
	pylab.xlabel("Central Frequency (Hz)")
	pylab.ylabel("Rate Density (Triggers/s per Hz)")

	pylab.savefig(desc.filename)


#
# Make a plot and send to client
#

description = Plot().parse_form()

makeplot(description, webplot.gettriggers(description))

webplot.SendImage(description)
