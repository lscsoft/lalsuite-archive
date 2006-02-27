#!/usr/bin/python

from matplotlib.patches import Patch
import pylab

from glue import lal
from glue import segments

import webplot


#
# Cummulative rate vs. snr plot description
#

class Plot(webplot.PlotDescription):
	pass


#
# How to make a cummulative rate vs. snr plot
#

def makeplot(desc, table):
	duration = float(segments.segmentlist.duration(desc.seglist & segments.segmentlist([desc.segment])))
	
	fig = pylab.figure(1)
	fig.set_figsize_inches(16,8)
	axes = pylab.gca()

	snr = pylab.sort(table.getColumnByName("snr").asarray())
	yvals = pylab.arrayrange(len(snr), 0.0, -1.0) / duration

	pylab.loglog(snr, yvals)

	pylab.grid(True)

	pylab.title(desc.instrument + " Excess Power Cummulative Trigger Rate vs. SNR\n(GPS Times %s ... %s, %d Triggers)" % (desc.segment[0], desc.segment[1], len(table)))
	pylab.xlabel("SNR")
	pylab.ylabel("Rate (Hz)")

	pylab.savefig(desc.filename)


#
# Make a plot and send to client
#

description = Plot().parse_form()

makeplot(description, webplot.gettriggers(description)[0])

webplot.SendImage(description)
