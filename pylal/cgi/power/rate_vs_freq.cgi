#!/usr/bin/python

import math
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

	central_freq = table.getColumnsByName("central_freq")[0].asarray()

	pylab.hist(central_freq, 64)

	pylab.set(axes, xlim = list(desc.band))
	pylab.grid(True)

	pylab.title(desc.instrument + " Excess Power Trigger Rate vs. Central Frequency\n(GPS Times %s ... %s)" % (desc.segment[0], desc.segment[1]))
	pylab.xlabel("Central Frequency (Hz)")
	pylab.ylabel("Count")

	pylab.savefig(desc.filename)


#
# Make a plot and send to client
#

class Row(object):
	__slots__ = ["central_freq"]

description = Plot().parse_form()

makeplot(description, webplot.gettriggers(description, Row))

webplot.SendOutput(description)
