#!/usr/bin/python

import math
from matplotlib.patches import Patch
import pylab

from glue import lal
from glue import segments

import webplot


#
# Cummulative rate vs. confidence plot description
#

class Plot(webplot.PlotDescription):
	pass


#
# How to make a cummulative rate vs. confidence plot
#

def makeplot(desc, table):
	duration = float(segments.segmentlist.duration(desc.seglist & segments.segmentlist([desc.segment])))
	
	fig = pylab.figure(1)
	fig.set_figsize_inches(16,8)
	axes = pylab.gca()

	confidence = pylab.sort(-table.getColumnsByName("confidence")[0].asarray())
	yvals = pylab.arrayrange(len(confidence), 0.0, -1.0) / duration

	pylab.loglog(confidence, yvals)

	pylab.grid(True)

	pylab.title(desc.instrument + " Cummulative Rate vs. Confidence\n(GPS Times %s ... %s)" % (desc.segment[0], desc.segment[1]))
	pylab.xlabel("|Confidence|")
	pylab.ylabel("Rate (Hz)")

	pylab.savefig(desc.filename)


#
# Make a plot and send to client
#

class Row(object):
	__slots__ = ["confidence"]

description = Plot().parse_form()

makeplot(description, webplot.gettriggers(description, Row))

webplot.SendOutput(description)
