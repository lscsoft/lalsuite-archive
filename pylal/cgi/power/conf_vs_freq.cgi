#!/usr/bin/python

import math
from matplotlib.patches import Patch
import pylab

from glue import lal
from glue import segments

import webplot


#
# Confidence vs. Frequency plot description
#

class Plot(webplot.PlotDescription):
	pass


#
# How to make a confidence vs. frequency plot
#

def makeplot(desc, table):
	fig = pylab.figure(1)
	fig.set_figsize_inches(16,8)
	axes = pylab.gca()

	confidence = -table.getColumnsByName("confidence")[0].asarray()
	central_freq = table.getColumnsByName("central_freq")[0].asarray()

	pylab.semilogy(central_freq, confidence, "b+")

	pylab.set(axes, xlim = list(desc.band))
	pylab.grid(True)

	pylab.title(desc.instrument + " Excess Power Confidence vs. Central Frequency")
	pylab.xlabel("Central Frequency (Hz)")
	pylab.ylabel("|Confidence|")

	pylab.savefig(desc.filename)


#
# Make a plot and send to client
#

class Row(object):
	__slots__ = ["central_freq", "confidence"]
	def get_peak(self):
		return lal.LIGOTimeGPS(self.peak_time, self.peak_time_ns)

description = Plot().parse_form()

makeplot(description, webplot.gettriggers(description, Row))

webplot.SendOutput(description)
