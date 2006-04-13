#!/usr/bin/python

import pylab

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

	confidence = -table.getColumnByName("confidence").asarray()
	central_freq = table.getColumnByName("central_freq").asarray()

	pylab.semilogy(central_freq, confidence, "b+")

	pylab.set(axes, xlim = list(desc.band))
	pylab.xticks(pylab.arange(desc.band[0], desc.band[1], 100))
	pylab.grid(True)

	pylab.title(desc.instrument + " Excess Power Trigger Confidence vs. Central Frequency\n(GPS Times %s ... %s, %d Triggers)" % (desc.segment[0], desc.segment[1], len(table)))
	pylab.xlabel("Central Frequency (Hz)")
	pylab.ylabel("|Confidence|")

	pylab.savefig(desc.filename)


#
# Make a plot and send to client
#

description = Plot().parse_form()

makeplot(description, webplot.gettriggers(description)[0])

webplot.SendImage(description)
