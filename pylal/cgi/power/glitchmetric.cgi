#!/usr/bin/python

import math
from matplotlib.patches import Patch
from numarray import nd_image
import pylab

from glue import lal
from glue import segments

import webplot


#
# Confidence rate ratio plot description
#

class Plot(webplot.PlotDescription):
	def __init__(self):
		webplot.PlotDescription.__init__(self)
		self.widthratio = 3.0
		self.glitchthreshold = 1.0 * math.sqrt(self.widthratio)

	def trig_segment(self):
		return self.segment.protract(5.0 * self.widthratio * self.ratewidth)


#
# How to make a confidence rate ratio plot
#

def glitchsegments(xvals, yvals, threshold):
	# find starts and ends of segments above threshold
	if yvals[0] >= threshold:
		starts = [xvals[0]]
	else:
		starts = []
	ends = []
	for i in range(1, len(yvals) - 1):
		if (yvals[i] < threshold) and (yvals[i + 1] >= threshold):
			starts.append(xvals[i])
		elif (yvals[i] < threshold) and (yvals[i - 1] >= threshold):
			ends.append(xvals[i])
	if yvals[-1] >= threshold:
		ends.append(xvals[-1])

	# found any?
	if (len(starts) < 1) or (len(ends) < 1):
		return segments.segmentlist([])

	# make sure we start with a start and end with an end
	if ends[0] < starts[0]:
		ends = ends[1:]
	if starts[-1] > ends[-1]:
		starts = starts[:-1]

	# turn start/end pairs into segments
	return segments.segmentlist(map(segments.segment, starts, ends))


def makeplot(desc, table):
	fig = pylab.figure(1)
	fig.set_figsize_inches(16,8)
	axes = pylab.gca()

	peaktime = [float(row.get_peak()) for row in table]
	confidence = -table.getColumnByName("confidence").asarray()

	xvals, yvals = webplot.smooth(peaktime, desc.trig_segment(), desc.ratewidth, weights = confidence)
	#pylab.plot(xvals, yvals)

	xvals2, yvals2 = webplot.smooth(peaktime, desc.trig_segment(), desc.widthratio * desc.ratewidth, weights = confidence)
	#pylab.plot(xvals2, yvals2, xvals2, desc.glitchthreshold * yvals2)

	yvals2 = nd_image.zoom(yvals2, desc.widthratio)
	yvals = pylab.where(yvals2 > 0.0, yvals, 0.0) / pylab.where(yvals2 > 0.0, yvals2, 1.0)
	pylab.semilogy(xvals, yvals)

	glitchsegs = glitchsegments(xvals, yvals, desc.glitchthreshold)

	pylab.axhline(desc.glitchthreshold, color = "r")

	pylab.set(axes, xlim = list(desc.segment))
	pylab.set(axes, ylim = [.1, 10])
	pylab.grid(True)

	for seg in ~desc.seglist & segments.segmentlist([desc.segment]):
		pylab.axvspan(seg[0], seg[1], facecolor = "k", alpha = 0.2)
	for seg in glitchsegs:
		pylab.axvspan(seg[0], seg[1], facecolor = "r", alpha = 0.2)

	pylab.title(desc.instrument + " Excess Power %g s Confidence Rate to %g s Confidence Rate Ratio vs. Time\n(%d Triggers)" % (desc.ratewidth, desc.widthratio * desc.ratewidth, len(table)))
	pylab.xlabel("GPS Time (s)")
	pylab.ylabel("Confidence Rate Ratio")

	pylab.savefig(desc.filename)


#
# Make a plot and send to client
#

description = Plot().parse_form()

makeplot(description, webplot.gettriggers(description))

webplot.SendImage(description)
