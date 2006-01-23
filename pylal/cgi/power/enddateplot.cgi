#!/usr/bin/python

import math
from matplotlib.patches import Patch
import os
import pylab
import shutil
import sys
import tempfile


from glue import lal
from glue import segments

import eventdisplay


#
# Rate plot description
#

handle, filename = tempfile.mkstemp("." + "png", "webplot_")
os.close(handle)


now = lal.LIGOTimeGPS(eventdisplay.runtconvert(eventdisplay.TconvertCommand("now")))
s5length = 1.0 * 365.25 * 24.0 * 60.0 * 60.0	# 1 year


trigsegs = eventdisplay.TrigSegs()
seglist = trigsegs.H1 & trigsegs.H2 & trigsegs.L1
del trigsegs


def end_date(t):
	l = seglist - segments.segmentlist([segments.segment(t, segments.infinity())])
	livetime = float(l.duration())
	if livetime > 0.0:
		return now + (s5length/livetime - 1) * float(l[-1][1] - l[0][0])
	else:
		return now


#
# How to make a plot of projected end date vs. time.
#

def makeplot():
	fig = pylab.figure(1)
	fig.set_figsize_inches(16,8)
	axes = pylab.gca()
	increment = 12 * 3600.0

	xvals = pylab.arange(float(seglist[0][0]) + increment, float(now), 12 * 3600.0)

	pylab.plot(map(float, xvals), map(float, map(end_date, xvals)))

	pylab.set(axes, xlim = [seglist[0][0], now])
	pylab.grid(True)

	pylab.title("Projected S5 End Time vs. Time")
	pylab.xlabel("GPS Time (s)")
	pylab.ylabel("Projected GPS End Time (s)")

	pylab.savefig(filename)


#
# Make a plot and send to client
#

makeplot()

#webplot.SendImage(filename)
print >>sys.stdout, "Content-Type: image/png\n"
shutil.copyfileobj(file(filename), sys.stdout)
os.remove(filename)
