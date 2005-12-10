#!/usr/bin/python

import os
import popen2
import math
from matplotlib.patches import Patch
import pylab
import tempfile
import urllib
from xml import sax

from glue import lal
from glue import segments
from glue import segmentsUtils
from glue.ligolw import ligolw
from glue.ligolw import lsctables
from glue.ligolw import docutils
from pylal import viz

#
# How to run tconvert
#

class TconvertCommand(object):
	def __init__(self, tspec = None):
		self._exec = "/home/kipp/local/bin/lalapps_tconvert"
		self.tspec = str(tspec)

	def __str__(self):
		s = self._exec
		if self.tspec:
			s += " " + self.tspec
		return s

def runtconvert(command):
	if type(command) != TconvertCommand:
		raise ValueError, "invalid argument to runtconvert(command): command must type TconvertCommand"
	child = popen2.Popen3(str(command), True)
	for line in child.childerr:
		pass
	for line in child.fromchild:
		result = lal.LIGOTimeGPS(line)
	status = child.wait()
	if not os.WIFEXITED(status) or os.WEXITSTATUS(status):
		raise Exception, "failure running \"" + str(command) + "\""
	return result


#
# How to get a table of triggers within a segment
#

class Row(object):
	__slots__ = ["start_time", "start_time_ns", "duration", "confidence", "central_freq", "bandwidth", "peak_time", "peak_time_ns"]
	def get_start(self):
		return lal.LIGOTimeGPS(self.start_time, self.start_time_ns)
	def get_peak(self):
		return lal.LIGOTimeGPS(self.peak_time, self.peak_time_ns)

lsctables.SnglBurstTable.RowType = Row

def CacheURLs(fname, seg):
	"""
	Return a list of CacheEntrys for files containing triggers of
	interest.
	"""
	urls = []
	for line in file(fname):
		c = lal.CacheEntry(line)
		if c.segment.intersects(seg):
			urls.append(c.url)
	return urls

def gettriggers(cachefile, segment):
	# load documents containing relevant triggers
	doc = ligolw.Document()
	handler = lsctables.LIGOLWContentHandler(doc)
	for url in CacheURLs(cachefile, segment):
		sax.parse(urllib.urlopen(url), handler)

	# if no files contain relevant triggers, return an empty table
	if not len(doc.childNodes):
		return lsctables.New(lsctables.SnglBurstTable)

	# extract trigger tables
	tables = []
	for e in doc.childNodes:
		tables += e.getChildrenByAttributes({"Name": lsctables.SnglBurstTable.tableName})

	# unlink tables to allow garbage collection
	for table in tables:
		table.parent.removeChild(table)
	del doc

	# merge trigger tables
	while len(tables) > 1:
		docutils.MergeTables(tables[0], tables[1])
		del tables[1]

	# FIXME: need to cluster!

	return tables[0]


#
# How to make a time-frequency plane plot
#

class TFPlotDescription(object):
	def __init__(self, instrument, segment, band, seglist):
		self.instrument = instrument
		self.segment = segment
		self.band = band
		self.seglist = seglist
		handle, self.filename = tempfile.mkstemp(".png", "tfplot_")
		os.close(handle)

	def trig_segment(self):
		return self.segment.protract(2)


def maketfplot(desc, table):
	fig = pylab.figure(1)
	fig.set_figsize_inches(8,8)
	axes = pylab.gca()

	bandwidth = table.getColumnsByName("bandwidth")[0].asarray()
	lo_freq = table.getColumnsByName("central_freq")[0].asarray() - 0.5 * bandwidth
	start_time = pylab.asarray([float(row.get_start()) for row in table.rows])

	if len(table.rows):
		viz.tfplot(start_time, table.getColumnsByName("duration")[0].asarray(), lo_freq, bandwidth, pylab.log(-table.getColumnsByName("confidence")[0].asarray()), alpha=0.3)

	pylab.set(axes, xlim = list(desc.segment), ylim = list(desc.band))

	for greyseg in ~desc.seglist & segments.segmentlist([desc.segment]):
		pylab.axvspan(greyseg[0], greyseg[1], facecolor = "k", alpha = 0.2)
	pylab.title(desc.instrument + " Excess Power Time-Frequency Plane")
	pylab.xlabel("GPS Time (s)")
	pylab.ylabel("Frequency (Hz)")

	pylab.savefig(desc.filename)
	pylab.clf()


#
# How to make a rate plot
#

class RatePlotDescription(object):
	def __init__(self, instrument, width, segment, seglist):
		self.instrument = instrument
		self.halfwidth = width / 2.0
		self.segment = segment
		self.seglist = seglist
		handle, self.filename = tempfile.mkstemp(".png", "rateplot_")
		os.close(handle)

	def trig_segment(self):
		return self.segment.protract(10 * self.halfwidth)


def makerateplot(desc, table):
	bins_per_second = 10.0 / desc.halfwidth
	fig = pylab.figure(1)
	fig.set_figsize_inches(8,8)
	axes = pylab.gca()

	def bin(x):
		return int(float(x - desc.trig_segment()[0]) * bins_per_second)

	window = pylab.exp(-pylab.arrayrange(-10.0 * desc.halfwidth, +10.0 * desc.halfwidth, 1.0/bins_per_second)**2.0 / (2.0 * desc.halfwidth**2.0)) / math.sqrt(2.0 * math.pi) / desc.halfwidth

	xvals = pylab.arrayrange(0.0, float(desc.trig_segment().duration()) + 1.0/bins_per_second, 1.0/bins_per_second)

	peaktimes = pylab.zeros(len(xvals), "Float32")
	#for t in pylab.arrayrange(float(desc.trig_segment()[0]), float(desc.trig_segment()[1]), 0.5):
	#	if desc.trig_segment()[0] <= t <= desc.trig_segment()[1]:
	#		peaktimes[bin(t)] += 1.0
	for t in [float(row.get_peak()) for row in table.rows]:
		if desc.trig_segment()[0] <= t <= desc.trig_segment()[1]:
			peaktimes[bin(t)] += 1.0

	pylab.plot(xvals + float(desc.trig_segment()[0]), pylab.convolve(peaktimes, window, mode=1))

	pylab.set(axes, xlim = list(desc.segment))

	for greyseg in ~desc.seglist & segments.segmentlist([desc.segment]):
		pylab.axvspan(greyseg[0], greyseg[1], facecolor = "k", alpha = 0.2)

	pylab.title(desc.instrument + " Excess Power Trigger Rate\n(%g s Average)" % (desc.halfwidth * 2.0))
	pylab.xlabel("GPS Time (s)")
	pylab.ylabel("Trigger Rate (Hz)")

	pylab.savefig(desc.filename)
	pylab.clf()
