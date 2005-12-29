import cgi
import cgitb ; cgitb.enable()
import math
from matplotlib.patches import Patch
import os
import pylab
import shutil
import sys
import tempfile
import urllib
from xml import sax

from glue import lal
from glue import segments
from glue import segmentsUtils
from glue.ligolw import ligolw
from glue.ligolw import lsctables
from glue.ligolw import docutils

from pylal import SnglBurstUtils

import eventdisplay


#
# CGI display request
#

class PlotDescription(object):
	def __init__(self):
		# set defaults
		now = eventdisplay.runtconvert(eventdisplay.TconvertCommand("now"))
		self.segment = segments.segment(now, now + (-1 * 3600))
		self.ratewidth = 60.0
		self.freqwidth = 16.0
		self.band = segments.segment(0.0, 2500.0)
		self.set_instrument("H1")
		self.seglist = segments.segmentlist([self.segment])
		handle, self.filename = tempfile.mkstemp(".png", "tfplot_")
		os.close(handle)

		return self

	def __del__(self):
		# remove temporary file
		os.remove(self.filename)

	def set_instrument(self, instrument):
		# update cache name along with instrument
		self.instrument = str(instrument)
		self.cache = "/home/kipp/cgi-bin/" + self.instrument + "/filelist.cache"

	def parse_form(self):
		# parse CGI form
		form = cgi.FieldStorage()

		start = lal.LIGOTimeGPS(form.getfirst("start", str(self.segment[0])))
		duration = lal.LIGOTimeGPS(form.getfirst("dur", str(self.segment.duration())))

		self.segment = segments.segment(start, start + duration)
		self.ratewidth = float(form.getfirst("ratewidth", str(self.ratewidth)))
		self.freqwidth = float(form.getfirst("freqwidth", str(self.freqwidth)))
		self.band = segments.segment(float(form.getfirst("lofreq", str(self.band[0]))), float(form.getfirst("hifreq", str(self.band[1]))))
		self.set_instrument(form.getfirst("inst", self.instrument))
		self.load_seglist()

		return self

	def load_seglist(self):
		# get a segmentlist from the cache file
		self.seglist = segmentsUtils.fromlalcache(file(self.cache), coltype=lal.LIGOTimeGPS).coalesce()

	def trig_segment(self):
		# interval in which triggers must be read in order to
		# produce a plot
		return self.segment


#
# How to get a table of triggers within a segment
#

def CacheURLs(cachename, seg):
	"""
	Return a list of CacheEntrys for files containing triggers of
	interest.
	"""
	return [c.url for c in map(lal.CacheEntry, file(cachename)) if c.segment.intersects(seg)]

def gettriggers(plotdesc):
	# load documents containing relevant triggers
	doc = ligolw.Document()
	handler = lsctables.LIGOLWContentHandler(doc)
	for url in CacheURLs(plotdesc.cache, plotdesc.segment):
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
		table.parentNode.removeChild(table)
	del doc

	# merge trigger tables
	while len(tables) > 1:
		docutils.MergeElements(tables[0], tables[1])
		del tables[1]

	# cluster
	SnglBurstUtils.ClusterSnglBurstTable(tables[0].rows, SnglBurstUtils.CompareSnglBurstByPeakTimeAndFreq, SnglBurstUtils.SnglBurstCluster, SnglBurstUtils.CompareSnglBurstByPeakTime)

	# remove triggers that lie outside the requested segment
	tables[0].filterRows(lambda row: row.get_peak() in plotdesc.segment)

	return tables[0]


#
# Convolve impulse events with a Gaussian window to produce a smooth
# function.
#

def smooth(impulses, segment, width):
	halfwidth = width / 2.0
	bins_per_unit = 10.0 / halfwidth

	window = pylab.exp(-pylab.arrayrange(-10.0 * halfwidth, +10.0 * halfwidth, 1.0/bins_per_unit)**2.0 / (2.0 * halfwidth**2.0)) / math.sqrt(2.0 * math.pi) / halfwidth

	xvals = pylab.arrayrange(0.0, float(segment.duration()) + 1.0/bins_per_unit, 1.0/bins_per_unit)

	def bin(x):
		return int(float(x - segment[0]) * bins_per_unit)

	yvals = pylab.zeros(len(xvals), "Float32")
	#for x in pylab.arrayrange(float(segment[0]), float(segment[1]), 0.5):
	#	if x in segment:
	#		yvals[bin(x)] += 1.0
	for x in impulses:
		if x in segment:
			yvals[bin(x)] += 1.0

	return (xvals + float(segment[0]), pylab.convolve(yvals, window, mode=1))


#
# How to send image to client
#

def SendPNG(plotdesc):
	print >>sys.stdout, "Content-Type: image/png\n"
	shutil.copyfileobj(file(plotdesc.filename), sys.stdout)
