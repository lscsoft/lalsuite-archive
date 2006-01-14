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
		self.filename = None
		self.format = None	# force update
		self.set_format("png")

		return self

	def __del__(self):
		self.remove_tmpfile()

	def remove_tmpfile(self):
		if self.filename:
			os.remove(self.filename)

	def set_format(self, format):
		if format not in ["eps", "png"]:
			raise Exception, "unrecognized format %s" % format
		if self.format == format:
			return
		self.format = format
		if format == "eps":
			self.set_extension("eps")
		elif format == "png":
			self.set_extension("png")

	def set_extension(self, extension):
		self.remove_tmpfile()
		handle, self.filename = tempfile.mkstemp("." + extension.strip(), "webplot_")
		os.close(handle)

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
		self.set_format(form.getfirst("format", self.format))

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

def SnglBurstOnlyHandler(doc):
	"""
	Construct a document handler that reads only SnglBurst tables.
	"""
	return docutils.PartialLIGOLWContentHandler(doc, lambda name, attrs: (name == ligolw.Table.tagName) and (attrs["Name"] == lsctables.SnglBurstTable.tableName))

def CacheURLs(cachename, seg):
	"""
	Open a trigger cache, and return a list of URLs for files
	intersecting seg.
	"""
	return [c.url for c in map(lal.CacheEntry, file(cachename)) if c.segment.intersects(seg)]

def gettriggers(plotdesc):
	# load SnglBurst tables containing relevant triggers
	doc = ligolw.Document()
	handler = SnglBurstOnlyHandler(doc)
	for url in CacheURLs(plotdesc.cache, plotdesc.segment):
		try:
			sax.parse(urllib.urlopen(url), handler)
		except ligolw.ElementError, e:
			raise Exception, "error parsing file %s: %s" % (url, str(e))

	# merge trigger tables
	docutils.MergeCompatibleTables(doc)
	if len(doc.childNodes) == 0:
		# no trigger tables found, return an empty table
		return lsctables.New(lsctables.SnglBurstTable)
	elif len(doc.childNodes) > 1:
		# couldn't merge trigger tables
		raise Exception, "files contain incompatible SnglBurst tables"
	table = doc.childNodes[0]

	# cluster
	SnglBurstUtils.ClusterSnglBurstTable(table.rows, SnglBurstUtils.CompareSnglBurstByPeakTimeAndFreq, SnglBurstUtils.SnglBurstCluster, SnglBurstUtils.CompareSnglBurstByPeakTime)

	# remove triggers that lie outside the required segment
	table.filterRows(lambda row: row.get_peak() in plotdesc.trig_segment())

	return table


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
	if True:
		for x in impulses:
			if x in segment:
				yvals[bin(x)] += 1.0
	else:
		# inject pulses at regular intervals to test normalization
		for x in pylab.arrayrange(float(segment[0]), float(segment[1]), 0.5):
			if x in segment:
				yvals[bin(x)] += 1.0

	return (xvals + float(segment[0]), pylab.convolve(yvals, window, mode=1))


#
# How to send image to client
#

def SendImage(plotdesc):
	if plotdesc.format == "png":
		print >>sys.stdout, "Content-Type: image/png\n"
	elif plotdesc.format == "eps":
		print >>sys.stdout, "Content-Type: application/postscript\n"
	shutil.copyfileobj(file(plotdesc.filename), sys.stdout)
