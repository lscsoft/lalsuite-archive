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
from glue.ligolw import metaio
from glue.ligolw import lsctables
from glue.ligolw import docutils

from pylal import SnglBurstUtils

import eventdisplay

#
# =============================================================================
#
#                             CGI display request
#
# =============================================================================
#

class PlotDescription(object):
	def __init__(self):
		# set defaults
		now = lal.LIGOTimeGPS(eventdisplay.runtconvert(eventdisplay.TconvertCommand("now")))
		self.segment = segments.segment(now, now + (-1 * 3600))
		self.ratewidth = 60.0
		self.freqwidth = 16.0
		self.band = segments.segment(0.0, 2500.0)
		self.set_instrument("H1")
		self.seglist = segments.segmentlist([self.segment])
		self.filename = None
		self.format = None	# force update
		self.set_format("png")
		self.cluster = 0

		return self

	def __del__(self):
		self.remove_tmpfile()

	def remove_tmpfile(self):
		if self.filename:
			os.remove(self.filename)

	def set_format(self, format):
		if format not in ["eps", "png", "xml"]:
			raise Exception, "unrecognized format %s" % format
		if self.format == format:
			return
		self.format = format
		if format == "eps":
			self.set_extension("eps")
		elif format == "png":
			self.set_extension("png")
		elif format == "xml":
			self.set_extension("xml")

	def set_extension(self, extension):
		self.remove_tmpfile()
		handle, self.filename = tempfile.mkstemp("." + extension.strip(), "webplot_")
		os.close(handle)

	def set_instrument(self, instrument):
		# update cache name along with instrument
		self.instrument = str(instrument)

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
		self.set_format(form.getfirst("format", self.format))
		self.cluster = int(form.getfirst("cluster", "0"))

		return self

	def trig_segment(self):
		# interval in which triggers must be read in order to
		# produce a plot
		return self.segment


#
# =============================================================================
#
#               How to get a table of triggers within a segment
#
# =============================================================================
#

def SnglBurstAndSearchSummOnlyHandler(doc):
	"""
	Construct a document handler that reads only sngl_burst and search
	summary tables.
	"""
	return docutils.PartialLIGOLWContentHandler(doc, lambda name, attrs: (name == ligolw.Table.tagName) and (metaio.StripTableName(attrs["Name"]) in map(metaio.StripTableName, [lsctables.SnglBurstTable.tableName, lsctables.SearchSummaryTable.tableName])))


def CacheURLs(cachename, seg):
	"""
	Open a trigger cache, and return a list of URLs for files
	intersecting seg.
	"""
	return [c.url for c in map(lal.CacheEntry, file(cachename)) if c.segment.intersects(seg)]

def gettriggers(plotdesc):
	# load SnglBurst tables containing relevant triggers
	doc = ligolw.Document()
	handler = SnglBurstAndSearchSummOnlyHandler(doc)
	for url in CacheURLs(eventdisplay.cache[plotdesc.instrument], plotdesc.segment):
		try:
			ligolw.make_parser(handler).parse(urllib.urlopen(url))
		except ligolw.ElementError, e:
			raise Exception, "error parsing file %s: %s" % (url, str(e))
	docutils.MergeCompatibleTables(doc)

	# get segment list from search summary table
	tables = lsctables.getTablesByType(doc, lsctables.SearchSummaryTable)
	if len(tables) == 0:
		# no search summary tables found, set an empty segment list
		plotdesc.seglist = segments.segmentlist([])
	elif len(tables) == 1:
		plotdesc.seglist = tables[0].get_inlist().coalesce()
	else:
		raise Exception, "files contain incompatible search summary tables"

	# get triggers from single_burst table
	tables = lsctables.getTablesByType(doc, lsctables.SnglBurstTable)
	if len(tables) == 0:
		# no trigger tables found, return an empty table
		return lsctables.New(lsctables.SnglBurstTable)
	elif len(tables) > 1:
		# couldn't merge trigger tables
		raise Exception, "files contain incompatible SnglBurst tables"

	# cluster
	if plotdesc.cluster:
		SnglBurstUtils.ClusterSnglBurstTable(tables[0].rows, SnglBurstUtils.CompareSnglBurstByPeakTimeAndFreq, SnglBurstUtils.SnglBurstCluster, SnglBurstUtils.CompareSnglBurstByPeakTime)

	# remove triggers that lie outside the required segment
	tables[0].filterRows(lambda row: row.get_peak() in plotdesc.trig_segment())

	return tables[0]


#
# =============================================================================
#
#     Convolve impulse events with a Gaussian window to compute the rate.
#
# =============================================================================
#

def Window(halfwidth):
	"""
	Generate a normalized Gaussian window (integral = 1).
	"""
	bins_per_unit = 10.0 / halfwidth
	return pylab.exp(-pylab.arrayrange(-10.0 * halfwidth, +10.0 * halfwidth, 1.0/bins_per_unit)**2.0 / (2.0 * halfwidth**2.0)) / math.sqrt(2.0 * math.pi) / halfwidth


class Rate(object):
	"""
	An object for binning and smoothing data to compute a moving average
	rate.
	"""
	def __init__(self, segment, halfwidth):
		"""
		Initialize the bins for the given segment and width.
		"""
		self.halfwidth = halfwidth
		self.bins_per_unit = 10.0 / halfwidth
		self.start = segment[0]
		self.xvals = pylab.arrayrange(0.0, float(segment.duration()), 1.0/self.bins_per_unit) + float(segment[0])
		self.yvals = pylab.zeros(len(self.xvals), "Float32")

	def bin(self, x):
		"""
		Return the index for the bin corresponding to x.
		"""
		return int(float(x - self.start) * self.bins_per_unit)

	def __getitem__(self, x):
		"""
		Retrieve the weight in bin corresponding to x.
		"""
		return self.yvals[self.bin(x)]

	def __setitem__(self, x, weight):
		"""
		Add weight to the bin corresponding to x.
		"""
		self.yvals[self.bin(x)] += weight

	def convolve(self):
		"""
		Generate a window, and convolve it with the binned weights
		to generate rate data.
		"""
		self.yvals = pylab.convolve(self.yvals, Window(self.halfwidth), mode=1)
		return self


def smooth(impulses, segment, width, weights = None):
	rate = Rate(segment, width / 2.0)
	if True:
		if weights == None:
			for x in impulses:
				if segment[0] <= x < segment[1]:
					rate[x] = 1.0
		else:
			for n, x in enumerate(impulses):
				if segment[0] <= x < segment[1]:
					rate[x] = weights[n]
	else:
		# inject pulses at regular intervals to test normalization
		for x in pylab.arrayrange(float(segment[0]), float(segment[1]), 0.5):
			rate[x] = 1.0
	rate.convolve()
	return rate.xvals, rate.yvals


#
# =============================================================================
#
#                                    Output
#
# =============================================================================
#

def SendImage(plotdesc):
	if plotdesc.format == "png":
		print >>sys.stdout, "Content-Type: image/png\n"
	elif plotdesc.format == "eps":
		print >>sys.stdout, "Content-Type: application/postscript\n"
	elif plotdesc.format == "xml":
		print >>sys.stdout, "Content-Type: text/xml\n"
	shutil.copyfileobj(file(plotdesc.filename), sys.stdout)
