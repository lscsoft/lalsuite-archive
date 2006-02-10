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
