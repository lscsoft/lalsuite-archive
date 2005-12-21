import cgi
import cgitb ; cgitb.enable()
import os
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

from eventdisplay import *


#
# CGI display request
#

class PlotDescription(object):
	def __init__(self):
		# set defaults
		now = runtconvert(TconvertCommand("now"))
		self.segment = segments.segment(now, now + (-1 * 3600))
		self.ratewidth = 60.0
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
		self.band = segments.segment(float(form.getfirst("lofreq", str(self.band[0]))), float(form.getfirst("hifreq", str(self.band[1]))))
		self.set_instrument(form.getfirst("inst", self.instrument))
		self.load_seglist()

		return self

	def load_seglist(self):
		# get a segmentlist from the cache file
		self.seglist = segmentsUtils.fromlalcache(file(self.cache), coltype=lal.LIGOTimeGPS).coalesce() & segments.segmentlist([self.segment])

	def trig_segment(self):
		# interval in which triggers must be read in order to
		# produce a plot
		return self.segment


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

def CacheURLs(cachename, seg):
	"""
	Return a list of CacheEntrys for files containing triggers of
	interest.
	"""
	#urls = [c.url for c in map(lal.CacheEntry, file(cachename)) if c.segment.intersects(seg)]
	urls = []
	for line in file(cachename):
		c = lal.CacheEntry(line)
		if c.segment.intersects(seg):
			urls.append(c.url)
	return urls

def gettriggers(plotdesc):
	# load documents containing relevant triggers
	doc = ligolw.Document()
	handler = lsctables.LIGOLWContentHandler(doc)
	for url in CacheURLs(plotdesc.cachefile, plotdesc.segment):
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
		docutils.MergeTables(tables[0], tables[1])
		del tables[1]

	# FIXME: need to cluster!

	return tables[0]


#
# How to send image to client
#

def SendOutput(plotdesc):
	print >>sys.stdout, "Content-Type: image/png\n"
	shutil.copyfileobj(file(plotdesc.filename), sys.stdout)
