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
		self.tspec = tspec

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
		table.parentNode.removeChild(table)
	del doc

	# merge trigger tables
	while len(tables) > 1:
		docutils.MergeTables(tables[0], tables[1])
		del tables[1]

	# FIXME: need to cluster!

	return tables[0]

