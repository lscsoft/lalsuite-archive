#!/usr/bin/python

import os
import popen2

from glue import lal
from glue import LSCsegFindClient
from glue import segments
from glue import segmentsUtils

#
# Some info
#

s5start = lal.LIGOTimeGPS(815155213)
cache = {
	"G1": os.getcwd() + "/G1/filelist.cache",
	"H1": os.getcwd() + "/H1/filelist.cache",
	"H2": os.getcwd() + "/H2/filelist.cache",
	"L1": os.getcwd() + "/L1/filelist.cache"
}


#
# How to run lalapps_lladd
#

class LLAddCommand(object):
	def __init__(self, urls, output = None):
		self._exec = "/home/kipp/local/bin/ligolw_add"
		self.urls = urls
		self.output = output

	def __str__(self):
		s = self._exec
		for url in self.urls:
			s += " \"" + url + "\""
		if self.output:
			s += " --output=\"" + self.output + "\""
		return s

def runlladd(command):
	if type(command) != LLAddCommand:
		raise ValueError, "invalid argument to runlladd(command): command must type LLAddCommand"
	child = popen2.Popen3(str(command), True)
	for line in child.childerr:
		pass
	result = reduce(str.__add__, child.fromchild, "")
	status = child.wait()
	if not os.WIFEXITED(status) or os.WEXITSTATUS(status):
		raise Exception, "failure running \"" + str(command) + "\""
	return result


#
# Trigger file segment lists
#

class TrigSegs(object):
	def __init__(self):
		self.G1 = segmentsUtils.fromlalcache(file(cache["G1"]), coltype = lal.LIGOTimeGPS).coalesce()
		self.H1 = segmentsUtils.fromlalcache(file(cache["H1"]), coltype = lal.LIGOTimeGPS).coalesce()
		self.H2 = segmentsUtils.fromlalcache(file(cache["H2"]), coltype = lal.LIGOTimeGPS).coalesce()
		self.L1 = segmentsUtils.fromlalcache(file(cache["L1"]), coltype = lal.LIGOTimeGPS).coalesce()


#
# Segment querying
#

class SegFindConfig(object):
	def __init__(self, host, port, instrument):
		self.host = host
		self.port = port
		self.instrument = instrument

SegFindConfigH1 = SegFindConfig("ldas.ligo-wa.caltech.edu", None, "H1")
SegFindConfigH2 = SegFindConfig("ldas.ligo-wa.caltech.edu", None, "H2")
SegFindConfigL1 = SegFindConfig("ldas.ligo-la.caltech.edu", None, "L1")

def getsegments(config, types, bounds):
	if config.port:
		client = LSCsegFindClient.LSCsegFind(config.host, config.port)
	else:
		client = LSCsegFindClient.LSCsegFind(config.host)
	list = client.findStateSegments({"interferometer" : config.instrument, "type" : types, "start" : str(int(bounds[0])), "end" : str(int(bounds[1])), "lfns" : False, "strict" : True})
	return segments.segmentlist([segments.segment(*map(lal.LIGOTimeGPS, seg)) for seg in list])

