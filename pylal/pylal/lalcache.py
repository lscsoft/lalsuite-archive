"""
LAL cache file.
"""

__author__ = "Kipp Cannon <kipp@gravity.phys.uwm.edu>"
__date__ = "$Date$"
__version__ = "$Revision$"


#
# Preamble
#

import re

from glue import lal
from glue import segments


#
# How to parse a line in a cache file.  Five white-space delimited columns.
#

CacheRegEx = re.compile(r"\A\s*(?P<observatory>\S+)\s+(?P<description>\S+)\s+(?P<start>\S+)\s+(?P<duration>\S+)\s+(?P<url>\S+)\s*\Z")


#
# Cache file entry
#

class CacheEntry(object):
	def __init__(self, string = None):
		if type(string) == type(None):
			self.observatory = None
			self.description = None
			self.segment = segments.segment(None, None)
			self.url = None
		else:
			match = CacheRegEx.search(string)
			if not match:
				raise ValueError, "could not convert \"%s\" to CacheEntry" % string
			self.observatory = match.group("observatory")
			self.description = match.group("description")
			self.segment = segments.segment(lal.LIGOTimeGPS(match.group("start")), lal.LIGOTimeGPS(match.group("start")) + lal.LIGOTimeGPS(match.group("duration")))
			self.url = match.group("url")

	def __str__(self):
		return "%s %s %s %s %s" % (str(self.observatory), str(self.description), str(self.segment[0]), str(self.segment.duration()), str(self.url))

	def __cmp__(self, other):
		if type(other)  != CacheEntry:
			raise TypeError, "can only compare CacheEntry to CacheEntry"
		return cmp((self.observatory, self.description, self.segment, self.url), (other.observatory, other.description, other.segment, other.url))
