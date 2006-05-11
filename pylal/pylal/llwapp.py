# $Id$
#
# Copyright (C) 2006  Kipp C. Cannon
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#

"""
A collection of utilities to assist in writing applications that manipulate
data in LIGO Light-Weight format.
"""

import bisect
import os
import socket
import stat
import sys
import time
import urllib

from glue import segments
from glue.ligolw import ligolw
from glue.ligolw import metaio
from glue.ligolw import lsctables
from glue.ligolw import docutils
from pylal.date import XLALUTCToGPS

__author__ = "Kipp Cannon <kipp@gravity.phys.uwm.edu>"
__version__ = "$Revision$"[11:-2]
__date__ = "$Date$"[7:-2]


#
# =============================================================================
#
#                                 Input/Output
#
# =============================================================================
#

def sort_files(filenames, verbose = False):
	"""
	Sort the files from largest to smallest.
	"""
	if verbose:
		print >>sys.stderr, "sorting files from largest to smallest..."
	l = [(os.stat(name)[stat.ST_SIZE], name) for name in filenames]
	l.sort()
	l.reverse()
	return [pair[1] for pair in l]


def load_filename(filename, verbose = False):
	if verbose:
		print >>sys.stderr, "reading %s..." % (filename or "stdin")
	doc = ligolw.Document()
	if filename:
		ligolw.make_parser(lsctables.LIGOLWContentHandler(doc)).parse(file(filename))
	else:
		ligolw.make_parser(lsctables.LIGOLWContentHandler(doc)).parse(sys.stdin)
	return doc


def load_url(url, verbose = False):
	if verbose:
		print >>sys.stderr, "reading %s..." % (url or "stdin")
	doc = ligolw.Document()
	if url:
		ligolw.make_parser(lsctables.LIGOLWContentHandler(doc)).parse(urllib.urlopen(url))
	else:
		ligolw.make_parser(lsctables.LIGOLWContentHandler(doc)).parse(sys.stdin)
	return doc


def write_filename(doc, filename, verbose = False):
	if verbose:
		print >>sys.stderr, "writing %s..." % (filename or "stdout")
	if filename:
		doc.write(file(filename, "w"))
	else:
		doc.write(sys.stdout)


#
# =============================================================================
#
#                                    Tables
#
# =============================================================================
#

def get_table(doc, name):
	"""
	Scan doc for a table named name.  Raises ValueError if not exactly
	1 such table is found.
	"""
	tables = metaio.getTablesByName(doc, name)
	if len(tables) != 1:
		raise ValueError, "document must contain exactly one %s table" % metaio.StripTableName(name)
	return tables[0]


class SearchSummarySegListDict(dict):
	"""
	A dictionary of segmentlist objects constructed from a search
	summary table, with the ability to apply offsets to the segment
	lists.
	"""
	def __new__(cls, searchsummtable, proc_ids):
		"""
		Construct a new SearchSummarySegListDict object from the
		contents of a search summary table.  The result is
		dictionary-like object whose keys are the instrument names
		appearing in the search summary table, and whose values are
		segmentlist representations of the times spanned by the
		rows.  Only those rows whose process_ids appear in the
		sorted proc_ids list are included.
		"""
		seglists = dict.__new__(cls, None)
		for row in searchsummtable:
			if not bisect_contains(proc_ids, row.process_id):
				continue
			if row.ifos in seglists:
				seglists[row.ifos].append(row.get_out())
			else:
				seglists[row.ifos] = segments.segmentlist([row.get_out()])
		for seglist in seglists.itervalues():
			seglist.coalesce()
		return seglists

	def __init__(self, *args):
		self.offsets = {}
		for instrument in self.iterkeys():
			self.offsets[instrument] = 0.0

	def set_offsets(self, offsetdict):
		"""
		Shift the segment lists according to the instrument/offset
		pairs in offsetdict.
		"""
		for instrument, offset in offsetdict.iteritems():
			delta = offset - self.offsets[instrument]
			if delta != 0.0:
				self[instrument] = self[instrument].shift(delta)
				self.offsets[instrument] += delta
		return self

	def remove_offsets(self):
		"""
		Remove the offsets from all segmentlists.
		"""
		for instrument, offset in self.offsets.iteritems():
			if offset:
				self[instrument] = self[instrument].shift(-offset)
				self.offsets[instrument] = 0.0
		return self

	def protract(self, x):
		"""
		Protract the segmentlists.
		"""
		for key in self.iterkeys():
			self[key] = self[key].protract(x)
		return self

	def intersect(self, instruments):
		"""
		Return the intersection of the segmentlists for the
		instruments in instruments.
		"""
		return reduce(segments.segmentlist.__and__, map(self.__getitem__, instruments))


def get_seglistdict(xmldoc, live_time_program):
	"""
	Convenience wrapper for the common case usage of the
	SearchSummarySegListDict:  searches the process table for
	occurances of a program named live_time_program, then scans the
	search summary table for matching process IDs and constructs a
	SearchSummarySegListDict object from those rows.
	"""
	return SearchSummarySegListDict(get_table(xmldoc, lsctables.SearchSummaryTable.tableName), get_process_ids_by_program(get_table(xmldoc, lsctables.ProcessTable.tableName), live_time_program))


#
# =============================================================================
#
#                               Process Metadata
#
# =============================================================================
#

def append_process(doc, program = "", version = "", cvs_repository = "", cvs_entry_time = "", comment = "", is_online = False, jobid = 0, domain = "", ifos = ""):
	"""
	Add an entry to the process table in doc.  program, version,
	cvs_repository, comment, domain, and ifos should all be strings.
	cvs_entry_time should be a string in the format "YYYY/MM/DD
	HH:MM:SS".  is_online should be a boolean, jobid an integer.
	"""
	proctable = get_table(doc, lsctables.ProcessTable.tableName)
	process = lsctables.Process()
	process.program = program
	process.version = version
	process.cvs_repository = cvs_repository
	process.cvs_entry_time = XLALUTCToGPS(time.strptime(cvs_entry_time, "%Y/%m/%d %H:%M:%S")).seconds
	process.comment = comment
	process.is_online = int(is_online)
	process.node = socket.gethostbyaddr(socket.gethostname())[0]
	process.username = os.environ["LOGNAME"]
	process.unix_procid = os.getpid()
	process.start_time = XLALUTCToGPS(time.gmtime()).seconds
	process.end_time = 0
	process.jobid = jobid
	process.domain = domain
	process.ifos = ifos
	process.process_id = docutils.NewILWDs(proctable, "process_id").next()
	proctable.append(process)
	return process


def set_process_end_time(process):
	"""
	Set the end time in a row in a process table to the current time.
	"""
	process.end_time = XLALUTCToGPS(time.gmtime()).seconds
	return process


def get_process_ids_by_program(processtable, program):
	"""
	From a process table, return a sorted list of the process IDs for
	the given program.
	"""
	ids = [row.process_id for row in processtable if row.program == program]
	ids.sort()
	return ids


def append_process_params(doc, process, params):
	"""
	doc is an XML document tree, process is the row in the process
	table for which these are the parameters, and params is a list of
	(name, type, value) tuples one for each parameter.
	"""
	paramtable = get_table(doc, lsctables.ProcessParamsTable.tableName)
	for name, type, value in params:
		param = lsctables.ProcessParams()
		param.program = process.program
		param.process_id = process.process_id
		param.param = str(name)
		param.type = str(type)
		param.value = str(value)
		paramtable.append(param)
	return process


def doc_includes_process(doc, program):
	"""
	Return True if the process table in doc includes entries for a
	program named program.
	"""
	return program in get_table(doc, lsctables.ProcessTable.tableName).getColumnByName("program")


#
# =============================================================================
#
#                               Search Metadata
#
# =============================================================================
#

def append_search_summary(doc, process, shared_object = "standalone", lalwrapper_cvs_tag = "", lal_cvs_tag = "", comment = None, ifos = None, inseg = None, outseg = None, nevents = 0, nnodes = 1):
	"""
	Append search summary information associated with the given process
	to the search summary table in doc.
	"""
	summary = lsctables.SearchSummary()
	summary.process_id = process.process_id
	summary.shared_object = shared_object
	summary.lalwrapper_cvs_tag = lalwrapper_cvs_tag
	summary.lal_cvs_tag = lal_cvs_tag
	summary.comment = comment or process.comment
	summary.ifos = ifos or process.ifos
	summary.set_in(inseg)
	summary.set_out(outseg)
	summary.nevents = nevents
	summary.nnodes = nnodes
	get_table(doc, lsctables.SearchSummaryTable.tableName).append(summary)
	return summary


#
# =============================================================================
#
#                                    Other
#
# =============================================================================
#

def bisect_contains(array, val):
	"""
	Uses a bisection search to determine if val is in array.  Returns
	True or False.
	"""
	try:
		return array[bisect.bisect_left(array, val)] == val
	except IndexError:
		return False

