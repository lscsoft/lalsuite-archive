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
import pickle
import socket
import stat
import sys
import time
import urllib

from glue import segments
from glue.ligolw import ligolw
from glue.ligolw import metaio
from glue.ligolw import param
from glue.ligolw import array
from glue.ligolw import lsctables
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

ContentHandler = ligolw.LIGOLWContentHandler

def measure_file_sizes(filenames, reverse = False):
	"""
	From a list of file names, return a list of (size, name) tuples
	sorted in ascending order by size.
	"""
	l = [(os.stat(name)[stat.ST_SIZE], name) for name in filenames if name]
	l.sort()
	if reverse:
		l.reverse()
	return l


def sort_files_by_size(filenames, verbose = False, reverse = False):
	"""
	Sort the file names from smallest file to largest file (or largest
	to smallest if reverse is set to True).
	"""
	if verbose:
		if reverse:
			print >>sys.stderr, "sorting files from largest to smallest..."
		else:
			print >>sys.stderr, "sorting files from smallest to largest..."
	return [pair[1] for pair in measure_file_sizes(filenames, reverse = reverse)]


def load_filename(filename, verbose = False):
	if verbose:
		print >>sys.stderr, "reading %s..." % (filename or "stdin")
	doc = ligolw.Document()
	if filename:
		ligolw.make_parser(ContentHandler(doc)).parse(file(filename))
	else:
		ligolw.make_parser(ContentHandler(doc)).parse(sys.stdin)
	return doc


def load_url(url, verbose = False):
	if verbose:
		print >>sys.stderr, "reading %s..." % (url or "stdin")
	doc = ligolw.Document()
	if url:
		ligolw.make_parser(ContentHandler(doc)).parse(urllib.urlopen(url))
	else:
		ligolw.make_parser(ContentHandler(doc)).parse(sys.stdin)
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


def segmentlistdict_fromsearchsummary(xmldoc, live_time_program = None):
	"""
	Convenience wrapper for a common case usage of the segmentlistdict
	class: searches the process table in xmldoc for occurances of a
	program named live_time_program, then scans the search summary
	table for matching process IDs and constructs a segmentlistdict
	object from those rows.
	"""
	procids = get_process_ids_by_program(xmldoc, live_time_program)
	seglistdict = segments.segmentlistdict()
	for row in get_table(xmldoc, lsctables.SearchSummaryTable.tableName):
		if (live_time_program == None) or bisect_contains(procids, row.process_id):
			if row.ifos in seglistdict:
				seglistdict[row.ifos].append(row.get_out())
			else:
				seglistdict[row.ifos] = segments.segmentlist([row.get_out()])
	return seglistdict.coalesce()


#
# =============================================================================
#
#                                    Arrays
#
# =============================================================================
#

def get_array(doc, name):
	"""
	Scan doc for an array named name.  Raises ValueError if not exactly
	1 such array is found.
	"""
	arrays = array.getArraysByName(doc, name)
	if len(arrays) != 1:
		raise ValueError, "document must contain exactly one %s array" % array.StripArrayName(name)
	return arrays[0]


#
# =============================================================================
#
#                                    Params
#
# =============================================================================
#

def get_param(doc, name):
	"""
	Scan doc for a param named name.  Raises ValueError if not exactly
	1 such param is found.
	"""
	params = param.getParamsByName(doc, name)
	if len(params) != 1:
		raise ValueError, "document must contain exactly one %s param" % param.StripParamName(name)
	return params[0]


def pickle_to_param(obj, name):
	"""
	Return the top-level element of a document sub-tree containing the
	pickled serialization of a Python object.
	"""
	# WARNING:  this is an abuse of most of the things involved;  but,
	# hey, it's a one-liner so it's hard to resist.  Do NOT use this in
	# files you care about.  If there proves to be a real desire for
	# this kind of thing, then we need to think of a proper way to do
	# it.  For example, there is at least one library that can
	# serialize a Python object into proper XML.
	return param.new_param("pickle:%s" % name, "lstring", urllib.quote(pickle.dumps(obj)))


def pickle_from_param(elem, name):
	"""
	Retrieve a pickled Python object from the document tree rooted at
	elem.
	"""
	return pickle.loads(urllib.unquote(get_param(elem, "pickle:%s:param" % name).pcdata))


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
	process.process_id = lsctables.NewILWDs(proctable, "process_id").next()
	proctable.append(process)
	return process


def set_process_end_time(process):
	"""
	Set the end time in a row in a process table to the current time.
	"""
	process.end_time = XLALUTCToGPS(time.gmtime()).seconds
	return process


def get_process_ids_by_program(xmldoc, program):
	"""
	Extract the process table, and return a sorted list of the process
	IDs for the given program.
	"""
	ids = [row.process_id for row in get_table(xmldoc, lsctables.ProcessTable.tableName) if row.program == program]
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
#                               Iteration Tools
#
# =============================================================================
#

class MultiIter(object):
	"""
	An iterator class for iterating over the elements of multiple lists
	simultaneously.  An instance of the class is initialized with a
	list of lists.  A call to next() returns a list of elements, one
	from each of the lists.  Subsequent calls to next() iterate over
	all combinations of elements from the lists.
	"""
	def __init__(self, lists):
		self.lists = tuple(lists)
		self.index = [0] * len(lists)
		self.length = tuple(map(len, lists))
		self.stop = 0 in self.length

	def __len__(self):
		return reduce(int.__mul__, self.length)

	def __iter__(self):
		return self

	def next(self):
		if self.stop:
			raise StopIteration
		l = map(lambda l, i: l[i], self.lists, self.index)
		for i in xrange(len(self.index)):
			self.index[i] += 1
			if self.index[i] < self.length[i]:
				break
			self.index[i] = 0
		else:
			self.stop = True
		return l


def choices(vals, n):
	"""
	Return a list of all choices of n elements from the list vals.
	Order is preserved.
	"""
	if n == len(vals):
		return [vals]
	if n == 1:
		return [[v] for v in vals]
	if n < 1:
		raise ValueError, n
	l = []
	for i in xrange(len(vals) - n + 1):
		for c in choices(vals[i+1:], n - 1):
			c.insert(0, vals[i])
			l.append(c)
	return l


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


def get_coincident_segmentlistdict(seglistdict, offsetdictlist):
	"""
	This function answers the question "Given a set of segment lists,
	and a set of time slides to apply to those segment lists, what
	segments do I need to keep in the original lists so that I have all
	the segments that will participate in a coincidence analysis done
	over all those time slides.

	This function constructs and returns a segmentlistdict object that,
	for each key in seglistdict, contains the segments from the
	corresponding list in seglistdict which are coincident under at
	least one of the time slides described by offsetdictlist.

	The elements of offsetlistdict are free to contain only subsets of
	the keys in seglistdict.  In those cases, the coincidence is
	computed only between the lists corresponding to the given keys.

	For example, let us say that "input" is a segmentlistdict object
	containing segment lists for three instruments, "H1", "H2" and
	"L1".  And let us say that "slides" is a list of dictionaries, and
	is equal to [{"H1":0, "H2":0, "L1":0}, {"H1":0, "H2":10}].  Then if

	output = get_coincident_segmentlistdict(input, slides)

	output will contain, for each instrument, the segments (or parts
	thereof) from the original lists that are required in order to
	perform a triple-coincident analysis at zero lag betwen the three
	instruments, *and* a double-coincident analysis between H1 and H2
	with H2 offset by 10 seconds.

	During the computations, the input segmentlistdict object will have
	offsets applied to it in place, but they will be restored to their
	original value upon exit.
	"""
	origoffsets = {}
	origoffsets.update(seglistdict.offsets)
	coincseglists = segments.segmentlistdict()
	for offsetdict in offsetdictlist:
		seglistdict.offsets.update(offsetdict)
		intersection = seglistdict.extract_common(offsetdict.keys())
		intersection.offsets.update(origoffsets)
		coincseglists |= intersection
	seglistdict.offsets.update(origoffsets)
	return coincseglists
