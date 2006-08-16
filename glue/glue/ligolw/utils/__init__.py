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
Collection of library back-ends for LIGO Light Weight XML utilities.
"""

import gzip
import os
import urllib
import signal
import stat
import sys

from glue.ligolw import ligolw

__author__ = "Kipp Cannon <kipp@gravity.phys.uwm.edu>"
__date__ = "$Date$"[7:-2]
__version__ = "$Revision$"[11:-2]

__all__ = []


#
# =============================================================================
#
#                                 Input/Output
#
# =============================================================================
#

class IOTrappedSignal(Exception):
	"""
	Raised by I/O functions upon completion if they trapped a signal
	during the operation
	"""
	def __init__(self, signum):
		self.signum = signum

	def __str__(self):
		return "trapped signal %d" % self.signum


ContentHandler = ligolw.LIGOLWContentHandler


def measure_file_sizes(filenames, reverse = False):
	"""
	From a list of file names, return a list of (size, name) tuples
	sorted in ascending order by size (or descending order if reverse
	is set to True).
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


def load_filename(filename, verbose = False, gz = False):
	if verbose:
		print >>sys.stderr, "reading %s ..." % (filename or "stdin")
	doc = ligolw.Document()
	if filename:
		fileobj = file(filename)
	else:
		fileobj = sys.stdin
	if gz:
		fileobj = gzip.GzipFile(mode = "rb", fileobj = fileobj)
	ligolw.make_parser(ContentHandler(doc)).parse(fileobj)
	return doc


def load_url(url, verbose = False, gz = False):
	if verbose:
		print >>sys.stderr, "reading %s ..." % (url or "stdin")
	doc = ligolw.Document()
	if url:
		fileobj = urllib.urlopen(url)
	else:
		fileobj = sys.stdin
	if gz:
		fileobj = gzip.GzipFile(mode = "rb", fileobj = fileobj)
	ligolw.make_parser(ContentHandler(doc)).parse(fileobj)
	return doc


def write_filename(doc, filename, verbose = False, gz = False):
	# trap SIGTERM to prevent Condor eviction while the document is
	# being written
	global __llwapp_write_filename_got_sigterm
	__llwapp_write_filename_got_sigterm = False
	def newsigterm(signum, frame):
		global __llwapp_write_filename_got_sigterm
		__llwapp_write_filename_got_sigterm = True
	oldsigterm = signal.getsignal(signal.SIGTERM)
	signal.signal(signal.SIGTERM, newsigterm)

	# write the document
	if verbose:
		print >>sys.stderr, "writing %s ..." % (filename or "stdout")
	if filename:
		fileobj = file(filename, "w")
	else:
		fileobj = sys.stdout
	if gz:
		fileobj = gzip.GzipFile(mode = "wb", fileobj = fileobj)
	doc.write(fileobj)

	# restore original signal handler, and report the signal if it was
	# received
	signal.signal(signal.SIGTERM, oldsigterm)
	if __llwapp_write_filename_got_sigterm:
		raise IOTrappedSignal(signal.SIGTERM)
