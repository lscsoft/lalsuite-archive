# $Id$
#
# Copyright (C) 2006  Kipp C. Cannon
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
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
Library of utility code for LIGO Light Weight XML applications.
"""


import codecs
import gzip
import md5
import os
import urllib2
import urlparse
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
	l = []
	for filename in filenames:
		if filename is None:
			# used internally by most ligolw codes to indicate
			# stdin
			l.append((0, None))
		else:
			l.append((os.stat(filename)[stat.ST_SIZE], filename))
	l.sort()
	if reverse:
		l.reverse()
	return l


def sort_files_by_size(filenames, verbose = False, reverse = False):
	"""
	Return a new list of the file names sorted in order of smallest
	file to largest file (or largest to smallest if reverse is set to
	True).
	"""
	if verbose:
		if reverse:
			print >>sys.stderr, "sorting files from largest to smallest ..."
		else:
			print >>sys.stderr, "sorting files from smallest to largest ..."
	return [pair[1] for pair in measure_file_sizes(filenames, reverse = reverse)]


class MD5File(file):
	def __init__(self, fileobj, md5obj = None):
		self.f = fileobj
		if md5obj is None:
			self.md5obj = md5.new()
		else:
			self.md5obj = md5obj

	def flush(self):
		self.f.flush()

	def next(self):
		buf = self.f.next()
		self.md5obj.update(buf)
		return buf

	def read(self, *args, **kwargs):
		buf = self.f.read(*args, **kwargs)
		self.md5obj.update(buf)
		return buf

	def readline(self, *args, **kwargs):
		buf = self.f.readline(*args, **kwargs)
		self.md5obj.update(buf)
		return buf

	def readlines(self, *args, **kwargs):
		raise NotImplementedError

	def xreadlines(self, *args, **kwargs):
		raise NotImplementedError

	def seek(self, *args, **kwargs):
		raise NotImplementedError

	def write(self, buf):
		self.md5obj.update(buf)
		return self.f.write(buf)

	def writelines(self, *args, **kwargs):
		raise NotImplementedError


def load_filename(filename, verbose = False, gz = False, xmldoc = None):
	"""
	Parse the contents of the file identified by filename, and return
	the contents as a LIGO Light Weight document tree.  Helpful
	verbosity messages are printed to stderr if verbose is True, and
	the file is gzip decompressed while reading if gz is set to True.
	If filename is None, then stdin is parsed.  If the optional xmldoc
	argument is provided and not None, the parsed XML tree will be
	appended to that document, otherwise a new document will be
	created.

	Example:

	>>> from glue.ligolw import utils
	>>> xmldoc = utils.load_filename(name, verbose = True, gz = (name or "stdin").endswidth(".gz"))
	"""
	if verbose:
		print >>sys.stderr, "reading %s ..." % (filename and ("'%s'" % filename) or "stdin")
	if filename is not None:
		fileobj = file(filename)
	else:
		fileobj = sys.stdin
	if gz:
		fileobj = gzip.GzipFile(mode = "rb", fileobj = fileobj)
	# FIXME: enable this when I know it works
	#fileobj = MD5File(fileobj)
	#md5obj = fileobj.md5obj
	if xmldoc is None:
		xmldoc = ligolw.Document()
	ligolw.make_parser(ContentHandler(xmldoc)).parse(fileobj)
	return xmldoc


def load_url(url, verbose = False, gz = False, xmldoc = None):
	"""
	This function has the same behaviour as load_filename() but accepts
	a URL instead of a filename.  Any source from which Python's
	urllib2 library can read data is acceptable.  stdin is parsed if
	the URL is None.  If the optional xmldoc argument is provided and
	is not None, the parsed XML tree will be appended to that document,
	otherwise a new document will be created.

	Example:

	>>> from glue.ligolw import utils
	>>> xmldoc = utils.load_url("file://localhost/tmp/data.xml")

	Bugs:
	  - Due to limitations in Python's gzip support and in the way its
	    URL library transfers data, it is not possible to read gzipped
	    XML files from remote locations.  Reading gzipped XML files
	    locally should work correctly.
	"""
	if verbose:
		print >>sys.stderr, "reading %s ..." % (url and ("'%s'" % url) or "stdin")
	if url is not None:
		# hack to detect local files:  urlopen() returns an object
		# that does not support seeking, which prevents GzipFile
		# from working correctly;  by opening local files as
		# regular files, this gets gzip support working for the
		# local case;  still can't read .xml.gz files from a remote
		# host, though;  what's needed is some kind of buffering
		# wrapper to provide seeking (I don't think GzipFile wants
		# to seek very far, so it wouldn't need a big buffer)
		(scheme, host, path, nul, nul, nul) = urlparse.urlparse(url)
		if scheme.lower() in ("", "file") and host.lower() in ("", "localhost"):
			fileobj = file(path)
		else:
			fileobj = urllib2.urlopen(url)
	else:
		fileobj = sys.stdin
	if gz:
		fileobj = gzip.GzipFile(mode = "rb", fileobj = fileobj)
	# FIXME: enable this when I know it works
	#fileobj = MD5File(fileobj)
	#md5obj = fileobj.md5obj
	if xmldoc is None:
		xmldoc = ligolw.Document()
	ligolw.make_parser(ContentHandler(xmldoc)).parse(fileobj)
	return xmldoc


def write_filename(xmldoc, filename, verbose = False, gz = False):
	"""
	Writes the LIGO Light Weight document tree rooted at xmldoc to the
	file name filename.  Friendly verbosity messages are printed while
	doing so if verbose is True.  The output data is gzip compressed on
	the fly if gz is True.
	
	This function traps SIGTERM and SIGTSTP during the write process,
	and it does this by temporarily installing its own signal handlers
	in place of the current handlers.  This is done to prevent Condor
	eviction during the write process.  If a signal is trapped, then
	when the write process has successfully concluded, the last thing
	this function does is raise IOTrappedSignal, with the most-recently
	trapped signal number as the argument.  This is the only condition
	in which this function will raise that exception, so calling code
	that wishes its own handler to be executed can arrange for that to
	happen by trapping the IOTrappedSignal exception.

	Example:

	>>> from glue.ligolw import utils
	>>> utils.write_filename(xmldoc, "data.xml")
	"""
	# initialize SIGTERM and SIGTSTP trap
	global __llwapp_write_filename_got_sig
	__llwapp_write_filename_got_sig = []
	def newsigterm(signum, frame):
		global __llwapp_write_filename_got_sig
		__llwapp_write_filename_got_sig.append(signum)
	oldhandlers = {}
	for sig in (signal.SIGTERM, signal.SIGTSTP):
		oldhandlers[sig] = signal.getsignal(sig)
		signal.signal(sig, newsigterm)

	# write the document
	if verbose:
		print >>sys.stderr, "writing %s ..." % (filename and ("'%s'" % filename) or "stdout")
	if filename is not None:
		fileobj = file(filename, "w")
	else:
		fileobj = sys.stdout
	if gz:
		fileobj = gzip.GzipFile(mode = "wb", fileobj = fileobj)
	# FIXME: enable this when I know it works
	#fileobj = MD5File(fileobj)
	#md5obj = fileobj.md5obj
	fileobj = codecs.EncodedFile(fileobj, "unicode_internal", "utf_8")
	xmldoc.write(fileobj)
	fileobj.flush()

	# restore original handlers, and report the most recently trapped
	# signal if any were
	for sig, oldhandler in oldhandlers.iteritems():
		signal.signal(sig, oldhandler)
	if __llwapp_write_filename_got_sig:
		raise IOTrappedSignal(__llwapp_write_filename_got_sig.pop())
