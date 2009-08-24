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
try:
	# >= 2.5.0
	from hashlib import md5
except ImportError:
	# < 2.5.0
	from md5 import new as md5
import os
import urllib2
import urlparse
import signal
import stat
import sys

try:
	os.SEEK_SET
except:
	# pre Python 2.5.x is missing these symbols
	os.SEEK_SET, os.SEEK_CUR, os.SEEK_END = range(3)


from glue.ligolw import ligolw


__author__ = "Kipp Cannon <kcannon@ligo.caltech.edu>"
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
	l.sort(reverse = reverse)
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
	return [filename for size, filename in measure_file_sizes(filenames, reverse = reverse)]


class RewindableInputFile(object):
	"""
	DON'T EVER USE THIS FOR ANYTHING!  I'M NOT EVEN KIDDING!
	"""
	# The GzipFile class in Python's standard library is, in my
	# opinion, somewhat weak.  Instead of relying on the return values
	# from the file object's .read() method, GzipFile checks for EOF
	# using calls to .seek().  Furthermore, it uses .seek() instead of
	# buffering data internally as required.  This makes GzipFile
	# gratuitously unable to work with pipes, urlfile objects, and
	# anything else that does not support seeking (including the
	# MD5File class in this module).  To hack around this, this class
	# provides the buffering needed by GzipFile.  It also does proper
	# EOF checking, and uses the results to emulate the results of
	# GzipFile's .seek() games.
	#
	# By wrapping your file object in this class before passing it to
	# GzipFile, you can use GzipFile to read from non-seekable files.
	#
	# How GzipFile checks for EOF == call .tell() to get current
	# position, seek to end of file with .seek(0, 2), call .tell()
	# again and check if the number has changed from before, if it has
	# then we weren't at EOF so call .seek() with original position.

	def __init__(self, fileobj, buffer_size = 16384):
		# the real source of data
		self.fileobj = fileobj
		# how many octets of the internal buffer to return before
		# getting more data
		self.reuse = 0
		# the internal buffer
		self.buf = buffer(" " * buffer_size)
		# flag indicating a .seek()-based EOF test is in progress
		self.gzip_hack_pretend_to_be_at_eof = False

	def __iter__(self):
		return self

	def next(self):
		if self.gzip_hack_pretend_to_be_at_eof:
			return buffer()
		if self.reuse:
			buf = self.buf[-self.reuse:]
			self.reuse = 0
		else:
			buf = self.fileobj.next()
			self.buf = (self.buf + buf)[-len(self.buf):]
		return buf

	def read(self, size = None):
		if self.gzip_hack_pretend_to_be_at_eof:
			return buffer()
		if self.reuse:
			if self.reuse < 0:
				buf = self.fileobj.read(size - self.reuse)
				self.buf = (self.buf + buf)[-len(self.buf):]
				buf = buf[-self.reuse:]
				self.reuse = 0
			elif 0 <= size < self.reuse:
				buf = self.buf[-self.reuse:-self.reuse + size]
				self.reuse -= size
			else:
				buf = self.buf[-self.reuse:]
				self.reuse = 0
				if len(buf) < size:
					buf += self.read(size - len(buf))
		else:
			buf = self.fileobj.read(size)
			self.buf = (self.buf + buf)[-len(self.buf):]
		return buf

	def seek(self, offset, whence = os.SEEK_SET):
		self.gzip_hack_pretend_to_be_at_eof = False
		if whence == os.SEEK_SET:
			pos = self.fileobj.tell()
			if offset >= 0 and pos - len(self.buf) <= offset <= pos:
				self.reuse = pos - offset
			else:
				raise IOError, "seek out of range"
		elif whence == os.SEEK_CUR:
			if self.reuse - len(self.buf) <= offset:
				self.reuse -= offset
			else:
				raise IOError, "seek out of range"
		elif whence == os.SEEK_END:
			if offset == 0:
				self.gzip_hack_pretend_to_be_at_eof = True
			else:
				raise IOError, "seek out of range"

	def tell(self):
		if self.gzip_hack_pretend_to_be_at_eof:
			# check to see if we are at EOF by seeing if we can
			# read 1 character.  save it in the internal buffer
			# to not loose it.
			c = self.fileobj.read(1)
			self.buf = (self.buf + c)[-len(self.buf):]
			self.reuse += len(c)
			if c:
				# since we have read a character, this will
				# not return the same answer as when
				# GzipFile called it
				return self.fileobj.tell()
		return self.fileobj.tell() - self.reuse

	def close(self):
		return self.fileobj.close()


class MD5File(object):
	def __init__(self, fileobj, md5obj = None):
		self.fileobj = fileobj
		if md5obj is None:
			self.md5obj = md5()
		else:
			self.md5obj = md5obj

	def __iter__(self):
		return self

	def next(self):
		buf = self.fileobj.next()
		self.md5obj.update(buf)
		return buf

	def read(self, size = None):
		buf = self.fileobj.read(size)
		self.md5obj.update(buf)
		return buf

	def write(self, buf):
		self.md5obj.update(buf)
		return self.fileobj.write(buf)

	def tell(self):
		return self.fileobj.tell()

	def flush(self):
		return self.fileobj.flush()

	def close(self):
		return self.fileobj.close()


def load_fileobj(fileobj, gz = False, xmldoc = None):
	"""
	Parse the contents of the file object fileobj, and return the
	contents as a LIGO Light Weight document tree.  The file object
	does not need to be seekable.  The file is gzip decompressed while
	reading if gz is set to True.  If the optional xmldoc argument is
	provided and not None, the parsed XML tree will be appended to that
	document, otherwise a new document will be created.  The return
	value is a tuple, the first element of the tuple is the XML
	document and the second is a string containing the MD5 digest in
	hex digits of the bytestream that was parsed.

	Example:

	>>> import sys
	>>> xmldoc, digest = utils.load_fileobj(sys.stdin, verbose = True, gz = True)
	"""
	fileobj = MD5File(fileobj)
	md5obj = fileobj.md5obj
	if gz:
		fileobj = gzip.GzipFile(mode = "rb", fileobj = RewindableInputFile(fileobj))
	if xmldoc is None:
		xmldoc = ligolw.Document()
	ligolw.make_parser(ContentHandler(xmldoc)).parse(fileobj)
	return xmldoc, md5obj.hexdigest()


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
	xmldoc, hexdigest = load_fileobj(fileobj, gz = gz, xmldoc = xmldoc)
	if verbose:
		print >>sys.stderr, "md5sum: %s  %s" % (hexdigest, filename or "")
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
	"""
	if verbose:
		print >>sys.stderr, "reading %s ..." % (url and ("'%s'" % url) or "stdin")
	if url is not None:
		scheme, host, path, nul, nul, nul = urlparse.urlparse(url)
		if scheme.lower() in ("", "file") and host.lower() in ("", "localhost"):
			fileobj = file(path)
		else:
			fileobj = urllib2.urlopen(url)
	else:
		fileobj = sys.stdin
	xmldoc, hexdigest = load_fileobj(fileobj, gz = gz, xmldoc = xmldoc)
	if verbose:
		print >>sys.stderr, "md5sum: %s  %s" % (hexdigest, url or "")
	return xmldoc


def write_fileobj(xmldoc, fileobj, gz = False, xsl_file = None):
	"""
	Writes the LIGO Light Weight document tree rooted at xmldoc to the
	given file object.  The file object need not be seekable.  The
	output data is gzip compressed on the fly if gz is True.  The
	return value is a string containing the hex digits of the MD5
	digest of the output bytestream.

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

	>>> import sys
	>>> utils.write_fileobj(xmldoc, sys.stdout)
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
	fileobj = MD5File(fileobj)
	md5obj = fileobj.md5obj
	if gz:
		fileobj = gzip.GzipFile(mode = "wb", fileobj = fileobj)
	fileobj = codecs.EncodedFile(fileobj, "unicode_internal", "utf_8")
	xmldoc.write(fileobj, xsl_file = xsl_file)
	fileobj.flush()
	del fileobj

	# restore original handlers, and report the most recently trapped
	# signal if any were
	for sig, oldhandler in oldhandlers.iteritems():
		signal.signal(sig, oldhandler)
	if __llwapp_write_filename_got_sig:
		raise IOTrappedSignal(__llwapp_write_filename_got_sig.pop())

	# return the hex digest of the bytestream that was written
	return md5obj.hexdigest()


def write_filename(xmldoc, filename, verbose = False, gz = False, xsl_file = None):
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
	if verbose:
		print >>sys.stderr, "writing %s ..." % (filename and ("'%s'" % filename) or "stdout")
	if filename is not None:
		fileobj = file(filename, "w")
	else:
		fileobj = sys.stdout
	hexdigest = write_fileobj(xmldoc, fileobj, gz = gz, xsl_file = xsl_file)
	fileobj.close()
	if verbose:
		print >>sys.stderr, "md5sum: %s  %s" % (hexdigest, filename or "")


def write_url(xmldoc, url, verbose = False, gz = False):
	"""
	Writes the LIGO Light Weight document tree rooted at xmldoc to the
	URL name url.  Friendly verbosity messages are printed while doing
	so if verbose is True.  The output data is gzip compressed on the
	fly if gz is True.

	See write_filename() for more information about signal trapping.

	NOTE:  only URLs that point to local files can be written to at
	this time.
	
	Example:

	>>> from glue.ligolw import utils
	>>> utils.write_url(xmldoc, "file:///data.xml")
	"""
	if url is None:
		scheme, host, path = "", "", None
	else:
		scheme, host, path, nul, nul, nul = urlparse.urlparse(url)
	if scheme.lower() in ("", "file") and host.lower() in ("", "localhost"):
		return write_filename(xmldoc, path, verbose = verbose, gz = gz)
	else:
		raise ValueError, "%s is not a local file" % repr(url)
