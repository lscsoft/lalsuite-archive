# Copyright (C) 2006  Kipp Cannon
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
import warnings
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


from glue import git_version
from glue.ligolw import ligolw


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date


__all__ = []


#
# =============================================================================
#
#                                 Input/Output
#
# =============================================================================
#


# FIXME:  remove, use parameter passed to load_*() functions instead
ContentHandler = ligolw.LIGOLWContentHandler
__orig_ContentHandler = ContentHandler	# to detect when ContentHandler symbol has been modified


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
	# then we weren't at EOF so call .seek() with original position and
	# keep going.  ?!

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
		self.pos = 0

	def __iter__(self):
		return self

	def next(self):
		buf = self.fileobj.next()
		self.md5obj.update(buf)
		self.pos += len(buf)
		return buf

	def read(self, size = None):
		buf = self.fileobj.read(size)
		self.md5obj.update(buf)
		self.pos += len(buf)
		return buf

	def write(self, buf):
		self.pos += len(buf)
		self.md5obj.update(buf)
		return self.fileobj.write(buf)

	def tell(self):
		try:
			return self.fileobj.tell()
		except (IOError, AttributeError):
			# some streams that don't support seeking, like
			# stdin, report IOError.  the things returned by
			# urllib don't have a .tell() method at all.  fake
			# it without our own count of bytes
			return self.pos

	def flush(self):
		return self.fileobj.flush()

	def close(self):
		return self.fileobj.close()


def load_fileobj(fileobj, gz = None, xmldoc = None, contenthandler = None):
	"""
	Parse the contents of the file object fileobj, and return the
	contents as a LIGO Light Weight document tree.  The file object
	does not need to be seekable.

	If the gz parameter is None (the default) then gzip compressed data
	will be automatically detected and decompressed, otherwise
	decompression can be forced on or off by setting gz to True or
	False respectively.

	If the optional xmldoc argument is provided and not None, the
	parsed XML tree will be appended to that document, otherwise a new
	document will be created.  The return value is a tuple, the first
	element of the tuple is the XML document and the second is a string
	containing the MD5 digest in hex digits of the bytestream that was
	parsed.

	Example:

	>>> import StringIO
	>>> f = StringIO.StringIO('<?xml version="1.0" encoding="utf-8" ?><!DOCTYPE LIGO_LW SYSTEM "http://ldas-sw.ligo.caltech.edu/doc/ligolwAPI/html/ligolw_dtd.txt"><LIGO_LW><Table Name="demo:table"><Column Name="name" Type="lstring"/><Column Name="value" Type="real8"/><Stream Name="demo:table" Type="Local" Delimiter=",">"mass",0.5,"velocity",34</Stream></Table></LIGO_LW>')
	>>> xmldoc, digest = load_fileobj(f)
	>>> digest
	'03d1f513120051f4dbf3e3bc58ddfaa6'

	The optional contenthandler argument allows the SAX content handler
	to be customized.  Previously, customization of the content handler
	was accomplished by replacing the ContentHandler symbol in this
	module with the custom handler, and although that technique is
	still supported a warning will be emitted if modification of that
	symbol is detected.  See
	glue.ligolw.ligolw.PartialLIGOLWContentHandler and
	glue.ligolw.ligolw.FilteringLIGOLWContentHandler for examples of
	custom content handlers used to load subsets of documents into
	memory.
	"""
	fileobj = MD5File(fileobj)
	md5obj = fileobj.md5obj
	if gz != False:
		fileobj = RewindableInputFile(fileobj)
		magic = fileobj.read(2)
		fileobj.seek(0, os.SEEK_SET)
		if gz == True or magic == '\037\213':
			fileobj = gzip.GzipFile(mode = "rb", fileobj = fileobj)
	if xmldoc is None:
		xmldoc = ligolw.Document()
	if contenthandler is None:
		if ContentHandler is not __orig_ContentHandler:
			warnings.warn("modification of glue.ligolw.utils.ContentHandler global variable for input customization is deprecated.  Use contenthandler keyword argument of glue.ligolw.utils.load_*() functions instead", DeprecationWarning)
		contenthandler = ContentHandler
	ligolw.make_parser(contenthandler(xmldoc)).parse(fileobj)
	return xmldoc, md5obj.hexdigest()


def load_filename(filename, verbose = False, gz = None, xmldoc = None, contenthandler = None):
	"""
	Parse the contents of the file identified by filename, and return
	the contents as a LIGO Light Weight document tree.  Helpful
	verbosity messages are printed to stderr if verbose is True.  All
	other parameters are passed verbatim to load_fileobj(), see that
	function for more information.

	Example:

	>>> xmldoc = load_filename(name, verbose = True)
	"""
	if verbose:
		print >>sys.stderr, "reading %s ..." % (filename and ("'%s'" % filename) or "stdin")
	if filename is not None:
		fileobj = open(filename, "rb")
	else:
		fileobj = sys.stdin
	xmldoc, hexdigest = load_fileobj(fileobj, gz = gz, xmldoc = xmldoc, contenthandler = contenthandler)
	if verbose:
		print >>sys.stderr, "md5sum: %s  %s" % (hexdigest, filename or "")
	return xmldoc


def load_url(url, verbose = False, gz = None, xmldoc = None, contenthandler = None):
	"""
	This function has the same behaviour as load_filename() but accepts
	a URL instead of a filename.  Any source from which Python's
	urllib2 library can read data is acceptable.  stdin is parsed if
	the URL is None.

	Example:

	>>> xmldoc = load_url("file://localhost/tmp/data.xml")
	"""
	if verbose:
		print >>sys.stderr, "reading %s ..." % (url and ("'%s'" % url) or "stdin")
	if url is not None:
		scheme, host, path, nul, nul, nul = urlparse.urlparse(url)
		if scheme.lower() in ("", "file") and host.lower() in ("", "localhost"):
			fileobj = open(path)
		else:
			fileobj = urllib2.urlopen(url)
	else:
		fileobj = sys.stdin
	xmldoc, hexdigest = load_fileobj(fileobj, gz = gz, xmldoc = xmldoc, contenthandler = contenthandler)
	if verbose:
		print >>sys.stderr, "md5sum: %s  %s" % (hexdigest, url or "")
	return xmldoc


def write_fileobj(xmldoc, fileobj, gz = False, xsl_file = None, trap_signals = (signal.SIGTERM, signal.SIGTSTP)):
	"""
	Writes the LIGO Light Weight document tree rooted at xmldoc to the
	given file object.  The file object need not be seekable.  The
	output data is gzip compressed on the fly if gz is True.  The
	return value is a string containing the hex digits of the MD5
	digest of the output bytestream.

	This function traps the signals in the trap_signals iterable during
	the write process (the default is signal.SIGTERM and
	signal.SIGTSTP), and it does this by temporarily installing its own
	signal handlers in place of the current handlers.  This is done to
	prevent Condor eviction during the write process.  When the file
	write is concluded the original signal handlers are restored.
	Then, if signals were trapped during the write process, the signals
	are then resent to the current process in the order in which they
	were received.  The signal.signal() system call cannot be invoked
	from threads, and trap_signals must be set to None or an empty
	sequence if this function is used from a thread.

	Example:

	>>> import sys
	>>> write_fileobj(xmldoc, sys.stdout)
	"""
	# initialize SIGTERM and SIGTSTP trap
	deferred_signals = []
	def newsigterm(signum, frame):
		deferred_signals.append(signum)
	oldhandlers = {}
	if trap_signals is not None:
		for sig in trap_signals:
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

	# restore original handlers, and send outselves any trapped signals
	# in order
	for sig, oldhandler in oldhandlers.iteritems():
		signal.signal(sig, oldhandler)
	while deferred_signals:
		os.kill(os.getpid(), deferred_signals.pop(0))

	# return the hex digest of the bytestream that was written
	return md5obj.hexdigest()


def write_filename(xmldoc, filename, verbose = False, gz = False, xsl_file = None, trap_signals = (signal.SIGTERM, signal.SIGTSTP)):
	"""
	Writes the LIGO Light Weight document tree rooted at xmldoc to the
	file name filename.  Friendly verbosity messages are printed while
	doing so if verbose is True.  The output data is gzip compressed on
	the fly if gz is True.

	See write_fileobj() for information about signal trapping during
	the write process.

	Example:

	>>> write_filename(xmldoc, "data.xml")
	"""
	if verbose:
		print >>sys.stderr, "writing %s ..." % (filename and ("'%s'" % filename) or "stdout")
	if filename is not None:
		if not gz and filename.endswith(".gz"):
			warnings.warn("filename '%s' ends in '.gz' but file is not being gzip-compressed" % filename, UserWarning)
		fileobj = open(filename, "w")
	else:
		fileobj = sys.stdout
	hexdigest = write_fileobj(xmldoc, fileobj, gz = gz, xsl_file = xsl_file, trap_signals = trap_signals)
	fileobj.close()
	if verbose:
		print >>sys.stderr, "md5sum: %s  %s" % (hexdigest, filename or "")


def write_url(xmldoc, url, verbose = False, gz = False, xsl_file = None, trap_signals = (signal.SIGTERM, signal.SIGTSTP)):
	"""
	Writes the LIGO Light Weight document tree rooted at xmldoc to the
	URL name url.  Friendly verbosity messages are printed while doing
	so if verbose is True.  The output data is gzip compressed on the
	fly if gz is True.

	See write_fileobj() for information about signal trapping during
	the write process.

	NOTE:  only URLs that point to local files can be written to at
	this time.
	
	Example:

	>>> write_url(xmldoc, "file:///data.xml")
	"""
	if url is None:
		scheme, host, path = "", "", None
	else:
		scheme, host, path, nul, nul, nul = urlparse.urlparse(url)
	if scheme.lower() not in ("", "file") or host.lower() not in ("", "localhost"):
		raise ValueError, "%s is not a local file" % repr(url)
	return write_filename(xmldoc, path, verbose = verbose, gz = gz, xsl_file = xsl_file, trap_signals = trap_signals)
