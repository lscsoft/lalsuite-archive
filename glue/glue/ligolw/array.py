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
While the ligolw module provides classes and parser support for reading and
writing LIGO Light Weight XML documents, this module supplements that code
with classes and parsers that add intelligence to the in-RAM document
representation.

In particular, the document tree associated with an Array element is
enhanced.  During parsing, the Stream element in this module converts the
character data contained within it into the elements of a numpy array
object.  The array has the appropriate dimensions and type.  When the
document is written out again, the Stream element serializes the array back
into character data.

The array is stored as an attribute of the Array element.
"""


import numpy
import re
import sys
from xml.sax.saxutils import escape as xmlescape
from xml.sax.xmlreader import AttributesImpl as Attributes


from glue import git_version
from glue import iterutils
from glue.ligolw import ligolw
from glue.ligolw import tokenizer
from glue.ligolw import types as ligolwtypes


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date


#
# =============================================================================
#
#                           Array Name Manipulation
#
# =============================================================================
#


#
# Regular expression used to extract the signifcant portion of an array or
# stream name, according to LIGO LW naming conventions.
#


ArrayPattern = re.compile(r"(?P<Name>[a-z0-9_:]+):array\Z")


def StripArrayName(name):
	"""
	Return the significant portion of an array name according to LIGO
	LW naming conventions.
	"""
	try:
		return ArrayPattern.search(name).group("Name")
	except AttributeError:
		return name


def CompareArrayNames(name1, name2):
	"""
	Convenience function to compare two array names according to LIGO
	LW naming conventions.
	"""
	return cmp(StripArrayName(name1), StripArrayName(name2))


def getArraysByName(elem, name):
	"""
	Return a list of arrays with name name under elem.
	"""
	return elem.getElements(lambda e: (e.tagName == ligolw.Array.tagName) and (CompareArrayNames(e.getAttribute("Name"), name) == 0))


#
# =============================================================================
#
#                                  Utilities
#
# =============================================================================
#


def from_array(name, array, dim_names = None):
	"""
	Construct a LIGO Light Weight XML Array document subtree from a
	numpy array object.
	"""
	doc = Array(Attributes({u"Name": u"%s:array" % name, u"Type": ligolwtypes.FromNumPyType[str(array.dtype)]}))
	s = list(array.shape)
	s.reverse()
	for n, dim in enumerate(s):
		attrs = {}
		if dim_names is not None:
			attrs[u"Name"] = dim_names[n]
		child = ligolw.Dim(Attributes(attrs))
		child.pcdata = unicode(dim)
		doc.appendChild(child)
	child = ArrayStream(Attributes({u"Type": u"Local", u"Delimiter": u" "}))
	doc.appendChild(child)
	doc.array = array
	return doc


def get_array(xmldoc, name):
	"""
	Scan xmldoc for an array named name.  Raises ValueError if not
	exactly 1 such array is found.
	"""
	arrays = getArraysByName(xmldoc, name)
	if len(arrays) != 1:
		raise ValueError, "document must contain exactly one %s array" % StripArrayName(name)
	return arrays[0]


#
# =============================================================================
#
#                               Element Classes
#
# =============================================================================
#


def IndexIter(shape):
	return iterutils.MultiIter(*map(range, shape))


class ArrayStream(ligolw.Stream):
	"""
	High-level Stream element for use inside Arrays.  This element
	knows how to parse the delimited character stream into the parent's
	array attribute, and knows how to turn the parent's array attribute
	back into a character stream.
	"""
	def __init__(self, attrs):
		ligolw.Stream.__init__(self, attrs)
		self._tokenizer = tokenizer.Tokenizer(self.getAttribute(u"Delimiter"))
		self._index = None

	def config(self, parentNode):
		# some initialization that can only be done once parentNode
		# has been set.
		self._tokenizer.set_types([parentNode.pytype])
		parentNode.array = numpy.zeros(parentNode.get_shape(), parentNode.arraytype)
		self._index = iter(IndexIter(parentNode.array.shape))
		return self

	def appendData(self, content):
		# tokenize buffer, and assign to array
		a = self.parentNode.array
		n = self._index.next
		for token in self._tokenizer.append(content):
			a[n()] = token

	def unlink(self):
		"""
		Break internal references within the document tree rooted
		on this element to promote garbage collection.
		"""
		self._index = None
		ligolw.Stream.unlink(self)

	def write(self, file = sys.stdout, indent = u""):
		# This is complicated because we need to not put a
		# delimiter after the last element.
		delim = self.getAttribute(u"Delimiter")
		format = ligolwtypes.FormatFunc[self.parentNode.getAttribute(u"Type")]
		a = self.parentNode.array
		n = iter(IndexIter(a.shape)).next
		file.write(self.start_tag(indent))
		try:
			indeces = n()
			while True:
				if not indeces[0]:
					file.write(u"\n" + indent + ligolw.Indent)
				file.write(xmlescape(format(a[indeces])))
				indeces = n()
				file.write(delim)
		except StopIteration:
			file.write(u"\n" + self.end_tag(indent) + u"\n")


class Array(ligolw.Array):
	"""
	High-level Array element.
	"""
	def __init__(self, *attrs):
		"""
		Initialize a new Array element.
		"""
		ligolw.Array.__init__(self, *attrs)
		self.pytype = ligolwtypes.ToPyType[self.getAttribute(u"Type")]
		self.arraytype = ligolwtypes.ToNumPyType[self.getAttribute(u"Type")]
		self.array = None

	def get_shape(self):
		"""
		Return a tuple of this array's dimensions.  This is done by
		querying the Dim children.  Note that once it has been
		created, it is also possible to examine an Array object's
		array attribute directly, and doing that is much faster.
		"""
		s = [int(c.pcdata) for c in self.getElementsByTagName(ligolw.Dim.tagName)]
		s.reverse()
		return tuple(s)


	#
	# Element methods
	#

	def unlink(self):
		"""
		Break internal references within the document tree rooted
		on this element to promote garbage collection.
		"""
		ligolw.Array.unlink(self)
		self.array = None


#
# =============================================================================
#
#                               Content Handler
#
# =============================================================================
#


#
# Override portions of ligolw.LIGOLWContentHandler class
#


__parent_startStream = ligolw.LIGOLWContentHandler.startStream
__parent_endStream = ligolw.LIGOLWContentHandler.endStream


def startStream(self, attrs):
	if self.current.tagName == ligolw.Array.tagName:
		return ArrayStream(attrs).config(self.current)
	return __parent_startStream(self, attrs)


def endStream(self):
	# stream tokenizer uses delimiter to identify end of each token, so
	# add a final delimiter to induce the last token to get parsed.
	if self.current.parentNode.tagName == ligolw.Array.tagName:
		self.current.appendData(self.current.getAttribute(u"Delimiter"))
	else:
		__parent_endStream(self)


def startArray(self, attrs):
	return Array(attrs)


ligolw.LIGOLWContentHandler.startStream = startStream
ligolw.LIGOLWContentHandler.endStream = endStream
ligolw.LIGOLWContentHandler.startArray = startArray
