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
High-level support for Param elements.
"""

__author__ = "Kipp Cannon <kipp@gravity.phys.uwm.edu>"
__date__ = "$Date$"[7:-2]
__version__ = "$Revision$"[11:-2]

import ligolw
import tokenizer
import types


#
# =============================================================================
#
#                           Param Name Manipulation
#
# =============================================================================
#

# Regular expression used to extract the signifcant portion of a param
# name, according to LIGO LW naming conventions.

ParamPattern = re.compile(r"(?:\A[a-z0-9_]+:|\A)(?P<Name>[a-z0-9_]+):param\Z")


def StripParamName(name):
	"""
	Return the significant portion of a param name according to LIGO LW
	naming conventions.
	"""
	try:
		return ParamPattern.search(name).group("Name")
	except AttributeError:
		return name


def CompareParamNames(name1, name2):
	"""
	Convenience function to compare two param names according to LIGO
	LW naming conventions.
	"""
	return cmp(StripParamName(name1), StripParamName(name2))


def getParamsByName(elem, name):
	"""
	Return a list of params with name name under elem.
	"""
	return elem.getElements(lambda e: (e.tagName == ligolw.Param.tagName) and (CompareParamNames(e.getAttribute("Name"), name) == 0))


#
# =============================================================================
#
#                                  Utilities
#
# =============================================================================
#

def new_param(name, type, value, comment = None):
	"""
	Construct a LIGO Light Weight XML Param document subtree.
	"""
	elem = Param({"Name": name, "Type": type})
	elem.pcdata = value
	if comment != None:
		elem.appendChild(ligolw.Comment())
		elem.childNodes[-1].pcdata = comment
	return elem


#
# =============================================================================
#
#                               Element Classes
#
# =============================================================================
#

class Param(ligolw.Param):
	"""
	High-level Param element.  The parsed value is stored in the pcdata
	attribute.
	"""
	def __init__(self, *attrs):
		"""
		Initialize a new Param element.
		"""
		ligolw.Param.__init__(self, *attrs)
		self.tokenizer = tokenizer.Tokenizer(" ")
		t = self.getAttribute("Type")
		if t in types.IntTypes:
			self.tokenizer.set_types([int])
		elif t in types.FloatTypes:
			self.tokenizer.set_types([float])
		elif t in types.StringTypes:
			self.tokenizer.set_types([str])
		else:
			raise TypeError, t
		self.pcdata = None

	def appendData(self, content):
		for token in self.tokenizer.add(content):
			if self.pcdata != None:
				raise ligolw.ElementError, "extra data %s in Param" % content
			self.pcdata = token

	def write(self, file = sys.stdout, indent = ""):
		print >>file, self.start_tag(indent)
		for c in self.childNodes:
			if c.tagName not in self.validchildren:
				raise ElementError, "invalid child %s for %s" % (c.tagName, self.tagName)
			c.write(file, indent + Indent)
		if self.pcdata:
			if self.getAttribute("Type") in StringTypes:
				print >>file, indent + Indent + "\"" + self.pcdata + "\""
			else:
				print >>file, indent + Indent + str(self.pcdata)
		print >>file, self.end_tag(indent)


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

def startParam(self, attrs):
	return Param(attrs)

def endParam(self):
	# param tokenizer uses delimiter to identify end of token, so add a
	# final delimiter to induce the last token to get parsed.
	self.current.appendData(" ")

ligolw.LIGOLWContentHandler.startParam = startParam
ligolw.LIGOLWContentHandler.endParam = endParam
