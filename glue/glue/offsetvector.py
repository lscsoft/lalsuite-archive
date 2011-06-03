# Copyright (C) 2010  Kipp Cannon
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


from glue import git_version


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date


#
# =============================================================================
#
#                                Offset Vector
#
# =============================================================================
#


class offsetvector(dict):
	"""
	Subclass of the dict built-in type for storing mappings of
	instrument to time offset.

	Examples:

	>>> x = offsetvector({"H1": 0, "L1": 10, "V1": 20})
	>>> x["H1"]
	0

	The motivation for introducing this class, instead of using
	dictionaries, is that it provides a number of tools for comparing
	offset vectors besides strict value-for-value equality.  For
	example the Python cmp() operation compares two offset vectors by
	the relative offsets between instruments rather than their absolute
	offsets.  There is also the ability to check if one offset vector
	is a subset of another one.
	"""

	@property
	def refkey(self):
		"""
		min(self.keys())

		Raises ValueError if the offsetvector is empty.
		"""
		if not self:
			raise ValueError, "offsetvector is empty"
		return min(self)

	@property
	def deltas(self):
		"""
		Dictionary of relative offsets.  The keys in the result are
		pairs of keys from the offset vector, (a, b), and the
		values are the relative offsets, (offset[b] - offset[a]).

		Example:

		>>> x = offsetvector({"H1": 0, "L1": 10, "V1": 20})
		>>> x.deltas
		{('H1', 'H1'): 0, ('H1', 'L1'): 10, ('H1', 'V1'): 20}
		>>> y = offsetvector({'H1': 100, 'L1': 110, 'V1': 120})
		>>> y.deltas == x.deltas
		True

		Note that the result always includes a "dummy" entry,
		giving the relative offset of self.refkey with respect to
		itself, which is always 0.

		See also .fromdeltas().

		BUGS:  I think the keys in each tuple should be reversed.
		I can't remember why I put them in the way they are.
		Expect them to change in the future.
		"""
		# NOTE:  the arithmetic used to construct the offsets
		# *must* match the arithmetic used by
		# time_slide_component_vectors() so that the results of the
		# two functions can be compared to each other without worry
		# of floating-point round off confusing things.
		refkey = self.refkey
		refoffset = self[refkey]
		return dict(((refkey, key), self[key] - refoffset) for key in self)

	def __str__(self, compact = False):
		"""
		Return a human-readable string representation of an offset
		vector.

		Example:

		>>> a = offsetvector({"H1": -10.1234567, "L1": 0.1})
		>>> str(a)
		'H1 = -10.1234567 s, L1 = +0.1 s'
		>>> a.__str__(compact = True)
		'H1=-10.123,L1=0.1'
		"""
		if compact:
			return ",".join(("%s=%.5g" % x) for x in sorted(self.items()))
		return ", ".join(("%s = %+.16g s" % x) for x in sorted(self.items()))

	def __repr__(self):
		"""
		Return a string representation of the offset vector.

		Example:

		>>> a = offsetvector({"H1": -10.1234567, "L1": 0.1})
		>>> repr(a)
		"offsetvector({'H1': -10.1234567, 'L1': 0.10000000000000001})"
		>>> b = eval(repr(a))
		>>> b
		offsetvector({'H1': -10.1234567, 'L1': 0.10000000000000001})
		>>> b == a
		True
		>>> b is a
		False
		"""
		return "%s(%s)" % (self.__class__.__name__, dict.__repr__(self))

	def __cmp__(self, other):
		"""
		Compare two offset vectors by their relative offsets.  The
		return value is 0 if the relative offsets are all equal,
		nonzero otherwise.

		Example:

		>>> a = offsetvector({"H1": 0.0, "H2": 0.0, "L1": 0.0})
		>>> b = offsetvector({"H1": 10.0, "H2": 10.0, "L1": 10.0})
		>>> cmp(a, b)
		0
		>>> a == b
		False

		Note that cmp() and testing for equality are different
		tests!  The equality test returns False because the offset
		vectors are not identical, however the cmp() function
		returns 0 because the relative offsets are all equal.
		"""
		return cmp(self.deltas, other.deltas)

	def contains(self, other):
		"""
		Returns True if offset vector @other can be found in @self,
		False otherwise.  An offset vector is "found in" another
		offset vector if the latter contains all of the former's
		instruments and the relative offsets among those
		instruments are equal (the absolute offsets need not be).

		Example:

		>>> a = offsetvector({"H1": 10, "L1": 20, "V1": 30})
		>>> b = offsetvector({"H1": 20, "V1": 40})
		>>> a.contains(b)
		True

		Note the distinction between this and the "in" operator:

		>>> b in a
		False
		>>> "H1" in a
		True
		"""
		return offsetvector((key, offset) for key, offset in self.items() if key in other).deltas == other.deltas

	def normalize(self, **kwargs):
		"""
		Adjust the offsetvector so that a particular instrument has
		the desired offset.  All other instruments have their
		offsets adjusted so that the relative offsets are
		preserved.  The instrument to noramlize, and the offset one
		wishes it to have, are provided as a key-word argument.
		The return value is the time slide dictionary, which is
		modified in place.

		If more than one key-word argument is provided the keys are
		sorted and considered in order until a key is found that is
		in the offset vector.  The offset vector is normalized to
		that value.  This function is a no-op if no key-word
		argument is found that applies.

		Example:

		>>> a = offsetvector({"H1": -10, "H2": -10, "L1": -10})
		>>> a.normalize(L1 = 0)
		{'H2': 0, 'H1': 0, 'L1': 0}
		>>> a = offsetvector({"H1": -10, "H2": -10})
		>>> a.normalize(L1 = 0, H2 = 5)
		{'H2': 5, 'H1': 5}
		"""
		# FIXME:  should it be performed in place?  if it should
		# be, the should there be no return value?
		for key, offset in sorted(kwargs.items()):
			if key in self:
				delta = offset - self[key]
				for key in self.keys():
					self[key] += delta
				break
		return self

	@classmethod
	def fromdeltas(cls, deltas):
		"""
		Construct an offsetvector from a dictionary of offset
		deltas as returned by the .deltas attribute.

		Example:

		>>> x = offsetvector({"H1": 0, "L1": 10, "V1": 20})
		>>> y = offsetvector.fromdeltas(x.deltas)
		>>> y
		offsetvector({'V1': 20, 'H1': 0, 'L1': 10})

		See also .deltas, .fromkeys()
		"""
		return cls((key, value) for (refkey, key), value in deltas.items())
