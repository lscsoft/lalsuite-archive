"""
This module contains bits and pieces of use when interacting with LAL and
LAL-derived code (eg. LALApps programs)
"""

__author__ = "Kipp Cannon <kipp@gravity.phys.uwm.edu>"
__date__ = "$Date$"
__version__ = "$Revision$"


#
# Preamble
#

import re

from glue import segments


#
# High precision time object
#

class LIGOTimeGPS(object):
	"""
	An object for storing times with nanosecond resolution.  LAL defines an
	equivalent object which is used through-out the search algorithms to
	represent times.  Many LALApps routines input and output times in a
	manner that meshes well with this object.

	Internally the time is represented as a signed integer "seconds" part
	and an unsigned integer "nanoseconds" part.  The actual time is always
	constructed by adding the nanoseconds to the seconds.  So -0.5 s is
	represented by setting seconds = -1, and nanoseconds to 500000000.
	That's the way LAL does it.
	"""

	# basic class methods

	def __atoparts(self, s):
		"""
		Internal routine for ASCII string conversion.  Users should
		call LIGOTimeGPS() to perform this function.
		"""
		parts = (s.strip() + ".0").split(".")[:2]
		parts[1] = parts[1][:9]
		parts[1] += "000000000"[:9-len(parts[1])]
		try:
			parts = map(int, parts)
		except:
			raise TypeError, "Cannot convert \"%s\" to a LIGOTimeGPS" % s
		if s[0] == "-":
			parts[0] -= 1
			parts[1] = 1000000000 - parts[1]
		return parts

	def __init__(self, seconds, nanoseconds = 0):
		"""
		Create a LIGOTimeGPS instance.  The first parameter is the
		count of seconds, and the second (optional) parameter is the
		count of nanoseconds.  If the nanoseconds parameter is not
		supplied, it is assumed to be 0.  Either parameter can be
		a numeric type or an ASCII string.

		Example use (all are equivalent):
			LIGOTimeGPS(100.5)
			LIGOTimeGPS("100.5")
			LIGOTimeGPS(100, 500000000)
			LIGOTimeGPS(0, 100500000000L)
			LIGOTimeGPS(100.2, 300000000)
		"""
		if type(nanoseconds) == str:
			nanoseconds = float(nanoseconds)
		elif not type(nanoseconds) in [float, int, long]:
			raise TypeError, "Cannot convert \"%s\" to LIGOTimeGPS" % repr(seconds)
		if type(seconds) == LIGOTimeGPS:
			nanoseconds += seconds.nanoseconds
			seconds = seconds.seconds
		elif type(seconds) == str:
			[seconds, self.nanoseconds] = self.__atoparts(seconds)
			nanoseconds += self.nanoseconds
		elif not type(seconds) in [float, int, long]:
			raise TypeError, "Cannot convert \"%s\" to LIGOTimeGPS" % repr(seconds)
		frac_seconds = round(seconds % 1 * 1000000000)
		if seconds < 0 and frac_seconds:
			seconds -= 1
		(self.seconds, self.nanoseconds) = map(int, divmod(frac_seconds + nanoseconds, 1000000000))
		self.seconds += int(seconds)

	def __repr__(self):
		return "LIGOTimeGPS(%d, %u)" % (self.seconds, self.nanoseconds)
	
	def __str__(self):
		"""
		Return an ASCII string representation of a LIGOTimeGPS.
		"""
		if (self.seconds >= 0) or (self.nanoseconds == 0):
			s = "%d.%09u" % (self.seconds, self.nanoseconds)
		elif (self.seconds < -1):
			s = "%d.%09u" % (self.seconds + 1, 1000000000 - self.nanoseconds)
		else:
			s = "-0.%09u" % (1000000000 - self.nanoseconds)
		return s.rstrip("0").rstrip(".")

	# type conversion

	def __float__(self):
		"""
		Convert a LIGOTimeGPS to seconds as a float.

		Example use:
			float(LIGOTimeGPS(100.5))
		"""
		return self.seconds + self.nanoseconds * 1e-9

	def __int__(self):
		"""
		Return the integer part (seconds) of a LIGOTimeGPS as an int.

		Example use:
			int(LIGOTimeGPS(100.5))
		"""
		return self.seconds

	def __long__(self):
		"""
		Return the integer part (seconds) of a LIGOTimeGPS as a long.

		Example use:
			long(LIGOTimeGPS(100.5))
		"""
		return long(self.seconds)

	def ns(self):
		"""
		Convert a LIGOTimeGPS to a count of nanoseconds as a long.

		Example use:
			LIGOTimeGPS(100.5).ns()
		"""
		return self.seconds * 1000000000L + self.nanoseconds

	# comparison

	def __cmp__(self, other):
		"""
		Compare a value to a LIGOTimeGPS.  If the value being compared
		to the LIGOTimeGPS is not also a LIGOTimeGPS, then an attempt
		is made to convert it to a LIGOTimeGPS.

		Example use (all are equivalent):
			LIGOTimeGPS(100.5) < LIGOTimeGPS(200)
			LIGOTimeGPS(100.5) < 200
			LIGOTimeGPS(100.5) < "200"
		"""
		if not type(other) == LIGOTimeGPS:
			try:
				other = LIGOTimeGPS(other)
			except:
				# if other can't be converted, then the two
				# args aren't equal.
				return 1
		return cmp((self.seconds, self.nanoseconds), (other.seconds, other.nanoseconds))

	def __nonzero__(self):
		"""
		Return True if the LIGOTimeGPS is nonzero.

		Example use:
			bool(LIGOTimeGPS(100.5))
		"""
		return self.seconds or self.nanoseconds

	# arithmetic

	def __add__(self, other):
		"""
		Add a value to a LIGOTimeGPS.  If the value being added to the
		LIGOTimeGPS is not also a LIGOTimeGPS, then an attempt is made
		to convert it to a LIGOTimeGPS.

		Example use (all are equivalent):
			LIGOTimeGPS(100.5) + LIGOTimeGPS(3)
			LIGOTimeGPS(100.5) + 3
			LIGOTimeGPS(100.5) + "3"
		"""
		if not type(other) == LIGOTimeGPS:
			other = LIGOTimeGPS(other)
		return LIGOTimeGPS(self.seconds + other.seconds, self.nanoseconds + other.nanoseconds)

	# addition is commutative.
	__radd__ = __add__

	def __sub__(self, other):
		"""
		Subtract a value from a LIGOTimeGPS.  If the value being
		subtracted from the LIGOTimeGPS is not also a LIGOTimeGPS, then
		an attempt is made to convert it to a LIGOTimeGPS.

		Example use (all are equivalent):
			LIGOTimeGPS(100.5) - LIGOTimeGPS(3)
			LIGOTimeGPS(100.5) - 3
			LIGOTimeGPS(100.5) - "3"
		"""
		if not type(other) == LIGOTimeGPS:
			other = LIGOTimeGPS(other)
		return LIGOTimeGPS(self.seconds - other.seconds, self.nanoseconds - other.nanoseconds)

	def __rsub__(self, other):
		"""
		Subtract a LIGOTimeGPS from a value.
		"""
		if not type(other) == LIGOTimeGPS:
			other = LIGOTimeGPS(other)
		return LIGOTimeGPS(other.seconds - self.seconds, other.nanoseconds - self.nanoseconds)

	def __mul__(self, other):
		"""
		Multiply a LIGOTimeGPS by a number.

		Example use:
			LIGOTimeGPS(100.5) * 2
		"""
		return LIGOTimeGPS(self.seconds * other, self.nanoseconds * other)

	# multiplication is commutative
	__rmul__ = __mul__

	def __div__(self, other):
		"""
		Divide a LIGOTimeGPS by a number.

		Example use:
			LIGOTimeGPS(100.5) / 2
		"""
		return LIGOTimeGPS(0, self.ns() / other)

	def __mod__(self, other):
		"""
		Compute the remainder when a LIGOTimeGPS is divided by a number.

		Example use:
			LIGOTimeGPS(100.5) % 3
		"""
		return LIGOTimeGPS(0, self.ns() % (other * 1000000000L))

	# unary arithmetic

	def __pos__(self):
		return self

	def __neg__(self):
		return LIGOTimeGPS(0, -self.ns())

	def __abs__(self):
		return LIGOTimeGPS(0, abs(self.ns()))


#
# LAL cache file manipulation
#

class CacheEntry(object):
	"""
	An object representing one line in a LAL cache file.
	"""
	# How to parse a line in a LAL cache file.  Five white-space
	# delimited columns.
	_regex = re.compile(r"\A\s*(?P<observatory>\S+)\s+(?P<description>\S+)\s+(?P<start>\S+)\s+(?P<duration>\S+)\s+(?P<url>\S+)\s*\Z")

	def __init__(self, line = None, coltype = LIGOTimeGPS):
		"""
		Create a CacheEntry object, parsing an optional string
		argument (one line from a LAL cache file).
		"""
		if line != None:
			match = self._regex.search(line)
			if not match:
				raise ValueError, "could not convert \"%s\" to CacheEntry" % args[0]
			self.observatory = match.group("observatory")
			self.description = match.group("description")
			self.segment = segments.segment(coltype(match.group("start")), coltype(match.group("start")) + coltype(match.group("duration")))
			self.url = match.group("url")

	def __str__(self):
		"""
		String representation is the line in the cache file.
		"""
		return "%s %s %s %s %s" % (self.observatory, self.description, self.segment[0], self.segment.duration(), self.url)

	def __cmp__(self, other):
		"""
		Sort by observatory, then description, then segment, then
		URL.
		"""
		if type(other)  != CacheEntry:
			raise TypeError, "can only compare CacheEntry to CacheEntry"
		return cmp((self.observatory, self.description, self.segment, self.url), (other.observatory, other.description, other.segment, other.url))
