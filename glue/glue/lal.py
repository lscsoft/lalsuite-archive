"""
This module contains bits and pieces of use when interacting with LAL and
LAL-derived code (eg. LALApps programs)
"""

__author__ = "Kipp Cannon <kipp@gravity.phys.uwm.edu>"
__date__ = "$Date$"
__version__ = "$Revision$"

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

	__slots__ = ["seconds", "nanoseconds"]

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
		supplied, it is assumed to be 0.  The seconds parameter can be
		a numeric type or an ASCII string.

		Example use (all are equivalent):
			LIGOTimeGPS(100.5)
			LIGOTimeGPS("100.5")
			LIGOTimeGPS(100, 500000000)
			LIGOTimeGPS(0, 100500000000L)
			LIGOTimeGPS(100.2, 300000000)

		Bugs:
			If the time in seconds is provided as an ASCII string,
			then the nanoseconds parameter is ignored.
		"""
		if type(seconds) == str:
			[self.seconds, self.nanoseconds] = self.__atoparts(seconds)
			return
		if not type(seconds) in [float, int, long]:
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
			return "%d.%09u" % (self.seconds, self.nanoseconds)
		elif (self.seconds < -1):
			return "%d.%09u" % (self.seconds + 1, 1000000000 - self.nanoseconds)
		else:
			return "-0.%09u" % (1000000000 - self.nanoseconds)

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
		Convert a LIGOTimeGPS to a cound of nanoseconds as a long.

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
			other = LIGOTimeGPS(other)
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

		Example use (all are equivalent:
			LIGOTimeGPS(100.5) + LIGOTimeGPS(3)
			LIGOTimeGPS(100.5) + 3
			LIGOTimeGPS(100.5) + "3"
		"""
		if not type(other) == LIGOTimeGPS:
			other = LIGOTimeGPS(other)
		return LIGOTimeGPS(self.seconds + other.seconds, self.nanoseconds + other.nanoseconds)

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

	def __mul__(self, other):
		"""
		Multiply a LIGOTimeGPS by a number.

		Example use:
			LIGOTimeGPS(100.5) * 2
		"""
		return LIGOTimeGPS(self.seconds * other, self.nanoseconds * other)

	def __div__(self, other):
		"""
		Divide a LIGOTimeGPS by a number.

		Example use:
			LIGOTimeGPS(100.5) / 2
		"""
		return LIGOTimeGPS(self.seconds / other, ((self.seconds % other * 1000000000) + self.nanoseconds) / other)
