import random
import sys
from glue.lal import *
import unittest


#
# Define the components of the test suite.
#

def maxLIGOTimeGPS():
	return LIGOTimeGPS(2**32 - 1, 999999999)

def randomLIGOTimeGPS():
	return LIGOTimeGPS(random.randint(-100000000, +100000000), random.randint(0, 999999999))


class test_LIGOTimeGPS(unittest.TestCase):
	def test__init__(self):
		correct = LIGOTimeGPS(100, 500000000)
		tests = [
			(100.5,),
			(100.500000000,),
			(100.50000000000000000000000,),
			(100, "500000000"),
			(100, "500000000.0000000000000"),
			(101, "-500000000"),
			(101, "-500000000.0000000000000"),
			("100.5",),
			("100.500000000",),
			("100.50000000000000000000000",),
			("100", 500000000),
			("100", 500000000.0000000000000),
			("101", -500000000),
			("101", -500000000.0000000000000),
			("100", "500000000"),
			("100", "500000000.0000000000000"),
			("101", "-500000000"),
			("101", "-500000000.0000000000000"),
			(0, 100500000000),
			(0, 100500000000.0000000000000),
			(99, 1500000000),
			(99.5, 1000000000),
			(-10, 110500000000),
			(-10.5, 111000000000)
		]
		for num, test in enumerate(tests):
			try:
				self.assertEqual(correct, LIGOTimeGPS(*test))
			except AssertionError, e:
				raise AssertionError, "Test %d failed: " % (num) + str(e)

	def test__float__(self):
		self.assertEqual(100.5, float(LIGOTimeGPS(100.5)))

	def test__int__(self):
		self.assertEqual(100, int(LIGOTimeGPS(100.1)))
		self.assertEqual(100, int(LIGOTimeGPS(100.9)))

	def testns(self):
		self.assertEqual(100500000000, LIGOTimeGPS(100.5).ns())

	def test__nonzero__(self):
		self.assertEqual(True, bool(LIGOTimeGPS(100.5)))
		self.assertEqual(False, bool(LIGOTimeGPS(0)))

	def test__add__(self):
		self.assertEqual(LIGOTimeGPS(110.5), LIGOTimeGPS(100.5) + 10)
		self.assertEqual(LIGOTimeGPS(110.5), LIGOTimeGPS(100.5) + LIGOTimeGPS(10))

	def test__mul__(self):
		self.assertEqual(LIGOTimeGPS(10), LIGOTimeGPS(5) * 2)
		self.assertEqual(LIGOTimeGPS(10), LIGOTimeGPS(20) * 0.5)

	def test__div__(self):
		self.assertEqual(LIGOTimeGPS(10), LIGOTimeGPS(20) / 2)
		self.assertEqual(LIGOTimeGPS(10), LIGOTimeGPS(5) / .5)

	def test__mod__(self):
		self.assertEqual(LIGOTimeGPS(3), LIGOTimeGPS(13) % 5.0)

	def test_pylal_comparison(self):
		try:
			from pylal.xlal.date import LIGOTimeGPS as pylalLIGOTimeGPS
		except ImportError:
			print >>sys.stderr, "pylal not available:  skipping test"
			return

		operators = {
			"add": (LIGOTimeGPS.__add__, pylalLIGOTimeGPS.__add__),
			"sub": (LIGOTimeGPS.__sub__, pylalLIGOTimeGPS.__sub__)
		}

		for i in xrange(1000000):
			key = random.choice(operators.keys())
			op, pylalop = operators[key]
			arg1 = randomLIGOTimeGPS() / 2
			arg2 = randomLIGOTimeGPS() / 2
			try:
				self.assertEqual(op(arg1, arg2), pylalop(pylalLIGOTimeGPS(arg1), pylalLIGOTimeGPS(arg2)))
			except AssertionError, s:
				raise AssertionError, "%s(%s, %s) comparison failed: %s" % (key, str(arg1), str(arg2), str(s))

		# FIXME:  mod tests fail, fix then enable
		operators = {
			"mul": (LIGOTimeGPS.__mul__, pylalLIGOTimeGPS.__mul__),
			"div": (LIGOTimeGPS.__div__, pylalLIGOTimeGPS.__div__)#,
			#"mod": (LIGOTimeGPS.__mod__, pylalLIGOTimeGPS.__mod__)
		}

		for i in xrange(1000000):
			key = random.choice(operators.keys())
			op, pylalop = operators[key]
			arg1 = randomLIGOTimeGPS() / 100
			arg2 = 100**(random.random() * 2 - 1)
			try:
				# FIXME:  max allowed discrepancy should be
				# smaller
				self.assertEqual(abs(op(arg1, arg2) - pylalop(pylalLIGOTimeGPS(arg1), arg2)) < 3e-8, True)
			except AssertionError, s:
				raise AssertionError, "%s(%s, %s) comparison failed: %s != %s" % (key, str(arg1), "%.17g" % arg2, str(op(arg1, arg2)), str(pylalop(pylalLIGOTimeGPS(arg1), arg2)))


#
# Construct and run the test suite.
#

suite = unittest.TestSuite()
suite.addTest(unittest.makeSuite(test_LIGOTimeGPS))

unittest.TextTestRunner(verbosity=2).run(suite)
