from glue.segments import *
import random
import unittest


#
# Some useful code.
#

def randomlist(n):
	"""
	Return a coalesced segmentlist of n elements with random boundaries.
	"""
	if n < 1:
		raise ValueError, "randomlist(n): n must be >= 1"
	x = random.random()
	list = segmentlist([segment(x, x + random.random())])
	for i in range(n - 1):
		x = list[-1][1] + random.random()
		list.append(segment(x, x + random.random()))
	return list

def iscoalesced(l):
	"""
	Return True if the segmentlist l is coalesced.
	"""
	for i in range(len(l)-1):
		if l[i][1] >= l[i+1][0]:
			return False
	return True


#
# Define the components of the test suite.
#

class test_infinity(unittest.TestCase):
	def test__cmp__(self):
		a = infinity()
		self.assertEqual( 0, cmp(-a, -a))
		self.assertEqual(-1, cmp(-a,  0))
		self.assertEqual(-1, cmp(-a,  a))
		self.assertEqual( 1, cmp( 0, -a))
		self.assertEqual(-1, cmp( 0,  a))
		self.assertEqual( 1, cmp( a, -a))
		self.assertEqual( 1, cmp( a,  0))
		self.assertEqual( 0, cmp( a,  a))

	def test__add__(self):
		a = infinity()
		b = infinity()
		self.assertEqual( b, (  a) + ( 10))
		self.assertEqual( b, (  a) + (-10))
		self.assertEqual(-b, ( -a) + ( 10))
		self.assertEqual(-b, ( -a) + (-10))
		self.assertEqual( b, ( 10) + (  a))
		self.assertEqual( b, (-10) + (  a))
		self.assertEqual(-b, ( 10) + ( -a))
		self.assertEqual(-b, (-10) + ( -a))
		self.assertEqual( b, (  a) + (  a))
		self.assertEqual(-b, ( -a) + ( -a))

	def test__sub__(self):
		a = infinity()
		b = infinity()
		self.assertEqual( b, (  a) - ( 10))
		self.assertEqual( b, (  a) - (-10))
		self.assertEqual(-b, ( -a) - ( 10))
		self.assertEqual(-b, ( -a) - (-10))
		self.assertEqual(-b, ( 10) - (  a))
		self.assertEqual(-b, (-10) - (  a))
		self.assertEqual( b, ( 10) - ( -a))
		self.assertEqual( b, (-10) - ( -a))
		self.assertEqual(None, ( a) - (  a))
		self.assertEqual(None, (-a) - ( -a))
		self.assertEqual( b, (  a) - ( -a))
		self.assertEqual(-b, ( -a) - (  a))


class test_segment(unittest.TestCase):
	set1 = [
		segment(-2, 2),
		segment(-2, 2),
		segment(-2, 2),
		segment(-2, 2),
		segment(-2, 2),
		segment(-2, 2),
		segment(-2, 2),
		segment(-2, 2),
		segment(-2, 2),
		segment(-2, 2),
		segment(-2, 2),
		segment(-2, 2),
		segment(-2, 2),
		segment(-2, 2)
	]
	set2 = [
		segment(-4, -3),
		segment(-4, -2),
		segment(-4,  0),
		segment(-4,  2),
		segment(-4,  4),
		segment(-2,  4),
		segment( 0,  4),
		segment( 2,  4),
		segment( 3,  4),
		segment(-2,  2),
		segment(-1,  1),
		segment(-infinity(), infinity()),
		segment(0, infinity()),
		segment(-infinity(), 0)
	]

	def test__new__(self):
		self.assertEqual((-2, 2), tuple(segment(-2, 2)))
		self.assertEqual((-2, 2), tuple(segment(2, -2)))
		self.assertEqual((-infinity(), 2), tuple(segment(-infinity(), 2)))
		self.assertEqual((-infinity(), 2), tuple(segment(2, -infinity())))
		self.assertEqual((2, infinity()), tuple(segment(infinity(), 2)))
		self.assertEqual((2, infinity()), tuple(segment(2, infinity())))
		self.assertEqual((-infinity(), infinity()), tuple(segment(-infinity(), infinity())))

	def testduration(self):
		results = [
			1,
			2,
			4,
			6,
			8,
			6,
			4,
			2,
			1,
			4,
			2,
			infinity(),
			infinity(),
			infinity()
		]
		map(lambda i, r, a: self.assertEqual((i, r), (i, a.duration())), range(len(results)), results, self.set2)

	def testintersects(self):
		results = [
			False,
			False,
			True,
			True,
			True,
			True,
			True,
			False,
			False,
			True,
			True,
			True,
			True,
			True
		]
		map(lambda i, r, a, b: self.assertEqual((i, r), (i, a.intersects(b))), range(len(results)), results, self.set1, self.set2)

	def test__contains__(self):
		results = [
			False,
			False,
			False,
			False,
			False,
			False,
			False,
			False,
			False,
			True,
			True,
			False,
			False,
			False
		]
		map(lambda i, r, a, b: self.assertEqual((i, r), (i, a.__contains__(b))), range(len(results)), results, self.set1, self.set2)

	def testcontinuous(self):
		results = [
			False,
			True,
			True,
			True,
			True,
			True,
			True,
			True,
			False,
			True,
			True,
			True,
			True,
			True
		]
		map(lambda i, r, a, b: self.assertEqual((i, r), (i, a.continuous(b))), range(len(results)), results, self.set1, self.set2)

	def testcontract(self):
		results = [
			segment(-5, -2),
			segment(-4, -2),
			segment(-2, -2),
			segment(-2,  0),
			segment(-2,  2),
			segment( 0,  2),
			segment( 2,  2),
			segment( 2,  4),
			segment( 2,  5),
			segment( 0,  0),
			segment(-1,  1),
			segment(-infinity(), infinity()),
			segment(2, infinity()),
			segment(-infinity(), -2)
		]
		map(lambda i, r, a: self.assertEqual((i, r), (i, a.contract(2))), range(len(results)), results, self.set2)


class test_segmentlist(unittest.TestCase):
	def test__sub__(self):
		self.assertEqual(segmentlist([]), segmentlist([]) - segmentlist([]))
		self.assertEqual(segmentlist([]), segmentlist([]) - segmentlist([segment(-1,1)]))
		self.assertEqual(segmentlist([segment(-1,1)]) - segmentlist([segment(-1,1)]), segmentlist([]))
		self.assertEqual(segmentlist([]), segmentlist([segment(-1,1)]) - segmentlist([segment(-1,1)]))

		self.assertEqual(segmentlist([segment(0,1)]), segmentlist([segment(0,1)]) - segmentlist([segment(2,3)]))
		self.assertEqual(segmentlist([segment(0,1)]), segmentlist([segment(0,1)]) - segmentlist([segment(2,3), segment(4,5)]))
		self.assertEqual(segmentlist([segment(0,1)]), segmentlist([segment(0,1), segment(2,3)]) - segmentlist([segment(2,3)]))
		self.assertEqual(segmentlist([segment(2,3)]), segmentlist([segment(0,1), segment(2,3)]) - segmentlist([segment(0,1)]))
		self.assertEqual(segmentlist([segment(0,1), segment(4,5)]), segmentlist([segment(0,1), segment(2,3), segment(4,5)]) - segmentlist([segment(2,3)]))

		self.assertEqual(segmentlist([segment(0,1)]), segmentlist([segment(0,2)]) - segmentlist([segment(1,2)]))
		self.assertEqual(segmentlist([segment(0.8, 0.9), segment(1.0, 1.8)]), segmentlist([segment(0, 2)]) - segmentlist([segment(0, 0.8), segment(0.9, 1.0), segment(1.8, 2)]))

		self.assertEqual(segmentlist([segment(-5, 10)]), segmentlist([segment(-10,10)]) - segmentlist([segment(-15,-5)]))
		self.assertEqual(segmentlist([segment(-10, -5), segment(5, 10)]), segmentlist([segment(-10,10)]) - segmentlist([segment(-5,5)]))
		self.assertEqual(segmentlist([segment(-10, 5)]), segmentlist([segment(-10,10)]) - segmentlist([segment(5,15)]))

		self.assertEqual(segmentlist([segment(0,5), segment(45,50)]), segmentlist([segment(0,10), segment(20,30), segment(40,50)]) - segmentlist([segment(5, 45)]))
		self.assertEqual(segmentlist([segment(-5,5)]), segmentlist([segment(-5,5)]) - segmentlist([segment(0,0)]))

	def test__invert__(self):
		self.assertEqual(segmentlist([segment(-infinity(), infinity())]), ~segmentlist([]))
		self.assertEqual(segmentlist([]), ~segmentlist([segment(-infinity(), infinity())]))
		self.assertEqual(segmentlist([segment(-infinity(), -5), segment(5, infinity())]), ~segmentlist([segment(-5,5)]))

	def test__or__(self):
		for i in range(10000):
			a = randomlist(random.randint(1, 50))
			b = randomlist(random.randint(1, 50))
			c = a | b
			try:
				# make sure c is coalesced
				self.assertEqual(True, iscoalesced(c))
				# make sure c contains all of a
				self.assertEqual(a, c & a)
				# make sure c contains all of b
				self.assertEqual(b, c & b)
				# make sure c contains nothing except a and b
				self.assertEqual(segmentlist([]), c - a - b)
			except AssertionError, e:
				raise AssertionError, str(e) + "\na = " + str(a) + "\nb = " + str(b)

	def testcontract(self):
		self.assertEqual(segmentlist([segment(0, 20)]), segmentlist([segment(3, 7), segment(13, 17)]).contract(-3))


#
# Construct and run the test suite.
#

suite = unittest.TestSuite()
suite.addTest(unittest.makeSuite(test_infinity))
suite.addTest(unittest.makeSuite(test_segment))
suite.addTest(unittest.makeSuite(test_segmentlist))

unittest.TextTestRunner(verbosity=2).run(suite)
