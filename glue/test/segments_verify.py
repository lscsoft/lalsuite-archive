from glue.segments import *
import unittest

class TestInfinity(unittest.TestCase):
	def testcmp(self):
		# __cmp__()
		a = infinity()
		self.assertEqual( 0, cmp(-a, -a))
		self.assertEqual(-1, cmp(-a,  0))
		self.assertEqual(-1, cmp(-a,  a))
		self.assertEqual( 1, cmp( 0, -a))
		self.assertEqual(-1, cmp( 0,  a))
		self.assertEqual( 1, cmp( a, -a))
		self.assertEqual( 1, cmp( a,  0))
		self.assertEqual( 0, cmp( a,  a))

	def testadd(self):
		# __add__()
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

	def testsub(self):
		# __sub__()
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


class TestSegment(unittest.TestCase):
	def testdefinition(self):
		# segment creation
		self.assertEqual((-2, 2), tuple(segment(-2, 2)))
		self.assertEqual((-2, 2), tuple(segment(2, -2)))
		self.assertEqual((-infinity(), 2), tuple(segment(-infinity(), 2)))
		self.assertEqual((-infinity(), 2), tuple(segment(2, -infinity())))
		self.assertEqual((2, infinity()), tuple(segment(infinity(), 2)))
		self.assertEqual((2, infinity()), tuple(segment(2, infinity())))
		self.assertEqual((-infinity(), infinity()), tuple(segment(-infinity(), infinity())))

	def testduration(self):
		# duration()
		self.assertEqual(4, segment(-2, 2).duration())
		self.assertEqual(infinity(), segment(0, infinity()).duration())
		self.assertEqual(infinity(), segment(-infinity(), 0).duration())
		self.assertEqual(infinity(), segment(-infinity(), infinity()).duration())

	def testintersects(self):
		self.assertEqual(False, segment(-2, 2).intersects(segment(-4, -3)))
		self.assertEqual(False, segment(-2, 2).intersects(segment(-4, -2)))
		self.assertEqual( True, segment(-2, 2).intersects(segment(-4,  0)))
		self.assertEqual( True, segment(-2, 2).intersects(segment(-4,  2)))
		self.assertEqual( True, segment(-2, 2).intersects(segment(-4,  4)))
		self.assertEqual( True, segment(-2, 2).intersects(segment(-2,  4)))
		self.assertEqual( True, segment(-2, 2).intersects(segment( 0,  4)))
		self.assertEqual(False, segment(-2, 2).intersects(segment( 2,  4)))
		self.assertEqual(False, segment(-2, 2).intersects(segment( 3,  4)))
		self.assertEqual( True, segment(-2, 2).intersects(segment(-2,  2)))
		self.assertEqual( True, segment(-2, 2).intersects(segment(-1,  1)))

	def testcontains(self):
		self.assertEqual(False, segment(-4, -3) in segment(-2, 2))
		self.assertEqual(False, segment(-4, -2) in segment(-2, 2))
		self.assertEqual(False, segment(-4, 0) in segment(-2, 2))
		self.assertEqual(False, segment(-4, 2) in segment(-2, 2))
		self.assertEqual(False, segment(-4, 4) in segment(-2, 2))
		self.assertEqual(False, segment(-2, 4) in segment(-2, 2))
		self.assertEqual(False, segment(0, 4) in segment(-2, 2))
		self.assertEqual(False, segment(2, 4) in segment(-2, 2))
		self.assertEqual(False, segment(3, 4) in segment(-2, 2))
		self.assertEqual( True, segment(-2, 2) in segment(-2, 2))
		self.assertEqual( True, segment(-1, 1) in segment(-2, 2))

	def testcontinuous(self):
		self.assertEqual(False, segment(-2, 2).continuous(segment(-4, -3)))
		self.assertEqual( True, segment(-2, 2).continuous(segment(-4, -2)))
		self.assertEqual( True, segment(-2, 2).continuous(segment(-4,  0)))
		self.assertEqual( True, segment(-2, 2).continuous(segment(-4,  2)))
		self.assertEqual( True, segment(-2, 2).continuous(segment(-4,  4)))
		self.assertEqual( True, segment(-2, 2).continuous(segment(-2,  4)))
		self.assertEqual( True, segment(-2, 2).continuous(segment( 0,  4)))
		self.assertEqual( True, segment(-2, 2).continuous(segment( 2,  4)))
		self.assertEqual(False, segment(-2, 2).continuous(segment( 3,  4)))
		self.assertEqual( True, segment(-2, 2).continuous(segment(-2,  2)))
		self.assertEqual( True, segment(-2, 2).continuous(segment(-1,  1)))

	def testcontraction(self):
		self.assertEqual(segment(5, 15), segment(0, 20).contract(5))


class TestSegmentList(unittest.TestCase):
	def testsub(self):
		# __sub__()
		self.assertEqual(segmentlist([]), segmentlist([]) - segmentlist([]))
		self.assertEqual(segmentlist([]), segmentlist([]) - segmentlist([segment(-1,1)]))
		self.assertEqual(segmentlist([segment(-1,1)]) - segmentlist([segment(-1,1)]), segmentlist([]))
		self.assertEqual(segmentlist([]), segmentlist([segment(-1,1)]) - segmentlist([segment(-1,1)]))

		self.assertEqual(segmentlist([segment(0,1)]), segmentlist([segment(0,1)]) - segmentlist([segment(2,3)]))
		self.assertEqual(segmentlist([segment(0,1)]), segmentlist([segment(0,1)]) - segmentlist([segment(2,3), segment(4,5)]))
		self.assertEqual(segmentlist([segment(0,1)]), segmentlist([segment(0,1), segment(2,3)]) - segmentlist([segment(2,3)]))
		self.assertEqual(segmentlist([segment(2,3)]), segmentlist([segment(0,1), segment(2,3)]) - segmentlist([segment(0,1)]))
		self.assertEqual(segmentlist([segment(0,1), segment(4,5)]), segmentlist([segment(0,1), segment(2,3), segment(4,5)]) - segmentlist([segment(2,3)]))

		self.assertEqual(segmentlist([segment(-5, 10)]), segmentlist([segment(-10,10)]) - segmentlist([segment(-15,-5)]))
		self.assertEqual(segmentlist([segment(-10, -5), segment(5, 10)]), segmentlist([segment(-10,10)]) - segmentlist([segment(-5,5)]))
		self.assertEqual(segmentlist([segment(-10, 5)]), segmentlist([segment(-10,10)]) - segmentlist([segment(5,15)]))

		self.assertEqual(segmentlist([segment(0,5), segment(45,50)]), segmentlist([segment(0,10), segment(20,30), segment(40,50)]) - segmentlist([segment(5, 45)]))
		self.assertEqual(segmentlist([segment(-5,5)]), segmentlist([segment(-5,5)]) - segmentlist([segment(0,0)]))

	def testinvert(self):
		# __invert__()
		self.assertEqual(segmentlist([segment(-infinity(), infinity())]), ~segmentlist([]))
		self.assertEqual(segmentlist([]), ~segmentlist([segment(-infinity(), infinity())]))
		self.assertEqual(segmentlist([segment(-infinity(), -5), segment(5, infinity())]), ~segmentlist([segment(-5,5)]))

	def testcontract(self):
		# contract()
		self.assertEqual(segmentlist([segment(0, 20)]), segmentlist([segment(3, 7), segment(13, 17)]).contract(-3))


if __name__ == "__main__":
	unittest.main()
