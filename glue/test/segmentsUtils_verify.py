import StringIO
import unittest

from glue.segments import *
from glue import segmentsUtils


#
# Define the components of the test suite.
#

class test_fromsegwizard(unittest.TestCase):
	def test(self):
		data = StringIO.StringIO("""# This is a comment
 # This is another comment
	# Again a comment
1  10 100 90
2 110 120 10# Here's a comment
3 125 130 5 # Another one

4   0 200 200""")
		correct = segmentlist([segment(10, 100), segment(110, 120), segment(125, 130), segment(0, 200)])
		self.assertEqual(correct, segmentsUtils.fromsegwizard(data, strict=True))


#
# Construct and run the test suite.
#

suite = unittest.TestSuite()
suite.addTest(unittest.makeSuite(test_fromsegwizard))

unittest.TextTestRunner(verbosity=2).run(suite)
