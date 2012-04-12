from ctypes import cast,create_string_buffer,POINTER,c_ssize_t,byref,addressof
import unittest

import numpy as np

from pylal.ctypes.datatypes.vector import REAL8Vector

class test_REAL8Vector(unittest.TestCase):
    
    def test__init__(self):
        
        zero=0.
        length=100
        
        vec=REAL8Vector(length)
        
        for i in range(100):
            self.assertEqual(zero, vec[i])
            
        self.assertRaises(IndexError,vec.__getitem__,length)
        
#
# Construct and run the test suite.
#

suite = unittest.TestSuite()
suite.addTest(unittest.makeSuite(test_REAL8Vector))

unittest.TextTestRunner(verbosity=2).run(suite)
