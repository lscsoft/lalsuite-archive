#!/usr/bin/env python

import unittest
import numpy

from pylal import Fr

class test_Fr(unittest.TestCase):
    def test_default_roundtrip(self):
        """ test call with default values """
        a = Fr.frgetvect1d("./test.dat","Adc1")
        Fr.frputvect('writetest.gwf', [{'name':'Adc1', 'data':a[0],
            'start':a[1], 'dx':a[3], 'kind':'ADC', 'x_unit':a[4],
            'y_unit':a[5]}])
        b = Fr.frgetvect1d("writetest.gwf", "Adc1")
        self.assert_(numpy.alltrue(a[0] == b[0]))
        self.assert_(numpy.alltrue(a[1:] == b[1:]))

    def test_keywords_roundtrip(self):
        """ test call with keyword arguments """
        a = Fr.frgetvect1d("./test.dat", "Adc1", span=1)
        Fr.frputvect('writetest.gwf', [{'name':'Adc1', 'data':a[0],
        'start':a[1], 'dx':a[3], 'kind':'ADC', 'x_unit':a[4],
        'y_unit':a[5]}])
        b = Fr.frgetvect1d("writetest.gwf", "Adc1")
        self.assert_(numpy.alltrue(a[0] == b[0]))
        self.assert_(numpy.alltrue(a[1:] == b[1:]))

    def test_two_channels_roundtrip(self):
        a = Fr.frgetvect1d("./test.dat","Adc1")
        Fr.frputvect('writetest.gwf', [{'name':'Adc1', 'data':a[0],
        'start':a[1], 'dx':a[3], 'kind':'ADC', 'x_unit':a[4],
        'y_unit':a[5]},{'name':'reverse', 'data':a[0][::-1], 'start':a[1],
        'dx':a[3], 'kind':'ADC', 'x_unit':a[4], 'y_unit': a[5]}])
        b = Fr.frgetvect1d("writetest.gwf", "Adc1")
        self.assert_(numpy.alltrue(a[0] == b[0]))
        self.assert_(numpy.alltrue(a[1:] == b[1:]))
        
        c = Fr.frgetvect1d("writetest.gwf", "reverse")
        self.assert_(numpy.alltrue(a[0][::-1] == c[0]))
        self.assert_(numpy.alltrue(a[1:] == c[1:]))

# construct and run the test suite.
suite = unittest.TestSuite()
suite.addTest(unittest.makeSuite(test_Fr))
unittest.TextTestRunner(verbosity=2).run(suite)
