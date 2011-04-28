#!/usr/bin/env python

import unittest
import numpy

from pylal import upper_limit_utils

class test_ulutils(unittest.TestCase):
    def test_upper_limit(self):
        '''
        Give the upper_limit function some known distributions with known 95% upper limits.
        '''
        # Gaussian
        mu = numpy.linspace(-5,5,1e6)
        post = numpy.exp(-(mu**2)/2)
        muhi = upper_limit_utils.compute_upper_limit(mu, post, alpha = 0.95)
        self.assertTrue( 1.6448 < muhi < 1.6449 ) # get upper limit to 4 sig figs

        # Exponential
        mu = numpy.linspace(0,15,1e6)
        post = numpy.exp(-mu)
        muhi = upper_limit_utils.compute_upper_limit(mu, post, alpha = 0.95)
        self.assertTrue( 2.9957 < muhi < 2.9958 ) # get upper limit to 4 sig figs

        # Uniform
        # Hmm... what should be the behavior in this case?

    def test_lamdba(self):
        '''
        Check the additivity property of posterior results.
        '''
        # lambda = 0
        print "lambda = 0"
        mu, post = upper_limit_utils.compute_posterior(1, 0, 0)
        muhi = upper_limit_utils.compute_upper_limit(mu, post, alpha = 0.90)
        print muhi
        #self.assertTrue( 2.303 < muhi < 2.304 )

        # lambda = 1
        print "lambda = 1"
        mu, post = upper_limit_utils.compute_posterior(1, 0, 1)
        muhi = upper_limit_utils.compute_upper_limit(mu, post, alpha = 0.90)
        print muhi
        #self.assertTrue( 3.272 < muhi < 3.273 )

        # lambda = infinity
        print "lambda = inf"
        mu, post = upper_limit_utils.compute_posterior(1, 0, 1e6)
        muhi = upper_limit_utils.compute_upper_limit(mu, post, alpha = 0.90)
        print muhi
        #self.assertTrue( 3.890 < muhi < 3.891 )



# construct and run the test suite.
suite = unittest.TestSuite()
suite.addTest(unittest.makeSuite(test_ulutils))
unittest.TextTestRunner(verbosity=2).run(suite)
