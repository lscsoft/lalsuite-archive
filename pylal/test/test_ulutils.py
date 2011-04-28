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

    def test_volume_lambda(self):
        '''
        Check the dependence of upper limits on volume and lambda.
        '''
        # volumes to test
        volumes = numpy.linspace(1e-3,1,100)

        for vol in volumes:
            # lambda = 0
            mu, post = upper_limit_utils.compute_posterior(vol, 0, 0)
            muhi = upper_limit_utils.compute_upper_limit(mu, post, alpha = 0.90)
            self.assertTrue( 2.30/vol < muhi < 2.31/vol )

            # lambda = 1
            mu, post = upper_limit_utils.compute_posterior(vol, 0, 1)
            muhi = upper_limit_utils.compute_upper_limit(mu, post, alpha = 0.90)
            self.assertTrue( 3.27/vol < muhi < 3.28/vol )

            # lambda = infinity
            mu, post = upper_limit_utils.compute_posterior(vol, 0, 1e6)
            muhi = upper_limit_utils.compute_upper_limit(mu, post, alpha = 0.90)
            self.assertTrue( 3.88/vol < muhi < 3.90/vol )

    def test_prior(self):
        '''
        Check that the code handles priors correctly.
        '''
        mu_in = numpy.linspace(0,10,1e6)
        prior = numpy.ones(mu_in.shape)

        # with a flat prior explicitly given
        mu, post = upper_limit_utils.compute_posterior(1, 0, 0, mu_in, prior)
        muhi_p = upper_limit_utils.compute_upper_limit(mu, post, alpha = 0.90)

        # no prior specified
        mu, post = upper_limit_utils.compute_posterior(1, 0, 0)
        muhi_np = upper_limit_utils.compute_upper_limit(mu, post, alpha = 0.90)
        self.assertTrue( abs(muhi_p-muhi_np) < 0.01 )

        # zero volume check
        mu, post = upper_limit_utils.compute_posterior(0, 0, 0, mu, post)
        muhi_zv = upper_limit_utils.compute_upper_limit(mu, post, alpha = 0.90)
        self.assertTrue( abs(muhi_zv - muhi_np) < 0.01 )

        # volume additivity check
        mu, post = upper_limit_utils.compute_posterior(1, 0, 0)
        mu, post = upper_limit_utils.compute_posterior(9, 0, 0, mu, post)
        muhi_1p9 = upper_limit_utils.compute_upper_limit(mu, post, alpha = 0.90)

        mu, post = upper_limit_utils.compute_posterior(10, 0, 0)
        muhi_10 = upper_limit_utils.compute_upper_limit(mu, post, alpha = 0.90)
        self.assertTrue( abs( muhi_1p9 - muhi_10 ) < 0.01 )

# construct and run the test suite.
suite = unittest.TestSuite()
suite.addTest(unittest.makeSuite(test_ulutils))
unittest.TextTestRunner(verbosity=2).run(suite)
