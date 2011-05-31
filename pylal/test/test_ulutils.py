#!/usr/bin/env python

import unittest
import numpy
from numpy import random

from pylal import upper_limit_utils

class test_ulutils(unittest.TestCase):
    def test_gaussian_upper_limit(self):
        '''
        Give the upper_limit function some known distributions with known 95% upper limits.
        '''
        # Gaussian
        mu = numpy.linspace(-5,5,1e6)
        post = numpy.exp(-(mu**2)/2)
        muhi = upper_limit_utils.compute_upper_limit(mu, post, alpha = 0.95)
        self.assertTrue( 1.6448 < muhi < 1.6449 ) # get upper limit to 4 sig figs

    def test_exponential_upper_limit(self):
        # Exponential
        mu = numpy.linspace(0,15,1e6)
        post = numpy.exp(-mu)
        muhi = upper_limit_utils.compute_upper_limit(mu, post, alpha = 0.95)
        self.assertTrue( 2.9957 < muhi < 2.9958 ) # get upper limit to 4 sig figs

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

    def test_uniform_prior(self):
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

    def test_zero_volume_search(self):
        '''
        Check that running a search with zero volume has
        no effect on upper limits.
        '''
        # no prior specified, unit volume
        mu, post = upper_limit_utils.compute_posterior(1, 0, 0)
        muhi_np = upper_limit_utils.compute_upper_limit(mu, post, alpha = 0.90)

        # posterior for prior, zero volume
        mu, post = upper_limit_utils.compute_posterior(0, 0, 0, mu, post)
        muhi_zv = upper_limit_utils.compute_upper_limit(mu, post, alpha = 0.90)
        self.assertTrue( abs(muhi_zv - muhi_np) < 0.01 )

    def test_volume_additivity(self):
        '''
        For a series of independent searches with 0 lambda,
        the search volumes add. Check that upper limits computed
        by summing volumes is same as upper limit obtained by
        iteratively applying the individual searches.
        '''
        # volume additivity check
        mu, post = upper_limit_utils.compute_posterior(1, 0, 0)
        mu, post = upper_limit_utils.compute_posterior(9, 0, 0, mu, post)
        muhi_1p9 = upper_limit_utils.compute_upper_limit(mu, post, alpha = 0.90)

        mu, post = upper_limit_utils.compute_posterior(10, 0, 0)
        muhi_10 = upper_limit_utils.compute_upper_limit(mu, post, alpha = 0.90)
        self.assertTrue( abs( muhi_1p9 - muhi_10 ) < 0.01 )

    def test_integrate_efficiency_lind(self):
        '''
        Check that the numerical accuracy of the integration
        methods are sufficient.
        '''
        xbins = numpy.linspace(0,1,50)
        centres = (xbins[:-1]+xbins[1:])/2
        mockeff = numpy.exp(-centres)/(4*numpy.pi*centres**2)
        v, verr = upper_limit_utils.integrate_efficiency(xbins, mockeff, logbins=False)
        vexpect = 1 - numpy.exp(-1)
        self.assertTrue(abs(v -vexpect ) < 0.01)

    def test_integrate_efficiency_logd(self):
        '''
        Check that the numerical accuracy of the integration
        methods are sufficient.
        '''
        xbins = numpy.logspace(-2,0,50)
        centres = numpy.exp((numpy.log(xbins[1:])+numpy.log(xbins[:-1]))/2) # log midpoint
        mockeff = numpy.exp(-centres)/(4*numpy.pi*centres**2)
        v, verr = upper_limit_utils.integrate_efficiency(xbins, mockeff, logbins=True)
        vexpect = 1 - numpy.exp(-1)
        self.assertTrue(abs(v -vexpect ) < 0.01)

    def test_integrate_realistic_efficiency(self):
        '''
        Check that the volume calculation gives the expected results in the
        ideal case of an analytically known efficiency curve.
        '''
        Rmax = 150
        def mockeff(r):
            # This efficiency curve has the general features
            # of a real efficiency curve seen empirically.
            # It also integrates nicely and so the results of
            # integrating this efficiency curve numerically can be
            # compared against an analytic reference.
            return (1-r/Rmax)**4

        rbins = numpy.linspace(0,Rmax,20)
        centres = (rbins[1:]+rbins[:-1])/2
        v, verr = upper_limit_utils.integrate_efficiency(rbins, mockeff(centres))
        vexpect = (1./35)*(4*numpy.pi/3)*(Rmax**3)
        self.assertTrue(abs(1-v/vexpect) < 0.01)

    def test_mean_efficiency(self):
        '''
        Check the mean efficiency calculation in a controlled way.
        '''
        # Assume a sigmoidal mock efficiency curve.
        def eff_model(inj, rchar = 25.0, order = 6.0):
            return 1./(1+(inj/rchar)**order)

        # think of these as log-d bins between 1 and 100 Mpc
        bins = numpy.logspace(0, 2, num=50)
        centres = 10**((numpy.log10(bins[1:])+numpy.log10(bins[:-1]))/2) # log midpoint

        # generate a bunch of injections uniform in log-d
        injections = random.uniform(1, 100, size=10000)
        class MiniInj(object):
            def __init__(self, distance):
                self.distance = distance

        found = []
        missed = []
        for inj in injections:
            if random.binomial(1, eff_model(inj) ):
                found.append( MiniInj(inj) )
            else:
                missed.append( MiniInj(inj) )

        eff, err = upper_limit_utils.mean_efficiency(found, missed, bins, bootnum=10)
        # the computed mean efficiency should agree with the efficiency
        # model to at least ~5% (though this can fluctuate)
        self.assertTrue( (eff - eff_model(centres)).sum()/len(eff) < 0.05 )

# construct and run the test suite.
suite = unittest.TestSuite()
suite.addTest(unittest.makeSuite(test_ulutils))
unittest.TextTestRunner(verbosity=2).run(suite)
