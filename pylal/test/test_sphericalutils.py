#!/usr/bin/env python

import math
import random
import unittest

import numpy as np
np.seterr(all="raise")

from pylal import sphericalutils as su
from pylal.xlal.constants import *

n_tests = 1000

#
# Utility functions
#

def random_vector():
    return math.acos(2 * random.uniform(0, 1) - 1), \
        random.uniform(0, 2 * LAL_PI)

def random_euler_angle():
    return random.uniform(0, LAL_2_PI), random.uniform(0, LAL_PI), \
        random.uniform(0, 2 * LAL_PI)

#
# Unit tests
#

class test_sphericalutils(unittest.TestCase):
    def test_identity(self):
        """
        Test that Euler angles of zero are the identity operation.
        """
        for i in xrange(n_tests):
            vector = np.array([random_vector()])
            err = su.angle_between_points(vector,
                su.rotate_euler(vector, 0., 0., 0.))
            self.assertAlmostEqual(err, 0.)

    def test_inverse(self):
        """
        Test that rotation * inverse = identity.
        """
        for i in xrange(n_tests):
            vector = np.array([random_vector()])
            a, b, c = random_euler_angle()
            err = su.angle_between_points(vector, \
                su.rotate_euler(su.rotate_euler(vector, a, b, c), -c, -b, -a))
            self.assertAlmostEqual(err, 0.)

    def test_angle_preserved(self):
        """
        Test that the angle between two vectors is preserved under rotation.
        """
        for i in xrange(n_tests):
            vectors_0 = np.array((random_vector(), random_vector()),
                dtype=float)
            d_0 = su.angle_between_points(*vectors_0)

            a, b, c = random_euler_angle()
            vectors_1 = su.rotate_euler(vectors_0, a, b, c)
            d_1 = su.angle_between_points(*vectors_1)
            self.assertAlmostEqual(d_0, d_1)

    def test_north_pole(self):
        """
        Test that new_z_to_euler() gives Euler angles that will rotate the
        original z axis to the requested z axis.
        """
        north_pole = np.array([[0., 0.]])
        for i in xrange(n_tests):
            new_dir = random_vector()
            a, b = su.new_z_to_euler(new_dir)
            new_north = su.rotate_euler(north_pole, a, b, 0)
            self.assertAlmostEqual(0., \
                su.angle_between_points(new_north, np.array([new_dir])))

# construct and run the test suite
suite = unittest.TestSuite()
suite.addTest(unittest.makeSuite(test_sphericalutils))
unittest.TextTestRunner(verbosity=2).run(suite)
