#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  test_lalinference_prior.py
#
#  Copyright 2011 Ben Aylott <benjamin.aylott@ligo.org>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
#

import unittest

from pylal import git_version

from pylal.ctypes.lalinference import LALInferencePrior
from pylal.ctypes.datatypes.primitives import REAL8,INT4,REAL4
from pylal.ctypes.utils import XLAL_Error
from pylal.ctypes.lalinference import LALInferenceVariables
from pylal.ctypes.lalinference import LALINFERENCE_PARAM_FIXED,LALINFERENCE_REAL8_t,LALINFERENCE_REAL4_t,LALINFERENCE_INT4_t,LALINFERENCE_PARAM_LINEAR


__author__="Ben Aylott <benjamin.aylott@ligo.org>"
__version__= "git id %s"%git_version.id
__date__= git_version.date

class test_LALInferencePrior(unittest.TestCase):
    def setUp(self):
        self.li_prior=LALInferencePrior()
        
    def test_add_min_max(self):
        self.li_prior.addMinMaxPrior("mass",0.,10.,LALINFERENCE_REAL8_t)
        self.assertTrue(self.li_prior.checkMinMaxPrior("mass"))
        self.assertFalse(self.li_prior.checkMinMaxPrior("mas"))
        
    def test_get_min_max(self):
        self.li_prior.addMinMaxPrior("mass",0.,10.,LALINFERENCE_REAL8_t)
        self.assertTrue(len(self.li_prior.getMinMaxPrior("mass"))==2)
        self.assertTrue(self.li_prior.getMinMaxPrior("mass")[0]==0.)
        self.assertFalse(self.li_prior.getMinMaxPrior("mass")[1]==0.)

    def test_remove_min_max(self):
        self.li_prior.addMinMaxPrior("mass",0.,10.,LALINFERENCE_REAL8_t)
        self.assertTrue(self.li_prior.checkMinMaxPrior("mass"))
        self.li_prior.removeMinMaxPrior("mass")
        self.assertFalse(self.li_prior.checkMinMaxPrior("mass"))
        self.li_prior.removeMinMaxPrior("mass")

    def test_add_gaussian(self):
        self.li_prior.addGaussianPrior("time",0.,10.,LALINFERENCE_REAL8_t)
        self.assertTrue(self.li_prior.checkGaussianPrior("time"))
        self.assertFalse(self.li_prior.checkGaussianPrior("mass"))
        
    def test_get_gaussian(self):
        self.li_prior.addGaussianPrior("time",0.,10.,LALINFERENCE_REAL8_t)
        self.assertTrue(len(self.li_prior.getGaussianPrior("time"))==2)
        self.assertTrue(self.li_prior.getGaussianPrior("time")[0]==0.)
        self.assertFalse(self.li_prior.getGaussianPrior("time")[1]==0.)
        
    def test_remove_gaussian(self):
        self.li_prior.addGaussianPrior("time",0.,10.,LALINFERENCE_REAL8_t)
        self.assertTrue(self.li_prior.checkGaussianPrior("time"))
        self.li_prior.removeGaussianPrior("time")
        self.assertFalse(self.li_prior.checkGaussianPrior("time"))
        self.li_prior.removeGaussianPrior("time")
        

suite = unittest.TestSuite()
suite.addTest(unittest.makeSuite(test_LALInferencePrior))

unittest.TextTestRunner(verbosity=2).run(suite)
