#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  test_lalinference_runstate.py
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

from pylal.ctypes.lalinference import LALInferenceRunState,BaseLALInferencePriorFunction,LALInferenceLikelihoodFunction

__author__="Ben Aylott <benjamin.aylott@ligo.org>"
__version__= "git id %s"%git_version.id
__date__= git_version.date

class test_LALInferenceRunState(unittest.TestCase):
    def setUp(self):
        self.li_run_state=LALInferenceRunState()

    def test_prior(self):
        lirs=self.li_run_state
        self.assertTrue(lirs.prior==None)
        lirs.prior=LALInferencePriorFunction()
    
    def test_likelihood(self):
        lirs=self.li_run_state
        
        self.assertTrue(lirs.likelihood==None)
        lirs.prior=LALInferenceLikelihoodFunction()
        
    

suite = unittest.TestSuite()
suite.addTest(unittest.makeSuite(test_LALInferenceRunState))

unittest.TextTestRunner(verbosity=2).run(suite)
