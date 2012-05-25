#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  test_lalinference_proposal.py
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

from pylal.ctypes.lalinference import LALInferenceProposalFunction

__author__="Ben Aylott <benjamin.aylott@ligo.org>"
__version__= "git id %s"%git_version.id
__date__= git_version.date

class test_LALInferenceProposalFunction(unittest.TestCase):
    def setUp(self):
        pass

    def test_init(self):
        self.li_prop_func=LALInferenceProposalFunction("Default")
    
    def test_get(self):
        self.li_prop_func=LALInferenceProposalFunction("TaylorF2")
        print self.li_prop_func._funcName
        self.li_prop_func._funcName="Default"
        print self.li_prop_func._funcName
        
    def test_set(self):
        self.li_prop_func=LALInferenceProposalFunction()
        
        print self.li_prop_func._funcName
        
        def test_func(o,x):
            o=x
                
        self.li_prop_func._funcName=None
        print self.li_prop_func._funcName

suite = unittest.TestSuite()
suite.addTest(unittest.makeSuite(test_LALInferenceProposalFunction))

unittest.TextTestRunner(verbosity=2).run(suite)
