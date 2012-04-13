#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  test_lalinference.py
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
from pylal.ctypes.datatypes.primitives import REAL8,INT4,REAL4
from pylal.ctypes.lalinference import LALInferenceVariables
from pylal.ctypes.lalinference import LALINFERENCE_PARAM_FIXED,LALINFERENCE_REAL8_t,LALINFERENCE_INT4_t,LALINFERENCE_PARAM_FIXED

__author__="Ben Aylott <benjamin.aylott@ligo.org>"
#__version__= "git id %s"%git_version.id
#__date__= git_version.date

class test_LALInferenceVariables(unittest.TestCase):
    """
    Unit tests for the LALInferenceVariables data structure(s).
    """
    def setUp(self):
        self.livars=LALInferenceVariables()

    def test_add(self):
        self.livars.addVariable("myvar1",REAL8(9.0),LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED)
        self.livars.addVariable("myvar2",REAL8(9.0),LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED)
        self.livars.addVariable("myvar3",REAL8(9.0),LALINFERENCE_INT4_t,LALINFERENCE_PARAM_FIXED)

    def test_rm(self):

        self.livars=LALInferenceVariables()

        self.livars.addVariable("myvar4",REAL8(9.0),LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED)
        self.livars.removeVariable("myvar1")
        
    def test_print(self):
        self.livars.addVariable("myvar4",REAL8(9.0),LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED)
        self.livars.addVariable("myvar5",REAL8(9.0),LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED)
        self.livars.addVariable("myvar6",INT4(2),LALINFERENCE_INT4_t,LALINFERENCE_PARAM_FIXED)

        self.livars.printVariables()
        #self.assertTrue(len(self.livars.printSample().split())==3)


    def test_print_nf(self):
        self.livars.addVariable("myvar1",REAL8(9.0),LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED)
        self.livars.addVariable("myvar2",INT4(9.0),LALINFERENCE_INT4_t,LALINFERENCE_PARAM_FIXED)
        self.livars.addVariable("myvar3",REAL4(9.0),LALINFERENCE_REAL4_t,LALINFERENCE_PARAM_FIXED)
        self.livars.addVariable("myvar4",REAL8(9.0),LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_NON_FIXED)
        self.livars.addVariable("myvar5",INT4(9.0),LALINFERENCE_INT4_t,LALINFERENCE_PARAM_NON_FIXED)
        self.livars.addVariable("myvar6",REAL8(9.0),LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED)

        self.assertFalse(self.livars.printSample()==self.livars.printSampleFixed())

    def test_get(self):
        self.livars.addVariable("myvar1",REAL8(9.0),LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED)
        self.livars.getVariable("myvar1")
        self.assertRaises(KeyError,self.livars.getVariable,"myvar2")

    def test_get_dim(self):
        self.livars.addVariable("myvar1",REAL8(9.0),LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED)
        self.assertTrue(self.livars.getVariableDimension()==1)
        self.livars.addVariable("myvar2",REAL8(9.0),LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED)
        self.livars.addVariable("myvar3",REAL8(9.0),LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_NON_FIXED)
        self.assertTrue(self.livars.getVariableDimension()==3)

    def test_get_dim_nf(self):
        self.livars.addVariable("myvar1",REAL8(9.0),LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED)
        self.assertTrue(self.livars.getVariableDimensionNonFixed()==0)
        self.livars.addVariable("myvar2",REAL8(2.0),LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_NON_FIXED)
        self.livars.addVariable("myvar3",REAL8(9.0),LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED)
        self.assertTrue(self.livars.getVariableDimensionNonFixed()==1)

    def test_get_type(self):
        self.livars.addVariable("myvar1",REAL4(9.0),LALINFERENCE_REAL4_t,LALINFERENCE_PARAM_FIXED)
        self.livars.addVariable("myvar1",REAL8(9.0),LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED)
        self.assertFalse(self.livars.getVariableType("myvar1")==self.livars.getVariableType("myvar2"))

    def test_get_type_by_index(self):
        self.livars.addVariable("myvar1",REAL8(9.0),LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED)
        self.livars.addVariable("myvar2",REAL8(9.0),LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED)
        self.assertTrue(self.livars.getVariableTypeByIndex(1)=="REAL8")

    def test_set(self):
        self.livars.addVariable("myvar1",REAL4(9.0),LALINFERENCE_REAL4_t,LALINFERENCE_PARAM_FIXED)
        self.livars.setVariable("myvar1",REAL4(5.0))
        self.assertTrue(self.livars.getVariable("myvar1")==5.0)

    def test_check(self):
        self.livars.addVariable("myvar1",REAL8(9.0),LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED)
        self.assertTrue(self.livars.checkVariable("myvar1"))
        self.assertFalse(self.livars.checkVariable("myvar2"))

    def test_copy(self):
        self.livars.addVariable("myvar1",REAL8(9.0),LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED)
        self.livars.addVariable("myvar2",REAL8(9.0),LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED)

        copy=self.livars.copyVariables()
        self.assertTrue(copy.checkVariable("myvar1"))
        self.assertTrue(copy.checkVariable("myvar2"))

    def test_cmp(self):
        pass

suite = unittest.TestSuite()
suite.addTest(unittest.makeSuite(test_LALInferenceVariables))

unittest.TextTestRunner(verbosity=2).run(suite)
