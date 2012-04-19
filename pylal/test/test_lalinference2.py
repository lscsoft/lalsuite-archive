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
from pylal.ctypes.lalinference import LALINFERENCE_PARAM_FIXED,LALINFERENCE_REAL8_t,LALINFERENCE_REAL4_t,LALINFERENCE_INT4_t,LALINFERENCE_PARAM_LINEAR


livars=LALInferenceVariables()
livars.addVariable("myvar1",9.0,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED)
livars.getVariable("myvar2")
