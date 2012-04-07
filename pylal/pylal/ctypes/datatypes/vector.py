#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#
#  vector.py
#  
#  Part of the ctypes wrapper library for LAL.
#
#  Copyright 2012 Ben Aylott <beaylott@gmail.com>
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

from ctypes import Structure,byref

import pylal.ctypes
from pylal.ctypes.utils import _set_types
from pylal.ctypes.datatypes.primitives import *
from pylal.ctypes.datatypes.complex import *
from pylal.ctypes.datatypes.ligotimegps import LIGOTimeGPS

class REAL8Vector(Structure):
    
    _fields_ = [("length",UINT4),("data",POINTER(REAL8))]
    
    def __init__(self,length):
        Structure.__init__(self)
    
    def __new__(cls,length):
        return XLALCreateREAL8Vector(UINT4(length)).contents
    
    def __getitem__(self,item_idx):
        return cast(byref(self.data.contents,item_idx),POINTER(REAL8)).contents
        
class COMPLEX16Vector(Structure):
    _fields_ = [("length",UINT4),("data",POINTER(COMPLEX16))]
    
    def __init__(self,length):
        Structure.__init__(self)
        
    def __new__(cls,length):
        return XLALCreateCOMPLEX16Vector(UINT4(length)).contents

REAL8Sequence=REAL8Vector
COMPLEX16Sequence=COMPLEX16Vector

XLALCreateREAL8Vector=_set_types(pylal.ctypes.liblal,"XLALCreateREAL8Vector",POINTER(REAL8Vector),[UINT4])
XLALCreateCOMPLEX16Vector=_set_types(pylal.ctypes.liblal,"XLALCreateCOMPLEX16Vector",POINTER(COMPLEX16Vector),[UINT4])

class LIGOTimeGPSVector(Structure):
    _fields_ = [("length",UINT4),("data",POINTER(LIGOTimeGPS)),("deltaT",REAL8)]
