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

from datatypes import *
from complex import *
from ctypes import *
from laltime import *

class REAL8Vector(Structure):
    
    _fields_ = [("length",UINT4),("data",POINTER(REAL8))]
    
    def __init__(self,np_vector):
        Structure.__init__(self)
        self.data=np_vector.ctypes.data_as(REAL8)
        
class COMPLEX16Vector(Structure):
    _fields_ = [("length",UINT4),("data",POINTER(COMPLEX16))]
    
    def __init__(self,np_vector):
        Structure.__init__(self)
        self.data=np_vector.ctypes.data_as(COMPLEX16)
    
class LIGOTimeGPSVector(Structure):
    _fields_ = [("length",UINT4),("data",POINTER(LIGOTimeGPS)),("deltaT",REAL8)]

    
REAL8Sequence=REAL8Vector
COMPLEX16Sequence=COMPLEX16Vector
