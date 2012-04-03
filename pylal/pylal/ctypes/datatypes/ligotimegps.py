#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#
#  time.py
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

from ctypes import *

import pylal.ctypes
from pylal.ctypes.datatypes.primitives import *
from pylal.ctypes.utils import __set_types

XLALINT8NSToGPS=_set_types(pylal.ctypes.liblal,"XLALINT8NSToGPS",POINTER(LIGOTimeGPS),[POINTER(LIGOTimeGPS),c_longlong])
XLALGPSToINT8NS=_set_types(pylal.ctypes.liblal,"XLALGPSToINT8NS",INT8,[POINTER(LIGOTimeGPS)])
XLALGPSAddGPS=_set_types(pylal.ctypes.liblal,"XLALGPSAddGPS",POINTER(LIGOTimeGPS),[POINTER(LIGOTimeGPS),POINTER(LIGOTimeGPS)])
XLALGPSDivide=_set_types(pylal.ctypes.liblal,"XLALGPSDivide",POINTER(LIGOTimeGPS),[POINTER(LIGOTimeGPS),REAL8])
XLALGPSGetREAL8=_set_types(pylal.ctypes.liblal,"XLALGPSGetREAL8",REAL8,[POINTER(LIGOTimeGPS)])

class LIGOTimeGPS(Structure):
    _fields_ = [("gpsSeconds",INT4),("gpsNanoSeconds",INT4)]

    def __init__(self,gpsSeconds,gpsNanoSeconds):
        Structure.__init__(self)
        
        self.gpsSeconds=INT4(gpsSeconds)
        self.gpsNanoSeconds=INT4(gpsNanoSeconds)
    
    def __copy__(self):
        new_gps=LIGOTimeGPS(self.gpsSeconds,self.gpsNanoSeconds)
    
    def __str__(self):
        print "LIGOTimeGPS(%s,%s)"%(self.gpsSeconds,self.gpsNanoSeconds)
        
    def __abs__(self):
        
        new_gps=copy(self)
        XLALINT8NSToGPS(pointer(new_gps),abs(XLALGPSToINT8NS(pointer(self))))
        
        return new_gps
        
    def __add__(self,other):
        
        return XLALGPSAddGPS(self,other)
        
    def __div__(self,other):
        
        new_gps=copy(self)
        return XLALGPSDivide(pointer(new_gps),other)
        
    def __float__(self)
    
        new_gps=copy(self)
        return XLALGPSGetREAL8(pointer(new_gps))
