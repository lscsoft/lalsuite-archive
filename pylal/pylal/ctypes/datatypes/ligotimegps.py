#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#
#  ligotimegps.py
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

from math import floor
from copy import copy

from ctypes import Structure,byref,POINTER,c_longlong

import pylal.ctypes
from pylal.ctypes.datatypes.primitives import INT4,INT8,REAL4,REAL8
from pylal.ctypes.utils import _set_types


class LIGOTimeGPS(Structure):
    _fields_ = [("gpsSeconds",INT4),("gpsNanoSeconds",INT4)]

    def __init__(self,*args):
        Structure.__init__(self)
        
        if len(args)==1:
            gpsTime=float(args[0])
            self.gpsSeconds=INT4(int(floor(gpsTime)))
            self.gpsNanoSeconds=INT4(int(round((gpsTime-floor(gpsTime))/10**(-9))))
            
        elif len(args)==2:
            gpsNanoSeconds=int(args[1])
            self.gpsSeconds=INT4(int(args[0])+int(floor(float(gpsNanoSeconds)/1e9)))
            self.gpsNanoSeconds=INT4(gpsNanoSeconds%1000000000)
        
        elif len(args)==0:
            self.gpsSeconds=INT4(0)
            self.gpsNanoSeconds=INT4(0)
            
        else:
            raise ValueError
    
    def __copy__(self):
        new_gps=LIGOTimeGPS()
        new_gps.gpsSeconds=self.gpsSeconds
        new_gps.gpsNanoSeconds=self.gpsNanoSeconds
        return new_gps
        
    def __str__(self):
        return "LIGOTimeGPS(%i,%i)"%(self.gpsSeconds.value,self.gpsNanoSeconds.value)
        
    def __abs__(self):
        
        new_gps=copy(self)
        XLALINT8NSToGPS(byref(new_gps),abs(XLALGPSToINT8NS(byref(self)).value))
        
        return new_gps
        
    def __add__(self,other):
        new_gps=copy(self)
        
        XLALGPSAddGPS(byref(new_gps),byref(other))
        
        return new_gps
        
    def __div__(self,other):
        
        new_gps=copy(self)
        return XLALGPSDivide(byref(new_gps),other).contents
        
    def __float__(self):
        new_gps=copy(self)
        return float(XLALGPSGetREAL8(byref(new_gps)).value)
        
    def __int__(self):
        return int(self.gpsSeconds.value)
        
    def __mod__(self,mod):
        new_gps=LIGOTimeGPS()
        return XLALINT8NSToGPS(byref(new_gps),INT8(int(round(XLALGPSToINT8NS(byref(self)).value%(1e9*mod))))).contents
        
    def __mul__(self,other):
        
        new_gps=copy(self)
        return XLALGPSMultiply(byref(new_gps),other).contents
        
    def __neg__(self):
        new_gps=copy(self)
        return XLALINT8NSToGPS(byref(new_gps), INT8(-XLALGPSToINT8NS(byref(new_gps)).value)).contents
        
        
    def __nonzero__(self):
        if bool(self.gpsSeconds) and bool(self.gpsNanoSeconds):
            return True
        else:
            return False
            
    def ns(self):
        new_gps=copy(self)
        return int(XLALGPSToINT8NS(byref(new_gps)).value)
        
    def __pos__(self):
        raise NotImplemented
        
    def __reduce__(self):
        return (LIGOTimeGPS,(self.gpsSeconds,self.gpsNanoSeconds))
        
    def __repr__(self):
        return "LIGOTimeGPS(%i,%i)"%(self.gpsSeconds.value,self.gpsNanoSeconds.value)

    def __cmp__(self,other):
        return self.__eq__(other)
            
    def __eq__(self,other):
        if int(self.gpsSeconds.value)==int(other.gpsSeconds.value) and int(self.gpsNanoSeconds.value)==int(other.gpsNanoSeconds.value):
            return True
        else:
            return False
    
    def __neq__(self,other):
        return not self.__eq__(other)
            
    def __hash__(self):
        raise NotImplemented
        
    def sub(self,other):
        return XLALINT8NSToGPS(XLALGPSToINT8NS(byref(self)))-XLALGPSToINT8NS(byref(other))

XLALINT8NSToGPS=_set_types(pylal.ctypes.liblal,"XLALINT8NSToGPS",POINTER(LIGOTimeGPS),[POINTER(LIGOTimeGPS),c_longlong])
XLALGPSToINT8NS=_set_types(pylal.ctypes.liblal,"XLALGPSToINT8NS",INT8,[POINTER(LIGOTimeGPS)])
XLALGPSAddGPS=_set_types(pylal.ctypes.liblal,"XLALGPSAddGPS",POINTER(LIGOTimeGPS),[POINTER(LIGOTimeGPS),POINTER(LIGOTimeGPS)])
XLALGPSDivide=_set_types(pylal.ctypes.liblal,"XLALGPSDivide",POINTER(LIGOTimeGPS),[POINTER(LIGOTimeGPS),REAL8])
XLALGPSGetREAL8=_set_types(pylal.ctypes.liblal,"XLALGPSGetREAL8",REAL8,[POINTER(LIGOTimeGPS)])
XLALGPSMultiply=_set_types(pylal.ctypes.liblal,"XLALGPSMultiply",POINTER(LIGOTimeGPS),[POINTER(LIGOTimeGPS),REAL8])
