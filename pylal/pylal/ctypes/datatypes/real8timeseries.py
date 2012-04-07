#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#
#  timeseries.py
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


from ctypes import Structure,POINTER,byref,create_string_buffer,cast,c_char_p,c_ssize_t

import pylal.ctypes
from pylal.ctypes.utils import _set_types
from pylal.ctypes.datatypes.primitives import CHAR,REAL8
from pylal.ctypes.datatypes.ligotimegps import LIGOTimeGPS
from pylal.ctypes.datatypes.lalunit import LALUnit
from pylal.ctypes.datatypes.vector import REAL8Sequence

class REAL8TimeSeries(Structure):
    
    _fields_ = [
        ("name",CHAR*100),
        ("epoch",LIGOTimeGPS),
        ("deltaT",REAL8),
        ("f0",REAL8),
        ("sampleUnits",LALUnit),
        ("data",POINTER(REAL8Sequence))
    ]
    
    def __init__(self,name,epoch,f0,deltaT,sampleUnits,length):
        Structure.__init__(self)
    
    def __new__(cls,name,epoch,f0,deltaT,sampleUnits,length):
        name_p=cast(create_string_buffer(name),POINTER(CHAR))
        epoch_p=byref(epoch)
        sampleUnits_p=byref(sampleUnits)
        return XLALCreateREAL8TimeSeries(name_p,epoch_p,REAL8(f0),REAL8(deltaT),sampleUnits_p,c_ssize_t(length)).contents
    
XLALCreateREAL8TimeSeries=_set_types(pylal.ctypes.liblal,"XLALCreateREAL8TimeSeries",POINTER(REAL8TimeSeries),[POINTER(CHAR),POINTER(LIGOTimeGPS),REAL8,REAL8,POINTER(LALUnit),c_ssize_t])
