#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#
#  real8frequencyseries.py
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


from ctypes import Structure,POINTER,pointer,byref,create_string_buffer,cast,c_char_p,c_size_t

import pylal.ctypes
from pylal.ctypes.utils import _set_types
from pylal.ctypes.datatypes.primitives import CHAR,REAL8
from pylal.ctypes.datatypes.ligotimegps import LIGOTimeGPS
from pylal.ctypes.datatypes.lalunit import LALUnit,lalDimensionlessUnit
from pylal.ctypes.datatypes.vector import REAL8Sequence

class REAL8FrequencySeries(Structure):
    
    _fields_ = [
        ("name",CHAR),
        ("epoch",LIGOTimeGPS),
        ("f0",REAL8),
        ("deltaF",REAL8),
        ("sampleUnits",LALUnit),
        ("data",POINTER(REAL8Sequence))
    ]
    
    def __init__(self,name="",epoch=LIGOTimeGPS(0.),f0=20.,deltaT=1e-9,sampleUnits=lalDimensionlessUnit,length=0):
        Structure.__init__(self)
        self.name=name
        self.epoch=epoch
        self.sampleUnits=sampleUnits
        self.f0=REAL8(f0)
        self.deltaT=REAL8(deltaT)
        self.data=pointer(REAL8Sequence(np.array(length*[0.],dtype=np.float64)))
    
#XLALCreateREAL8TimeSeries=_set_types(pylal.ctypes.liblal,"XLALCreateREAL8TimeSeries",POINTER(REAL8TimeSeries),[POINTER(CHAR),POINTER(LIGOTimeGPS),REAL8,REAL8,POINTER(LALUnit),c_size_t])
