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

from vector import *
from ctypes import Structure

from pylal.ctypes.ligotimegps import LIGOTimeGPS
from pylal.ctypes.lalunit import LALUnit

class COMPLEX16TimeSeries(Structure):
    
    _fields_ = [
        ("name",CHAR),
        ("epoch",LIGOTimeGPS),
        ("deltaT",REAL8),
        ("sampleUnits",LALUnit),
        ("data",POINTER(COMPLEX16Sequence))
    ]
