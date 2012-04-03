#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#
#  laldetector.py
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

from ctypes import Structure,POINTER

from datatypes import REAL4,REAL8

class LALDetectorType(Structure): pass

class LALFrDetector(Structure):
    _fields_= [
        ("name",CHAR*LALNameLength),
        ("prefix",CHAR*3),
        ("vertexLongitudeRadians",REAL8),
        ("vertexLatitudeRadians",REAL8),
        ("vertexElevation",REAL4),
        ("xArmAltitudeRadians",REAL4),
        ("xArmAzimuthRadians",REAL4),
        ("yArmAltitudeRadians",REAL4),
        ("yArmAzimuthRadians",REAL4),
        ("xArmMidpoint",REAL4),
        ("yArmMidpoint",REAL4)
    ]
    
class LALDetector(Structure):
    _fields_= [
        ("location",REAL8*3),
        ("response",(REAL4*3)*3),
        ("type",LALDetectorType),
        ("frDetector",LALFrDetector)   
    ]

