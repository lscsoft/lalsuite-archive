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

from ctypes import Structure,POINTER

from datatypes import INT2,UINT2
from utils import make_li_enum_typedef     

make_li_enum_typedef("LALUnit",
    [
        "LALUnitIndexMeter",
        "LALUnitIndexKiloGram",
        "LALUnitIndexSecond"
        "LALUnitIndexAmpere",
        "LALUnitIndexKelvin",
        "LALUnitIndexStrain",
        "LALUnitIndexADCCount",
        "LALNumUnits"
    ]
)

class LALUnit(Structure):
    _fields_=[
        ("powerOfTen",INT2),
        ("unitNumerator",POINTER(INT2)),
        ("unitDenominatorMinusOne",POINTER(UINT2))
    ]
