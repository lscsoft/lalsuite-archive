#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#
#  lalunit.py
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
from pylal.ctypes.utils import make_enum_anon,_set_types

globals().update(
    make_enum_anon(
        [
            "LALUnitIndexMeter",
            "LALUnitIndexKiloGram",
            "LALUnitIndexSecond",
            "LALUnitIndexAmpere",
            "LALUnitIndexKelvin",
            "LALUnitIndexStrain",
            "LALUnitIndexADCCount",
            "LALNumUnits"
        ]
    )
)

class LALUnit(Structure):
    _fields_=[
        ("powerOfTen",INT2),
        ("unitNumerator",INT2*int(LALNumUnits.value)),
        ("unitDenominatorMinusOne",UINT2*int(LALNumUnits.value))
    ]
    
    def __set(self,powerOfTen,unitNumerator,unitDenominatorMinusOne):
        self.powerOfTen=INT2(powerOfTen)
        
        for i in range(len(self.unitNumerator)):
            self.unitNumerator[i]=INT2(unitNumerator[i])
    
        for i in range(len(self.unitDenominatorMinusOne)):
            self.unitDenominatorMinusOne[i]=UINT2(unitDenominatorMinusOne[i])
    
    def __init__(self,powerOfTen=0,unitNumerator=[0]*int(LALNumUnits.value),unitDenominatorMinusOne=[0]*int(LALNumUnits.value)):
        
        self.__set(powerOfTen,unitNumerator,unitDenominatorMinusOne)
    
    def __repr__(self):
        rep_str=create_string_buffer(50)
        XLALUnitAsString(rep_str,sizeof(rep_str),pointer(self))
        return str(rep_str)
        
    def __str__(self):
        return self.__repr__()
        
    def __cmp__(self,other):
        return bool(XLALUnitCompare(pointer(self),pointer(other)).value)
        
    def __hash__(self): 
        raise NotImplemented
    
    def __mul__(self,other):
        return XLALUnitMultiply(pointer(self),pointer(other)).contents
        
    def __div__(self,other):
        return XLALUnitDivide(pointer(self),pointer(other)).contents
        
    def __pow__(self):
        raise NotImplemented
        
    def __float__(self):
        if not bool(XLALUnitIsDimensionless(self)):
            raise ValueError
        
        return pow(10.,int(self.powerOfTen))
        
    def __invert__(self):
        new=LALUnit()
        return XLALUnitInvert(pointer(new),pointer(self)).contents

XLALUnitAsString=_set_types(pylal.ctypes.liblal,"XLALUnitAsString",char_p,[char_p,UINT4,POINTER(LALUnit)])
XLALUnitCompare=_set_types(pylal.ctypes.liblal,"XLALUnitCompare",c_int,[POINTER(LALUnit),POINTER(LALUnit)])
XLALUnitMultiply=_set_types(pylal.ctypes.liblal,"XLALUnitMultiply",POINTER(LALUnit),[POINTER(LALUnit),POINTER(LALUnit)])
XLALUnitDivide=_set_types(pylal.ctypes.liblal,"XLALUnitDivide",POINTER(LALUnit),[POINTER(LALUnit),POINTER(LALUnit)])
XLALUnitIsDimensionless=_set_types(pylal.ctypes.liblal,"XLALUnitIsDimensionless",c_int,[POINTER(LALUnit)])
XLALUnitInvert=_set_types(pylal.ctypes.liblal,"XLALUnitInvert",POINTER(LALUnit),[POINTER(LALUnit),POINTER(LALUnit)])

lalDimensionlessUnit=cast(pylal.ctypes.liblal.lalDimensionlessUnit,POINTER(LALUnit))
lalMeterUnit=cast(pylal.ctypes.liblal.lalMeterUnit,POINTER(LALUnit))
lalKiloGramUnit=cast(pylal.ctypes.liblal.lalKiloGramUnit,POINTER(LALUnit))
lalSecondUnit=cast(pylal.ctypes.liblal.lalSecondUnit,POINTER(LALUnit))
lalAmpereUnit=cast(pylal.ctypes.liblal.lalAmpereUnit,POINTER(LALUnit))
lalKelvinUnit=cast(pylal.ctypes.liblal.lalKelvinUnit,POINTER(LALUnit))
lalStrainUnit=cast(pylal.ctypes.liblal.lalStrainUnit,POINTER(LALUnit))
lalADCCountUnit=cast(pylal.ctypes.liblal.lalADCCountUnit,POINTER(LALUnit))
lalYottaUnit=cast(pylal.ctypes.liblal.lalYottaUnit,POINTER(LALUnit))
