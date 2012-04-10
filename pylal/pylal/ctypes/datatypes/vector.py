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

from ctypes import Structure,byref,memset,sizeof,cast

import numpy as np

import pylal.ctypes
from pylal.ctypes.utils import _set_types,ptr_add
from pylal.ctypes.datatypes.primitives import *
from pylal.ctypes.datatypes.complex import *
from pylal.ctypes.datatypes.ligotimegps import LIGOTimeGPS

#= AVFactory API =#

vectorTypes=[
    "REAL8",
    "REAL4",
    "COMPLEX8",
    "COMPLEX16",
    "INT4",
    "INT8",
    "UINT4"
]

vectorFunctionTable=[
    ["Create","POINTER(%(vtype)sVector)","[UINT4]"],
    ["Resize","POINTER(%(vtype)sVector)","[POINTER(%(vtype)sVector),UINT4]"],
    ["Destroy","None","[POINTER(%(vtype)sVector)]"]
]

def VectorFactory(vtype):
    class Vector(Structure):
        
        def __init__(self,length):
            Structure.__init__(self)
            self.data=cast(pointer((self.vtype*length)(0.)),POINTER(self.vtype))
            self.length=UINT4(length)
        
        def __getitem__(self,item_idx):
            
            if item_idx>=self.length.value or item_idx<0:
                raise IndexError
                
            elif item_idx==0:
                return self.data.contents.value
                
            else:
                a=ptr_add(self.data,item_idx*sizeof(self.vtype))
                return a.contents.value

        def as_array(self):
            return np.ctypeslib.as_array(self.data,(int(self.length.value),))
            
        def resize(self,size):
            pass
            
        def delete(self):
            pass

    Vector.vtype=vtype
    Vector._fields_ = [("length",UINT4),("data",POINTER(vtype))]
    
    return Vector

def __create_vector_and_sequence_types(vectorTypes):
    for vectorType in vectorTypes:
        globals()[vectorType+"Vector"]=VectorFactory(eval(vectorType))
        globals()[vectorType+"Sequence"]=eval(vectorType+"Vector")
        
def __create_vector_functions(vectorTypes,vectorFunctionTable):
    
    for vtype in vectorTypes:
        for vfname,vfrest,vfargt in vectorFunctionTable:
            function_name="XLAL"+vfname+vtype+"Vector"
            globals()[function_name]=_set_types(pylal.ctypes.liblal,function_name,eval(vfrest%{"vtype":vtype}),eval(vfargt%{"vtype":vtype}))

def _generate_avfactory_api(vectorTypes,vectorFunctionTable):
    __create_vector_and_sequence_types(vectorTypes)
    __create_vector_functions(vectorTypes,vectorFunctionTable)


_generate_avfactory_api(vectorTypes,vectorFunctionTable)

#= LIGOTimeGPSVector =#

class LIGOTimeGPSVector(Structure):
    _fields_ = [("length",UINT4),("data",POINTER(LIGOTimeGPS)),("deltaT",REAL8)]
