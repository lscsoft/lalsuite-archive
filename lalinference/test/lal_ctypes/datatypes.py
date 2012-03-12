#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#
#  datatypes.py
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

#INT2
class INT2(c_int16): pass
#UINT2    
class UINT2(c_uint16): pass
#INT4
class INT4(c_int32): pass
#UINT4
class UINT4(c_uint32): pass
#INT8
class INT8(c_int64): pass
#UINT8
class UINT8(c_uint64): pass
#REAL4
class REAL4(c_float): pass
#REAL8
class REAL8(c_double): pass
#CHAR
class CHAR(c_char): pass
#UCHAR
class UCHAR(c_ubyte): pass
#BOOLEAN
class BOOLEAN(c_ubyte): pass
