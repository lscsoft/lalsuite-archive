from ctypes import cast,create_string_buffer,POINTER,c_ssize_t,byref,addressof

import numpy as np

from pylal.ctypes.datatypes.primitives import CHAR,REAL8
from pylal.ctypes.datatypes.ligotimegps import LIGOTimeGPS
from pylal.ctypes.datatypes.lalunit import lalDimensionlessUnit
from pylal.ctypes.datatypes.real8timeseries import REAL8TimeSeries,XLALCreateREAL8TimeSeries
from pylal.ctypes.datatypes.vector import REAL8Vector,REAL4Vector

vec=REAL4Vector(100)

veca=vec.as_array()
print veca

for i in range(0,4):
    print i,vec[i]

