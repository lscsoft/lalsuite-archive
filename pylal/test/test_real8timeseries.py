from ctypes import cast,create_string_buffer,POINTER,c_ssize_t,byref

from pylal.ctypes.datatypes.primitives import CHAR,REAL8
from pylal.ctypes.datatypes.ligotimegps import LIGOTimeGPS
from pylal.ctypes.datatypes.lalunit import lalDimensionlessUnit
from pylal.ctypes.datatypes.real8timeseries import REAL8TimeSeries,XLALCreateREAL8TimeSeries
from pylal.ctypes.datatypes.vector import REAL8Vector

vec=REAL8Vector(1000)
print vec

print vec[90].value

name_p=cast(create_string_buffer("ASD"),POINTER(CHAR))
epoch_p=byref(LIGOTimeGPS(0.))
sampleUnits_p=byref(lalDimensionlessUnit)

vec=XLALCreateREAL8TimeSeries(name_p,epoch_p,REAL8(30.),REAL8(0.0001),sampleUnits_p,c_ssize_t(10000))
print vec
print vec.contents.name

