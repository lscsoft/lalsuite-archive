from ctypes import *

class gsl_rng_type(Structure):pass

class gsl_rng(Structure):
    
    _fields_=[
        ("type",POINTER(gsl_rng_type)),
        ("state",c_void_p)
    ]

  
