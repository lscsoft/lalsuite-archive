#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#
#   utils.py
#  
#   Part of the ctypes wrapper library for LAL.
#
#   Adapted from part of ctypesGSL by Szymon Jaroszewicz. 
#
#   Copyright 2012 Szymon Jaroszewicz , Ben Aylott <beaylott@gmail.com>
#  
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.

import os
import fnmatch
from functools import partial 
from types import MethodType
from ctypes import CDLL,c_uint,cast,c_void_p
from ctypes.util import find_library
from ctypes import *
import pylal.ctypes
from pylal.ctypes.datatypes.primitives import *

def _parse_xlal_error_codes(path_to_xlalerror_h):
    import re
    
    xlal_error_h=open(path_to_xlalerror_h,'r')
    xlal_error_h_str=''.join(xlal_error_h.readlines()).replace('\n','')
    
    re_match_enum=re.compile('enum*\sXLALErrorValue.*?{.*?}')
    re_match_key_value_pairs=re.compile('XLAL_(.+?),')
    
    enum_str=re_match_enum.search(xlal_error_h_str)
    key_value_str=re_match_key_value_pairs.findall(enum_str.group())
    
    xlal_error_codes_dict={}
    for key_value in key_value_str:
        pair=key_value.replace(' ','').split('=')
        xlal_error_codes_dict[int(pair[1])]='XLAL_'+pair[0]
    
    return xlal_error_codes_dict
            

def makeglobal(vdict):
    
    globals().update(vdict)

class enum_typedef(c_uint): pass

def make_enum_typedef(enum_typedef_name,enum_names):
    vdict={}
    new_enum_typedef=type(enum_typedef_name,(enum_typedef,),{})
    vdict[enum_typedef_name]=new_enum_typedef
    i=0
    for enum_name in enum_names:
        vdict[enum_name]=new_enum_typedef(i)
        i+=1
    return vdict

def make_enum_anon(enum_names):
    vdict={}
    i=0
    for enum_name in enum_names:
        vdict[enum_name]=enum_typedef(i)
        i+=1
    return vdict
        
class PkgConfig(object):
    def __init__(self, names):
        def stripfirsttwo(string):
            return string[2:]
        self.libs = map(stripfirsttwo, os.popen("pkg-config --libs-only-l %s" % names).read().split())
        self.libdirs = map(stripfirsttwo, os.popen("pkg-config --libs-only-L %s" % names).read().split())
        self.incdirs = map(stripfirsttwo, os.popen("pkg-config --cflags-only-I %s" % names).read().split())
        self.extra_cflags = os.popen("pkg-config --cflags-only-other %s" % names).read().split()

def _xlal_check_status(xlal_call):
    restype=xlal_call.restype
    
    if restype is c_int:
        def error_code_test_int(code):
            if code!=0:
                return True
            else:  
                return False
        error_code_test=error_code_test_int
    
    elif restype is c_void_p:
        def error_code_test_void_p(code):
            
            if code is None:
                return True
            else:  
                return False
        error_code_test=error_code_test_void_p
                
    elif restype is REAL4:
        def error_code_test_real4(code):
            error_code=XLALIsREAL4FailNaN(code)
            return bool(error_code.value)
        error_code_test=error_code_test_real4
        
    elif restype is REAL8:
        def error_code_test_real8(code):
            error_code=XLALIsREAL8FailNaN(code)
            return bool(error_code.value)
        error_code_test=error_code_test_real8
    
    else:
        def error_code_test_none(code):
            return False
        error_code_test=error_code_test_none
        print "WARNING: XLAL return type not recognised, XLAL errors will not be propagated as exceptions."+str(restype)
    
    def wrapper(*args, **kwargs):
        status_code=xlal_call(*args, **kwargs)
        if error_code_test(status_code):
            return xlal_status_handler(status_code)
        else:
            return status_code
        
    return wrapper

def _set_xlal_types(lib,fname, restype, argtypes):
    return _xlal_check_status(_set_types(lib,fname, restype, argtypes))
    
def _set_types(lib,fname, restype, argtypes):
    f = lib.__getattr__(fname)
    f.restype = restype
    f.argtypes = argtypes
    return f

XLALIsREAL4FailNaN=None#_set_types(pylal.ctypes.liblal,"XLALIsREAL4FailNaN",c_int,[REAL8])
XLALIsREAL8FailNaN=None#_set_types(pylal.ctypes.liblal,"XLALIsREAL8FailNaN",c_int,[REAL8])


def _add_function(fname, restype, argtypes, globs = globals()):
    new_m = _set_types(fname, restype, argtypes)
    globs[fname] = new_m

def _new_method(func,*args, **kwargs):
    """Helper for dynamically adding methods to classes."""
    return func(*args, **kwargs)
    
def _add_method(lib,fname,cls,method_name,restype,argtypes):
    
    new_m = _set_types(lib,fname,restype,argtypes)
    method = new_m
    setattr(cls, "_"+method_name, method)

def make_class(my_class,class_table): 

    for lib,c_lib_fname,py_cls_fname,restype,argtypes in class_table:
        _add_method(lib,c_lib_fname,my_class,py_cls_fname,restype,argtypes)

def _int_to_bool(i):
    return bool(i)
    
    from ctypes import *

#ptr_add
#From: http://permalink.gmane.org/gmane.comp.python.ctypes/4343
def ptr_add(ptr, other):
    
    try:
        ofs = other.__index__()
    except AttributeError:
        raise TypeError("Can only add integer to pointer")
    p = cast(ptr, c_void_p)
    
    p.value += ofs
    
    return cast(p, type(ptr))

#ptr_sub
#From: http://permalink.gmane.org/gmane.comp.python.ctypes/4343
def ptr_sub(ptr, other):
    
    if type(ptr) == type(other):
        return cast(ptr, c_void_p).value - cast(other, c_void_p).value
    try:
        ofs = other.__index__()
    except AttributeError:
        raise TypeError("Can only substract pointer or integer from pointer")
    p = cast(ptr, c_void_p)
    p.value -= ofs
    return cast(p, type(ptr))
    


### RETURN CODES
# XLAL error codes go here?

#   globals().update(__parse_xlal_error_codes('./XLALError.h'))


### ERROR HANDLING

#_add_function("strerror", c_char_p, [c_int])

class XLAL_Error(RuntimeError):
    def __init__(self, xlal_err_code, fl = None, line = None, reason = None, result = None):
        msg = "xlal: "
        if fl is not None:
            msg += fl+":"
        if line is not None:
            msg += str(line)+":"
        if msg[-1] != " ":
            msg += " "
        msg += "ERROR: "
        if reason is not None:
            msg += reason
        else:
            #TODO: convert error codes to messages!
            msg+=str(xlal_err_code)
            
        RuntimeError.__init__(self, msg)
        self.xlal_err_code = xlal_err_code
        self.fl = fl # file in which error occurred
        self.line = line
        self.reason = reason
        self.result = result # result to be returned by the function

# helper function to test allocation
#def _xlal_check_null_pointer(p):
    #"""Raises an exception if a allocation failed in a XLAL function."""
    #p = cast(p, c_void_p)
    #if p == 0:
        #raise XLAL_Error(XLAL_ENOMEM)
    #return p

### XLAL error handling
# handling of XLAL status returned by functions
def _xlal_status_exception(status_code):
    raise XLAL_Error(status_code)
def _xlal_status_warning(status_code, result):
    print "WARNING: " + str(XLAL_Error(status_code, result = result))
    return status_code
def _xlal_status_off(status_code, result):
    return status_code

# current status handler function
xlal_status_handler = _xlal_status_exception

#def _xlal_check_status_void_p(xlal_call):
    #def wrapper(*args, **kwargs):
        #rvalue=xlal_call(*args, **kwargs)
        #if rvalue!=0:
            #xlal_status_handler(status_code, result)
        #return rvalue
    #return wrapper
    
#def _xlal_check_status_int(xlal_call):
    #def wrapper(*args, **kwargs):
        #rvalue=xlal_call(*args, **kwargs)
        #if rvalue!=0:
            #return xlal_status_handler(status_code, result)
        #else:
            #return rvalue
    #return wrapper
    

    
def set_status_handler(h):
    global xlal_status_handler
    old_handler = xlal_status_handler
    xlal_status_handler = h
    return old_handler
    
def set_status_handler_off():    
    return set_status_handler(_xlal_status_off)
    
def set_status_handler_warning():
    return set_status_handler(_xlal_status_warning)
    
def set_status_handler_exception():
    return set_status_handler(_xlal_status_exception)

#### internal gsl error handling
XLAL_ERROR_HANDLER_T=CFUNCTYPE(None,c_char_p,c_char_p,c_int,c_int);
XLALSetErrorHandler=_set_types(pylal.ctypes.liblal,"XLALSetErrorHandler", XLAL_ERROR_HANDLER_T, [XLAL_ERROR_HANDLER_T])
XLALSetDefaultErrorHandler=_set_types(pylal.ctypes.liblal,"XLALSetDefaultErrorHandler", XLAL_ERROR_HANDLER_T, [])

# Python function callback definitions which convert XLAL errors into exceptions
def __xlal_error_handler_warning(reason, fl, line, xlal_errno):
    print "WARNING: " + str(XLAL_Error(xlal_errno, fl, line, reason))
    # !!! currently, due to limitations in ctypes, this exception is
    # !!! not thrown to the program, just printed
    # raise XLAL_Error(gsl_errno, fl, line, reason)
def __xlal_error_handler_exception(reason, fl, line, xlal_errno):
    # !!! currently, due to limitations in ctypes, this exception is
    # !!! not thrown to the program, just printed
    raise XLAL_Error(xlal_errno, fl, line, reason)

# Convert above Python callback functions into ctypes callback functions
_xlal_error_handler_warning   = XLAL_ERROR_HANDLER_T(__xlal_error_handler_warning)
_xlal_error_handler_exception = XLAL_ERROR_HANDLER_T(__xlal_error_handler_exception)


def set_error_handler_warning():
    return XLALSetErrorHandler(_xlal_error_handler_warning)
def set_error_handler_exception():
    return XLALSetErrorHandler(_xlal_error_handler_exception)

# Set the default handler
xlal_default_error_handler = set_error_handler_warning()

