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

### CONVENIENCE FUNCTIONS
#   Various convenience functions, some from ctypesGSL and may not be neeeded.

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

def _set_xlal_types(lib,fname, restype, argtypes):
    """
    Convenience function for setting the types on ctypes calls. This version
    decorates the calls with 
    
    """
    return _xlal_check_status(_set_types(lib,fname, restype, argtypes))
    
def _set_types(lib,fname, restype, argtypes):
    f = lib.__getattr__(fname)
    f.restype = restype
    f.argtypes = argtypes
    return f

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

### XLAL RETURN CODES
# XLAL error codes go here?

#   globals().update(__parse_xlal_error_codes('./XLALError.h'))

### EXTERNAL XLAL ERROR HANDLING

XLAL_REAL4_FAIL_NAN_INT=0x7fc001a1
XLAL_REAL8_FAIL_NAN_INT=0x7ff80000000001a1

def _xlal_check_status(xlal_call):
    """
    This function decorates XLAL calls with the appropriate (external) 
    error handler function.
    """
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
            if code.value is XLAL_REAL4_FAIL_NAN_INT:
                return True
            else: 
                return False
        error_code_test=error_code_test_real4
        
    elif restype is REAL8:
        def error_code_test_real8(code):
            if code.value is XLAL_REAL8_FAIL_NAN_INT:
                return True
            else:
                return False
        error_code_test=error_code_test_real8
    
    else:
        def error_code_test_none(code):
            return False
        error_code_test=error_code_test_none
        print "WARNING: XLAL return type not recognised ( %s ). XLAL errors will not be propagated as exceptions."%str(restype)
    
    def wrapper(*args, **kwargs):
        status_code=xlal_call(*args, **kwargs)
        if error_code_test(status_code):
            return xlal_status_handler(status_code)
        else:
            return status_code
        
    return wrapper

class XLAL_Error(RuntimeError):
    """
    Exception convenience class for propagating XLAL errors as Python RuntimeErrors.
    """
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

#### INTERNAL GSL ERROR HANDLING WRAPPER
# The purpose of these wrapper classes/functions is to implement a ctypes 
# call-back function to replace the internal error handlers. At the moment
# ctypes does not propagate a Python exception from within such a call-back
# (which would be the most simple/convenient mechanism for replacing SIGABRT/SIGKILL
# signals emanating from within xlal code. At the moment the ctypes wrapper is setup to
# replace the internal error handlers and print the information to the console instead.
# It is possible to conceive of a small modification to the library (such as checking a global 
# flag or implementing a Python error handler using the C API) which would allow the xlal library code to 
# exit gracefully into the Python runtime and generate an exception.
#
# NOTE that this mechanism for handling XLAL errors in Python is independent from the above routines which test
# the output value after it has been returned to Python. However, by using them in conjunction you can get information
# about internal errors which do not propagate if the errors get handled within the library. 


XLAL_ERROR_HANDLER_T=CFUNCTYPE(None,c_char_p,c_char_p,c_int,c_int);
XLALSetErrorHandler=_set_types(pylal.ctypes.liblal,"XLALSetErrorHandler", XLAL_ERROR_HANDLER_T, [XLAL_ERROR_HANDLER_T])
XLALSetDefaultErrorHandler=_set_types(pylal.ctypes.liblal,"XLALSetDefaultErrorHandler", XLAL_ERROR_HANDLER_T, [])

# Python function callback definitions which convert XLAL errors into exceptions
def __xlal_error_handler_warning(reason, fl, line, xlal_errno):
    """
    This callback just prints the information found in the internal error structures.
    """    
    print "WARNING: " + str(XLAL_Error(xlal_errno, fl, line, reason))

def __xlal_error_handler_exception(reason, fl, line, xlal_errno):
    """ 
    This function does not work as intended at the moment as ctypes callbacks do not
    propagate exceptions.
    """
    raise XLAL_Error(xlal_errno, fl, line, reason) 

# Convert above Python callback functions into ctypes callback functions...
_xlal_error_handler_warning   = XLAL_ERROR_HANDLER_T(__xlal_error_handler_warning)
_xlal_error_handler_exception = XLAL_ERROR_HANDLER_T(__xlal_error_handler_exception)

#Generate a pointer to the appropriate error handler
def set_error_handler_warning():
    return XLALSetErrorHandler(_xlal_error_handler_warning)
def set_error_handler_exception():
    return XLALSetErrorHandler(_xlal_error_handler_exception)

# Set the default handler (the warning handler atm as the exception one doesnt have an effect)
xlal_default_error_handler = set_error_handler_warning()

