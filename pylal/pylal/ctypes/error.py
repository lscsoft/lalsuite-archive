#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#
#  error.py
#  
#  Copyright 2012 Ben Aylott <ben@amber>
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
from pylal.ctypes.utils import __parse_xlal_error_codes,_set_types

### RETURN CODES
# XLAL error codes go here?

#   globals().update(__parse_xlal_error_codes('./XLALError.h'))

#XLALIsREAL4FailNaN=_set_types(pylal.ctypes.liblal,"XLALIsREAL4FailNaN",c_int,[REAL8])
#XLALIsREAL8FailNaN=_set_types(pylal.ctypes.liblal,"XLALIsREAL8FailNaN",c_int,[REAL8])

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
            msg += strerror(gsl_err_code)
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
def _xlal_status_exception(status_code, result):
    raise XLAL_Error(status_code, result = result)
def _xlal_status_warning(status_code, result):
    print "WARNING: " + str(XLAL_Error(status_code, result = result))
    return status_code
def _xlal_status_off(status_code, result):
    return status_code

# current status handler function
xlal_status_handler = _xlal_status_exception

def _xlal_check_status_void_p(xlal_call):
    def wrapper(*args, **kwargs):
        rvalue=xlal_call(*args, **kwargs)
        if rvalue!=0x0:
            xlal_status_handler(status_code, result)
        return rvalue

def _xlal_check_status_int(xlal_call):
    def wrapper(*args, **kwargs):
        rvalue=xlal_call(*args, **kwargs)
        if rvalue!=0:
            return xlal_status_handler(status_code, result)
        else:
            return rvalue

def _xlal_check_status_real4(xlal_call):
    def wrapper(*args, **kwargs):
        status_code=xlal_call(*args, **kwargs)
        if rvalue!=0:
            return xlal_status_handler(status_code, result)
        else:
            return status_code

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

# method for testing return values
def _xlal_check_status(status_code, result = None):
    """Raises an exception if a XLAL function returns an error
    condition."""
    if status_code != XLAL_SUCCESS:
        xlal_status_handler(status_code, result)
    return status_code

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
