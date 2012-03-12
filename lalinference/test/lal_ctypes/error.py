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

def _gsl_return_to_bool(i):
    """Convert return value to Bool.

    No exceptions are raised.  True return value means successful
    execution (the meaning depends on specific function)."""
    return i == 1 # XLAL_SUCCESS

### RETURN CODES
# XLAL error codes go here?

### ERROR HANDLING

#_add_function("strerror", c_char_p, [c_int])

#class GSL_Error(RuntimeError):
    #def __init__(self, gsl_err_code, fl = None, line = None, reason = None, result = None):
        #msg = "gsl: "
        #if fl is not None:
            #msg += fl+":"
        #if line is not None:
            #msg += str(line)+":"
        #if msg[-1] != " ":
            #msg += " "
        #msg += "ERROR: "
        #if reason is not None:
            #msg += reason
        #else:
            #msg += strerror(gsl_err_code)
        #RuntimeError.__init__(self, msg)
        #self.gsl_err_code = gsl_err_code
        #self.fl = fl # file in which error occurred
        #self.line = line
        #self.reason = reason
        #self.result = result # result to be returned by the function

## helper function to test allocation
#def _gsl_check_null_pointer(p):
    #"""Raises an exception if a allocation failed in a GSL function."""
    #p = cast(p, c_void_p)
    #if p == 0:
        #raise GSL_Error(GSL_ENOMEM)
    #return p

#### ctypesGsl error handling
## handling of gsl status returned by functions
#def _ctypesGsl_status_exception(status_code, result):
    #raise GSL_Error(status_code, result = result)
#def _ctypesGsl_status_warning(status_code, result):
    #print "WARNING: " + str(GSL_Error(status_code, result = result))
    #return status_code
#def _ctypesGsl_status_off(status_code, result):
    #return status_code

## current status handler function
#ctypesGsl_status_handler = _ctypesGsl_status_exception

#def set_status_handler(h):
    #global ctypesGsl_status_handler
    #old_handler = ctypesGsl_status_handler
    #ctypesGsl_status_handler = h
    #return old_handler
#def set_status_handler_off():
    #return set_status_handler(_ctypesGsl_status_off)
#def set_status_handler_warning():
    #return set_status_handler(_ctypesGsl_status_warning)
#def set_status_handler_exception():
    #return set_status_handler(_ctypesGsl_status_exception)

## method for testing return values
#def _gsl_check_status(status_code, result = None):
    #"""Raises an exception if a GSL function returns an error
    #condition."""
    #if status_code != GSL_SUCCESS:
        #ctypesGsl_status_handler(status_code, result)
    #return status_code

#### internal gsl error handling
#GSL_ERROR_HANDLER_T = CFUNCTYPE(None, c_char_p, c_char_p, c_int, c_int)
#_add_function("set_error_handler", GSL_ERROR_HANDLER_T, [GSL_ERROR_HANDLER_T])
#_add_function("set_error_handler_off", GSL_ERROR_HANDLER_T, [])

## create our own error handler to raise exceptions instead of aborts
#def __ctypesGsl_error_handler_warning(reason, fl, line, gsl_errno):
    #print "WARNING: " + str(GSL_Error(gsl_errno, fl, line, reason))
    ## !!! currently, due to limitations in ctypes, this exception is
    ## !!! not thrown to the program, just printed
    ## raise GSL_Error(gsl_errno, fl, line, reason)
#def __ctypesGsl_error_handler_exception(reason, fl, line, gsl_errno):
    ## !!! currently, due to limitations in ctypes, this exception is
    ## !!! not thrown to the program, just printed
    #raise GSL_Error(gsl_errno, fl, line, reason)
#_ctypesGsl_error_handler_warning   = GSL_ERROR_HANDLER_T(__ctypesGsl_error_handler_warning)
#_ctypesGsl_error_handler_exception = GSL_ERROR_HANDLER_T(__ctypesGsl_error_handler_exception)

#def set_error_handler_warning():
    #return set_error_handler(_ctypesGsl_error_handler_warning)
#def set_error_handler_exception():
    #return set_error_handler(_ctypesGsl_error_handler_exception)

## set the default handler
#gsl_default_error_handler = set_error_handler_warning()
