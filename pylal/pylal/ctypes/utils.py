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
    

# load the libraries
def __load_lib(libname, mode = None):
    lib_path = find_library(libname)
    if lib_path is None:
        pkg_config=PkgConfig(libname)
        targets=[]
        
        if pkg_config.libdirs:
            lib_path=pkg_config.libdirs[0]
            for f in os.listdir(lib_path):
                if fnmatch.fnmatch(f,'lib'+libname+'.so.*') or fnmatch.fnmatch(f,'lib'+libname+'.dylib*'):
                    targets.append(f)
        
        if not targets:
            lib_path=None
        else:
            targets.sort()
            lib_path+='/'+targets[0]
            
        #print pkg_config.libs
        if lib_path is None:
            raise RuntimeError("Could not find the " + libname + " library")
    
    if mode is None:
        
        
        lib = CDLL(lib_path)
    else:
        lib = CDLL(lib_path, mode)
    return lib

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
