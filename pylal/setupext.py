"""
Some helper functions for building the C extensions

you may need to edit basedir to point to the default location of your
required libs, eg, png, z, freetype

"""

import os

basedir = {
    'win32'  : ['win32_static',],
    'linux2' : ['/usr/local', '/usr',],
    'linux'  : ['/usr/local', '/usr',],
    'darwin' : ['/usr/local', '/usr', '/sw', '/usr/X11R6'],
    'freebsd4' : ['/usr/local', '/usr'],
    'sunos5' : [os.getenv('MPLIB_BASE') or '/usr/local',],
}

import sys, os, stat
from distutils.core import Extension
import glob

major, minor1, minor2, s, tmp = sys.version_info
if major<2 or (major==2 and minor1<3):
    True = 1
    False = 0
else:
    True = True
    False = False

try:
  lal_location = os.environ['LAL_LOCATION']
except:
  print >>sys.stderr, "Error getting environment LAL_LOCATION"
  print >>sys.stderr, "Please set LAL_LOCATION environment variable"
  sys.exit(1)

lal_include = lal_location + '/include'
lal_lib = lal_location + '/lib'

BUILT_METAIO   = False

class CleanUpFile:
    """CleanUpFile deletes the specified filename when self is destroyed."""
    def __init__(self, name):
        self.name = name
    def __del__(self):
        os.remove(self.name)

def temp_copy(_from, _to):
    """temp_copy copies a named file into a named temporary file.
    The temporary will be deleted when the setupext module is destructed.
    """
    # Copy the file data from _from to _to
    s = open(_from).read()
    open(_to,"w+").write(s)
    # Suppress object rebuild by preserving time stamps.
    stats = os.stat(_from)
    os.utime(_to, (stats.st_atime, stats.st_mtime))
    # Make an object to eliminate the temporary file at exit time.
    globals()["_cleanup_"+_to] = CleanUpFile(_to)

def add_base_flags(module):
    incdirs = [os.path.join(p, 'include') for p in basedir[sys.platform]
               if os.path.exists(p)]
    libdirs = [os.path.join(p, 'lib')     for p in basedir[sys.platform]
               if os.path.exists(p)]
    module.include_dirs.extend(incdirs)
    module.library_dirs.extend(libdirs)


def getoutput(s):
    'get the output of a system command'

    ret =  os.popen(s).read().strip()
    return ret


def build_metaio(ext_modules, packages):
    global BUILT_METAIO
    if BUILT_METAIO: return # only build it if you you haven't already


    module = Extension( 'lgen.metaio', ['src/metaio.c'], 
        libraries = ['stdc++', 'lal', 'lalmetaio', 'lalsupport',
          'metaio'], 
        include_dirs = [lal_include, 
          '/opt/lscsoft/libmetaio/include'], 
        library_dirs = [lal_lib,
        '/opt/lscsoft/libmetaio/lib'])    
    ext_modules.append(module)    

    BUILT_METAIO = True
