# make Python extension from SWIG wrapping code
# Author: Karl Wette, 2011

import sys, os
from distutils import sysconfig
from distutils.core import setup, Extension
from distutils.command.build import build as _build
from distutils.command.clean import clean as _clean

# setup build directories
def setup_build_dir(self):
    self.build_base    = 'build'
    self.build_lib     = 'build/lib'
    self.build_scripts = 'build/scripts'
    self.build_temp    = 'build/temp'
    self.bdist_base    = 'build/bdist'

# subclass distutils build class
class build(_build):
    def initialize_options(self):
        _build.initialize_options(self)
        setup_build_dir(self)

# subclass distutils clean class
class clean(_clean):
    def initialize_options(self):
        _clean.initialize_options(self)
        setup_build_dir(self)

# remove these flags from default Python compiler flags:
noopts = set([
              '-DNDEBUG',             # handled elsewhere
              '-Wstrict-prototypes'   # not valid for C++
             ])
opts = sysconfig.get_config_var('OPT').split()
opts = [opt for opt in opts if not opt in noopts]
sysconfig.get_config_vars()['OPT'] = ' '.join(opts)

# module name and details
modname = os.environ['swig_wrapname']
modversion = os.environ['PACKAGE_VERSION']
moddesc = 'SWIG-generated wrapping of the ' + os.environ['PACKAGE_NAME'] + ' library'
modurl = 'https://www.lsc-group.phys.uwm.edu/daswg/projects/lalsuite.html'
modkeywords = {
}

# macro definitions
defines = map(lambda s : (s, None), os.environ['swig_defines'].split())

# include directories
incldirs = os.environ['swig_inclpath'].split()

# compile-time library directories
libdirs = os.environ['swig_libpath'].split()

# run-time library directory
rtlibdir = os.environ['swig_libdir']

# libraries to link against
libs = os.environ['swig_libs'].split()

# source file
srcfile = os.environ['swig_wrapfile']

# create extension class for SWIG module
swigmodule = Extension('_' + modname,
                       define_macros = defines,
                       include_dirs = incldirs,
                       library_dirs = libdirs,
                       runtime_library_dirs = [rtlibdir],
                       libraries = libs,
                       sources = [srcfile])

# run setup commands for SWIG module
setup(cmdclass = {'build':build, 'clean':clean},
      name = modname,
      version = modversion,
      description = moddesc,
      url = modurl,
      keywords = modkeywords,
      ext_modules = [swigmodule],
      py_modules = [modname])
