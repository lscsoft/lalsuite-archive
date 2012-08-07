import os
import fnmatch
from ctypes import CDLL
from ctypes.util import find_library

class PkgConfig(object):
    def __init__(self, names):
        def stripfirsttwo(string):
            return string[2:]
        self.libs = map(stripfirsttwo, os.popen("pkg-config --libs-only-l %s" % names).read().split())
        self.libdirs = map(stripfirsttwo, os.popen("pkg-config --libs-only-L %s" % names).read().split())
        self.incdirs = map(stripfirsttwo, os.popen("pkg-config --cflags-only-I %s" % names).read().split())
        self.extra_cflags = os.popen("pkg-config --cflags-only-other %s" % names).read().split()


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

libc = __load_lib("c",mode=RTLD_GLOBAL)
liblal = __load_lib("lal",mode=RTLD_GLOBAL)
liblalinference = __load_lib("lalinference",mode=RTLD_GLOBAL)
