# check SWIG Python module wrapping the LALXML library
# Author: Karl Wette, 2011

import sys, os

def msg(str):
    print(os.path.basename(__file__) + ": " + str)

# check that the module loads
try:
    sys.path.insert(0, os.getcwd())
    from swiglalxml import *
except:
    msg("FAILED module load")
    exit(1)

# passed all tests!
msg("================")
msg("PASSED all tests")
msg("================")
exit(0)
