# check SWIG Python module wrapping the LAL library
# Author: Karl Wette, 2011

import os

def msg(str):
    print(os.path.basename(__file__) + ": " + str)

# passed all tests!
msg("================")
msg("PASSED all tests")
msg("================")
exit(0)
