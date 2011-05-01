# check SWIG Python module wrapping the LAL library
# Author: Karl Wette, 2011

import os

def msg(str):
    print(os.path.basename(__file__) + ": " + str)

# check memory allocation
if not cvar.swiglal_debug:
    msg("skipping memory allocation")
else:
    try:
        LALCheckMemoryLeaks()   # should NOT throw an exception
        mem1 = LALDetector()
        mem2 = LALStringVector()
        mem3 = COMPLEX8Vector()
        mem4 = XLALCreateREAL8Vector(3)
        try:
            msg("*** below should be an error message from LALCheckMemoryLeaks() ***");
            LALCheckMemoryLeaks()   # should throw an exception
            msg("FAILED memory allocation #2")
            exit(1)
        except:
            msg("*** above should be an error message from LALCheckMemoryLeaks() ***");
            del mem1, mem2, mem3
            XLALDestroyREAL8Vector(mem4)
            try:
                LALCheckMemoryLeaks()   # should NOT throw an exception
            except:
                msg("FAILED memory allocation #3")
                exit(1)
    except:
        msg("FAILED memory allocation #1")
        exit(1)

# passed all tests!
msg("================")
msg("PASSED all tests")
msg("================")
exit(0)
