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

# check complex number conversions
try:
    args = [complex(1, 3), complex(2, -5), complex(3, -2)]
    assert(XLALCOMPLEX8Add(args[0], args[1]) == args[2])
    args = [complex(4, 2), complex(10, 5), complex(-6, -3)]
    assert(XLALCOMPLEX8Sub(args[0], args[1]) == args[2])
    args = [complex(10, -7), complex(30, 17), complex(419, -40)]
    assert(XLALCOMPLEX16Mul(args[0], args[1]) == args[2])
    args = [complex(111.75, -120.50), complex(5, 12), complex(-5.25, -11.5)]
    assert(XLALCOMPLEX16Div(args[0], args[1]) == args[2])
    LALCheckMemoryLeaks()
except:
    msg("FAILED complex number conversions")
    exit(1)

# check string conversions
try:
    strs = ["a", "bc", "def"]
    sv = XLALCreateStringVector(*strs)
    assert(sv.length == 3)
    assert(all(sv.data == strs))
    strs[0] = "ghijk"
    sv.data_setel(0, strs[0])
    strs.append("lmnopq")
    XLALAppendString2Vector(sv, strs[3])
    for i in range(0, sv.length):
        assert(sv.data_getel(i) == strs[i])
    XLALDestroyStringVector(sv)
    LALCheckMemoryLeaks()
except:
    msg("FAILED string conversions")
    exit(1)

# check static vector/matrix conversions
if not cvar.swiglal_debug:
    msg("skipping static vector/matrix conversions")
else:
    try:
        sts = swiglal_static_test_struct()
        assert(len(sts.vector) == 3)
        assert(len(sts.enum_vector) == 3)
        assert(sts.matrix.shape == (2, 3))
        assert(sts.enum_matrix.shape == (2, 3))
        sts.vector = [3, 2, 1]
        assert(all(sts.vector == [3, 2, 1]))
        sts.matrix = [[4, 5, 6], (9, 8, 7)]
        try:
            sts.matrix = [[1.1, 2.3, 4.5], [6.5, 4.3, 2.1]]
            msg("FAILED static vector/matrix conversions #2")
            exit(1)
        except:
            pass
        assert((sts.matrix == [[4, 5, 6], [9, 8, 7]]).all())
        for i in range(0, 3):
            sts.enum_vector_setel(i, 2*i + 3)
            assert(sts.enum_vector_getel(i) == (2*i + 3))
        del sts
        assert(not any(cvar.swiglal_static_test_vector))
        assert(not cvar.swiglal_static_test_matrix.any())
        assert(not any(cvar.swiglal_static_test_enum_vector))
        assert(not cvar.swiglal_static_test_enum_matrix.any())
        cvar.swiglal_static_test_vector = cvar.swiglal_static_test_const_vector
        assert(all(cvar.swiglal_static_test_vector == [1, 2, 4]))
        assert(swiglal_static_test_const_vector_getel(2) == 4)
        try:
            swiglal_static_test_const_vector_getel(20)
            msg("FAILED static vector/matrix conversions #3")
            exit(1)
        except:
            pass
    except:
        msg("FAILED static vector/matrix conversions #1")
        exit(1)

# passed all tests!
msg("================")
msg("PASSED all tests")
msg("================")
exit(0)
