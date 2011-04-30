## check SWIG Octave module wrapping the LAL library
## Author: Karl Wette, 2011

function msg(varargin)
  disp([program_name, ": ", sprintf(varargin{:})]);
endfunction

## check memory allocation
if !cvar.swiglal_debug
  msg("skipping memory allocation");
else
  try
    LALCheckMemoryLeaks();   # should NOT throw an exception
    mem1 = new_LALDetector();
    mem2 = new_LALStringVector();
    mem3 = new_COMPLEX8Vector();
    mem4 = XLALCreateREAL8Vector(3);
    try
      msg("*** below should be an error message from LALCheckMemoryLeaks() ***");
      LALCheckMemoryLeaks();   # should throw an exception
      msg("FAILED memory allocation #2");
      exit(1);
    catch
      msg("*** above should be an error message from LALCheckMemoryLeaks() ***");
      clear mem1 mem2 mem3;
      XLALDestroyREAL8Vector(mem4);
      try
        LALCheckMemoryLeaks();   # should NOT throw an exception
      catch
        msg("FAILED memory allocation #3");
        exit(1);
      end_try_catch
    end_try_catch
  catch
    msg("FAILED memory allocation #1");
    exit(1);
  end_try_catch
endif

## check complex number conversions
try
  args = [complex(1, 3), complex(2, -5), complex(3, -2)];
  assert(XLALCOMPLEX8Add(args(1), args(2)) == args(3));
  args = [complex(4, 2), complex(10, 5), complex(-6, -3)];
  assert(XLALCOMPLEX8Sub(args(1), args(2)) == args(3));
  args = [complex(10, -7), complex(30, 17), complex(419, -40)];
  assert(XLALCOMPLEX16Mul(args(1), args(2)) == args(3));
  args = [complex(111.75, -120.50), complex(5, 12), complex(-5.25, -11.5)];
  assert(XLALCOMPLEX16Div(args(1), args(2)) == args(3));
  LALCheckMemoryLeaks();
catch
  msg("FAILED complex number conversions");
  exit(1);
end_try_catch

## passed all tests!
msg("================");
msg("PASSED all tests");
msg("================");
exit(0)
