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

## passed all tests!
msg("================");
msg("PASSED all tests");
msg("================");
exit(0)
