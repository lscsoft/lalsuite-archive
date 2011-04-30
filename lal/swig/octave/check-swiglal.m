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

## check string conversions
try
  strs = {"a", "bc", "def"};
  sv = XLALCreateStringVector(strs{:});
  assert(sv.length == 3);
  assert(all(strcmp(sv.data, strs)));
  strs{1} = "ghijk";
  sv.data_setel(0, strs{1});
  strs{end+1} = "lmnopq";
  XLALAppendString2Vector(sv, strs{4});
  for i = 1:sv.length
    assert(strcmp(sv.data_getel(i-1), strs{i}));
  endfor
  XLALDestroyStringVector(sv);
  LALCheckMemoryLeaks();
catch
  msg("FAILED string conversions");
  exit(1);
end_try_catch

## check static vector/matrix conversions
if !cvar.swiglal_debug
  msg("skipping static vector/matrix conversions");
else
  try
    sts = new_swiglal_static_test_struct();
    assert(length(sts.vector) == 3);
    assert(length(sts.enum_vector) == 3);
    assert(all(size(sts.matrix) == [2, 3]));
    assert(all(size(sts.enum_matrix) == [2, 3]));
    sts.vector = [3, 2, 1];
    assert(all(sts.vector == [3, 2, 1]));
    sts.matrix = [4, 5, 6; 9, 8, 7];
    try
      sts.matrix = [1.1, 2.3, 4.5; 6.5, 4.3, 2.1];
      msg("FAILED static vector/matrix conversions #2");
      exit(1);
    end_try_catch
    assert(all(sts.matrix == [4, 5, 6; 9, 8, 7]));
    for i = 0:2
      sts.enum_vector_setel(i, 2*i + 3);
      assert(sts.enum_vector_getel(i) == (2*i + 3));
    endfor
    clear sts;
    assert(!any(cvar.swiglal_static_test_vector));
    assert(!any(cvar.swiglal_static_test_matrix(:)));
    assert(!any(cvar.swiglal_static_test_enum_vector));
    assert(!any(cvar.swiglal_static_test_enum_matrix(:)));
    cvar.swiglal_static_test_vector = cvar.swiglal_static_test_const_vector;
    assert(all(cvar.swiglal_static_test_vector == [1, 2, 4]));
    assert(swiglal_static_test_const_vector_getel(2) == 4);
    try
      swiglal_static_test_const_vector_getel(20);
      msg("FAILED static vector/matrix conversions #3");
      exit(1);
    end_try_catch
  catch
    msg("FAILED static vector/matrix conversions #1");
    exit(1);
  end_try_catch
endif

## passed all tests!
msg("================");
msg("PASSED all tests");
msg("================");
exit(0)
