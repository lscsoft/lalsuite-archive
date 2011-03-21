## check SWIG Octave module wrapping the LAL library
## Author: Karl Wette, 2011

function msg(varargin)
  disp([program_name, ": ", sprintf(varargin{:})]);
endfunction

## check that the module loads
try
  addpath(pwd(), 0);
  swiglal;
catch
  msg("FAILED module load");
  exit(1);
end_try_catch
swiglal.cvar.lalDebugLevel = 1;

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

## check dynamic vector/matrix conversions
function check_dynamic_vector_matrix(id, iv, ivl, rv, rvl, cm, cms1, cms2)
  try
    assert(ivl == 5);
    iv.data = [1, 3, 2, 4, 3];
    assert(all(iv.data == [1, 3, 2, 4, 3]));
    iv.data_setel(3, 7);
    assert(iv.data_getel(3) == 7);
    assert(rvl == 5);
    rv.data = [1.2, 3.4, 2.6, 4.8, 3.5];
    assert(all(rv.data == [1.2, 3.4, 2.6, 4.8, 3.5]));
    rv.data_setel(rvl - 1, 7.5);
    assert(rv.data_getel(rvl - 1) == 7.5);
    try
      rv.data_setel(rvl, 99.9);
      msg("FAILED dynamic vector/matrix conversions #2 (%s)", id);
      exit(1);
    end_try_catch
    try
      iv.data = rv.data;
      msg("FAILED dynamic vector/matrix conversions #3 (%s)", id);
      exit(1);
    end_try_catch
    try
      rv.data = iv.data;
      assert(all(rv.data == iv.data));
    catch
      msg("FAILED dynamic vector/matrix conversions #4 (%s)", id);
      exit(1);
    end_try_catch
    assert(cms1 == 4);
    assert(cms2 == 6);
    for i = 0:cms1 - 1
      for j = 0:cms2 - 1
        cm.data_setel(i, j, complex(i / 4.0, j / 2.0));
      endfor
    endfor
    assert(cm.data_getel(2, 3) == complex(0.5, 1.5));
    assert(cm.data_getel(3, 2) == complex(0.75, 1.0));
    try
      iv.data_setel(0, cm.data_getel(2, 3));
      msg("FAILED dynamic vector/matrix conversions #5 (%s)", id);
      exit(1);
    end_try_catch
    try
      rv.data_setel(0, cm.data_getel(3, 2));
      msg("FAILED dynamic vector/matrix conversions #6 (%s)", id);
      exit(1);
    end_try_catch
  catch
    msg("FAILED dynamic vector/matrix conversions #1 (%s)", id);
    exit(1);
  end_try_catch
endfunction
try
  # check LAL vector and matrix datatypes
  iv = XLALCreateINT4Vector(5);
  rv = XLALCreateREAL8Vector(5);
  cm = XLALCreateCOMPLEX8VectorSequence(4, 6);
  check_dynamic_vector_matrix("LAL", iv, iv.length, rv, rv.length,
                              cm, cm.length, cm.vectorLength);
  XLALDestroyINT4Vector(iv);
  XLALDestroyREAL8Vector(rv);
  XLALDestroyCOMPLEX8VectorSequence(cm);
  LALCheckMemoryLeaks();
  ## check GSL vectors and matrices
  iv = gsl_vector_int_calloc(5);
  rv = gsl_vector_calloc(5);
  cm = gsl_matrix_complex_float_calloc(4, 6);
  check_dynamic_vector_matrix("GSL", iv, iv.size, rv, rv.size,
                              cm, cm.size1, cm.size2);
  gsl_vector_int_free(iv);
  gsl_vector_free(rv);
  gsl_matrix_complex_float_free(cm);
catch
  msg("FAILED dynamic vector/matrix conversions");
  exit(1);
end_try_catch

## check 'tm' struct conversions
try
  gps = 989168284;
  utc = [2011, 5, 11, 16, 57, 49, 4, 131, 0];
  assert(all(XLALGPSToUTC([], gps) == utc));
  assert(XLALUTCToGPS(utc) == gps);
  assert(XLALUTCToGPS(utc(1:6)) == gps);
  utc(7) = utc(8) = 0;
  for f = [-1, 0, 1]
    utc(9) = f;
    assert(XLALUTCToGPS(utc) == gps);
  endfor
  utcd = utc;
  for d = 0:10
    utcd(3) = utc(3) + d;
    utcd = XLALGPSToUTC([], XLALUTCToGPS(utcd));
    assert(utcd(7) == weekday(datenum(utcd(1:6))));
  endfor
catch
  msg("FAILED 'tm' struct conversions");
  exit(1);
end_try_catch

## passed all tests!
msg("================");
msg("PASSED all tests");
msg("================");
exit(0)
