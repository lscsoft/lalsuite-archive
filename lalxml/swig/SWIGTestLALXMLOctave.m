## Check SWIG Octave bindings for LALXML
## Author: Karl Wette, 2011--2014

page_screen_output(0);
crash_dumps_octave_core(0);

## check module load
disp("checking module load ...");
lalxml;
assert(exist("lalxml", "var"));
lal;
assert(exist("lal", "var"));
disp("PASSED module load");

## check object parent tracking
disp("checking object parent tracking ...");
a = lalxml.new_swig_lalxml_test_parent_map_struct();
for i = 1:7
  b = a.s;
  c = lalxml.swig_lalxml_test_parent_map.s;
  lalxml.swig_lalxml_test_parent_map.s = lal.swig_lal_test_struct_const;
endfor
clear c;
clear b;
clear a;
clear ans;
LALCheckMemoryLeaks();
disp("PASSED object parent tracking");

## passed all tests!
disp("PASSED all tests");
