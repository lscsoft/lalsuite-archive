## check SWIG Octave module wrapping the LALXML library
## Author: Karl Wette, 2011

function msg(varargin)
  disp([program_name, ": ", sprintf(varargin{:})]);
endfunction

## check that the module loads
try
  addpath(pwd(), 0);
  swiglalxml;
catch
  msg("FAILED module load");
  exit(1);
end_try_catch

## passed all tests!
msg("================");
msg("PASSED all tests");
msg("================");
exit(0)
