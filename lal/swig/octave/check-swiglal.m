## check SWIG Octave module wrapping the LAL library
## Author: Karl Wette, 2011

function msg(varargin)
  disp([program_name, ": ", sprintf(varargin{:})]);
endfunction

## passed all tests!
msg("================");
msg("PASSED all tests");
msg("================");
exit(0)
