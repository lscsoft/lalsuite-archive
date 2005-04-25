# Glue environment setup file
#
# $Id$
#
# source this file to set up your environment for glue

if ( ! $?LGT_LOCATION ) then
    echo "ERROR: environment variable LGT_LOCATION not defined"
    exit 1
endif

if ( ! $?LD_LIBRARY_PATH ) then
    setenv LD_LIBRARY_PATH ''
endif

if ( $?LGT_PATH ) then
    setenv PATH `echo "${PATH}" | sed -e "s%:${LGT_PATH}[^:]*%%g" -e "s%^${LGT_PATH}[^:]*:\{0,1\}%%"`
    setenv PYTHONPATH `echo "${PYTHONPATH}" | sed -e "s%:${LGT_PATH}[^:]*%%g" -e "s%^${LGT_PATH}[^:]*:\{0,1\}%%"`
    setenv LD_LIBRARY_PATH `echo "${LD_LIBRARY_PATH}" | sed -e "s%:${LGT_PATH}[^:]*%%g" -e "s%^${LGT_PATH}[^:]*:\{0,1\}%%"`
    if ( $?MANPATH ) then
        setenv MANPATH `echo "${MANPATH}" | sed -e "s%:${LGT_PATH}[^:]*%%g" -e "s%^${LGT_PATH}[^:]*:\{0,1\}%%"`
    endif
endif

setenv PATH `echo "${PATH}" | sed -e "s%:${LGT_LOCATION}[^:]*%%g" -e "s%^${LGT_LOCATION}[^:]*:\{0,1\}%%"`
setenv PYTHONPATH `echo "${PYTHONPATH}" | sed -e "s%:${LGT_LOCATION}[^:]*%%g" -e "s%^${LGT_LOCATION}[^:]*:\{0,1\}%%"`
setenv LD_LIBRARY_PATH `echo "${LD_LIBRARY_PATH}" | sed -e "s%:${LGT_LOCATION}[^:]*%%g" -e "s%^${LGT_LOCATION}[^:]*:\{0,1\}%%"`
if ( $?MANPATH ) then
    setenv MANPATH `echo "${MANPATH}" | sed -e "s%:${LGT_LOCATION}[^:]*%%g" -e "s%^${LGT_LOCATION}[^:]*:\{0,1\}%%"`
endif

setenv LGT_PATH "${LGT_LOCATION}"
setenv PATH "${LGT_LOCATION}/bin:${LGT_LOCATION}/sbin:${PATH}";
setenv PYTHONPATH "${LGT_LOCATION}/lib/python:${PYTHONPATH}";

if ( $?MANPATH ) then
    set DELIM
    if ( "X${MANPATH}" != "X" ) then
        set DELIM=:
    endif
    setenv MANPATH "${LGT_LOCATION}/man${DELIM}${MANPATH}"
endif

set DELIM=
if ( "X${LD_LIBRARY_PATH}" != "X" ) then
    set DELIM=:
endif
setenv LD_LIBRARY_PATH "${LGT_LOCATION}/lib${DELIM}${LD_LIBRARY_PATH}"
