# Glue environment setup file
#
# $Id$
#
# source this file to set up your environment for glue

if ( ! $?GLUE_LOCATION ) then
    echo "ERROR: environment variable GLUE_LOCATION not defined"
    exit 1
endif

if ( ! $?LD_LIBRARY_PATH ) then
    setenv LD_LIBRARY_PATH ''
endif

if ( $?GLUE_PATH ) then
    setenv PATH `echo "${PATH}" | sed -e "s%:${GLUE_PATH}[^:]*%%g" -e "s%^${GLUE_PATH}[^:]*:\{0,1\}%%"`
    setenv PYTHONPATH `echo "${PYTHONPATH}" | sed -e "s%:${GLUE_PATH}[^:]*%%g" -e "s%^${GLUE_PATH}[^:]*:\{0,1\}%%"`
    setenv LD_LIBRARY_PATH `echo "${LD_LIBRARY_PATH}" | sed -e "s%:${GLUE_PATH}[^:]*%%g" -e "s%^${GLUE_PATH}[^:]*:\{0,1\}%%"`
    if ( $?MANPATH ) then
        setenv MANPATH `echo "${MANPATH}" | sed -e "s%:${GLUE_PATH}[^:]*%%g" -e "s%^${GLUE_PATH}[^:]*:\{0,1\}%%"`
    endif
endif

setenv PATH `echo "${PATH}" | sed -e "s%:${GLUE_LOCATION}[^:]*%%g" -e "s%^${GLUE_LOCATION}[^:]*:\{0,1\}%%"`
setenv PYTHONPATH `echo "${PYTHONPATH}" | sed -e "s%:${GLUE_LOCATION}[^:]*%%g" -e "s%^${GLUE_LOCATION}[^:]*:\{0,1\}%%"`
setenv LD_LIBRARY_PATH `echo "${LD_LIBRARY_PATH}" | sed -e "s%:${GLUE_LOCATION}[^:]*%%g" -e "s%^${GLUE_LOCATION}[^:]*:\{0,1\}%%"`
if ( $?MANPATH ) then
    setenv MANPATH `echo "${MANPATH}" | sed -e "s%:${GLUE_LOCATION}[^:]*%%g" -e "s%^${GLUE_LOCATION}[^:]*:\{0,1\}%%"`
endif

setenv GLUE_PATH "${GLUE_LOCATION}"
setenv PATH "${GLUE_LOCATION}/bin:${GLUE_LOCATION}/sbin:${PATH}";
setenv PYTHONPATH "${GLUE_LOCATION}/lib/python:${PYTHONPATH}";

if ( $?MANPATH ) then
    set DELIM
    if ( "X${MANPATH}" != "X" ) then
        set DELIM=:
    endif
    setenv MANPATH "${GLUE_LOCATION}/man${DELIM}${MANPATH}"
endif

set DELIM=
if ( "X${LD_LIBRARY_PATH}" != "X" ) then
    set DELIM=:
endif
setenv LD_LIBRARY_PATH "${GLUE_LOCATION}/lib${DELIM}${LD_LIBRARY_PATH}"
