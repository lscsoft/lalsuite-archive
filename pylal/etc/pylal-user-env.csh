# Glue environment setup file
#
# $Id$
#
# source this file to set up your environment for glue

if ( ! $?PYLAL_LOCATION ) then
    echo "ERROR: environment variable PYLAL_LOCATION not defined"
    exit 1
endif

#python -V outputs to stderr, so have to redirect to catch it
setenv PYTHON_VERSION `python -V 2>&1 | cut -d" " -f2 | cut -d. -f1-2`

setenv mylibdir "lib"
if(`uname -p` == "x86_64") then
  setenv mylibdir "lib64"
endif

if ( ! $?LD_LIBRARY_PATH ) then
    setenv LD_LIBRARY_PATH ''
endif

if ( $?PYLAL_PATH ) then
    setenv PATH `echo "${PATH}" | sed -e "s%:${PYLAL_PATH}[^:]*%%g" -e "s%^${PYLAL_PATH}[^:]*:\{0,1\}%%"`
    setenv PYTHONPATH `echo "${PYTHONPATH}" | sed -e "s%:${PYLAL_PATH}[^:]*%%g" -e "s%^${PYLAL_PATH}[^:]*:\{0,1\}%%"`
    setenv LD_LIBRARY_PATH `echo "${LD_LIBRARY_PATH}" | sed -e "s%:${PYLAL_PATH}[^:]*%%g" -e "s%^${PYLAL_PATH}[^:]*:\{0,1\}%%"`
    if ( $?MANPATH ) then
        setenv MANPATH `echo "${MANPATH}" | sed -e "s%:${PYLAL_PATH}[^:]*%%g" -e "s%^${PYLAL_PATH}[^:]*:\{0,1\}%%"`
    endif
endif

setenv PATH `echo "${PATH}" | sed -e "s%:${PYLAL_LOCATION}[^:]*%%g" -e "s%^${PYLAL_LOCATION}[^:]*:\{0,1\}%%"`
setenv PYTHONPATH `echo "${PYTHONPATH}" | sed -e "s%:${PYLAL_LOCATION}[^:]*%%g" -e "s%^${PYLAL_LOCATION}[^:]*:\{0,1\}%%"`
setenv LD_LIBRARY_PATH `echo "${LD_LIBRARY_PATH}" | sed -e "s%:${PYLAL_LOCATION}[^:]*%%g" -e "s%^${PYLAL_LOCATION}[^:]*:\{0,1\}%%"`
if ( $?MANPATH ) then
    setenv MANPATH `echo "${MANPATH}" | sed -e "s%:${PYLAL_LOCATION}[^:]*%%g" -e "s%^${PYLAL_LOCATION}[^:]*:\{0,1\}%%"`
endif

setenv PYLAL_PATH "${PYLAL_LOCATION}"
setenv PATH "${PYLAL_LOCATION}/bin:${PYLAL_LOCATION}/sbin:${PATH}";
setenv PYTHONPATH "${PYLAL_LOCATION}/${mylibdir}/python$PYTHON_VERSION/site-packages:${PYTHONPATH}";

if ( $?MANPATH ) then
    set DELIM
    if ( "X${MANPATH}" != "X" ) then
        set DELIM=:
    endif
    setenv MANPATH "${PYLAL_LOCATION}/man${DELIM}${MANPATH}"
endif

set DELIM=
if ( "X${LD_LIBRARY_PATH}" != "X" ) then
    set DELIM=:
endif
setenv LD_LIBRARY_PATH "${PYLAL_LOCATION}/${mylibdir}${DELIM}${LD_LIBRARY_PATH}"

unsetenv PYTHON_VERSION
unsetenv mylibdir
