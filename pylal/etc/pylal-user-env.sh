# Glue environment setup file
#
# $Id$
#
# source this file to set up your environment for glue

if [ -z "${PYLAL_LOCATION}" ]; then
    echo "ERROR: environment variable PYLAL_LOCATION not defined"  1>&2
    return 1
fi

PATH=`echo "${PATH}" | sed -e "s%:${PYLAL_LOCATION}[^:]*%%g" -e "s%^${PYLAL_LOCATION}[^:]*:\{0,1\}%%"`
PYTHONPATH=`echo "${PYTHONPATH}" | sed -e "s%:${PYLAL_LOCATION}[^:]*%%g" -e "s%^${PYLAL_LOCATION}[^:]*:\{0,1\}%%"`
LD_LIBRARY_PATH=`echo "${LD_LIBRARY_PATH}" | sed -e "s%:${PYLAL_LOCATION}[^:]*%%g" -e "s%^${PYLAL_LOCATION}[^:]*:\{0,1\}%%"`

if [ -n "${MANPATH}" ]; then
    MANPATH=`echo "${MANPATH}" | sed -e "s%:${PYLAL_LOCATION}[^:]*%%g" -e "s%^${PYLAL_LOCATION}[^:]*:\{0,1\}%%"`
fi

LGT_PATH=${PYLAL_LOCATION}
PATH="${PYLAL_LOCATION}/bin:${PYLAL_LOCATION}/sbin:${PATH}";
PYTHONPATH="${PYLAL_LOCATION}/lib/python:${PYTHONPATH}"

if [ -n "${MANPATH}" ]; then
    MANPATH="${PYLAL_LOCATION}/man:${MANPATH}"
fi

DELIM=
if [ -n "${LD_LIBRARY_PATH}" ]; then
    DELIM=:
fi
LD_LIBRARY_PATH="${PYLAL_LOCATION}/lib${DELIM}${LD_LIBRARY_PATH}"

export LGT_PATH PATH MANPATH LD_LIBRARY_PATH PYTHONPATH
