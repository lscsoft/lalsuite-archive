# Glue environment setup file
#
# $Id$
#
# source this file to set up your environment for glue

if [ -z "${GLUE_LOCATION}" ]; then
    echo "ERROR: environment variable GLUE_LOCATION not defined"  1>&2
    return 1
fi

PATH=`echo "${PATH}" | sed -e "s%:${GLUE_LOCATION}[^:]*%%g" -e "s%^${GLUE_LOCATION}[^:]*:\{0,1\}%%"`
LD_LIBRARY_PATH=`echo "${LD_LIBRARY_PATH}" | sed -e "s%:${GLUE_LOCATION}[^:]*%%g" -e "s%^${GLUE_LOCATION}[^:]*:\{0,1\}%%"`
if [ -n "${MANPATH}" ]; then
if [ -n "${MANPATH}" ]; then
    MANPATH=`echo "${MANPATH}" | sed -e "s%:${GLUE_LOCATION}[^:]*%%g" -e "s%^${GLUE_LOCATION}[^:]*:\{0,1\}%%"`
fi
PYTHONPATH=`echo "${PYTHONPATH}" | sed -e "s%:${GLUE_LOCATION}[^:]*%%g" -e "s%^${GLUE_LOCATION}[^:]*:\{0,1\}%%"`

GLUE_PATH=${GLUE_LOCATION}
PATH="${GLUE_LOCATION}/bin:${GLUE_LOCATION}/sbin:${PATH}";
PYTHONPATH="${GLUE_LOCATION}/lib/python:${PYTHONPATH}"

if [ -n "${MANPATH}" ]; then
    MANPATH="${GLUE_LOCATION}/man:${MANPATH}"
fi

DELIM=
if [ -n "${LD_LIBRARY_PATH}" ]; then
    DELIM=:
fi
LD_LIBRARY_PATH="${GLUE_LOCATION}/lib${DELIM}${LD_LIBRARY_PATH}"

export GLUE_PATH PATH MANPATH LD_LIBRARY_PATH PYTHONPATH
