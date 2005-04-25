# Glue environment setup file
#
# $Id$
#
# source this file to set up your environment for glue

if [ -z "${LGT_LOCATION}" ]; then
    echo "ERROR: environment variable LGT_LOCATION not defined"  1>&2
    return 1
fi

PATH=`echo "${PATH}" | sed -e "s%:${LGT_LOCATION}[^:]*%%g" -e "s%^${LGT_LOCATION}[^:]*:\{0,1\}%%"`
PYTHONPATH=`echo "${PYTHONPATH}" | sed -e "s%:${LGT_LOCATION}[^:]*%%g" -e "s%^${LGT_LOCATION}[^:]*:\{0,1\}%%"`
LD_LIBRARY_PATH=`echo "${LD_LIBRARY_PATH}" | sed -e "s%:${LGT_LOCATION}[^:]*%%g" -e "s%^${LGT_LOCATION}[^:]*:\{0,1\}%%"`

if [ -n "${MANPATH}" ]; then
    MANPATH=`echo "${MANPATH}" | sed -e "s%:${LGT_LOCATION}[^:]*%%g" -e "s%^${LGT_LOCATION}[^:]*:\{0,1\}%%"`
fi

LGT_PATH=${LGT_LOCATION}
PATH="${LGT_LOCATION}/bin:${LGT_LOCATION}/sbin:${PATH}";
PYTHONPATH="${LGT_LOCATION}/lib/python:${PYTHONPATH}"

if [ -n "${MANPATH}" ]; then
    MANPATH="${LGT_LOCATION}/man:${MANPATH}"
fi

DELIM=
if [ -n "${LD_LIBRARY_PATH}" ]; then
    DELIM=:
fi
LD_LIBRARY_PATH="${LGT_LOCATION}/lib${DELIM}${LD_LIBRARY_PATH}"

export LGT_PATH PATH MANPATH LD_LIBRARY_PATH PYTHONPATH
