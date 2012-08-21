# Source this file to access GLUE
#GLUE_PREFIX=/opt/lscsoft/glue
GLUE_PREFIX=/usr
export GLUE_PREFIX
# dynamically get python version (Carsten Aulbert)
PYSITE_PATH=python`python -V 2>&1 | cut -d' ' -f2 | cut -d. -f-2`
PATH=${GLUE_PREFIX}/bin:${PATH}
PYTHONPATH=${GLUE_PREFIX}/lib/${PYSITE_PATH}/site-packages:${PYTHONPATH}
LD_LIBRARY_PATH=${GLUE_PREFIX}/lib/${PYSITE_PATH}/site-packages:${LD_LIBRARY_PATH}
DYLD_LIBRARY_PATH=${GLUE_PREFIX}/lib/${PYSITE_PATH}/site-packages:${DYLD_LIBRARY_PATH}
export PATH PYTHONPATH LD_LIBRARY_PATH DYLD_LIBRARY_PATH
