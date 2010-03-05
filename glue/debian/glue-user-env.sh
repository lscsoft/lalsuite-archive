# Source this file to access GLUE
GLUE_PREFIX=/opt/lscsoft/glue
export GLUE_PREFIX
# dynamically get python version (Carsten Aulbert)
PYSITE_PATH=python`python -V 2>&1 | cut -d' ' -f2 | cut -d. -f-2`
PATH=/opt/lscsoft/glue/bin:${PATH}
PYTHONPATH=/opt/lscsoft/glue/lib/${PYSITE_PATH}/site-packages:${PYTHONPATH}
LD_LIBRARY_PATH=/opt/lscsoft/glue/lib/${PYSITE_PATH}/site-packages:${LD_LIBRARY_PATH}
DYLD_LIBRARY_PATH=/opt/lscsoft/glue/lib/${PYSITE_PATH}/site-packages:${DYLD_LIBRARY_PATH}
export PATH PYTHONPATH LD_LIBRARY_PATH DYLD_LIBRARY_PATH
# Source this file to access GLUE
GLUE_PREFIX=/opt/lscsoft/glue
export GLUE_PREFIX
# dynamically get python version (Carsten Aulbert)
PYSITE_PATH=python`python -V 2>&1 | cut -d' ' -f2 | cut -d. -f-2`
PATH=/opt/lscsoft/glue/bin:${PATH}
PYTHONPATH=/opt/lscsoft/glue/lib/${PYSITE_PATH}/site-packages:${PYTHONPATH}
LD_LIBRARY_PATH=/opt/lscsoft/glue/lib/${PYSITE_PATH}/site-packages:${LD_LIBRARY_PATH}
DYLD_LIBRARY_PATH=/opt/lscsoft/glue/lib/${PYSITE_PATH}/site-packages:${DYLD_LIBRARY_PATH}
export PATH PYTHONPATH LD_LIBRARY_PATH DYLD_LIBRARY_PATH
