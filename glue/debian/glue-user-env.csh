# Source this file to access GLUE
setenv GLUE_PREFIX /opt/lscsoft/glue
setenv PATH /opt/lscsoft/glue/bin:${PATH}
# dynamically get python version (Carsten Aulbert)
setenv PYSITE_PATH python`python -V 2>&1 | cut -d' ' -f2 | cut -d. -f-2`
if ( $?PYTHONPATH ) then
  setenv PYTHONPATH /opt/lscsoft/glue/lib/${PYSITE_PATH}/site-packages:${PYTHONPATH}
else
  setenv PYTHONPATH /opt/lscsoft/glue/lib/${PYSITE_PATH}/site-packages
endif
if ( $?LD_LIBRARY_PATH ) then
  setenv LD_LIBRARY_PATH /opt/lscsoft/glue/lib/${PYSITE_PATH}/site-packages:${LD_LIBRARY_PATH}
else
  setenv LD_LIBRARY_PATH /opt/lscsoft/glue/lib/${PYSITE_PATH}/site-packages
endif
if ( $?DYLD_LIBRARY_PATH ) then
  setenv DYLD_LIBRARY_PATH /opt/lscsoft/glue/lib/${PYSITE_PATH}/site-packages:${DYLD_LIBRARY_PATH}
else
  setenv DYLD_LIBRARY_PATH /opt/lscsoft/glue/lib/${PYSITE_PATH}/site-packages
endif
# Source this file to access GLUE
setenv GLUE_PREFIX /opt/lscsoft/glue
setenv PATH /opt/lscsoft/glue/bin:${PATH}
# dynamically get python version (Carsten Aulbert)
setenv PYSITE_PATH python`python -V 2>&1 | cut -d' ' -f2 | cut -d. -f-2`
if ( $?PYTHONPATH ) then
  setenv PYTHONPATH /opt/lscsoft/glue/lib/${PYSITE_PATH}/site-packages:${PYTHONPATH}
else
  setenv PYTHONPATH /opt/lscsoft/glue/lib/${PYSITE_PATH}/site-packages
endif
if ( $?LD_LIBRARY_PATH ) then
  setenv LD_LIBRARY_PATH /opt/lscsoft/glue/lib/${PYSITE_PATH}/site-packages:${LD_LIBRARY_PATH}
else
  setenv LD_LIBRARY_PATH /opt/lscsoft/glue/lib/${PYSITE_PATH}/site-packages
endif
if ( $?DYLD_LIBRARY_PATH ) then
  setenv DYLD_LIBRARY_PATH /opt/lscsoft/glue/lib/${PYSITE_PATH}/site-packages:${DYLD_LIBRARY_PATH}
else
  setenv DYLD_LIBRARY_PATH /opt/lscsoft/glue/lib/${PYSITE_PATH}/site-packages
endif
