# Source this file to access GLUE
#setenv GLUE_PREFIX /opt/lscsoft/glue
setenv GLUE_PREFIX /usr
setenv PATH ${GLUE_PREFIX}/bin:${PATH}
# dynamically get python version (Carsten Aulbert)
setenv PYSITE_PATH python`python -V |& cut -d' ' -f2 | cut -d. -f-2`
if ( $?PYTHONPATH ) then
  setenv PYTHONPATH ${GLUE_PREFIX}/lib/${PYSITE_PATH}/site-packages:${PYTHONPATH}
else
  setenv PYTHONPATH ${GLUE_PREFIX}/lib/${PYSITE_PATH}/site-packages
endif
if ( $?LD_LIBRARY_PATH ) then
  setenv LD_LIBRARY_PATH ${GLUE_PREFIX}/lib/${PYSITE_PATH}/site-packages:${LD_LIBRARY_PATH}
else
  setenv LD_LIBRARY_PATH ${GLUE_PREFIX}/lib/${PYSITE_PATH}/site-packages
endif
if ( $?DYLD_LIBRARY_PATH ) then
  setenv DYLD_LIBRARY_PATH ${GLUE_PREFIX}/lib/${PYSITE_PATH}/site-packages:${DYLD_LIBRARY_PATH}
else
  setenv DYLD_LIBRARY_PATH ${GLUE_PREFIX}/lib/${PYSITE_PATH}/site-packages
endif
