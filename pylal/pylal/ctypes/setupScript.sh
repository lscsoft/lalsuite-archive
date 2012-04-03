PYTHON=python2.7
PYTHON_INSTALL_PATH=/usr/lib/$PYTHON

$PYTHON $PYTHON_INSTALL_PATH/Tools/scripts/h2py.py $LAL_PREFIX/include/lal/LALConstants.h
mv LALCONSTANTS.py constants.py
$PYTHON $PYTHON_INSTALL_PATH/Tools/scripts/h2py.py $LAL_PREFIX/include/lal/LALDetectors.h
mv LALDETECTORS.py laldetectors.py
$PYTHON $PYTHON_INSTALL_PATH/Tools/scripts/h2py.py $LAL_PREFIX/include/lal/LALDatatypes.h
mv LALDATATYPES.py laldatatypes.py
$PYTHON $PYTHON_INSTALL_PATH/Tools/scripts/h2py.py $LAL_PREFIX/include/lal/LALAtomicDatatypes.h
mv LALATOMICDATATYPES.py lalatomicdatatypes.py


