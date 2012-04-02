PYTHON=python2.7
PYTHON_INSTALL_PATH=/usr/lib/$PYTHON

$PYTHON $PYTHON_INSTALL_PATH/Tools/scripts/h2py.py ~/reps/lalsuite/lal/include/lal/LALConstants.h
mv LALCONSTANTS.py constants.py
$PYTHON $PYTHON_INSTALL_PATH/Tools/scripts/h2py.py ~/reps/lalsuite/lal/include/lal/LALDetectors.h
mv LALDETECTORS.py laldetectors.py
