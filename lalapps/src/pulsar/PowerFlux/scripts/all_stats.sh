#!/bin/bash

UTILS=~/PowerFlux/scripts/

ROOT=`$UTILS/get_param.tcl ROOT_DIR`

rm -f *.dat
find $ROOT/output -name powerflux.log -exec $UTILS/split.tcl {} \;

