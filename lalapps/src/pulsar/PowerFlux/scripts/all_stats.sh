#!/bin/bash

UTILS=~/PowerFlux/scripts/

ROOT=`$UTILS/get_param.tcl ROOT_DIR`
STATS_DIR=`$UTILS/get_param.tcl STATS_DIR`

rm -f $STATS_DIR/*.dat
find $ROOT/output -name powerflux.log -exec $UTILS/split.tcl {} \;

