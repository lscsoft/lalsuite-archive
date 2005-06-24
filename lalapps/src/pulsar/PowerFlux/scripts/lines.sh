#!/bin/bash

UTILS=~/PowerFlux/scripts/

ROOT=`$UTILS/get_param.tcl ROOT_DIR`

find $ROOT -name powerflux.log -exec $UTILS/get_lines.tcl {} \;
 
