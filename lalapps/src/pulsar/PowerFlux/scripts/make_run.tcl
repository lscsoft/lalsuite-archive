#!/usr/bin/env tclsh
#
# This script file generates condor and PowerFlux configuration files
# for PowerFlux analysis run
#
# Do not forget to change the parameters in params.tcl
#
#
source params.tcl

foreach $PARAMS_FORMAT $PARAMS {
	set $var [subst -nocommands -nobackslashes $value]
	}

file mkdir $ROOT_DIR
file mkdir $CONF_DIR $OUTPUT_DIR $ERR_DIR

	
set CONDOR_FILE {
universe=standard
executable=/home/volodya/PowerFlux/powerflux
input=/dev/null
output=$ERR_DIR/out.\$(PID)
error=$ERR_DIR/err.\$(PID)
arguments=--config=$CONF_DIR/\$(PID)
log=$DAG_LOG
queue
}

proc sample { limits } {
set a [lindex $limits 0]
set b [lindex $limits 1]
return [expr $a+rand()*($b-$a)]
}

set i 0	
for { set f $FREQ_START } { $f < $FREQ_END } { set f [expr $f+$FREQ_STEP] } {
	set firstbin [expr round($f*1800)]

	set FILE [open "$CONF_DIR/$i" "w"]
	puts $FILE [subst -nocommands -nobackslashes $POWERFLUX_CONF_FILE]
	close $FILE
	incr i
	}
	
set FILE [open "$ROOT_DIR/condor" "w"]
puts $FILE [subst -nocommands $CONDOR_FILE]
close $FILE

set FILE [open "$ROOT_DIR/dag" "w"]
for { set k 0 } { $k < $i } { incr k } {
	puts $FILE "JOB A$k $ROOT_DIR/condor"
	puts $FILE "VARS A$k PID=\"$k\""
	}
close $FILE

exec cp params.tcl $ROOT_DIR
