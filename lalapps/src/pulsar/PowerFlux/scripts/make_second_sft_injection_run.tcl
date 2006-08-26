#!/usr/bin/env tclsh

#
# This is the same as make_sft_injection_runs.tcl but
# we reuse the injections made previously
#
source "second_sft_injection_params.tcl"

foreach $PARAMS_FORMAT $PARAMS {
	global $var
	set $var [subst -nocommands $value]
	}

global DELETE
set DELETE 0

foreach {var value} $argv {
	global $var
	set $var $value
	}

puts stderr "ROOT_DIR=$ROOT_DIR"

if { $DELETE } {
	file delete -force $ROOT_DIR
	file delete -force $CONF_DIR
	file delete -force $OUTPUT_DIR
	file delete -force $ERR_DIR
	file delete -force $DST_DIR
	file delete -force $LOG_FILE
	}
file mkdir $ROOT_DIR
file mkdir $CONF_DIR $OUTPUT_DIR $ERR_DIR $DST_DIR
	
proc sample { limits } {
set a [lindex $limits 0]
set b [lindex $limits 1]
return [expr $a+rand()*($b-$a)]
}

file copy -force $ORIGINAL_PARAMS_FILE "$ROOT_DIR/params.txt"

set first_bin [expr round($FREQ_START*1800)]	
set k 0
set PARAMS_FILE [open "$ROOT_DIR/params.txt" "r"]
gets $PARAMS_FILE header
set PARAMS_LIST [split $header "\t"]

set DAG_FILE [open "$ROOT_DIR/dag" "w"]
set node_idx 0
while { ! [eof $PARAMS_FILE] }  {
	gets $PARAMS_FILE line
	if { $line == "" } { continue }

	foreach var $PARAMS_LIST value [split $line "\t"] {
		set $var $value
		}

        puts $DAG_FILE "JOB A$k $ROOT_DIR/condor"
        puts $DAG_FILE "VARS A$k PID=\"$k\""

	set FILE [open "$CONF_DIR/$k" "w"]
	puts $FILE [subst -nocommands -nobackslashes $POWERFLUX_CONF_FILE]
	close $FILE

	set FILE [open "$DST_DIR/$k" "w"]
	puts $FILE [subst -nocommands -nobackslashes $DATASET_CONF_FILE]
	close $FILE

	incr k
	if { $k >= $MAX_COUNT } { break }
	}
close $PARAMS_FILE
close $DAG_FILE

set CONDOR_FILE {
universe=standard
executable=$ANALYSIS_PROGRAM
input=/dev/null
output=$ERR_DIR/out.\$(PID)
error=$ERR_DIR/err.\$(PID)
arguments=--config=$CONF_DIR/\$(PID)
log=$LOG_FILE
queue
}

set FILE [open "$ROOT_DIR/condor" "w"]
puts $FILE [subst -nocommands $CONDOR_FILE]
close $FILE

file copy -force second_sft_injection_params.tcl $ROOT_DIR
puts stderr "Prepared input for $k jobs"
