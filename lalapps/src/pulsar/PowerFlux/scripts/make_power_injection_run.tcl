#!/usr/bin/env tclsh

source "power_injection_params.tcl"

foreach $PARAMS_FORMAT $PARAMS {
	set $var [subst -nocommands $value]
	}

puts stderr "ROOT_DIR=$ROOT_DIR"

file mkdir $ROOT_DIR
file mkdir $CONF_DIR $OUTPUT_DIR $ERR_DIR
	
proc sample { limits } {
set a [lindex $limits 0]
set b [lindex $limits 1]
return [expr $a+rand()*($b-$a)]
}

set first_bin [expr round($FREQ_START*1800)]	
set i 0
set k 0
while { 1 }  {
	set dec [sample $DEC_RANGE]
	set ra [sample $RA_RANGE]
	set orientation [sample $ORIENT_RANGE]
	set freq [sample [list [expr $first_bin/1800.0] [expr ($first_bin+450)/1800.0]]]
	set power [expr exp([sample $POWER_LOG10_RANGE] * log(10.0)) * $POWER_MAX]
        set spindown [expr exp([sample $SPINDOWN_LOG10_RANGE] * log(10.0)) * $SPINDOWN_MAX]

	#
	# And old, but proven method of sampling arbitrary distribution
	# Note that density_max should better be precise, or very close to precise
	#
	set density_x [sample [list 0 $INJECTION_DENSITY_MAX]]
	set density_cut [eval expr $INJECTION_DENSITY]
	if { $density_x > $density_cut } continue

	set FILE [open "$CONF_DIR/$i" "w"]
	puts $FILE [subst -nocommands -nobackslashes $POWERFLUX_CONF_FILE]
	close $FILE

	incr i
	incr k

	if { $k >= $INJECTIONS_PER_BAND } {
		puts stderr "Band [expr $first_bin/1800.0] complete"
		set first_bin [expr round($first_bin+$FREQ_STEP*1800)]
		if { $first_bin >= $FREQ_STOP*1800 } { break }
		set k 0
		}
	}

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

set FILE [open "$ROOT_DIR/dag" "w"]
for { set k 0 } { $k < $i } { incr k } {
	puts $FILE "JOB A$k $ROOT_DIR/condor"
	puts $FILE "VARS A$k PID=\"$k\""
	}
close $FILE

exec cp power_injection_params.tcl $ROOT_DIR
puts stderr "Prepared input for $i jobs"
