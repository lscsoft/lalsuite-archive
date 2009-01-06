#!/usr/bin/env tclsh

source "sft_injection_params.tcl"

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

puts stderr "ROOT_DIR= $ROOT_DIR"

if { $DELETE } {
	file delete -force $ROOT_DIR
	file delete -force $CONF_DIR
	file delete -force $OUTPUT_DIR
	file delete -force $ERR_DIR
	file delete -force $LOG_FILE
	}
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
set PARAMS_FILE [open "$ROOT_DIR/params.txt" "w"]
set PARAMS_LIST {i h0 ra dec psi phi iota f0 spindown aPlus aCross band_start}
puts $PARAMS_FILE [join $PARAMS_LIST "\t"]
set DAG_FILE [open "$ROOT_DIR/dag" "w"]

while { 1 }  {
	set dec [sample $DEC_RANGE]
	set ra [sample $RA_RANGE]
	set psi [sample $PSI_RANGE]
	set phi [sample $PHI_RANGE]
	set iota [sample $IOTA_RANGE]
	set freq [sample [list [expr $first_bin/1800.0] [expr ($first_bin+450)/1800.0]]]
	set h0  [expr exp([sample $POWER_LOG10_RANGE] * log(10.0)) * $POWER_MAX]
        set spindown [expr exp([sample $SPINDOWN_LOG10_RANGE] * log(10.0)) * $SPINDOWN_MAX]
	set aPlus [expr $h0 * (1.0+cos($iota)*cos($iota))/2.0]
	set aCross [expr $h0 * cos($iota)]
	set f0 $freq
	set band_start [expr $first_bin/1800.0]
	set seed [expr round(rand()*1e10)]

	#
	# And old, but proven method of sampling arbitrary distribution
	# Note that density_max should better be precise, or very close to precise
	#
	set density_x [sample [list 0 $INJECTION_DENSITY_MAX]]
	set density_cut [eval expr $INJECTION_DENSITY]
	if { $density_x > $density_cut } continue

	set L [list ]
	foreach var $PARAMS_LIST {
		lappend L [set $var]
		}
	puts $PARAMS_FILE [join $L "\t"]
	flush $PARAMS_FILE

        puts $DAG_FILE "JOB A$i $ROOT_DIR/condor"
        puts $DAG_FILE "VARS A$i PID=\"$i\""

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
notification=Error
notify_user=volodya@mindspring.com
want_remote_io=False
queue
}

set FILE [open "$ROOT_DIR/condor" "w"]
puts $FILE [subst -nocommands $CONDOR_FILE]
close $FILE

exec cp sft_injection_params.tcl $ROOT_DIR
puts stderr "Prepared input for $i jobs"
