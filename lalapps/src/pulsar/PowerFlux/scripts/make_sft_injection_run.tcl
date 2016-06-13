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

puts stderr "ROOT_DIR=$ROOT_DIR"
global N_NODES 
set N_NODES [llength $STORAGE_NODE_LIST]
puts stderr "$N_NODES storage nodes found"

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
set PARAMS_LIST {i h0 ra dec psi phi iota f0 spindown aPlus aCross band_start INJ_SFT_DIR}
puts $PARAMS_FILE [join $PARAMS_LIST "\t"]
set DAG_FILE [open "$ROOT_DIR/dag" "w"]
set node_idx 0
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

	#
	# And old, but proven method of sampling arbitrary distribution
	# Note that density_max should better be precise, or very close to precise
	#
	set density_x [sample [list 0 $INJECTION_DENSITY_MAX]]
	set density_cut [eval expr $INJECTION_DENSITY]
	if { $density_x > $density_cut } continue

	set nodenum [lindex $STORAGE_NODE_LIST $node_idx]
	set INJ_SFT_DIR [subst -nocommands "$SFT_OUTPUT_DIR_TEMPLATE/$i"]
	if { $DELETE } {
		file delete -force ${INJ_SFT_DIR}
		}
	file mkdir ${INJ_SFT_DIR}

	set L [list ]
	foreach var $PARAMS_LIST {
		lappend L [set $var]
		}
	puts $PARAMS_FILE [join $L "\t"]
	flush $PARAMS_FILE

        puts $DAG_FILE "JOB I$i $ROOT_DIR/icondor"
        puts -nonewline $DAG_FILE "VARS I$i PID=\"$i\" INJ_SFT_DIR=\"$INJ_SFT_DIR\"" 
	foreach var $PARAMS_LIST {
		puts -nonewline $DAG_FILE " $var=\"[set $var]\""
		}
	puts $DAG_FILE ""
        puts $DAG_FILE "JOB A$i $ROOT_DIR/condor"
        puts $DAG_FILE "VARS A$i PID=\"$i\""
        puts $DAG_FILE "PARENT I$i CHILD A$i"

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
	incr node_idx
	if { $node_idx >= $N_NODES } { 
		set node_idx 0
		}
	}
close $PARAMS_FILE
close $DAG_FILE

set CONDOR_INJECTION_FILE {
universe=vanilla
executable=$INJECTION_SCRIPT
input=/dev/null
output=$ERR_DIR/iout.\$(PID)
error=$ERR_DIR/ierr.\$(PID)
arguments= detector $IFO ra \$(ra) dec \$(dec) psi \$(psi) phi \$(phi) f0 \$(f0) spindown \$(spindown) aPlus \$(aPlus) aCross \$(aCross) band_start \$(band_start) noiseGlob ${SFT_INPUT}* outputDir \$(INJ_SFT_DIR)
log=$LOG_FILE
queue
}

set FILE [open "$ROOT_DIR/icondor" "w"]
puts $FILE [subst -nocommands $CONDOR_INJECTION_FILE]
close $FILE

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

exec cp sft_injection_params.tcl $ROOT_DIR
puts stderr "Prepared input for $i jobs"
