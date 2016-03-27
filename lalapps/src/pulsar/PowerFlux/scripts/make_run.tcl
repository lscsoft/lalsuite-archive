#!/usr/bin/env tclsh
#
# This script file generates condor and PowerFlux configuration files
# for PowerFlux analysis run
#
# Do not forget to change the parameters in params.tcl
#
#
# Assign defaults
foreach {var value} {
        SFT_LENGTH 1800
        } {
        global $var
        set $var $value
        }


source params.tcl

foreach $PARAMS_FORMAT $PARAMS {
	set $var [subst -nocommands $value]
	}

puts stderr "ROOT_DIR= $ROOT_DIR"

file mkdir $ROOT_DIR
file mkdir $CONF_DIR $OUTPUT_DIR $ERR_DIR

	
proc sample { limits } {
set a [lindex $limits 0]
set b [lindex $limits 1]
return [expr $a+rand()*($b-$a)]
}

set i 0	
set last_i 0
for { set band $FREQ_START } { $band < $FREQ_END } { set band [expr $band+$FREQ_STEP] } {
	set firstbin [expr round($band*$SFT_LENGTH)]

	set DATASET [GET_DATASET $firstbin]

	set max_spindown_count [expr round($MAX_SPINDOWN_COUNT)]
	for { set spindown_start $FIRST_SPINDOWN } { $spindown_start <= $LAST_SPINDOWN } {} {
		set spindown_count [expr round(1.0+($LAST_SPINDOWN-$spindown_start)/$SPINDOWN_STEP)]
		if { $spindown_count > $max_spindown_count } { set spindown_count $max_spindown_count }

		#
		# Loop over skymarks file
		#
		foreach skymarks $SKYMARKS {
			set FILE [open "$CONF_DIR/$i" "w"]
			puts $FILE [subst -nocommands -nobackslashes $POWERFLUX_CONF_FILE]
			close $FILE
			incr i
			}

		set spindown_start [expr $spindown_start + $SPINDOWN_STEP*$spindown_count]
		}
	if { $i > $last_i+500 } {
		puts stderr "$i jobs generated (band=$band max_spindown_count=$max_spindown_count)"
		set last_i $i
		}
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
