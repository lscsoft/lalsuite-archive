#!/usr/bin/env tclsh

source sft_params.tcl

file mkdir $control_info_dir
file mkdir $config_dir
file mkdir $err_dir
file mkdir $sfts_dir

# Open log file
global LOG
set LOG [open $generation_log_file "w"]

puts $LOG "Loading frame library from \"$frame_library\""
puts stderr "Loading frame library from \"$frame_library\""
set FILE [open $frame_library "r"]
set groups 0
while { ! [eof $FILE] } {
        gets $FILE line
	if { $line == "" } { continue }
        regexp -- $filename_regexp $line {} filename
        regexp -- $frame_start_length_regexp $filename {} frame_start frame_length
	set group [expr $frame_start / $seg_step]
        global frames_$group
        if { [info exists frames_$group] } {
                lappend frames_$group $filename
                } {
                set frames_$group $filename
                puts $LOG "New group $group"
		incr groups
                }
        }
close $FILE

puts $LOG "Total groups: $groups"

#
# Find frame files covering segment [start, end]
#
proc find_data { start end } {
global seg_start seg_step frame_start_length_regexp LOG
set START1 [expr ($start/$seg_step-1)*$seg_step+$seg_start]
set END1 [expr (($end+$seg_step-1)/$seg_step)*$seg_step+$seg_start]
set L [list ]
for { set i $START1 } { $i < $END1 } { incr i $seg_step } {
	set group [expr $i/$seg_step]
        global frames_$group
        if { [info exists frames_$group] } {
                foreach filename [set frames_$group] {
			set frame_start 0
                       	regexp -- $frame_start_length_regexp $filename {} frame_start frame_length
                        if { (($frame_start+$frame_length) >= $start) && ($frame_start < $end) } {
                                lappend L $filename
                                }
                        }
                }
        }
return [lsort $L]
}

#
# Find earliest time within [start, end] for which the data exists
#
proc bump_start { start end } {
global seg_start seg_step frame_start_length_regexp LOG
set START1 [expr ($start/$seg_step-1)*$seg_step+$seg_start]
set END1 [expr (($end+$seg_step-1)/$seg_step)*$seg_step+$seg_start]
set t 0
for { set i $START1 } { $i < $END1 } { incr i $seg_step } {
	set group [expr $i/$seg_step]
        global frames_$group
        if { [info exists frames_$group] } {
                foreach filename [set frames_$group] {
			set frame_start 0
			set frame_length 0
                       	regexp -- $frame_start_length_regexp $filename {} frame_start frame_length

                        if { (($frame_start+$frame_length) >= $start) && ($frame_start < $end) } {
				if { $t == 0 } { 
					set t $frame_start
					}
				if { $frame_start < $t } { 
					set t $frame_start
					}
                                }
                        }
                }
	if { $t > 0 } {
		# Found something - other segments will be later
		if { $t <= $start } { 
			set t $start
			} {
			puts $LOG "Moved start $start to $t"
			}
		return $t
		}
        }
return $end
}


set cwd [pwd]

puts -nonewline "\nInitializing flags:"
set veto(ALL) 0

set SEGMENT_FLAGS [open $segment_flags "r"]
while { ! [eof $SEGMENT_FLAGS] } {
	gets $SEGMENT_FLAGS line
	if { $line == "" } { continue }
	if { [regexp {^#} $line ] } { continue }

	set flag [lindex $line 0]

	puts -nonewline " $flag"

	# Initialize veto and non_veto sets
	set veto($flag)   0
	set non_veto($flag) 0
	}
puts ""

foreach value $veto_set {
	set veto($value) 1
	}
	
foreach value $non_veto_set {
	set non_veto($value) 1
	}

set SEGMENTS_FILE [open $segments_file "r"]
set i 0
set i_possible 0
set END_PREVIOUS 0
set START_PREVIOUS 0
while { ! [eof $SEGMENTS_FILE] } {
        gets $SEGMENTS_FILE line
        if { $line == "" } { continue }
	if { [regexp {^#} $line ] } { continue }
	#
	# Use symbolic flags
	#
	set flags [lrange $line 5 end]
	set do_veto 0
	foreach value  $flags {
		if { !$non_veto($value) && ($veto(ALL) || $veto($value)) } {
			puts $LOG "Skipping segment [lrange $line 0 2] due to flag $value"
			puts $LOG "\tFlags: $flags"
			set do_veto 1
			break
			}
		}
	if { $do_veto } { continue }

        set START [lindex $line 1]
        set END [lindex $line 2]
        puts $LOG "processing i=$i segment $START-$END"
        puts stderr "processing i=$i segment $START-$END"
        #
        # Make sure to have some valid data in case we don't use all data from the segment
        # This way we can be certain that all data made it to RDS frames (which are 16 sec long)
        #
	if { $END_PREVIOUS == $START } {
		# we can merge these two segments
		set START $START_PREVIOUS
		} {
		# Brand new chunk of data
		
		# Move START so there is a frame file covering it
		set START [bump_start $START $END]
		}
        for { set k 0 } { 1 } { incr k } {
                set start [expr $START+$k*$overlap]
                set end [expr $START+$k*$overlap+$timebase]
                if { $end > $END } { 
			set END_PREVIOUS $END
			set START_PREVIOUS $start
			break 
			}

		incr i_possible

                set FRAME_FILES [find_data $start $end]
		if { $FRAME_FILES == "" } { 
			puts $LOG "** Did not find any frame files for segment $start $end"
			continue 
			}

                set FILE [open "$config_dir/in.$i" "w"]
                set num [expr 1+$i]
                puts $FILE [subst -nocommands -nobackslashes $make_sft_config]
                puts $FILE [join [lsort -unique -dictionary $FRAME_FILES] "\n"]
                close $FILE
                incr i
                }
        }
close $SEGMENTS_FILE
puts $LOG "Generated $i input files (out of possible $i_possible)"
puts stderr "Generated $i input files (out of possible $i_possible)"

set SUBMIT_FILE [open $submit_file "w"]
puts $SUBMIT_FILE [subst -nocommands {
universe=vanilla
executable=$sft_program
input=$config_dir/in.\$(PID)
output=$err_dir/out.\$(PID)
error=$err_dir/err.\$(PID)
arguments=
log=$log_file
queue
} ]
close $SUBMIT_FILE

set DAG_FILE [open $dag_file "w"]
for { set k 0 } { $k < $i } { incr k } {
	puts $DAG_FILE "JOB A$k submit"
	puts $DAG_FILE "VARS A$k PID=\"$k\""
	}
close $DAG_FILE

