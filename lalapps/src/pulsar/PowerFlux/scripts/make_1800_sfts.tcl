#!/usr/bin/env tclsh

source sft_params.tcl

file mkdir $control_info_dir
file mkdir $config_dir
file mkdir $err_dir
file mkdir $sfts_dir


puts stderr "Loading frame library from \"$frame_library\""
set FILE [open $frame_library "r"]
set groups 0
while { ! [eof $FILE] } {
        gets $FILE line
	if { $line == "" } { continue }
        regexp -- $group_regexp $line {} group
        regexp -- $filename_regexp $line {} filename
	set group [expr $group / $seg_step]
        global frames_$group
        if { [info exists frames_$group] } {
                lappend frames_$group $filename
                } {
                set frames_$group $filename
                puts stderr "New group $group"
		incr groups
                }
        }
close $FILE

puts stderr "Total groups: $groups"

proc find_data { start end } {
global seg_start seg_step epoch_regexp
set START1 [expr ($start/$seg_step-1)*$seg_step+$seg_start]
set END1 [expr (($end+$seg_step-1)/$seg_step)*$seg_step+$seg_start]
set L [list ]
for { set i $START1 } { $i < $END1 } { incr i $seg_step } {
	set group [expr $i/$seg_step]
        global frames_$group
        if { [info exists frames_$group] } {
                foreach line [set frames_$group] {
			set epoch 0
                       	regexp -- $epoch_regexp $line {} epoch
                        if { (($epoch+16) >= $start) && ($epoch < $end) } {
                                lappend L $line
                                }
                        }
                }
        }
return $L
}


set cwd [pwd]

puts "\nProcessing flags"
set SEGMENT_FLAGS [open $segment_flags "r"]
while { ! [eof $SEGMENT_FLAGS] } {
	gets $SEGMENT_FLAGS line
	if { $line == "" } { continue }
	if { [regexp {^#} $line ] } { continue }
	regsub {\^} [lindex $line 1] "," args
	set FLAG_[lindex $line 0] [expr round(pow($args))]
	puts "FLAG_[lindex $line 0]=[set FLAG_[lindex $line 0]]"
	}

set VETO_VALUE [expr $veto_expr]
set NON_VETO_VALUE [expr $non_veto_expr]

set SEGMENTS_FILE [open $segments_file "r"]
set i 0
set i_possible 0
set END_PREVIOUS 0
set START_PREVIOUS 0
while { ! [eof $SEGMENTS_FILE] } {
        gets $SEGMENTS_FILE line
        if { $line == "" } { continue }
	if { [regexp {^#} $line ] } { continue }
        # skip segments with non-zero flags
        if { (([lindex $line 4] & ~$NON_VETO_VALUE) & $VETO_VALUE) != 0 } { 
		puts stderr "Skipping segment [lrange $line 0 2]"
		puts stderr "\tFlags: [lrange $line 5 end]"
		continue 
		}

        set START [lindex $line 1]
        set END [lindex $line 2]
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
	        if { (($END-$START) % 1800)>31 } { incr START 15 }
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
		if { $FRAME_FILES == "" } { continue }

                set FILE [open "$config_dir/in.$i" "w"]
                set num [expr 1+$i]
                puts $FILE [subst -nocommands -nobackslashes $make_sft_config]
                puts $FILE [join [lsort -unique -dictionary $FRAME_FILES] "\n"]
                close $FILE
                incr i
                }
        }
close $SEGMENTS_FILE
puts stderr "Generated $i input files (out of possible $i_possible)"

set SUBMIT_FILE [open $submit_file "w"]
puts $SUBMIT_FILE [subst -nocommands {
universe=vanilla
executable=/home/volodya/SFT-3/make_sft_plain
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
	puts $DAG_FILE "JOB A$k dag_sub.$interferometer"
	puts $DAG_FILE "VARS A$k PID=\"$k\""
	}
close $DAG_FILE

