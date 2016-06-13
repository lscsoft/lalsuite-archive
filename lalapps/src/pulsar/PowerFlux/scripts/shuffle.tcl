#!/usr/bin/env tclsh

foreach {filename batch step} $argv {
	set FILE [open $filename "r"]
	set LINES [split [read $FILE] '\n']
	close $FILE

	set FOUT [open "${filename}.shuffled" "w"]
	for { set i 0 } { $i < $step } { incr i } {
		for { set j [expr $batch*$i] } { $j < [llength $LINES] } { incr j [expr $step*$batch] } {
			for { set k 0 } { $k < $batch } { incr k } {
				puts $FOUT [lindex $LINES [expr $j+$k]]
				}
			}
		}
	}	
