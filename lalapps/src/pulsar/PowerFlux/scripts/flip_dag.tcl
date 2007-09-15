#!/usr/bin/env tclsh

set F [open $argv "r"]
set lines {}
while { ![eof $F] } {
	gets $F line
	if { $line == "" } { continue }
	lappend lines $line
	}
close $F

set L {}
foreach {line1 line2} $lines {
	set L [linsert $L 0 $line2]
	set L [linsert $L 0 $line1]
	}
puts [join $L "\n"]
	
