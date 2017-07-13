#!/usr/bin/env tclsh

set i 1
while { ![eof stdin] } {
	gets stdin line

	if { $line == "" } { continue }

	puts "JOB X$i condorX"
	puts "VARS X$i PID=\"$i\" INPUT=\"$line\" OUTPUT=\"${line}.output\""
	puts ""
	incr i
	}

