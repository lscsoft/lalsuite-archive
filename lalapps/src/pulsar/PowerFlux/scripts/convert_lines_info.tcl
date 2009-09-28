#!/usr/bin/env tclsh

puts "instance\tfirst_bin\tband\tdataset\tbin\tz\tstrength\tflag"
while { ! [eof stdin] } {
	gets stdin line
	if { $line == "" } { continue }

	if { ! [regexp {output/([^/]*)/([^/-]*)-([^/-]*)/powerflux.log:(.*) background line detected: bin=(.*) z=(.*) strength=(.*) flag=(.*)$} $line {} instance first_bin frequency dataset bin z strength flag] } {
		puts stderr "Could not parse line \"$line\""
		continue
		}
	puts "$instance\t$first_bin\t$frequency\t\"$dataset\"\t$bin\t$z\t$strength\t[expr $flag]"
	}
 
