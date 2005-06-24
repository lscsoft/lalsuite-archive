#!/usr/bin/env tclsh

source params.tcl

foreach $PARAMS_FORMAT $PARAMS {
        set $var [subst -nocommands -nobackslashes $value]
        }
	
foreach {flag value} {
	LINE_CANDIDATE 	(1<<0)
	LINE_HIGH	(1<<1)
	LINE_VERY_HIGH 	(1<<2)
	LINE_CLUSTERED	(1<<3)
	LINE_ISOLATED	(1<<4)
	LINE_YES	(1<<7)
	} {
	set $flag [expr $value]
	}	
	
set FILENAME $argv
set FILE [open $FILENAME "r"]

puts stderr "Reading $FILENAME"
while { ! [eof $FILE] } {
	gets $FILE s
	if { [regexp {^firstbin.*:(.*)$} $s {} FIRSTBIN ] } {
		continue
		}
	if { [regexp {^line detected: bin=(.*) z=(.*) strength=(.*) flag=(.*)$} $s {} fbin z strength flag ] } {
		set comment ""
		if { $strength > 2 } {
			append comment " 6-sigma outlier"
			}
		if { $flag & $LINE_CLUSTERED } {
			append comment " part of broad feature"
			}
		if { $flag & $LINE_ISOLATED } {
			append comment " bin-centered"
			}
		if { $comment != "" } {
			puts "[expr ($FIRSTBIN+$fbin)] [expr ($FIRSTBIN+$fbin)/1800.0] $z $strength $comment"
			}
		}
	}
close $FILE
