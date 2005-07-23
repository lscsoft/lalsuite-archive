#!/bin/env tclsh

#
# This script parses PowerFlux log files from a run 
# and creates files that can be loaded by R (or other analysis program)
#

source power_injection_params.tcl

foreach $PARAMS_FORMAT $PARAMS {
        set $var [subst -nocommands -nobackslashes $value]
        }

set FIELDS {
	"^useful band start"	band 0
	"^seconds elapsed"	cputime 0
	"^fine_grid npoints"	npoints 0
	"^max_high_ul:"	max_high_ul 3
	"^max_circ_ul:"	max_circ_ul 3
	"\\(TMedian\\):"	TMedian 0
	"median:"	median 0
	"qlines:"	qlines 0
	"qmost :"	qmost 0
	"^side_cut"	side_cut 0
	"^Maximum bin shift"	max_shift 3
	"^Minimum bin shift"	min_shift 3
	"^hist_residuals:"	hist_residuals 0
        "^max_dx:"       max_dx 1
        "^largest:"  largest 1
        "^masked:"   masked 1
	"^spindown  :" spindown 3
	"^fake ra :"  i_ra 0
	"^fake dec:"  i_dec 0
	"^fake orientation" i_orientation 0
	"^fake spindown" i_spindown 0
	"^fake strain" i_strain 0
	"^fake frequency" i_frequency 0
	"^i_closest"	i_closest 1
	"^i_largest"	i_largest 1
	}

set FIELDS_LAYOUT {exp var pol}

foreach $FIELDS_LAYOUT $FIELDS {
	set $var none
	}

# Comment this out to check for completed jobs
set cputime [list "seconds elapsed: NA"]

set spindown_count 1

set FILENAME $argv
set FILE [open $FILENAME "r"]

puts stderr "Reading $FILENAME"

while { ! [eof $FILE] } {
	gets $FILE s
	foreach $FIELDS_LAYOUT $FIELDS {
		if { [regexp $exp $s] } {
			if { [set $var] == "none" } {
				set $var [list $s]
				} {
				lappend $var $s
				}
			break
			}
		}
#	if { [regexp {^spindown count:} $s ] } {
#		set spindown_count [lindex $s 2]
#		} 
	}

close $FILE

set NPOL 5
set spindown_count [llength $spindown]

foreach $FIELDS_LAYOUT $FIELDS {
	if { [set $var] == "none" } {
		puts stderr "Field $var in file $FILENAME is missing !"
		exit -1
		}
	}

file mkdir $STATS_DIR

foreach $FIELDS_LAYOUT $FIELDS {
	set FOUT [open $STATS_DIR/$var${STATS_SUFFIX}.dat "a"]
	switch -exact $pol {
		1  	{
			foreach value [set $var] {
				puts $FOUT $value
				}
			} 
		3 	{
			foreach value [set $var ] {
				for { set i 0 } { $i < $NPOL } { incr i } {
					puts $FOUT $value
					}
				}
			}
		0	{
                        for { set i 0 } { $i < $NPOL*$spindown_count } { incr i } {
                                puts $FOUT [lindex [set $var] end]
                                }
			}
		}
	close $FOUT
	}

