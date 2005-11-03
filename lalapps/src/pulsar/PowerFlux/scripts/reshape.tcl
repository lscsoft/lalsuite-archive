#!/bin/env tclsh

#
# This script parses PowerFlux log files from a run 
# and creates files that can be loaded by R (or other analysis program)
#

source params.tcl

foreach $PARAMS_FORMAT $PARAMS {
        set $var [subst -nocommands $value]
        }

set FIELDS {
	"^useful band start"	band 0		{3 3}
	"^seconds elapsed"	cputime 0  	{2 2}
	"^fine_grid npoints"	npoints 0 	{3 3}
	"^max_high_ul:"	max_high_ul 3  		{1 end}
	"^max_circ_ul:"	max_circ_ul 3  		{1 end}
	"\\(TMedian\\):"	TMedian 0 	{4 4}
	"median:"	median 0 		{2 2}
	"qlines:"	qlines 0 		{2 2}
	"qmost :"	qmost 0			{3 3}
	"^side_cut"	side_cut 0		{2 2}
	"^Maximum bin shift"	max_shift 3	{3 3}
	"^Minimum bin shift"	min_shift 3	{3 3}
	"^hist_residuals:"	hist_residuals 0	{1 4}
        "^max_dx:"       max_dx 1	{1 end}
        "^largest:"  largest 1	{1 end}
        "^masked:"   masked 1	{1 end}
	"^spindown  :" spindown 3	{1 end}
	}

set FIELDS_LAYOUT {exp var pol fields}

for { set band 0 } { $band < $NBANDS } { incr band } {
	lappend FIELDS "^max_high_ul_band: $band " "max_high_ul_band.$band" 3	{1 end}
	lappend FIELDS "^max_circ_ul_band: $band " "max_circ_ul_band.$band" 3	{1 end}
        lappend FIELDS "^max_band: $band " "max_band.$band" 1	{1 end}
        lappend FIELDS "^masked_max_band: $band " "masked_max_band.$band" 1	{1 end}
	lappend FIELDS "^max_ratio: $band " "max_ratio.$band" 1	{1 end}
	lappend FIELDS "hist_.*_ks_test: $band " "ks_hist.$band" 1		{1 4}
        }


foreach $FIELDS_LAYOUT $FIELDS {
	set $var none
	}

# Return subsetted list
proc list_subset { L subset } {
set LL {}
foreach {first last} $subset {
	lappend LL [lrange $L $first $last]
	}
return [join [concat $LL] "\t"]
}

proc list_header { root L subset } {
set LL {}
foreach {first last} $subset {
	if { $last == "end" } { set last [expr [llength $L]-1] }
	for { set i $first } { $i <= $last } { incr i } {
		lappend LL "${root}.$i"
		}
	}
return [join $LL "\t"]
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

set FOUT [open $STATS_DIR/stat${STATS_SUFFIX}.dat "a"]

puts -nonewline $FOUT "LogFile"
# Write out header
foreach $FIELDS_LAYOUT $FIELDS {
	set header [list_header $var [lindex [set $var] end] $fields]
	puts -nonewline $FOUT "\t$header"
	}
puts $FOUT ""

for { set i 0 } { $i < $NPOL*$spindown_count } { incr i } {
	puts -nonewline $FOUT $argv
	foreach $FIELDS_LAYOUT $FIELDS {
		puts -nonewline $FOUT "\t"
		switch -exact $pol {
			1  	{
				set value [list_subset [lindex [set $var] [expr $i]] $fields]
				puts -nonewline $FOUT [join $value "\t"]
				} 
			3 	{
				set value [list_subset [lindex [set $var] [expr $i/$NPOL]] $fields]
				puts -nonewline $FOUT [join $value "\t"]
				}
			0	{
				set value [list_subset [lindex [set $var] end] $fields]
				puts -nonewline $FOUT [join $value "\t"]
				}
			}
		}
	
	puts $FOUT ""
	}

close $FOUT

