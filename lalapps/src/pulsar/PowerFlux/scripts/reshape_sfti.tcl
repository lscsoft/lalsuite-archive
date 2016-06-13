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
	"^label"		label 0  {1 1}  1
	"^useful band start"	band 0		{3 3}		1
	"^seconds elapsed"	cputime 0  	{2 2}		1
	"^fine_grid npoints"	npoints 0 	{3 3}		1
	"^max_high_ul:"	max_high_ul 3  		{1 4}		4
	"^max_circ_ul:"	max_circ_ul 3  		{1 4}		4
	"\\(TMedian\\):"	TMedian 0 	{4 4}		1
	"median:"	median 0 		{2 3}		2
	"qlines:"	qlines 0 		{2 2}		1
	"qmost :"	qmost 0			{4 4}		1
	"^side_cut"	side_cut 0		{2 2}		1
	"^Maximum bin shift"	max_shift 3	{3 3}		1
	"^Minimum bin shift"	min_shift 3	{3 3}		1
	"^hist_residuals:"	hist_residuals 0	{1 4}	4
        "^max_dx:"       max_dx 1	{1 7}			7
        "^largest:"  largest 1	{1 9}				9
        "^masked:"   masked 1	{1 2}				2
	"^spindown  :" spindown 3	{2 2}			1
	"^noise_floor:"  noise_floor  1  {1 1}			1
	"^Second pass processing time:" second_pass_time 1  {4 4} 	1
	"^optimized_candidates_count:" optimized_candidates_count 1  {1 1} 	1
	"^high_candidates_count:" high_candidates_count 1  {1 1} 	1
	"^candidates_count:" candidates_count 1  {1 1} 	1
	"^grid_points: 0" grid_points 1 {1 3} 3
	"^fake ra :" fake_ra 0 {3 3} 1
	"^fake dec:" fake_dec 0 {2 2} 1
	"^fake psi :" fake_psi 0 {3 3} 1
	"^fake phi :" fake_phi 0 {3 3} 1
	"^fake iota :" fake_iota 0 {3 3} 1
	"^fake spindown:" fake_spindown 0 {2 2} 1
	"^fake strain:" fake_strain 0 {2 2} 1
	"^fake frequency:" fake_frequency 0 {2 2} 1
	}

set FIELDS_LAYOUT {exp var pol fields EC}

for { set band 0 } { $band < $NBANDS } { incr band } {
	#lappend FIELDS "^max_high_ul_band: $band " "max_high_ul_band.$band" 3	{1 3}		3
	#lappend FIELDS "^max_circ_ul_band: $band " "max_circ_ul_band.$band" 3	{1 3}		3
        #lappend FIELDS "^max_band: $band " "max_band.$band" 1	{1 9}				9
        #lappend FIELDS "^masked_max_band: $band " "masked_max_band.$band" 1	{1 9}		9
	#lappend FIELDS "^max_ratio: $band " "max_ratio.$band" 1	{1 3}				3
	lappend FIELDS "hist_.*_ks_test: $band " "ks_hist.$band" 1		{1 4}		4
        }


foreach $FIELDS_LAYOUT $FIELDS {
	set $var none
	}

# Return subsetted list
proc list_subset { L subset EC } {
set LL {}
foreach {first last} $subset {
	lappend LL [lrange $L $first $last]
	}
set S [join $LL " "]
if { $EC != "X" } {
	set L [split $S " "]
	if { [llength $L ] != $EC } {
		set L  {}
		for { set i 0 } { $i < $EC } { incr i } { lappend L NaN }
		set S [join $L " "]
		}
	}
return $S
}

proc list_header { root L subset } {
set LL {}
foreach {first last} $subset {
	if { $last == "end" } { set last [expr [llength $L]-1] }
	for { set i $first } { $i <= $last } { incr i } {
		lappend LL "${root}.$i"
		}
	}
return [join $LL " "]
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

set OFILENAME "$STATS_DIR/stat${STATS_SUFFIX}.dat"

set appending [file exists $OFILENAME]

set FOUT [open $OFILENAME "a"]

if { !$appending } {
	puts -nonewline $FOUT "LogFile"
	# Write out header
	foreach $FIELDS_LAYOUT $FIELDS {
		set header [list_header $var [lindex [set $var] end] $fields]
		puts -nonewline $FOUT " $header"
		}
	puts $FOUT ""
	}

for { set i 0 } { $i < $NPOL*$spindown_count } { incr i } {
	puts -nonewline $FOUT $argv
	foreach $FIELDS_LAYOUT $FIELDS {
		puts -nonewline $FOUT " "
		switch -exact $pol {
			1  	{
				set value [list_subset [lindex [set $var] [expr $i]] $fields $EC]
				puts -nonewline $FOUT $value
				} 
			3 	{
				set value [list_subset [lindex [set $var] [expr $i/$NPOL]] $fields $EC]
				puts -nonewline $FOUT $value
				}
			0	{
				set value [list_subset [lindex [set $var] end] $fields $EC]
				puts -nonewline $FOUT $value
				}
			}
		}
	
	puts $FOUT ""
	}

close $FOUT

