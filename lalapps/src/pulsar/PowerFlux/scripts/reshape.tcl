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
	"^useful band start"	band 0		{3 3}		1	0
	"^seconds elapsed"	cputime 0  	{2 2}		1	0
	"^fine_grid npoints"	npoints 0 	{3 3}		1	0
	"^max_high_ul:"	max_high_ul 3  		{1 4}		4	0
	"^max_circ_ul:"	max_circ_ul 3  		{1 4}		4	0
	"\\(TMedian\\):"	TMedian 0 	{4 4}		1	0
	"median:"	median 0 		{2 3}		2	0
	"qlines:"	qlines 0 		{2 2}		1	0
	"qmost :"	qmost 0			{3 3}		1	0
	"^side_cut"	side_cut 0		{2 2}		1	0
	"^Maximum bin shift"	max_shift 3	{3 3}		1	0
	"^Minimum bin shift"	min_shift 3	{3 3}		1	0
	"^hist_residuals:"	hist_residuals 0	{1 4}	4	0
        "^max_dx:"       max_dx 1	{1 7}			7	0
        "^largest:"  largest 1	{1 9}				9	0
        "^masked:"   masked 1	{1 2}				2	0
	"^spindown  :" spindown 3	{2 2}			1	0
	"^noise_floor:"  noise_floor  0  {2 2}			1	0
	"^optimized_candidates_count:" optimized_candidates_count 3  {1 1} 	1	0
	"^high_candidates_count:" high_candidates_count 3  {1 1} 	1	0
	"^candidates_count:" candidates_count 3  {1 1} 	1	0

	"grid_points: ([^ ]*) " "grid_points" 3 {1 3} 3	1
	"^max_high_ul_band: ([^ ]*) " "max_high_ul_band" 3	{1 3}		3	1
        "^max_band: ([^ ]*) " "max_band" 1	{1 9}				9	1
	"^max_dx_band: ([^ ]*) " "max_dx_band" 3	{1 8}		8	1
	"^max_ratio: ([^ ]*) " "max_ratio" 1	{1 3}				3	1
	"hist_.*_ks_test: ([^ ]*) " "ks_hist" 1		{1 4}		4	1
	}

set OMITTED_FIELDS {
	"^max_circ_ul_band: ([^ ]*) " "max_circ_ul_band" 3	{1 3}		3	1
        "^masked_max_band: ([^ ]*) " "masked_max_band" 1	{1 9}		9	1
	}

set FIELDS_LAYOUT {exp var pol fields EC param}



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
set cputime [list "seconds elapsed: NaN"]

set spindown_count 1

set FILENAME $argv
set FILE [open $FILENAME "r"]

puts stderr "Reading $FILENAME"

set input {}

while { ! [eof $FILE] } {
	gets $FILE s
	foreach $FIELDS_LAYOUT $FIELDS {
		if { [regexp $exp $s {} band] } {
			if { $param } {
				set var "${var}.$band"
				if { $band >= $NBANDS } {
					set NBANDS [expr $band+1]
					puts stderr "WARNING: Increased band count to $NBANDS"
					}
				}
			lappend input $var $s
			break
			}
		}
#	if { [regexp {^spindown count:} $s ] } {
#		set spindown_count [lindex $s 2]
#		} 
	}

close $FILE

set EXPANDED_FIELDS {}
foreach $FIELDS_LAYOUT $FIELDS {
	if { $param } {
		set var_root $var
		for { set i 0 } { $i < $NBANDS } { incr i } {
			set var "${var_root}.$i"
			foreach v $FIELDS_LAYOUT { lappend EXPANDED_FIELDS [set $v] }
			}
		} {
		foreach v $FIELDS_LAYOUT { lappend EXPANDED_FIELDS [set $v] }
		}
	}

foreach $FIELDS_LAYOUT $EXPANDED_FIELDS {
	set $var {}
	}

foreach {var line } $input {
	lappend $var $line
	}

set NPOL 5
set spindown_count [llength $spindown]

foreach $FIELDS_LAYOUT $EXPANDED_FIELDS {
	if { [set $var] == "" } {
		puts stderr "Field $var in file $FILENAME is missing !"
		exit -1
		}
	}

file mkdir $STATS_DIR

set OFILENAME "$STATS_DIR/stat${STATS_SUFFIX}.dat"

#set appending [file exists $OFILENAME]

set FOUT [open $OFILENAME "a"]
fconfigure $FOUT -buffering full -buffersize 10000000

if { [tell $FOUT] == 0 } {
	puts stderr "\tWriting header"
	puts -nonewline $FOUT "LogFile"
	# Write out header
	foreach $FIELDS_LAYOUT $EXPANDED_FIELDS {
		set header [list_header $var [lindex [set $var] end] $fields]
		puts -nonewline $FOUT " $header"
		}
	puts $FOUT ""
	}

for { set i 0 } { $i < $NPOL*$spindown_count } { incr i } {
	puts -nonewline $FOUT $argv
	foreach $FIELDS_LAYOUT $EXPANDED_FIELDS {
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

