#!/bin/env tclsh

#
# This script parses PowerFlux log files from a run 
# and creates files that can be loaded by R (or other analysis program)
#

foreach {var value} {
	instance_label {""}
	} {
	global $var
	set $var $value
	}

set FIELDS {
	"^useful band start: "	{}	band 0		{0 0}		1
	"^seconds elapsed: "	{}	cputime 0  	{0 0}		1
	"^fine_grid npoints  : "	{}	npoints 0 	{0 0}		1
	"^max_high_ul: "	"^max_high_ul legend: "	max_high_ul 3  		{0 3}		4
	"^max_circ_ul: "	"^max_circ_ul legend: "	max_circ_ul 3  		{0 3}		4
	"\\(TMedian\\): "	{}	TMedian 0 	{0 0}		1
	"median: "	{}	median 0 		{0 1}		2
	"qlines: "	{}	qlines 0 		{0 0}		1
	"qmost : "	{}	qmost 0			{0 0}		1
	"^side_cut  : "	{}	side_cut 0		{0 0}		1
	"^Maximum bin shift: "	{}	max_shift 3	{0 0}		1
	"^Minimum bin shift: "	{}	min_shift 3	{0 0}		1
	"^hist_residuals: "	{}	hist_residuals 0	{0 3}	4
        "^max_dx: "       "^strongest signal: "	max_dx 1	{0 6}			7
        "^largest: "  "^largest signal: "	largest 1	{0 8}				9
        "^masked: "   {}	masked 1	{0 1}				2
	"^spindown  :" {}	spindown 3	{0 0}			1
	"^noise_floor: "  {}	noise_floor  0  {0 0}			1
	"^optimized_candidates_count: " {}	optimized_candidates_count 3  {0 0} 	1
	"^high_candidates_count: " {}	high_candidates_count 3  {0 0} 	1
	"^candidates_count: " {}	candidates_count 3  {0 0} 	1

	"grid_points: " {}	"grid_points" 3 {0 2} 3
	"^max_high_ul_band: " {}	"max_high_ul_band" 3	{0 2}		3
        "^masked_max_band: " "^max/masked band format: "	"masked_max_band" 1	{0 8}		9
        "^max_band: " "^max/masked band format: "	"max_band" 1	{0 8}				9
	"^max_dx_band: " {}	"max_dx_band" 3	{0 7}		8
	"^max_ratio: " {}	"max_ratio" 1	{0 2}				3
	"hist_.*_ks_test: " {}	"ks_hist" 1		{0 3}		4
	"^max_circ_ul_band: " {}	"max_circ_ul_band" 3	{0 2}		3
	}

set FIELDS_LAYOUT {exp header_exp var pol fields EC}

set IGNORE_LINES {^(candidate_(initial|optimized|final)|Marking disk)}

source params.tcl

foreach $PARAMS_FORMAT $PARAMS {
        set $var [subst -nocommands $value]
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

set match_list {}
set i 0

set L {}
set group_exp ""
foreach $FIELDS_LAYOUT $FIELDS {
	if { $i > 0 } { append group_exp "|" }
	append group_exp "($exp)"
	if { $header_exp!="" } { append group_exp "|($header_exp)" }
	foreach element $FIELDS_LAYOUT { lappend L [set $element] }
	incr i

	if { $i > 5 } {
		lappend match_list "$group_exp" $L
		set group_exp ""
		set L {}
		set i 0
		}
	}
if { [llength $L] > 0 } {
	lappend match_list "$group_exp" $L
	}

file mkdir $STATS_DIR
	
foreach FILENAME $argv {
	set FILE [open $FILENAME "r"]
	
	puts stderr "Reading $FILENAME"
	
	foreach $FIELDS_LAYOUT $FIELDS {
		set $var {}
		}

	# Comment this out to check for completed jobs
	set cputime [list "seconds elapsed: NaN"]

	while { ! [eof $FILE] } {
		gets $FILE s
		if { [regexp "^label: (.*)$" $s {} instance_label] } { continue }
		if { [regexp $IGNORE_LINES $s ] } { continue }
		foreach {group_exp match_sublist} $match_list {
			if { ! [regexp $group_exp $s] } { continue }
			foreach $FIELDS_LAYOUT $match_sublist {
				if { $header_exp!="" && [regexp "${header_exp}(.*)$" $s {} header] } {
					set ${var}_header $header
					}
				if { [regexp "${exp}(.*)$" $s {} data] } {
					lappend $var $data
					break
					}
				}
			}
		}
	
	close $FILE
	
	set next 0
	foreach $FIELDS_LAYOUT $FIELDS {
		if { [set $var] == "" } {
			puts stderr "Field $var in file $FILENAME is missing !"
			set next 1
			break
			}
		}

	if { $next } { continue }
	
	foreach $FIELDS_LAYOUT $FIELDS {
	
		puts -nonewline stderr " $var"
	
		
		#set appending [file exists $OFILENAME]

		set output_header 0

		if { ! [info exists OF_$var] } {
			set OFILENAME "$STATS_DIR/stat_${var}${STATS_SUFFIX}.dat"
			set FOUT [open $OFILENAME "a"]
			fconfigure $FOUT -buffering full -buffersize 10000000
			set OF_$var $FOUT

			if { [tell $FOUT] == 0 } {
				set output_header 1
				}

			} {
			set FOUT [set OF_$var]
			#puts -nonewline stderr "+"
			}

		#set FOUT [open $OFILENAME "a"]
		#fconfigure $FOUT -buffering full -buffersize 10000000
		
		if { $output_header } {
			puts -nonewline stderr "(*)"
			# Write out header
			if { $header_exp != {} } {
				set header [list_subset [set ${var}_header] $fields $EC]
				} {
				set header [list_header $var [lindex [set $var] end] $fields]
				}
			puts $FOUT "LogFile Label $header"
			}
		
		foreach line [set $var] {
			puts $FOUT "$FILENAME $instance_label [list_subset $line $fields $EC]"
			}
		
		}	
	
	puts stderr ""
	}

puts stderr "Closing handles"

foreach $FIELDS_LAYOUT $FIELDS {
	close [set OF_$var]
	}
