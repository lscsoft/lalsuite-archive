#!/usr/bin/env tclsh

#
# Usage: ./pick_job.tcl JOBID < condor_log_file
#

set data [read stdin]
set L [split [string map {... `} $data] `]

foreach par $L {
	if { [regexp $argv $par] } {

		puts $par
		}
	}



