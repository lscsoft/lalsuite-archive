#!/usr/bin/env tclsh

set cwd [pwd]

set keep 0

set PIPE ""

foreach dir $argv {
	if { $dir == "" } {
		puts stderr "Skipping empty directory"
		continue
		}

	if { $dir == "--keep" } {
		set keep 1
		continue
		}
	if { $dir == "--remove" } {
		set keep 0
		continue
		}

	if { $PIPE == "" } {
		set PIPE [open "|at now+1min >& /dev/null" "w"]
		puts $PIPE "cd '$cwd'"
		}
	if { $keep } {
		puts -nonewline $PIPE "tar zcf '$dir'.tgz '$dir'/ && tar tf '$dir'.tgz > '$dir'.lst"
		puts -nonewline $PIPE " && echo 'Done archiving directory \"$dir\". KEPT original files.'"
		puts $PIPE " || echo 'FAILED archiving directory \"$dir\"'"
		} {
		puts -nonewline $PIPE "tar zcf '$dir'.tgz '$dir'/ && tar tf '$dir'.tgz > '$dir'.lst && rm -rf '$dir'/"
		puts -nonewline  $PIPE " && echo 'Done archiving directory \"$dir\", original files REMOVED'"
		puts $PIPE " || echo 'FAILED archiving directory \"$dir\"'"
		}
	puts stderr "Scheduled archiving of '$dir'"
	}

if { $PIPE != "" } { 
	close $PIPE
	} {
	puts stderr "\nUsage: archive.tcl \[options\] directory \[ \[options\] directory ..\]"
	puts stderr "Options:"
	puts stderr "\t--keep     keep original files"
	puts stderr "\t--remove   remove original files (default)\n"
	}
