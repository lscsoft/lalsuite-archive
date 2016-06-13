#!/usr/bin/env tclsh

foreach { type name } $argv {
	set FILE [open $name "r"]
	fconfigure $FILE -translation binary
	switch -exact $type {
		"shorts"   {
			binary scan [read $FILE] "s*" data
			}
		"ints"   {
			binary scan [read $FILE] "i*" data
			}
		"floats" {
			binary scan "[read $FILE]" "f*" data
			}
		"doubles" {
			binary scan [read $FILE] "d*" data
			}
		}
	puts [join $data "\n"]
	close $FILE
	}
