#!/usr/bin/env tclsh

source params.tcl

foreach $PARAMS_FORMAT $PARAMS {
        set $var [subst -nocommands $value]
        }

puts [set $argv]

