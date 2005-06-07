#!/usr/bin/env tclsh

source params.tcl

foreach $PARAMS_FORMAT $PARAMS {
        set $var [subst -nocommands -nobackslashes $value]
        }

puts [set $argv]

