#!/usr/bin/env tclsh

foreach {var value} {
	SITE	"H"
	IFO	"H1"
	DIRECTORY	"/home/volodya/LIGO/S4/Injections/sfts/108/"
	PREV_SFT_TEMPLATE	"*sft*"
	NEW_SFT_TEMPLATE {${SITE}-1_${IFO}_1800SFT_v1C02hoftHann-${gps}-1800.sft}
	} {
	global $var
	set $var $value
	}

foreach {var value} $argv {
	global $var
	set $var $value
	}

proc get_header { name } {
set FILE [open $name "r"]
fconfigure $FILE -encoding binary -translation binary
set header [read $FILE 32]
binary scan $header diidii  key gps nsec timebase bin_start0 nbins0
close $FILE
return [list gps $gps timebase $timebase nbins $nbins0]
}

set files [glob -directory $DIRECTORY $PREV_SFT_TEMPLATE]
foreach F $files {
	foreach {var value} [get_header $F] { set $var $value }
	set NEW_NAME [subst $NEW_SFT_TEMPLATE]
	puts "$NEW_NAME"
	file rename $F [file join $DIRECTORY [subst $NEW_NAME]]
	}
