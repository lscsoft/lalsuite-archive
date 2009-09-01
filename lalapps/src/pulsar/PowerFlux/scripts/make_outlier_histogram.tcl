#!/usr/bin/env tclsh

foreach {var value} {
	FILENAME ""
	OUTPUT_HEADER_ONLY 0
	KIND "snr"
	FBIN_LEFT_CUTOFF 10
	FBIN_RIGHT_CUTOFF 491

	SNR_LEFT 5.0
	SNR_RIGHT 7.0
	SNR_STEP 0.25	
	} {
	global $var
	set $var $value
	}

foreach {var value} $argv {
	global $var
	set $var $value
	}

if { $OUTPUT_HEADER_ONLY } {
	puts "label\tset\tskyband\tindex\tcount"
	exit 0
	}

set HEADER {kind label index set pi pps_count template_count first_bin min_gps max_gps skyband frequency spindown ra dec iota psi snr ul ll M S ks_value ks_count m1_neg m3_neg m4 frequency_bin max_weight weight_loss_fraction max_ks_value max_m1_neg min_m1_neg max_m3_neg min_m3_neg max_m4 min_m4 max_weight_loss_fraction}

set KIND_INDEX [lsearch -exact $HEADER "kind"]
set SET_INDEX  [lsearch -exact $HEADER "set"]
set SNR_INDEX  [lsearch -exact $HEADER "snr"]
set FREQUENCY_INDEX [lsearch -exact $HEADER "frequency"]
set FREQUENCY_BIN_INDEX [lsearch -exact $HEADER "frequency_bin"]
set SPINDOWN_INDEX [lsearch -exact $HEADER "spindown"]
set RA_INDEX [lsearch -exact $HEADER "ra"]
set DEC_INDEX [lsearch -exact $HEADER "dec"]
set SKYBAND_INDEX [lsearch -exact $HEADER "skyband"]
set LABEL_INDEX [lsearch -exact $HEADER "label"]


set FILE [open $FILENAME "r"]



set SNR_MAX_INDEX [expr round(ceil(($SNR_RIGHT-$SNR_LEFT)/$SNR_STEP)+1)]

while { ![eof $FILE] } {
	gets $FILE s
	if { [lindex $s $KIND_INDEX]!=$KIND } { continue }
	set fbin [lindex $s $FREQUENCY_BIN_INDEX]
	if { ($fbin <= $FBIN_LEFT_CUTOFF) || ($fbin >= $FBIN_RIGHT_CUTOFF) } { continue }
	
	set SET [lindex $s $SET_INDEX]
        set SKYBAND [lindex $s $SKYBAND_INDEX]
        set SNR [lindex $s $SNR_INDEX]
        set LABEL [lindex $s $LABEL_INDEX]

	set index [expr round(ceil(($SNR-$SNR_LEFT)/$SNR_STEP))]
	if { $index < 0 } { set index 0 }
	if { $index > $SNR_MAX_INDEX } { set index $SNR_MAX_INDEX }

	set id [list $LABEL $SET $SKYBAND $index]

	if { [info exists DATA($id)] } {
		set a $DATA($id)
		set DATA($id) [expr $a+1]
		} {
		set DATA($id) 1
		}

	}

close $FILE

foreach {var value} [array get DATA ] {
	lappend var $value
	puts [join $var "\t"]
	}

