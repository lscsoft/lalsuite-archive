#!/usr/bin/env tclsh

foreach {var value} {
	FILENAME ""
	KIND "snr"
	FREQ_TOLERANCE 0.005
	DIST_TOLERANCE 0.07
	SPINDOWN_TOLERANCE 5e-10
	MIN_ALL_SNR 7
	FBIN_LEFT_CUTOFF 10
	FBIN_RIGHT_CUTOFF 491
	} {
	global $var
	set $var $value
	}

foreach {var value} $argv {
	global $var
	set $var $value
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

set FILE [open $FILENAME "r"]

#set DATA {}
set SETS {}

set COSDIST_TOLERANCE [expr cos($DIST_TOLERANCE)]

while { ![eof $FILE] } {
	gets $FILE s
	if { [lindex $s $KIND_INDEX]!=$KIND } { continue }
	set fbin [lindex $s $FREQUENCY_BIN_INDEX]
	if { ($fbin <= $FBIN_LEFT_CUTOFF) || ($fbin >= $FBIN_RIGHT_CUTOFF) } { continue }
	
	set SET [lindex $s $SET_INDEX]

	if { [regexp "_all$" $SET] } {
		if { [lindex $s $SNR_INDEX] < $MIN_ALL_SNR } { continue }
		} {
		lappend SETS $SET
		}

	set idx [expr round([lindex $s $FREQUENCY_INDEX]/$FREQ_TOLERANCE)]
	lappend DATA(${SET}_$idx) $s
	}

close $FILE

set SETS [lsort -unique $SETS]

#puts $SETS
#puts [array get DATA]

set DET_FOUND {}

foreach tag [array names DATA -regexp "_all_"] {
	regexp {^(.*)_all_(.*)} $tag {} ROOT BIN
	#puts "$ROOT $BIN"
	set DET_ROOTS [lsearch -all -inline -regexp $SETS "^${ROOT}_"]

	foreach line $DATA($tag) {

		set frequency [lindex $line $FREQUENCY_INDEX]
		set spindown [lindex $line $SPINDOWN_INDEX]
		set ra [lindex $line $RA_INDEX]
		set dec [lindex $line $DEC_INDEX]

		set LOCAL_DET_FOUND {}

		set found_roots 0
		foreach root $DET_ROOTS {
			set found 0
			#set DET_FOUND($root) {}
			foreach bin [list $BIN [expr $BIN-1] [expr $BIN+1]] {
				if { ! [info exists DATA(${root}_${bin})] } { continue }
				foreach s $DATA(${root}_${bin}) {
					if { abs($frequency-[lindex $s $FREQUENCY_INDEX])>$FREQ_TOLERANCE} { continue }
					if { abs($spindown-[lindex $s $SPINDOWN_INDEX])>$SPINDOWN_TOLERANCE} { continue }

					set ra2 [lindex $s $RA_INDEX]
					set dec2 [lindex $s $DEC_INDEX]

					if { sin($dec)*sin($dec2)+cos($dec)*cos($dec2)*cos($ra-$ra2) < $COSDIST_TOLERANCE } { continue }
					
					set found 1
					lappend LOCAL_DET_FOUND $s
					}
				}
			if { $found } { incr found_roots }
			}
		if { $found_roots == [llength $DET_ROOTS] } {
			puts $line
			foreach s $LOCAL_DET_FOUND { lappend DET_FOUND $s }
			}
		}
	}

if { [llength $DET_FOUND] > 0 } {
	puts [join [lsort -unique $DET_FOUND] "\n"]
	}
