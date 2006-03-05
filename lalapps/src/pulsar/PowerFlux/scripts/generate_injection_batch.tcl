#!/usr/bin/env tclsh

foreach {var value} {
	gps_start	793154935
	gps_stop	793157935
	gps_stop	795679077
	gps_step	900
	spindown_gps_start 793154935
	band_start	149.0
	band		0.25
	band_extra	0.25
	detector	"H1"
	noise		"2e-24"
	ephemDir 	"/home/volodya/LIGO/LAL/lal/packages/pulsar/test/"
	ephemYear	"05-09"
	noiseGlob	"/home/volodya/LIGO/SFT-DATA/S4/S4.H1.ht.geo/sft.*"
	outputDir	"./"
	DEC 0.0
	RA 0.0
	psi 0.0
	phi 0.0
	f0 149.125
	aPlus 3e-24
	aCross 3e-24
	spindown 0
	} {
	global $var
	set $var $value
	}

#
# Allow to override globals from command line
#	
foreach {var value} $argv {
	global $var
	set $var $value	
	}

proc get_geo_range { SFT bin_start nbins } {
set FILE [open $SFT "r"]
fconfigure $FILE -encoding binary -translation binary
set header [read $FILE 32]
binary scan $header diidii  key gps nsec timebase bin_start0 nbins0
puts stderr "$SFT header: gps=$gps timebase=$timebase bin_start=$bin_start nbins=$nbins"
if { $bin_start0 > $bin_start } {
	puts stderr "Could not find range: $bin_start0 vs $bin_start"
	exit -1
	}
seek $FILE [expr ($bin_start-$bin_start0)*8] current
set data_raw [read $FILE [expr $nbins*8]]
binary scan $data_raw f* data
close $FILE
return [list gps $gps timebase $timebase data $data]
}

proc apply_hann_noise { INPUT OUTPUT noise_data} {
set FILE [open	$INPUT "r"]
fconfigure $FILE -encoding binary -translation binary
set header [read $FILE 32]
set body [read $FILE [expr [llength $noise_data]*4]]
close $FILE

binary scan $header diidii  key gps nsec timebase bin_start nbins
binary scan $body f*  bins

puts stderr "makefakedata header: gps=$gps timebase=$timebase bin_start=$bin_start nbins=$nbins"

# Output input bins for testing
#puts [join $bins "\n"]

#
# Hann window is equivalent to applying averaging filter
# with coefficients of -1/4 1/2 -1/4
# 

#
# We clobber first and last bins - or well..
#
set bins_out {0 0}

set x0	0
set y0 	0
set x1	[lindex $bins 0]
set y1	[lindex $bins 1]
foreach {x y} [lrange $bins 2 end] {xn yn} [lrange $noise_data 1 end-1] {
	lappend bins_out [expr -$x0/4.0+$x1/2.0-$x/4.0 + $xn]
	#puts "$x $y"
	lappend bins_out [expr -$y0/4.0+$y1/2.0-$y/4.0 + $yn]
	set x0 $x1
	set x1 $x
	set y0 $y1
	set y1 $y
	}
lappend bins_out 0 0

# Output filtered bins for testing
#puts [join $bins_out "\n"]



set FILE [open $OUTPUT "w"]
fconfigure $FILE -encoding binary -translation binary
puts -nonewline $FILE $header
puts -nonewline $FILE [binary format f* $bins_out]
close $FILE
}
	

foreach filename [glob $noiseGlob] {
	set bin_start [expr round(floor($band_start-$band_extra)*1800)]
	set nbins [expr round(floor($band_start+$band+2*$band_extra)*1800)-$bin_start]
	set sft_name ${outputDir}/inj.[file tail ${filename}]

	foreach {var value} [get_geo_range $filename $bin_start $nbins] {
		set "noise_$var" $value
		}

	exec ./lalapps_Makefakedata.static \
		--fmin=[expr $bin_start/1800.0] \
		--Band=[expr $nbins/1800.0] \
		--duration=1800 \
		--startTime=$noise_gps \
		--outSFTbname=$sft_name \
		--detector=$detector \
		--ephemDir=$ephemDir \
		--refTime=$spindown_gps_start \
		--Tsft=1800 \
		--longitude=$DEC \
		--latitude=$RA \
		--psi=$psi \
		--phi0=$phi \
		--f1dot=$spindown \
		--f0=$f0 \
		--aPlus=$aPlus \
		--aCross=$aCross \
		--outSFTv1=1 \
		--ephemYear=$ephemYear
	apply_hann_noise "${sft_name}.00000" $sft_name	$noise_data
	break
	}
