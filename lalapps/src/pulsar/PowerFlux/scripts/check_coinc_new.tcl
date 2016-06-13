#!/usr/bin/env tclshexe

# Check for coincidence between pairs of IFO's from Vladimir's 
# single-ifo candidate lists. 

# Rewritten extensively October 2006 - K.Riles

# Constants

set pi 3.1415927
set twopi [expr 2.0*$pi]

# Set filetype list to scan

set filetypes [list AllPoints]
set filetypes [list RedPoints RedDiamonds BluePoints]

# Set pair of interferometers 

set ifolist [list H1 L1]

# Set file format (1 - original for S4)

set fileformat 1
    
# Set maximum allowed frequency difference

set deltaf0 0.010

# Set minimum max-SNR to consider a candidate

set dxthresh 7.0

# Set low / high frequencies for violin fundamental and next 3 harmonics

set flo_violin(1) 343     
set fhi_violin(1) 348
set nviolin_harmonics 4
for {set i 2} {$i <= $nviolin_harmonics} {incr i} {
    set flo_violin($i) [expr $i*$flo_violin(1)]
    set fhi_violin($i) [expr $i*$fhi_violin(1)]
}

# Set half-width of region to ignore around calibration lines (sideband issues)

set calib_halfwidth 0.75

# Set calibration line frequencies for each interferometer

set ncal(H1) 4
set fcalib(H1,1)   46.7
set fcalib(H1,2)  393.1
set fcalib(H1,3) 1144.3
set fcalib(H1,4)  973.3

set ncal(H2) 3
set fcalib(H2,1)   54.1
set fcalib(H2,2)  407.3
set fcalib(H2,3) 1159.7

set ncal(L1) 3
set fcalib(L1,1)   54.7
set fcalib(L1,2)  396.7
set fcalib(L1,3) 1151.5

# Compute resulting low/high ranges:

foreach ifo $ifolist {
    for {set i 1} {$i <= $ncal($ifo)} {incr i} {
	set flo_calib($ifo,$i) [expr $fcalib($ifo,$i) - $calib_halfwidth]
	set fhi_calib($ifo,$i) [expr $fcalib($ifo,$i) + $calib_halfwidth]
    }
}

# Set pulsar injection frequencies and halfwidths

set npulsar 12
set fpulsar(0) 265.577
set fpulsar(1) 849.070
set fpulsar(2) 575.164
set fpulsar(3) 108.857
set fpulsar(4) 1402.07
set fpulsar(5) 52.808
set fpulsar(6) 148.431
set fpulsar(7) 1220.93
set fpulsar(8) 193.938
set fpulsar(9) 763.847
set fpulsar(10) 501.239
set fpulsar(11) 376.070

# Set nominal half-widths to exclude near pulsars (maximum Doppler modulation)

for {set i 0} {$i < $npulsar} {incr i} {
    set pulsar_halfwidth($i) [expr $fpulsar($i)*1.1e-4]
}

# Adjust half-widths of pulsars 3 and 8 which are so strong that spectral leakage leads to
# candidates as far away as 0.02 and 0.08 Hz, respectively

set pulsar_halfwidth(3) 0.025
set pulsar_halfwidth(8) 0.085

# Adjust half-widths for two binary pulsars which have much larger modulations

set pulsar_halfwidth(10) 0.250
set pulsar_halfwidth(11) 0.350

# Compute resulting low/high ranges

for {set i 0} {$i < $npulsar} {incr i} {
    set flo_pulsar($i) [expr $fpulsar($i) - $pulsar_halfwidth($i)]
    set fhi_pulsar($i) [expr $fpulsar($i) + $pulsar_halfwidth($i)]
}

# Consolidate all frequency exclusion regions together

set nbad 0

for {set i 1} {$i <= $nviolin_harmonics} {incr i} {
    set index [expr $nbad+$i]
    set flo_bad($index) $flo_violin($i)
    set fhi_bad($index) $fhi_violin($i)
    set source_bad($index) "Violin harmonic "
    append source_bad($index) $i
}
set nbad [expr $nbad + $nviolin_harmonics]

foreach ifo $ifolist {
    for {set i 1} {$i <= $ncal($ifo)} {incr i} {
	set index [expr $nbad+$i]
	set flo_bad($index) $flo_calib($ifo,$i)
	set fhi_bad($index) $fhi_calib($ifo,$i)
	set source_bad($index) "Calibration line ("
	append source_bad($index) $ifo ")"
    }
    set nbad [expr $nbad + $ncal($ifo)]
}

for {set i 0} {$i < $npulsar} {incr i} {
    set index [expr $nbad+$i+1]
    set flo_bad($index) $flo_pulsar($i)
    set fhi_bad($index) $fhi_pulsar($i)
    set source_bad($index) "Injected pulsar "
    append source_bad($index) $i
}
set nbad [expr $nbad + $npulsar]

puts "Defined $nbad excluded bands:"
puts "  violin harmonics: $nviolin_harmonics"
set string "  calibration lines: "
foreach ifo $ifolist {
    append string $ifo "(" $ncal($ifo) ") "
}
puts $string
puts "  injected pulsars: $npulsar\n"

puts "Excluded bands:"
for {set i 1} {$i <= $nbad} {incr i} {
    puts "  $flo_bad($i) $fhi_bad($i)"
}

# Set frequency range of coincidence search 
# (low = 50-225 Hz with 51 spindowns; high = 200-1000 Hz with 11 spindowns)

set freqrange low
# set freqrange high

# Open log file for good matches

set tag "_2006_10_25"
set matchfilename "./matchfile_"
append matchfilename $freqrange $tag ".txt"
set matchfile [open $matchfilename w]
puts "Opened matches output file: $matchfilename"

# Set directories to look in for H1 and L1 files:

if { [string equal $freqrange "low"] } {
    
    set ifodir(1) "/home/volodya/public_html/Restricted/S4/UL.run.2/S4-H1-short-2-report/"
    set ifodir(2) "/home/volodya/public_html/Restricted/S4/UL.run.2/S4-L1-short-1d-report/"
    
} elseif { [string equal $freqrange "high"] } {
    
    set ifodir(1) "/home/volodya/public_html/Restricted/S4/UL.run.2/S4-H1-long-2-report/"
    set ifodir(2) "/home/volodya/public_html/Restricted/S4/UL.run.2/S4-L1-long-2-report/"
    
} else {
    puts "Oops, invalid freqency range: $freqrange"
    exit 1
}

# Set skyband to check for coincidence

set skyband 0

# Read in and save records for IFO 1

foreach type1 $filetypes {
 
    # Open IFO 1 input file
   
    set ifo1fname $ifodir(1)
    append ifo1fname $type1 "_" $skyband ".txt"
    set ifo1fid [open $ifo1fname]
    puts "Opened IFO 1 input file: $ifo1fname"
    puts $matchfile "Opened IFO 1 input file: $ifo1fname"

    # Initialize IFO 1 list and counters 

    set list1 ""
    set i1 0
    set i1all 0
    gets $ifo1fid line

    # Parse next IFO 1 file record

    while { [gets $ifo1fid line] >= 0 } {
	if { [regexp {^\"\d+\" +([\de.-]+) +([\de.-]+) +([\de.-]+) +([\de.-]+) +([\de.-]+) +([\de.-]+) +([\de.-]+) +([\de.-]+) +([\de.-]+) +([\de.-]+) +\"(.*)\" +(.*)$} $line match fstart1 linul1 circul1 resmax1 fdot1 maxdx1 ra1 dec1 f01 ul1 pol1 thetacos1] } {
	    incr i1all
	    if { [expr $maxdx1>$dxthresh] } {
		incr i1
		lappend list1 [list $fstart1 $linul1 $circul1 $resmax1 $fdot1 $maxdx1 $ra1 $dec1 $f01 $ul1 $pol1 $thetacos1]
	    }
	} elseif {[regexp {^.*?NA NA NA NA NA NA NA NA NA NA NA NA.*$} $line match] } {
	    puts "Skipping NA line"
	} elseif {[regexp {^.*?NA.*$} $line match] } {
	    puts "Skipping partial NA line"
	} else {
	    puts "Unrecognized record: /$line/"
	    exit 1
	}
    }

    # Print out summary
    
    puts "Read in $i1 records with maxdx>$dxthresh out of $i1all total records"
    puts $matchfile "Read in $i1 records with maxdx>$dxthresh out of $i1all total records"
    set list1 [lsort -index 8 -real $list1]
    puts "Sorted the $i1 records"

    # Read in and save records for IFO 2
    
    foreach type2 $filetypes {
 
	# Open IFO 2 input file
	
	set ifo2fname $ifodir(2)
	append ifo2fname $type2 "_" $skyband ".txt"
	set ifo2fid [open $ifo2fname]
	puts "Opened IFO 2 input file: $ifo2fname"
	puts $matchfile "Opened IFO 2 input file: $ifo2fname"

	# Initialize IFO 2 list and counters 

	set list2 ""
	set i2 0
	set i2all 0
	gets $ifo2fid line

	# Parse next IFO 2 file record
	
	while { [gets $ifo2fid line] >= 0 } {
	    if { [regexp {^\"\d+\" +([\de.-]+) +([\de.-]+) +([\de.-]+) +([\de.-]+) +([\de.-]+) +([\de.-]+) +([\de.-]+) +([\de.-]+) +([\de.-]+) +([\de.-]+) +\"(.*)\" +(.*)$} $line match fstart2 linul2 circul2 resmax2 fdot2 maxdx2 ra2 dec2 f02 ul2 pol2 thetacos2] } {
		incr i2all
		if { [expr $maxdx2>$dxthresh] } { 
		    incr i2
		    lappend list2 [list $fstart2 $linul2 $circul2 $resmax2 $fdot2 $maxdx2 $ra2 $dec2 $f02 $ul2 $pol2 $thetacos2]
		}
	    } elseif {[regexp {^.*?NA NA NA NA NA NA NA NA NA NA NA NA.*$} $line match] } {
		puts "Skipping NA line"
	    } elseif {[regexp {^.*?NA.*$} $line match] } {
		puts "Skipping partial NA line"
	    } else {
		puts "Unrecognized record: /$line/"
		exit 2
	    }
	}

	# Print out summary
	
	puts "Read in $i2 records with maxdx>$dxthresh out of $i2all total records"
	puts $matchfile "Read in $i2 records with maxdx>$dxthresh out of $i2all total records"
	set list2 [lsort -index 8 -real $list2]
	puts "Sorted the $i2 records"
	
	
	# Start search for coincidence
	
	set matchcnt 0
	set ii2start 0
	
	# Loop over sorted IFO 1 records
	
	for {set ii1 0} {$ii1<$i1} {incr ii1} {
	    set done 0
	    
	    # Get frequency of next IFO 1 candidate
	    
	    set record1 [lindex $list1 $ii1]
	    set f01 [lindex $record1 8]
	    set fdot1 [lindex $record1 4]
	    
	    # Set start frequency in IFO 2 and go to nearest frequency in IFO 2

	    set f02start [expr $f01-1.01*$deltaf0]
	    set foundstart 0
	    while { [expr !$foundstart && $ii2start>0] } {
		incr ii2start -1
		if { [expr $ii2start<0] } {
		    set ii2start 0
		}
		set record2 [lindex $list2 $ii2start]
		set f02 [lindex $record2 8]
		if { [expr $f02<$f02start] } {
		    set foundstart 1
		}
	    }

	    for {set ii2 $ii2start} {$ii2<$i2 && !$done} {incr ii2} {
		set record2 [lindex $list2 $ii2]
		set f02 [lindex $record2 8]

	    # Look for match in frequency

		if { [expr $f02-$f01>$deltaf0] } {
		    set done 1
		    set ii2start $ii2
		} elseif { [expr abs($f01-$f02) <= $deltaf0] } {
		    puts "Match: IFO 1 $type1 $ii1 at f0=$f01 matches with IFO 2 $type2 $ii2 at f0=$f02"
		    puts "  IFO1 record = $record1"
		    puts "  IFO2 record = $record2"
		    set fdot2 [lindex $record2 4]
		    set df0 [expr $f02-$f01]
		    set ra1 [lindex $record1 6]
		    set dec1 [lindex $record1 7]
		    set ra2 [lindex $record2 6]
		    set dec2 [lindex $record2 7]
		    set dra [expr $ra2-$ra1]
		    set absdra abs($dra)
		    if { [expr $absdra>$twopi] } {
#			puts "Resetting absdra from old=$absdra (>2pi) to new:"
			set absdra [expr $absdra-$twopi]
#			puts "$absdra where dra=$dra and ra1=$ra1, ra2=$ra2"
		    }
		    if { [expr $absdra>$pi] } {
#			puts "Resetting absdra from old=$absdra (>pi) to new:"
			set absdra [expr $twopi - $absdra]
#			puts "$absdra where dra=$dra and ra1=$ra1, ra2=$ra2"
		    }
		    set ddec [expr $dec2-$dec1]
		    set absddec [expr abs($ddec)]
		    set dfdot [expr $fdot2-$fdot1]
		    puts "  f0 diff = $df0 Hz; RA diff = $dra rad; decl diff = $ddec rad fdot diff = $dfdot"
		    incr matchcnt

		    # Check if in excluded band

		    set found 0
		    for {set i 1} {$i <= $nbad} {incr i} {
			if { [expr ($f01>=$flo_bad($i) && $f01<=$fhi_bad($i)) || ($f02>=$flo_bad($i) && $f02<=$fhi_bad($i))] } {
			    set source $source_bad($i)
			    puts "Likely source: $source"
			    set found 1
			}
		    }

		    # Check sky location match

		    if { [expr $absdra<0.5 && $absddec<0.5] } {
			set skypass 1
		    } else {
			set skypass 0
		    }

		    # Check spindown match

		    if { [expr abs($dfdot)<1.1E-9] } {
			set fdotpass 1
		    } else {
			set fdotpass 0
		    }

		    # Print out "golden" candidatest

		    if { [expr $skypass>0 && $fdotpass>0] } {
			puts "Passes BOTH sky and spindown match requirement\n"
			puts $matchfile "Match: IFO 1 $type1 $ii1 at f0=$f01 matches with IFO 2 $type2 $ii2 at f0=$f02"
			puts $matchfile "  IFO1 record = $record1"
			puts $matchfile "  IFO2 record = $record2"
			puts $matchfile "  f0 diff = $df0 Hz; RA diff = $dra rad; decl diff = $ddec rad fdot diff = $dfdot"
			if { [expr $found > 0] } {
			    puts $matchfile "Likely source: $source"
			}
			puts $matchfile "Passes BOTH sky and spindown match requirement\n"
			if { [expr $found==0] } {
			    puts "WHAT IS THIS?\n"
			    puts $matchfile "WHAT IS THIS?\n"
			}
		    } elseif { [expr $skypass>0] } {
			puts "Passes sky but fails spindown match requirment\n"
		    } elseif { [expr $fdotpass>0] } {
			puts "Passes spindown but Fails sky match requirment\n"
		    } else {
			puts "Fails both sky and spindown match requirement\n"
		    }
		    if { [expr $found==0] } {
			puts "\nWHAT IS THIS?\n"
		    }
		}
	    }
	}
	puts "Found $matchcnt matches"
	puts "****************************************************************************\n"
	puts $matchfile "Found $matchcnt matches"
	puts $matchfile "****************************************************************************\n"
    }
}
    
