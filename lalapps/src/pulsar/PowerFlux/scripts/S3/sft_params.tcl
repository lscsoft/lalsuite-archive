#!/usr/bin/env tclsh

# This is Tcl/Tk file meant to be loaded with make_1800_sfts script
#
foreach {var value} {
	sft_kind "ht"
	data_set "S3"
        interferometer "H1"
        instrument "H"
	channel "${interferometer}:Calibrated-Strain"
        storage_dir "/scratch4/volodya"
        frame_library "/scratch4/volodya/${data_set}.${interferometer}.${sft_kind}.txt"
        sc_files  "/home/volodya/${data_set}/sc.tcl"
        control_info_dir   "${storage_dir}/sfts_control/${data_set}.${interferometer}.${sft_kind}.geo/"
	config_dir "$control_info_dir/in/"
	err_dir "$control_info_dir/err/"
	dag_file "$control_info_dir/dag"
	submit_file "$control_info_dir/submit"
	log_file "/people/volodya/${data_set}.${interferometer}.${sft_kind}.log"
	sfts_dir "$storage_dir/SFT-3/${data_set}.${interferometer}.${sft_kind}.geo/"
	group_regexp {/([0-9]*)-([0-9]*)/.-.._RDS_C0._LX}
	filename_regexp {(/netdat./.*)$}
	epoch_regexp {-([0-9]*)-16.gwf}
        timebase 1800
        overlap 900
	seg_start 0
	seg_step    20499
	frame_length  16
	make_sft_prog "make_sft_plain"
	veto_expr {~0}
	non_veto_expr {$FLAG_NO_CALIB_LINE_V01 | $FLAG_AIRPLANE | $FLAG_DUST | $FLAG_SEISMIC_TRANSIENT | $FLAG_INJECTION_STOCHASTIC | $FLAG_HUMAN_INTRUSION | $FLAG_SEISMIC_ELEVATED | $FLAG_SEISMIC_HIGH}
        } {
        global $var
        set $var $value
        }

# Process command line arguments
foreach {var value} $argv {
        global $var
        set $var $value
        }

# Expand variables that depend on other variables
foreach {var} {frame_library control_info_dir 
	config_dir err_dir dag_file submit_file
	sfts_dir sc_files log_file channel } {
        global $var
        set $var [subst -nocommands [set $var]]
        }


source $sc_files

set make_sft_config {
#
# Channel to process
#
CHANNEL "$channel"
SAMPLES_PER_SECOND 16384

#
# Specify duration, times are inclusive
#
DURATION $start $end

#
#
#  Calibration
#
#  We need three files: R file, C file and file with alpha beta constants
#  If necessary I can add the ability to specify alpha beta constants directly
#  within this file.
#
#  It is possible to specify more than one file of alpha beta constants -
#  they would be merged in this case. Not sure whether this is useful or not.
#
#

#DONT_FAIL_ON_MISSING_CALIBRATION
BYPASS_LINE_SUPPRESSION
#BYPASS_FIRST_FFT

#  Do not perform high-pass Butterworth filter
BYPASS_HIGHPASS_FILTER

#  Do not perform low-pass Butterworth filter
BYPASS_LOWPASS_FILTER


#R_FILE "$Response_file"
#C_FILE "$CavGain_file"
#ALPHA_BETA_FILE "$AlphaBeta_file"

#
#
# Specify highpass filter parameters: ORDER F A
# (see LAL documentation for Butterworth filter)
# Default parameters are 5 100.0 0.5
#
HIGHPASS_FILTER 5 20 0.5

#
# Specify lowpass filter parameters: ORDER F A
# (see LAL documentation for Butterworth filter)
# Default parameters are -1 -1 -1 - which means the filter is not applied
#
LOWPASS_FILTER  5 2500 0.5

#
#
# Windowing:
#
#  Data is windowed three times:
#    * once in time domain on input to calibration code
#    * once in frequency domain before converting calibration result to time domain
#    * once in time domain before applying large SFT to the calibrated time domain data
#
#  Window type: one of Hann, Welch (other types from LAL library can be added easily).
#
#  Window tail size - the program creates Hann window of twice the size
#  then splits it into two equal tails and applies them appropriately.
# (in the ends of time-domain data and in the ends of valid frequency-domain data)
#
WINDOW1_TAIL_SIZE 491520
WINDOW1_TYPE "Hann"

#WINDOW2_TAIL_SIZE 491520
#WINDOW2_TYPE "Hann"

#
# Which frequency band to write into output file
# OUTPUT_BAND
#
OUTPUT_BAND     0       2000.0

#
# File to write results into
#

#OUTPUT_POWER "$sfts_dir/power.${interferometer}.$num"
OUTPUT_SFT "$sfts_dir/sft.${interferometer}.$num"


#
# Output mode.
#  OUTPUT_MODE_TEXT writes data as ascii file
#  OUTPUT_MODE_BINARY writes comments as text, but actual data as binary
#
# Debug output is always using OUTPUT_MODE_TEXT
# Default is OUTPUT_MODE_TEXT
#

OUTPUT_MODE_GEO

#
#
# Precision of decimal representation of floating point numbers.
# Used in debug output and when OUTPUT_MODE_TEXT is used.
#
# If not specified precision is 6
#

PRECISION 6

#
# Valid segments
#
# Segment specification is inclusive
#
# It can be either second granularity (and thus the entire last second is marked)
# or subsecond granularity (and thus the last sample is marked)
#
#
SEGMENTS
$start $end

#
#
# Frame files
#
FRAME_FILES
}
