#!/usr/bin/env tclsh

set PARAMS {
	IFO "H1"
	DETECTOR "H"
	RUN_NAME "fake_random.1"
	ROOT_DIR "/archive/home/volodya/runs/$RUN_NAME/"
	CONF_DIR "$ROOT_DIR/in/"
	OUTPUT_DIR "$ROOT_DIR/output/"
	ERR_DIR "$ROOT_DIR/err/"
	STATS_DIR "$ROOT_DIR/stats"
	STATS_SUFFIX ".$IFO"
	LOG_FILE "/usr1/volodya/$RUN_NAME.log"
	FREQ_START  150.0
	FREQ_STOP   250.0
	FREQ_STEP    0.25
	INJECTIONS_PER_BAND  20
	DEC_RANGE  { -1.57 1.57 }
	RA_RANGE { 0 6.283 }
	ORIENT_RANGE { 0 0.785 }
	ORIENT_RANGE { 0 3.1415 }
	SPINDOWN_LOG10_RANGE { -1.5 0}
	SPINDOWN_MAX 2e-9
	POWER_LOG10_RANGE { -1.5 0}
	POWER_MAX 4e-23
	ANALYSIS_PROGRAM "/archive/home/volodya/PowerFlux/powerflux"
	}

set PARAMS_FORMAT { var value }	

set POWERFLUX_CONF_FILE {
input /archive/home/volodya/SFT-3/S4.${IFO}.ht.geo/sft.${IFO}.
input-format GEO
detector L${DETECTOR}O
earth-ephemeris /archive/home/volodya/detresponse/earth05.dat
sun-ephemeris /archive/home/volodya/detresponse/sun05.dat
output $OUTPUT_DIR/$i
first-bin $first_bin
nbins 501
do-cutoff 1
three-bins 0
filter-lines 1
fake-linear
fake-ra	$ra
fake-dec $dec
fake-strain $power 
fake-freq $freq
fake-orientation $orientation
fake-spindown $spindown
write-dat NONE
write-png NONE
focus-ra $ra
focus-dec $dec
focus-radius 0.1
subtract-background 1
ks-test 1
	}
	
	
