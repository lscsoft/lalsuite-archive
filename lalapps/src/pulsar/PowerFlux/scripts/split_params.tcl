#!/usr/bin/env tclsh

set PARAMS { 
	RUN_NAME "S5.all.spindown0.run.1"
	ROOT_DIR "/home/volodya/runs/$RUN_NAME/"
	CONF_DIR "$ROOT_DIR/in/"
	OUTPUT_DIR "$ROOT_DIR/output/"
	ERR_DIR "$ROOT_DIR/err/"
	STATS_DIR "$ROOT_DIR/stats"
	STATS_SUFFIX ""
	DAG_LOG "/local/volodya/$RUN_NAME.condor.log"
	NBANDS	74

	FREQ_START	50
	FREQ_STEP	0.25
	FREQ_END	800

        FIRST_SPINDOWN  -1.00e-9
	FIRST_SPINDOWN 0
	LAST_SPINDOWN  0
        SPINDOWN_STEP   2e-11
        MAX_SPINDOWN_COUNT {1000*500*500/(\$band*\$band)}

        SKYMARKS {
                /home/volodya/S5.targeted/all_sky_marks.txt
                }


	ANALYSIS_PROGRAM	"/home/volodya/PowerFlux/powerflux-condor.1.4.26-64"
	}

set PARAMS_FORMAT { var value }	

proc GET_DATASET { firstbin } {
set fstart [expr round($firstbin/45000)*45000-2700]
return "/home/volodya/S5.split/dst/S5.all.${fstart}.dst"
}	

set POWERFLUX_CONF_FILE {
label "$RUN_NAME"
dataset $DATASET
input-format GEO
earth-ephemeris /home/volodya/detresponse/earth05-09.dat
sun-ephemeris /home/volodya/detresponse/sun05-09.dat
sky-marks $skymarks
first-bin $firstbin
nbins 501
do-cutoff 1
filter-lines 1
averaging-mode one 
#band-axis explicit(0.979306,-0.185684,-0.0805064)
#band-axis-norm 1.97319e-11
#large-S 3.08e-09
#band-axis explicit(0.98,-0.19,-0.08)
#band-axis-norm 1.97e-11
#large-S 1.32e-09
spindown-start $spindown_start
spindown-step  $SPINDOWN_STEP
spindown-count $spindown_count
spindown-start-time 846885755 
write-dat NONE
write-png NONE
subtract-background 1
ks-test 1
compute-betas 0
no-secondary-skymaps 1
output-initial 1
max-candidates 0
output $OUTPUT_DIR/$i
        }

set CONDOR_FILE {
universe=standard
executable=$ANALYSIS_PROGRAM
input=/dev/null
output=$ERR_DIR/out.\$(PID)
error=$ERR_DIR/err.\$(PID)
arguments=--config=$CONF_DIR/\$(PID)
log=$DAG_LOG
queue
}

