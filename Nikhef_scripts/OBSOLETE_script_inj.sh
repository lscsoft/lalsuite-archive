#!/bin/bash

test_par=dphi7
shift=1 ####### THIS IS IN PERCENTS!!
seed=7000
outname=injections_$test_par'_'$shift'pc_'$seed.xml

lalapps_inspinj \
--output $outname \
--f-lower 20.0 \
--gps-start-time 932170000 \
--gps-end-time 932765900 \
--seed $seed \
--waveform IMRPhenomFBTestthreePointFivePN \
--min-distance 3.00e+05 \
--max-distance 1.25e+06 \
--d-distr volume \
--l-distr random \
--i-distr uniform \
--min-mass1 5.0 \
--max-mass1 15.0 \
--min-mass2 5.0 \
--max-mass2 15.0 \
--m-distr componentMass \
--min-mtotal 10.0 \
--max-mtotal 30.0 \
--disable-spin \
--time-step 1.0e+03 \
--amp-order 0 \
--enable-dphi \
--$test_par `echo $(($shift/100.0))`

seed_line=`cat $outname | grep seed`
temp=`echo ${seed_line%\"*}`
temp2=`echo ${temp##*\"}`



echo "The file was created with seed $temp2 "

