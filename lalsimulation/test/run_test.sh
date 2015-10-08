#!/bin/bash

HASH=$(git rev-parse HEAD)
HASH_ORIGINAL=4f360267440a71d315fd4625f2b5e7097a1e0112

rm GenerateSimulation
make GenerateSimulation

/usr/bin/time -p --output=timing-$HASH.txt ./GenerateSimulation --approximant SEOBNRv3 --spin1x 0 --spin1y 0 --spin1z 0.8 --spin2x 0.5 --spin2y 0 --spin2z 0.0 --m1 40 --m2 30 --outname simV3-$HASH.out --f-min 10 --distance 1 --inclination 0 --sample-rate 2048 --domain TD
$HOME/spec/specMaster/Support/bin/ComputeComplexPhase -A 2 <simV3-$HASH.out >simV3-$HASH-ampphase.out

# compute maximum phase difference
python ./max_phase_difference.py "simV3-$HASH_ORIGINAL-ampphase.out" "simV3-$HASH-ampphase.out" >>"simV3-$HASH-ampphase.out"
