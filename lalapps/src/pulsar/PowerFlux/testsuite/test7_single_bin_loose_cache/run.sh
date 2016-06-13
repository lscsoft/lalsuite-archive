#!/bin/bash

PF=../../powerflux-${PFSUFFIX:=mt2}
PFVERSION=`$PF --version 2> /dev/null | cut -f 2 -d \  `
OUTPUT=output-$PFSUFFIX-$PFVERSION
REFERENCE=golden-$PFSUFFIX
PRECIOUS='(max_.*(ul|dx).*:)|initial|band_info'

echo $OUTPUT

rm -rf $OUTPUT
mkdir $OUTPUT

rm -rf $OUTPUT.bypassed
mkdir $OUTPUT.bypassed

echo "Running with cache bypassed"
$PF --num-threads=1 --config=config --bypass-powersum-cache=1 --output $OUTPUT.bypassed >& $OUTPUT.bypassed/err.txt

echo "Running with cache enabled"
$PF --num-threads=1 --config=config --output $OUTPUT >& $OUTPUT/err.txt

grep -E "$PRECIOUS" $REFERENCE/*/powerflux.log > golden_precious.txt
grep -E "$PRECIOUS" $OUTPUT.bypassed/*/powerflux.log > output_bypassed_precious.txt
grep -E "$PRECIOUS" $OUTPUT/*/powerflux.log > output_precious.txt 

diff golden_precious.txt output_precious.txt > changes.diff
diff output_bypassed_precious.txt output_precious.txt > changes_bypassed.diff

if [ "`wc -l changes.diff`" == "0 changes.diff" ] &&  [ "`wc -l changes_bypassed.diff`" == "0 changes_bypassed.diff" ] ; then
	exit 0
else
	echo "--------------------- changes.diff --------------------------"
	cat changes.diff
	echo "------------------ changes_bypassed.diff --------------------"
	cat changes_bypassed.diff
	exit -1
fi
