#!/bin/bash

PF=../../powerflux-${PFSUFFIX:=mt}
PFVERSION=`$PF --version 2> /dev/null | cut -f 2 -d \  `
OUTPUT=output-$PFSUFFIX-$PFVERSION
REFERENCE=golden-$PFSUFFIX
PRECIOUS='max_.*(ul|dx).*:'

echo $OUTPUT

rm -rf $OUTPUT
mkdir $OUTPUT

$PF --config=config --output $OUTPUT >& $OUTPUT/err.txt

grep -E "$PRECIOUS" $REFERENCE/*/powerflux.log > golden_precious.txt
grep -E "$PRECIOUS" $OUTPUT/*/powerflux.log > output_precious.txt 

diff golden_precious.txt output_precious.txt > changes.diff

if [ "`wc -l changes.diff`" == "0 changes.diff" ] ; then
	exit 0
else
	cat changes.diff
	exit -1
fi
