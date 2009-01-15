#!/bin/bash

for dir in test* ; do
echo Running $dir
(cd $dir ; ./run.sh)
done
