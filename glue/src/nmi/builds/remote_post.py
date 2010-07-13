#!/usr/bin/env python
#
# This script is responsible for packing up a completed build's
# installation dir or, if the build failed, any build artifacts
# needed for debugging.
#
# "Packing up" means simply creating a results.tar.gz file, which
# Metronome will subsequently look for and transfer back to the submit
# host automatically.

import os
import tarfile
import nmiOpts

# set up basic -v -q options
parser = nmiOpts.OptionParserInit()
(options, args) = parser.parse_args()

# if any part of the build failed, we want to pack up src/ for debugging
debugging_files = []
if os.getenv("_NMI_STEP_FAILED") is not None:
    debugging_files = [ "src" ]

tar = tarfile.open("results.tar.gz", "w:gz")
for name in [ "head" ] + debugging_files:
    tar.add(name)

if options.verbose:
    tar.list(verbose=False)

tar.close()
