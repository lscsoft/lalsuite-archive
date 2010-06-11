#!/bin/sh
# generates the metronome submit files to build lalsuite for a given git tag

TMPDIR=${TMPDIR:-/tmp}

# exit immediately if any command exits with a non-zero status.
set -e
# treat unset variables as an error when performing parameter expansion.
set -u

if [[ $# -lt 2 ]]; then
    echo usage: $0 git_id git_branch_name
    exit 1
fi

# env vars of the form _NMI_* can be used in the Metronome submit files
export _NMI_GIT_ID=$1
export _NMI_GIT_BRANCH=$2

# create a temp dir for the submit files
SUBMIT_DIR=$TMPDIR/$USER/$(basename $0).${_NMI_GIT_ID}.$$.$(date +%s)
mkdir -p $SUBMIT_DIR
cd $SUBMIT_DIR

# make a copy of Metronome run-specification (aka submit) and
# input-specification files
# TODO: these should probably come from the current lscsoft install
# (e.g., $GLUE_PREFIX/lib/nmi) instead
cp $LSCSOFT_SRCDIR/lalsuite/glue/src/nmi/builds/{cmdfile,scripts.git} .

# print some handy debugging info
echo SUBMIT_DIR=$SUBMIT_DIR
echo LSCSOFT_SRCDIR=$LSCSOFT_SRCDIR
env | grep ^_NMI_

# TODO: add option to create the submit files and stop here, but not
# actually submit

# submit the build to Metronome
nmi_submit --no-wait cmdfile

# nasty kludge to work around condor submit throttle edge case
(sleep 21; condor_reschedule) > /dev/null 2>&1 &
disown %
