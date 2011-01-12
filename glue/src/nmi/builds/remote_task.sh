#!/bin/sh
# code that runs on the execute node to build lalsuite

# this script expects:
#
# git installed in default system PATH or /usr/local/bin
# condor installed in PATH, and CONDOR_CONFIG defined if necessary
# libframe and libmetaio installed in /opt/lscsoft
# lalapps installed in ./head/opt/lscsoft/lalapps
#
# $NMI_ligo_reference_xml defined (by Metronome, from our submit file)
# $NMI_component_version (git id) defined (by Metronome, from our submit file)
# $NMI_ligo_add_params defined (by Metronome, from our submit file)

# wrap this whole script in a block combining stdout & stderr to ease debugging
{
# TODO: these should be replaced with proper Metronome prereqs
export PKG_CONFIG_PATH=/usr/lib64/pkgconfig
export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:/opt/lscsoft/libframe/lib64/pkgconfig
export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:/opt/lscsoft/libmetaio/lib64/pkgconfig

# exit immediately if any command exits with a non-zero status.
set -e
# treat unset variables as an error when performing parameter expansion.
set -u
# print (unexpanded) shell input lines as they are read
set -v
# print (expanded) commands before they're executed
set -x

if [ -d head ]; then
    echo Moving old install out of the way...
    mv -v head head.old.$$
fi

# most of what's below is adapted from
# <https://www.lsc-group.phys.uwm.edu/daswg/docs/howto/lal-install.html>
# (circa 2010-02-22)

BASE_DIR=$PWD
export LSCSOFT_SRCDIR=$BASE_DIR/src
export LSCSOFT_ROOTDIR=$BASE_DIR/opt/lscsoft
#export LAL_PREFIX=$LSCSOFT_ROOTDIR/lal
#export LALAPPS_PREFIX=$LSCSOFT_ROOTDIR/lalapps

mkdir -p ${LSCSOFT_SRCDIR}
cd ${LSCSOFT_SRCDIR}

# if local repo doesn't already exist, create it
if [[ ! -d ${LSCSOFT_SRCDIR}/lalsuite/.git ]]; then
    cd ${LSCSOFT_SRCDIR}
#    git clone /home/cbc/nmi/src/lalsuite/.git
    git clone $NMI_git_repo
fi

cd ${LSCSOFT_SRCDIR}/lalsuite
git checkout -f $NMI_git_id
#git clean -dqfx

#mkdir -p ${LSCSOFT_ROOTDIR}/etc
#echo "export LSCSOFT_LOCATION=${LSCSOFT_ROOTDIR}/opt/lscsoft" > ${LSCSOFT_ROOTDIR}/etc/lscsoftrc

# temporary debugging line
echo "condor_config_val LIB =" $(condor_config_val LIB)

./00boot
./configure --prefix=$LSCSOFT_ROOTDIR --enable-lalxml
make -j4 install

cd glue
# nasty kludge to work around glue install bug known to exist in
# s6abc_lowmass tag 85f56e4a1555a60fe2ee98dde0b6e22afface3ad
if [[ -L glue/misc/generate_vcs_info.py ]]; then
    rm glue/misc/generate_vcs_info.py
    cp lal/lib/generate_vcs_info.py glue/misc/
    touch glue/misc/__init__.py
fi
python ./setup.py install --prefix=$LSCSOFT_ROOTDIR
cd -

cd pylal
# we have to disable error-on-unbound-var before sourcing lal-user-env.sh (or else it fails)
set +u
. $LSCSOFT_ROOTDIR/etc/lal-user-env.sh
set -u
python ./setup.py build
python ./setup.py install --prefix=$LSCSOFT_ROOTDIR
cd -

# useful info to have logged for later debugging
env
${LSCSOFT_ROOTDIR}/bin/lal-version

# end of stdout/stderr-combining block
} 2>&1
