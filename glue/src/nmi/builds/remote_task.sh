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
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/lscsoft/libframe/lib64
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
export LSCSOFT_SRCDIR=$BASE_DIR/src/lscsoft
export LSCSOFT_ROOTDIR=$BASE_DIR/head
export LAL_PREFIX=$LSCSOFT_ROOTDIR/opt/lscsoft/lal
export LALAPPS_PREFIX=$LSCSOFT_ROOTDIR/opt/lscsoft/lalapps

mkdir -p ${LSCSOFT_SRCDIR}
cd ${LSCSOFT_SRCDIR}

# if local repo doesn't already exist, create it
if [[ ! -d ${LSCSOFT_SRCDIR}/lalsuite/.git ]]; then
    cd ${LSCSOFT_SRCDIR}
#    git clone git://ligo-vcs.phys.uwm.edu/lalsuite.git
    git clone /home/cbc/nmi/src/lalsuite/.git
fi

cd ${LSCSOFT_SRCDIR}/lalsuite
git checkout -f $NMI_git_id
git clean -dqfx

mkdir -p ${LSCSOFT_ROOTDIR}/etc
echo "export LSCSOFT_LOCATION=${LSCSOFT_ROOTDIR}/opt/lscsoft" > ${LSCSOFT_ROOTDIR}/etc/lscsoftrc

### LAL ###

# temporary debugging line
echo checking "checking condor_config_val LIB"...
condor_config_val LIB

for subsys in lal lalframe lalmetaio lalburst lalstochastic lalpulsar lalxml lalinspiral lalapps pylal glue; do
    if [[ ! -d ${LSCSOFT_SRCDIR}/lalsuite/$subsys ]]; then
	echo "Hmm, no $subsys source dir found; moving on..."
	continue;
    fi
    SUBSYS=$(echo $subsys | tr '[:lower:]' '[:upper:]')
    SUBSYS_PREFIX=$LSCSOFT_ROOTDIR/opt/lscsoft/$subsys
    eval export ${SUBSYS}_PREFIX=$LSCSOFT_ROOTDIR/opt/lscsoft/$subsys
    cd ${LSCSOFT_SRCDIR}/lalsuite/$subsys
	
    set +u
    source ${LSCSOFT_ROOTDIR}/etc/lscsoftrc
    set -u

    if [[ $subsys == "pylal" || $subsys == "glue" ]]; then
	python setup.py install --prefix=$SUBSYS_PREFIX
    else
	./00boot
# TODO: should the --enable-condor flag be passed to the new lalapps_inpsiral subsys instead of (or in addition to) lalapps?
	if [[ $subsys == "lalapps" ]]; then
#	    ./configure --prefix=$SUBSYS_PREFIX --enable-condor --disable-gcc-flags --disable-debug
#	    ./configure --prefix=$SUBSYS_PREFIX --enable-condor --with-extra-ldflags="-static -all-static"
	    ./configure --prefix=$SUBSYS_PREFIX --enable-condor --with-extra-ldflags=-static
	else
#	    ./configure --prefix=$SUBSYS_PREFIX --disable-gcc-flags --disable-debug --with-extra-ldflags="-static"
	    ./configure --prefix=$SUBSYS_PREFIX
	fi
        # using -j is "cheating" for a condor job, but I've found that
        # builds spend so much time on disk i/o that -j2 typically
        # tends to consume <=1 CPU
#	make -j2 install
	make -j4 install
    fi

    echo "# setup $SUBSYS for development:  " >> ${LSCSOFT_ROOTDIR}/etc/lscsoftrc
    echo "export ${SUBSYS}_LOCATION=\$LSCSOFT_LOCATION/$subsys" >> ${LSCSOFT_ROOTDIR}/etc/lscsoftrc
    echo "if [ -f "\$${SUBSYS}_LOCATION/etc/${subsys}-user-env.sh" ]; then" >> ${LSCSOFT_ROOTDIR}/etc/lscsoftrc
    echo "  source \$${SUBSYS}_LOCATION/etc/${subsys}-user-env.sh" >> ${LSCSOFT_ROOTDIR}/etc/lscsoftrc
    echo "fi" >> ${LSCSOFT_ROOTDIR}/etc/lscsoftrc
done

# useful info to have logged for later debugging
env
${LSCSOFT_ROOTDIR}/opt/lscsoft/lal/bin/lal-version

# end of stdout/stderr-combining block
} 2>&1
