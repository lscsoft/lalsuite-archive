#!/bin/sh
#PBS -l platform=LINUX     # Plateforme d'execution
#PBS -l T=30000            # Nombre d'unites normalisees (consommation cpu)
#PBS -l M=2048MB           # Memoire en MB
#PBS -l scratch=50MB       # Taille du scratch en MB
#PBS -l spool=50KB         # Taille du spool en KB

# CCIN2P3 specific script calling omega scan to analyze events in VSR1
# Call with
#      wscan_ccin2p3.sh event_time config_file output_directory
#
# where event time is the GPS time of the event in seconds
#       config_file is the qscan configuration file
#       output directory is the directory for outputing the html result
#
#     example : qscan-ccin2p3.sh 863559000 V1-raw-gravitational.txt test
#

if [ $# -lt 3 ]; then
  echo "Some of the arguments are missing, arguments needed are :"
  echo " - event GPS time"
  echo " - configuration file"
  echo " - output directory"
  exit
fi

eventTime=$1
configFile=$2
outputDirectory=$3

omegaDirectory="$THRONG_DIR/pro/omegadev/omega_r2062"
FFLFile="/afs/in2p3.fr/group/virgo/BKDB/VSR1/VSR1_raw.ffl"

# Set path for omega
testpath=`echo $PATH | grep -i 'omegadev/omega_r2062/bin'`

if [ -z $testpath ]; then
  export PATH=$omegaDirectory/bin:$PATH
fi

if [ -d $outputDirectory/$eventTime ]; then
  echo ""
  echo "Directory $outputDirectory/$eventTime exists already."
  echo "**********************"
  echo "*** Cleaning it... ***"
  echo "**********************"
  rm -rf $outputDirectory/$eventTime
fi

# Execute the wscanlite

OMEGASCAN="$omegaDirectory/bin/wpipeline scan -c $configFile -f $FFLFile -o $outputDirectory/$eventTime $eventTime"

echo "execute : $OMEGASCAN"
export LD_LIBRARY_PATH_SAV=${LD_LIBRARY_PATH}
source /usr/local/shared/bin/xrootd_env.sh
$OMEGASCAN
unset LD_PRELOAD
unset XROOTD_VMP
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH_SAV}

exit 0
