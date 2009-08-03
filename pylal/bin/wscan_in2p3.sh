#!/bin/sh
#PBS -l platform=LINUX     # Plateforme d'execution
#PBS -l T=30000            # Nombre d'unites normalisees (consommation cpu)
#PBS -l M=2048MB           # Memoire en MB

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

webDirectory="$GROUP_DIR/www/followups"
if [ $# -lt 3 ]; then
  echo "Some of the arguments are missing, arguments needed are :"
  echo " - event GPS time"
  echo " - configuration file"
  echo " - output directory"
  echo " - optional : web subdirectory (subdirectory of $webDirectory )"
  exit
fi

eventTime=$1
configFile=$2
outputDirectory=$3
if [ $# -eq 4 ]; then
  webSubDirectory=$4
else
  webSubDirectory=""
fi

omegaDirectory="$THRONG_DIR/pro/omegadev/omega-nightly/test/install"
FFLFile="/afs/in2p3.fr/group/virgo/BKDB/VSR1/VSR1_raw.ffl"
#webDirectory="buskulic@olserver14.virgo.infn.it:/opt/w3/MonitoringWeb/OmegaEvents/"

# Set path for omega
testpath=`echo $PATH | grep -i 'omegadev/omega-nightly/test/install/bin'`

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

# Execute the wscan

OMEGASCAN="$omegaDirectory/bin/wpipeline scan -c $configFile -f $FFLFile -o $outputDirectory $eventTime"

echo "execute : $OMEGASCAN"
export LD_LIBRARY_PATH_SAV=${LD_LIBRARY_PATH}
source /usr/local/shared/bin/xrootd_env.sh
$OMEGASCAN
unset LD_PRELOAD
unset XROOTD_VMP
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH_SAV}

echo "convert thumbnails...."

tempConvert="tmpConvert$QSUB_FILEID.sh"

for i in `ls -1 $outputDirectory/$eventTime` ; do
  pngEnd=`echo $i | sed "s/.*.png$/.png/"`
  fileName=$outputDirectory/$eventTime/$i
  if [ $pngEnd = ".png" ]; then
     echo $fileName | awk '{tmp = substr($1,1,length($1)-4);print "convert -resize 300x " $1 "  -strip -depth 8 -colors 256 " tmp"_thumbnail.png" }' >> $tempConvert
  fi
done
chmod u+x $tempConvert;
./$tempConvert; rm $tempConvert

echo "transfer files to web directory"
#scp -i ~/.ssh/id_rsa -r $outputDirectory/$eventTime $webDirectory/$webSubDirectory
mkdir -p $webDirectory/$webSubDirectory
cp -r $outputDirectory/$eventTime $webDirectory/$webSubDirectory

exit 0
