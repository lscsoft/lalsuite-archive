#!/usr/bin/python

import os, glob, optparse, shutil, warnings
from collections import namedtuple
from subprocess import Popen, PIPE, STDOUT
from pylal.xlal.datatypes.ligotimegps import LIGOTimeGPS
from pylal.xlal.date import XLALGPSToUTC
import seismon_psd, seismon_html

__author__ = "Michael Coughlin <michael.coughlin@ligo.org>"
__date__ = "2012/8/26"
__version__ = "0.1"

# =============================================================================
#
#                               DEFINITIONS
#
# =============================================================================

def parse_commandline():
    """
    Parse the options given on the command-line.
    """
    parser = optparse.OptionParser(usage=__doc__,version=__version__)

    parser.add_option("-p", "--paramsFile", help="Seismon params file.", default ="/home/mcoughlin/matapps/packages/seismic/trunk/seismon/python/input/seismon_params_H1PEM.txt")

    parser.add_option("-g", "--gps", help="GPS Time.", default="1030044950")

    parser.add_option("-v", "--verbose", action="store_true", default=False,
                      help="Run verbosely. (Default: False)")

    opts, args = parser.parse_args()

    # show parameters
    if opts.verbose:
        print >> sys.stderr, ""
        print >> sys.stderr, "running network_eqmon..."
        print >> sys.stderr, "version: %s"%__version__
        print >> sys.stderr, ""
        print >> sys.stderr, "***************** PARAMETERS ********************"
        for o in opts.__dict__.items():
          print >> sys.stderr, o[0]+":"
          print >> sys.stderr, o[1]
        print >> sys.stderr, ""

    return opts

def readParamsFromFile(file):
        
    params = {}
    if os.path.isfile(file):
        with open(file,'r') as f:
            for line in f:
                line_without_return = line.split("\n")
                line_split = line_without_return[0].split(" ")
                params[line_split[0]] = line_split[1]
    return params

def channel_struct(channelList):
    # Create channel structure
    structproxy_channel = namedtuple( "structproxy_channel", "station station_underscore samplef calibration latitude longitude" )

    channel = []

    with open(channelList,'r') as f:

       for line in f:

           line_without_return = line.split("\n")
           line_split = line_without_return[0].split(" ")

           station = line_split[0]
           station_underscore = station.replace(":","_")

           samplef = float(line_split[1])
           calibration = float(line_split[2])

           if station[0] == "H":
               latitude = 46.6475
               longitude = -119.5986;
           elif station[0] == "L":
               latitude = 30.4986
               longitude = -90.7483

           channel.append( structproxy_channel(station,station_underscore,samplef,calibration,latitude,longitude))
    return channel

def frame_struct(params):

    p = Popen("ligo_data_find -o %s -t %s -s %.0f -e %.0f -u file"%(params["ifo"][0],params["frameType"],params["gps"],params["gps"]),shell=True,stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    output = p.stdout.read()
    frameList = output.split("\n")
    frames = []
    frameGPS = []
    frameDur = []
    for frame in frameList:
        thisFrame = frame.replace("file://localhost","")
        if thisFrame is not "":
            frames.append(thisFrame)
            thisFrameSplit = thisFrame.split("-")
            frameGPS = float(thisFrameSplit[-2])
            frameDur = float(thisFrameSplit[-1].replace(".gwf",""))
    params["frame"] = frames
    params["frameGPS"] = frameGPS
    params["frameDur"] = frameDur
    return params

# =============================================================================
#
#                                    MAIN
#
# =============================================================================

warnings.filterwarnings("ignore")

# Parse command line
opts = parse_commandline()
params = readParamsFromFile(opts.paramsFile)
params["gps"] = float(opts.gps)
params["date"] = XLALGPSToUTC(LIGOTimeGPS(params["gps"]))

channelList = params["matappsPath"] + "/seismon/python/input/seismon-" + params["ifo"] + "-" + params["frameType"] + "-channel_list.txt"
channels = channel_struct(channelList)

params = frame_struct(params)

# Output path for run
params["path"] = params["dirPath"] + "/" + params["ifo"] + "/" + params["runName"] + "-" + str(params["frameGPS"])
if not os.path.isdir(params["path"]):
    os.makedirs(params["path"])

params["fmin"] = 1/float(64)
params["fmax"] = 10
params["doPlots"] = 1

for channel in channels:
    seismon_psd.mat(params,channel);
    seismon_psd.analysis(params,channel);

htmlPage = seismon_html.summary_page(params,channels)
if htmlPage is not None:
    f = open(os.path.join(params["path"],"summary.html"),"w")
    f.write(htmlPage)
    f.close()

# Public HTML output path
params["outputPath"] = os.path.join(params["publicPath"],params["ifo"]);
if not os.path.isdir(params["outputPath"]):
    os.makedirs(params["outputPath"])

if os.path.isdir(os.path.join(params["outputPath"],params["runName"])):
    p = Popen("rm -r %s/%s"%(params["outputPath"],params["runName"]),shell=True,stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)

# Softlink HTML to output path
p = Popen("cp -r %s %s/%s"%(params["path"],params["outputPath"],params["runName"]),shell=True,stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
#f = open(os.path.join(params["path"],"softlink.sh"),"w")
#f.write("#! /usr/bin/env bash\n")
#f.close()
#p = Popen("source %s"%(os.path.join(params["path"],"softlink.sh")),shell=True,stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)

