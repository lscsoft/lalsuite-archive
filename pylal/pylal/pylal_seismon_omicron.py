#!/usr/bin/python

import os, glob, optparse, shutil, warnings, re, matplotlib
import matplotlib.pyplot as plt
import numpy as np
from collections import namedtuple
from subprocess import Popen, PIPE, STDOUT
from pylal.xlal.datatypes.ligotimegps import LIGOTimeGPS
from pylal.xlal.date import XLALGPSToUTC

from glue.ligolw import table as ligolw_table
from glue.ligolw import utils as ligolw_utils
from glue.ligolw import lsctables

__author__ = "Michael Coughlin <michael.coughlin@ligo.org>"
__date__ = "2012/8/26"
__version__ = "0.1"

# =============================================================================
#
#               DEFINITIONS
#
# =============================================================================

# set up useful regular expressions
_re_xml = re.compile('(xml|xml.gz)\Z')

def read_triggers(fp, format_, columns=None):
    """@returns a LIGOLw xml table of events loaded from the given filepath

    @param fp
        path to file on disk
    @param format
        identifier of trigger type, or table name
    @param columns
        list of column name strings to load
    """
    format_ = format_.lower()
    isxml = _re_xml.search(fp)

    # verify columns
    if columns is not None:
        columns = map(str.lower, columns)
        if "time" in columns:
            columns.pop(columns.index("time"))
            if re.match("sim_inspiral", format_, re.I):
                columns.extend(["geocent_end_time", "geocent_end_time_ns"])
            elif re.match("sim_burst", format_, re.I):
                columns.extend(["time_geocent_gps", "time_geocent_gps_ns"])
            elif re.match("sim_ringdown", format_, re.I):
                columns.extend(["geocent_start_time", "geocent_start_time_ns"])
            elif re.search("inspiral", format_, re.I):
                columns.extend(["end_time", "end_time_ns"])
            elif re.search("burst", format_, re.I):
                columns.extend(["peak_time", "peak_time_ns"])
            elif re.search("(ringdown|stochastic)", format_, re.I):
                columns.extend(["start_time", "start_time_ns"])

    # read XML file
    if isxml:
        # get table class
        Table = lsctables.TableByName[format_]
        # verify columns are actually in the table. If not, load all of them
        for c in columns:
            if c.lower() not in Table.validcolumns.keys():
                columns = None
                break
        # set columns
        if columns is not None:
            Table.loadcolumns = columns
        # load table
        xmldoc = ligolw_utils.load_filename(fp, gz=fp.endswith("gz"))
        return ligolw_table.get_table(xmldoc, format_)
    # read ASCII format_ file
    else:
        if format_ == "omega":
            return omegautils.fromfile(fp, columns=columns)
        elif format_ == "kw":
            return kwutils.fromfile(fp, columns=columns)
        else:
            raise ValueError("No read function defined for ASCII format "
             "\"%s\"" % format_)

def plot_triggers(params,channel):

    columns = ["peak_time","peak_time_ns","start_time","start_time_ns","stop_time","stop_time_ns","duration","central_freq","flow","fhigh","bandwidth","amplitude","snr"]

    trigger_threshold = 5

    omicronDirectory = os.path.join(params["path"],"omicron")
    omicronPath = os.path.join(omicronDirectory,channel.station)
    omicronXMLs = glob.glob(os.path.join(omicronPath,"*.xml"))

    triggers = []

    for omicronXML in omicronXMLs:
        triggers_xml = read_triggers(omicronXML, "sngl_burst", columns)
        for trigger in triggers_xml:
            if trigger.snr < trigger_threshold:
                continue
            t = trigger.peak_time + trigger.peak_time_ns * 10**(-9)
            triggers.append([t,trigger.central_freq,trigger.snr])
      
    textLocation = params["path"] + "/" + channel.station_underscore
    if not os.path.isdir(textLocation):
        os.makedirs(textLocation)

    triggers = np.array(triggers)
    if len(triggers) > 0:
        triggers_t = triggers[:,0]
        triggers_f = triggers[:,1]
        triggers_snr = triggers[:,2]

    f = open(os.path.join(textLocation,"triggers.txt"),"w")
    for i in xrange(len(triggers)):
        f.write("%.1f %e %e\n"%(triggers_t[i],triggers_f[i],triggers_snr[i]))
    f.close()

    earthquakesDirectory = os.path.join(params["path"],"earthquakes")
    earthquakesFile = os.path.join(earthquakesDirectory,"earthquakes.txt")
    try:
        earthquakes = np.loadtxt(earthquakesFile)
    except:
        earthquakes = []

    if params["doPlots"]:

        plotLocation = params["path"] + "/" + channel.station_underscore
        if not os.path.isdir(plotLocation):
            os.makedirs(plotLocation)

        startTime = params["gpsStart"]
        timeDur = params["gpsEnd"] - startTime

        ax = plt.subplot(111)
        if len(triggers) > 0:
            triggers_t = triggers_t - startTime
            plt.scatter(triggers_t,triggers_f, c=triggers_snr,vmin=min(triggers_snr),vmax=max(triggers_snr))
        cbar = plt.colorbar(orientation='horizontal')  
        #cbar = plt.colorbar(orientation='vertical') 
        cbar.set_label('SNR')
        ax.set_yscale('log')

        if len(earthquakes) > 0:
            if len(earthquakes.shape) == 1:
                shape_x = 1
            else:
                [shape_x,shape_y] = earthquakes.shape
            for i in xrange(shape_x):
                if earthquakes[i,1] < 4.0:
                    continue

                Ptime = earthquakes[i,2] - startTime
                Stime = earthquakes[i,3] - startTime
                Rtime = earthquakes[i,4] - startTime

                plt.text(Ptime, 15, 'P', fontsize=18, ha='center', va='top')
                plt.text(Stime, 15, 'S', fontsize=18, ha='center', va='top')
                plt.text(Rtime, 15, 'R', fontsize=18, ha='center', va='top')

                plt.axvline(x=Ptime,color='r',linewidth=2,zorder = 0,clip_on=False)
                plt.axvline(x=Stime,color='b',linewidth=2,zorder = 0,clip_on=False)
                plt.axvline(x=Rtime,color='g',linewidth=2,zorder = 0,clip_on=False)

        plt.ylim([0.1,10])
        #plt.xlim([params["gpsStart"], params["gpsEnd"]])
        plt.xlim([0,timeDur])
        plt.ylabel("Frequency [Hz]")
        plt.xlabel("GPS [s] [%d]"%startTime)
        plt.grid
        plt.show()
        plt.savefig(os.path.join(plotLocation,"omicron.png"),dpi=200)
        plt.savefig(os.path.join(plotLocation,"omicron.eps"),dpi=200)
        plt.close('all')

def generate_triggers(params,channels):

    omicronDirectory = os.path.join(params["path"],"omicron")
    if not os.path.isdir(omicronDirectory):
        os.makedirs(omicronDirectory)

    gpsStart = 1e20
    gpsEnd = -1e20
    f = open(os.path.join(omicronDirectory,"frames.ffl"),"w")
    for frame in params["frame"]:
        f.write("%s %d %d 0 0\n"%(frame.path, frame.segment[0], frame.segment[1]-frame.segment[0]))
        gpsStart = min(gpsStart,frame.segment[0])
        gpsEnd = max(gpsEnd,frame.segment[1])
    f.close()

    paramsFile = omicron_params(params,channels)
    f = open(os.path.join(omicronDirectory,"params.txt"),"w")
    f.write("%s"%(paramsFile))
    f.close()

    f = open(os.path.join(omicronDirectory,"segments.txt"),"w")
    f.write("%d %d\n"%(gpsStart,gpsEnd))
    f.close()

    omicron = "/home/detchar/opt/virgosoft/Omicron/v0r3/Linux-x86_64/omicron.exe"
    environmentSetup = "CMTPATH=/home/detchar/opt/virgosoft; export CMTPATH; source /home/detchar/opt/virgosoft/Omicron/v0r3/cmt/setup.sh"
    omicronCommand = "%s; %s %s %s"%(environmentSetup, omicron, os.path.join(omicronDirectory,"segments.txt"),os.path.join(omicronDirectory,"params.txt"))

    p = Popen(omicronCommand,shell=True,stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    output = p.stdout.read()

def omicron_params(params,channels):

    omicronDirectory = os.path.join(params["path"],"omicron")

    channelList = ""
    samplerateList = ""
    for channel in channels:
        channelList = "%s %s"%(channelList,channel.station)
        samplerateList = "%s %d"%(samplerateList,channel.samplef)

    paramsFile = """

    DATA    FFL     %s/frames.ffl
    
    //** list of channels you want to process
    DATA    CHANNELS %s
    
    //** native sampling frequency (Hz) of working channels (as many as
    //the number of input channels)
    DATA    NATIVEFREQUENCY %s
    
    //** working sampling (one value for all channels)
    DATA    SAMPLEFREQUENCY 32
    
    //*************************************************************************************
    //************************        SEARCH PARAMETERS
    //*****************************
    //*************************************************************************************
    
    //** chunk duration in seconds (must be an integer)
    PARAMETER       CHUNKDURATION   864
    
    //** segment duration in seconds (must be an integer)
    PARAMETER       BLOCKDURATION   512
    
    //** overlap duration between segments in seconds (must be an integer)
    PARAMETER       OVERLAPDURATION  160
    
    //** search frequency range
    PARAMETER       FREQUENCYRANGE  0.1      10
    
    //** search Q range
    PARAMETER       QRANGE          3.3166  141
    
    //** maximal mismatch between 2 consecutive tiles (0<MM<1)
    PARAMETER       MISMATCHMAX     0.2
    
    //*************************************************************************************
    //************************            TRIGGERS
    //*****************************
    //*************************************************************************************
    
    //** tile SNR threshold
    TRIGGER         SNRTHRESHOLD    5
    
    //** maximum number of triggers per file
    TRIGGER         NMAX            500000
    
    //*************************************************************************************
    //************************             OUTPUT
    //*****************************
    //*************************************************************************************
    
    //** full path to output directory
    OUTPUT  DIRECTORY       %s/
    
    OUTPUT  FORMAT   xml
    
    """%(omicronDirectory,channelList,samplerateList,omicronDirectory)

    return paramsFile

