#!/usr/bin/python

import os, glob, optparse, shutil, warnings, matplotlib, pickle, math, copy, pickle
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal, scipy.stats, scipy.fftpack
from collections import namedtuple
from lxml import etree
from pylal.xlal.datatypes.ligotimegps import LIGOTimeGPS
from pylal.xlal.date import XLALGPSToUTC
from pylal import Fr
import glue.datafind, glue.segments, glue.segmentsUtils, glue.lal
import pylal.pylal_seismon_NLNM, pylal.pylal_seismon_html
import pylal.pylal_seismon_eqmon

import gwpy.time, gwpy.timeseries, gwpy.spectrum, gwpy.plotter

__author__ = "Michael Coughlin <michael.coughlin@ligo.org>"
__date__ = "2012/8/26"
__version__ = "0.1"

# =============================================================================
#
#                               DEFINITIONS
#
# =============================================================================

def mkdir(path):

    pathSplit = path.split("/")
    pathAppend = "/"
    for piece in pathSplit:
        if piece == "":
            continue
        pathAppend = os.path.join(pathAppend,piece)
        if not os.path.isdir(pathAppend):
            os.mkdir(pathAppend)   

def read_frames(start_time,end_time,channel,cache):

    time = []
    data = []

    #== loop over frames in cache
    for frame in cache:

        if end_time < frame.segment[0]:
            continue
        if start_time > frame.segment[1]:
            continue

        frame_data,data_start,_,dt,_,_ = Fr.frgetvect1d(frame.path,channel.station)
        frame_length = float(dt)*len(frame_data)
        frame_time = data_start+dt*np.arange(len(frame_data))

        for i in range(len(frame_data)):
            if frame_time[i] <= start_time:  continue
            if frame_time[i] >= end_time:  continue
            time.append(frame_time[i])
            data.append(frame_data[i])
    data = [e/channel.calibration for e in data]

    indexes = np.where(np.isnan(data))[0]
    meanSamples = np.mean(np.ma.masked_array(data,np.isnan(data)))
    for index in indexes:
        data[index] = meanSamples

    return time,data

def read_nds(start_time,end_time,channel,conn):

    try:
        buffers = conn.fetch(start_time, end_time,[channel.station])
    except:
        time = []
        data = []
        return time,data

    data = buffers[0].data
    data_start = buffers[0].gps_seconds + buffers[0].gps_nanoseconds
    dt = 1.0 / buffers[0].channel.sample_rate 
    time = data_start+dt*np.arange(len(data))
    data = [e/channel.calibration for e in data]

    indexes = np.where(np.isnan(data))[0]
    meanSamples = np.mean(np.ma.masked_array(data,np.isnan(data)))
    for index in indexes:
        data[index] = meanSamples

    return time,data

def normalize_timeseries(data):  

    dataSort = np.sort(data)
    index10 = np.floor(len(data) * 0.1)
    index90 = np.floor(len(data) * 0.9)
    dataMin = dataSort[index10]
    dataMax = dataSort[index90] 

    dataNorm = (1/(dataMax - dataMin)) * (data - dataMin)

    indexes = (np.absolute(dataNorm) >= 2).nonzero()[0]
    dataNorm[indexes] = 0

    dataNorm = dataNorm / 2
 
    return dataNorm

def envelope(data):

    hilb = scipy.fftpack.hilbert(data)
    data = (data ** 2 + hilb ** 2) ** 0.5
    return data

def read_eqmons(file):

    attributeDics = []

    if not os.path.isfile(file):
        print "Missing eqmon file: %s"%file
        return attributeDics

    tree = etree.parse(file)
    baseroot = tree.getroot()       # get the document root
    for root in baseroot.iterchildren():
        attributeDic = {}
        for element in root.iterchildren(): # now iter through it and print the text
            if element.tag == "traveltimes":
                attributeDic[element.tag] = {}
                for subelement in element.iterchildren():
                    attributeDic[element.tag][subelement.tag] = {}
                    for subsubelement in subelement.iterchildren():
                        textlist = subsubelement.text.replace("\n","").split(" ")
                        floatlist = [float(x) for x in textlist]
                        attributeDic[element.tag][subelement.tag][subsubelement.tag] = floatlist
            else:
                try:
                    attributeDic[element.tag] = float(element.text)
                except:
                    attributeDic[element.tag] = element.text

        magThreshold = 0
        if not "Magnitude" in attributeDic or attributeDic["Magnitude"] < magThreshold:
            return attributeDic

        attributeDic["doPlots"] = 0
        for ifoName, traveltimes in attributeDic["traveltimes"].items():
            arrivalMin = min([max(traveltimes["Rtimes"]),max(traveltimes["Stimes"]),max(traveltimes["Ptimes"])])
            arrivalMax = max([max(traveltimes["Rtimes"]),max(traveltimes["Stimes"]),max(traveltimes["Ptimes"])])
            attributeDic["traveltimes"][ifoName]["arrivalMin"] = arrivalMin
            attributeDic["traveltimes"][ifoName]["arrivalMax"] = arrivalMax
            #if params["gps"] <= attributeDic["traveltimes"][ifoName]["arrivalMax"]:
            #    attributeDic["doPlots"] = 1

        attributeDics.append(attributeDic)
    return attributeDics

def spectral_histogram(specgram,bins=None,lowBin=None,highBin=None,nbins=None):

    # Define bins for the spectral variation histogram
    if lowBin == None:
        lowBin = np.log10(np.min(specgram)/2)
    if highBin == None:
        highBin = np.log10(np.max(specgram)*2)
    if nbins == None:
        nbins = 500   
    if bins == None:
        bins = np.logspace(lowBin,highBin,num=nbins)

    # Ensure we work with numpy array data
    data = np.array(specgram)

    spectral_variation_norm = []
    rows, columns = data.shape

    # Loop over frequencies
    for i in xrange(columns):
        # calculate histogram for this frequency bin
        this_spectral_variation, bin_edges = np.histogram(data[:,i],bins)
        this_spectral_variation = np.array(this_spectral_variation)
        # Calculate weights for bins (to normalize)
        weight = (100/float(sum(this_spectral_variation))) + np.zeros(this_spectral_variation.shape)
        # stack output array
        if spectral_variation_norm == []:
            spectral_variation_norm = this_spectral_variation * weight
        else:
            spectral_variation_norm = np.vstack([spectral_variation_norm,this_spectral_variation * weight])
    spectral_variation_norm = np.transpose(spectral_variation_norm)

    return bins,spectral_variation_norm

def spectral_percentiles(specvar,bins,percentile):

    # Ensure we work with numpy array data
    data = np.array(specvar)

    percentiles = []
    rows, columns = specvar.shape

    # Loop over frequencies
    for i in xrange(columns):
        # Calculate cumulative sum for array
        cumsumvals = np.cumsum(data[:,i])

        # Find value nearest requested percentile
        abs_cumsumvals_minus_percentile = abs(cumsumvals - percentile)
        minindex = abs_cumsumvals_minus_percentile.argmin()
        val = bins[minindex]

        percentiles.append(val)

    return percentiles

def html_bgcolor(snr,data):

    data = np.append(data,snr)

    # Number of colors in array
    N = 256

    colormap = []
    for i in xrange(N):
        r,g,b,a = matplotlib.pyplot.cm.jet(i)
        r = int(round((r * 255),0))
        g = int(round((g * 255),0))
        b = int(round((b* 255),0))
        colormap.append((r,g,b))

    data = np.sort(data)
    itemIndex = np.where(data==snr)

    # Determine significance of snr (between 0 and 1)
    snrSig = itemIndex[0][0] / float(len(data)+1)

    # Determine color index of this significance
    index = int(np.floor(N * snrSig))

    # Return colors of this index
    thisColor = colormap[index]
    # Return rgb string containing these colors
    bgcolor = "rgb(%d,%d,%d)"%(thisColor[0],thisColor[1],thisColor[2])

    return snrSig, bgcolor

def segment_struct(params):

    if params["doSegmentsDatabase"]:
        segmentlist, segmentlistValid = pylal.dq.dqSegmentUtils.grab_segments(
                                               params["gpsStart"],params["gpsEnd"],
                                               params["segmentFlag"],params["segmentDatabase"],
                                               segment_summary=True)
    elif params["doSegmentsTextFile"]:
        segmentlist = glue.segments.segmentlist()
        segs = np.loadtxt(params["segmentsTextFile"])
        for seg in segs:
            segmentlist.append(glue.segments.segment(seg[0],seg[1]))
    else:
        segmentlist = [glue.segments.segment(params["gpsStart"],params["gpsEnd"])]

    params["segments"] = segmentlist

    return params

def frame_struct(params):

    gpsStart = params["gpsStart"]-1000
    gpsEnd = params["gpsEnd"]

    if params["ifo"] == "XG":
        frameDir = "/archive/frames/MBH/"
        frameList = [os.path.join(root, name)
            for root, dirs, files in os.walk(frameDir)
            for name in files]

        datacache = []
        for frame in frameList:
            thisFrame = frame.replace("file://localhost","")
            if thisFrame == "":
                continue

            thisFrameSplit = thisFrame.split(".")
            if thisFrameSplit[-1] == "log":
                continue

            thisFrameSplit = thisFrame.split("-")
            gps = float(thisFrameSplit[-2])
            dur = float(thisFrameSplit[-1].replace(".gwf",""))

            if gps+dur < gpsStart:
                continue
            if gps > gpsEnd:
                continue

            cacheFile = glue.lal.CacheEntry("%s %s %d %d %s"%("XG","Homestake",gps,dur,frame))
            datacache.append(cacheFile)

    else:
        if params["frameType"] == "nds":
            conn = nds2.connection(params["ndsServer"])
            #y = conn.find_channels('*',nds2.channel.CHANNEL_TYPE_RAW,\
            #    nds2.channel.DATA_TYPE_FLOAT32, 128, 16384)

            params["ndsConnection"] = conn

        else:
            connection = glue.datafind.GWDataFindHTTPConnection()
            datacache = connection.find_frame_urls(params["ifo"][0], params["frameType"],
                                                   gpsStart, gpsEnd,
                                                   urltype="file",
                                                   on_gaps="warn")
            connection.close()

            params["frame"] = datacache

    return params

def channel_struct(params,channelList):
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
           elif station[0] == "G":
               latitude = 52.246944
               longitude = 9.808333
           elif station[0] == "V":
               latitude = 43.631389
               longitude = 10.505
           elif station[0] == "C":
               latitude = 34.1391
               longitude = -118.1238
           elif station[0] == "M":
               latitude = 44.3465
               longitude = -103.7574

           if not params["channel"] == None:
               if not station in params["channel"]:
                   continue

           channel.append( structproxy_channel(station,station_underscore,samplef,calibration,latitude,longitude))
    return channel

def readParamsFromFile(file):

    params = {}
    if os.path.isfile(file):
        with open(file,'r') as f:
            for line in f:
                line_without_return = line.split("\n")
                line_split = line_without_return[0].split(" ")
                params[line_split[0]] = line_split[1]
    return params
