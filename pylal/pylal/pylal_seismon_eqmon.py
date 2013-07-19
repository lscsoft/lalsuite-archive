#!/usr/bin/python

import os, sys, time, glob, math, matplotlib, random, string
import numpy as np
from datetime import datetime
from operator import itemgetter
import xml.dom.minidom
import glue.GWDataFindClient, glue.segments, glue.segmentsUtils
from pylal.xlal.datatypes.ligotimegps import LIGOTimeGPS
from pylal.xlal.date import XLALUTCToGPS, XLALGPSToUTC
from lxml import etree
import scipy.spatial

import pylal.pylal_seismon_eqmon_plot

def run_earthquakes(params):

    timeseriesDirectory = os.path.join(params["path"],"timeseries")
    if not os.path.isdir(timeseriesDirectory):
        os.makedirs(timeseriesDirectory)

    earthquakesDirectory = os.path.join(params["path"],"earthquakes")
    if not os.path.isdir(earthquakesDirectory):
        os.makedirs(earthquakesDirectory)

    if params["ifo"] == "H1":
        ifo = "LHO"
    elif params["ifo"] == "L1":
        ifo = "LLO"
    elif params["ifo"] == "G1":
        ifo = "GEO"
    elif params["ifo"] == "V1":
        ifo = "VIRGO"
    elif params["ifo"] == "C1":
        ifo = "FortyMeter"

    if params["doEarthquakesAnalysis"]:
       params["earthquakesMinMag"] = 5
    else:
       params["earthquakesMinMag"] = 0

    attributeDics = retrieve_earthquakes(params)
    attributeDics = sorted(attributeDics, key=itemgetter("Magnitude"), reverse=True)

    if params["doEarthquakesMonitor"]:
        earthquakesFile = os.path.join(earthquakesDirectory,"%d-%d.txt"%(params["gpsStart"],params["gpsEnd"]))
        timeseriesFile = os.path.join(timeseriesDirectory,"%d-%d.txt"%(params["gpsStart"],params["gpsEnd"]))
    else:
        earthquakesFile = os.path.join(earthquakesDirectory,"earthquakes.txt")
        timeseriesFile = os.path.join(timeseriesDirectory,"amp.txt")
 
    f = open(earthquakesFile,"w+")

    amp = 0
    segmentlist = glue.segments.segmentlist()
 
    for attributeDic in attributeDics:

        #attributeDic = calculate_traveltimes(attributeDic)

        traveltimes = attributeDic["traveltimes"][ifo]

        gpsStart = max(traveltimes["Rtimes"]) - 200
        gpsEnd = max(traveltimes["Rtimes"]) + 200

        check_intersect = (gpsEnd >= params["gpsStart"]) and (params["gpsEnd"] >= gpsStart)

        if check_intersect:
            amp += traveltimes["Rfamp"][0]

            f.write("%.1f %.1f %.1f %.1f %.1f %.5e %d %d %.1f %.1f\n"%(attributeDic["GPS"],attributeDic["Magnitude"],max(traveltimes["Ptimes"]),max(traveltimes["Stimes"]),max(traveltimes["Rtimes"]),traveltimes["Rfamp"][0],gpsStart,gpsEnd,attributeDic["Latitude"],attributeDic["Longitude"]))

        segmentlist.append(glue.segments.segment(gpsStart,gpsEnd))

    f.close()

    f = open(timeseriesFile,"w+")
    f.write("%e\n"%(amp))
    f.close()

    if not params["doPlots"]:
        return segmentlist

    plotsDirectory = os.path.join(params["path"],"plots")
    if not os.path.isdir(plotsDirectory):
        os.makedirs(plotsDirectory)

    try:
        earthquakes = np.loadtxt(earthquakesFile)
    except:
        return segmentlist

    if earthquakes.ndim == 1:
        gpsStart = earthquakes[6]
        gpsEnd = earthquakes[7]
        latitude = earthquakes[8]
        longitude = earthquakes[9]
        magnitude = earthquakes[1]
    else:
        gpsStart = earthquakes[:,6]
        gpsEnd = earthquakes[:,7]
        latitude = earthquakes[:,8]
        longitude = earthquakes[:,9]
        magnitude = earthquakes[:,1]
    #for i in xrange(len(earthquakes)):
    #    print gpsStart[i], gpsEnd[i]

    files = glob.glob(os.path.join(timeseriesDirectory,"*.txt"))
    files = sorted(files)

    ttStart = []
    ttEnd = []
    amp = []

    for file in files:

        fileSplit = file.split("/")

        if fileSplit[-1] == "amp.txt":
            continue

        txtFile = fileSplit[-1].replace(".txt","")
        txtFileSplit = txtFile.split("-")
        thisTTStart = int(txtFileSplit[0])
        thisTTEnd = int(txtFileSplit[1])

        if (thisTTStart < params["gpsStart"]) or (thisTTEnd > params["gpsEnd"]):
            continue

        ttStart.append(thisTTStart)
        ttEnd.append(thisTTEnd)

        data_out = np.loadtxt(file)
        thisAmp = data_out

        amp.append(thisAmp)

    ttStart = np.array(ttStart)
    ttEnd = np.array(ttEnd)
    amp = np.array(amp)

    data = {}
    data["prediction"] = {}
    data["prediction"]["tt"] = np.array(ttStart)
    data["prediction"]["data"] = np.array(amp)

    data["channels"] = {}

    # Break up entire frequency band into 6 segments
    ff_ave = [1/float(128), 1/float(64),  0.1, 1, 3, 5, 10]

    for channel in params["channels"]:

        psdLocation = params["dirPath"] + "/Text_Files/PSD/" + channel.station_underscore
        if not os.path.isdir(psdLocation):
            os.makedirs(psdLocation)
        psdLocation = os.path.join(psdLocation,str(params["fftDuration"]))
        if not os.path.isdir(psdLocation):
            os.makedirs(psdLocation)

        files = glob.glob(os.path.join(psdLocation,"*.txt"))
        files = sorted(files)

        ttStart = []
        ttEnd = []
        amp = []
    
        for file in files:
    
            fileSplit = file.split("/")
            txtFile = fileSplit[-1].replace(".txt","")
            txtFileSplit = txtFile.split("-")
            thisTTStart = int(txtFileSplit[0])
            thisTTEnd = int(txtFileSplit[1])
   
            if (thisTTStart < params["gpsStart"]) or (thisTTEnd > params["gpsEnd"]):
                continue
 
            ttStart.append(thisTTStart)
            ttEnd.append(thisTTEnd)

            data_out = np.loadtxt(file)
            thisSpectra_out = data_out[:,1]
            thisFreq_out = data_out[:,0]

            freqAmps = []

            for i in xrange(len(ff_ave)-1):
                newSpectraNow = []
                for j in xrange(len(thisFreq_out)):
                    if ff_ave[i] <= thisFreq_out[j] and thisFreq_out[j] <= ff_ave[i+1]:
                        newSpectraNow.append(thisSpectra_out[j])
                freqAmps.append(np.mean(newSpectraNow)) 
               
            thisAmp = freqAmps[1]
            amp.append(thisAmp)
    
        ttStart = np.array(ttStart)
        ttEnd = np.array(ttEnd)
        amp = np.array(amp)

        data["channels"][channel.station_underscore] = {}
        data["channels"][channel.station_underscore]["tt"] = np.array(ttStart)
        data["channels"][channel.station_underscore]["data"] = np.array(amp)

    plotName = os.path.join(plotsDirectory,"%d-%d.png"%(params["gpsStart"],params["gpsEnd"]))
    pylal.pylal_seismon_eqmon_plot.prediction(data,plotName)
    plotName = os.path.join(plotsDirectory,"%d-%d-residual.png"%(params["gpsStart"],params["gpsEnd"]))
    pylal.pylal_seismon_eqmon_plot.residual(data,plotName)
    print plotName
    print penis

    for attributeDic in attributeDics:

        if not ifo in attributeDic["traveltimes"]:
            continue

        traveltimes = attributeDic["traveltimes"][ifo]

        gpsStart = max(traveltimes["Rtimes"]) - 200
        gpsEnd = max(traveltimes["Rtimes"]) + 200

        for channel in params["channels"]:

            envelopeLocation = params["dirPath"] + "/Text_Files/Envelope/" + channel.station_underscore
            if not os.path.isdir(envelopeLocation):
                os.makedirs(envelopeLocation)
            envelopeLocation = os.path.join(envelopeLocation,str(params["fftDuration"]))
            if not os.path.isdir(envelopeLocation):
                os.makedirs(envelopeLocation)

            envelopeFiles = glob.glob(os.path.join(envelopeLocation,"*"))
            envelopeFiles = sorted(envelopeFiles)

            time_envelope,data_envelope = get_envelope(gpsStart,gpsEnd,envelopeFiles)
            time_envelope = np.array(time_envelope)
            data_envelope = np.array(data_envelope)

            if len(time_envelope) > 0:

                data_envelope_argmax = data_envelope.argmax()
                time_envelope_argmax = time_envelope[data_envelope_argmax]

                timeEstimate = time_envelope_argmax
                attributeDic["traveltimes"][ifo]["Restimate"] = timeEstimate
                plotName = os.path.join(earthquakesDirectory,"%s-%d-%d.png"%(channel.station_underscore,\
                    gpsStart,gpsEnd))
                pylal.pylal_seismon_eqmon_plot.plot_envelope(params,time_envelope,data_envelope,\
                    attributeDic["traveltimes"][ifo],plotName)

    if params["doPlots"]:

        if params["doEarthquakesAnalysis"]:
            plotName = os.path.join(earthquakesDirectory,"worldmap_magnitudes.png")
            pylal.pylal_seismon_eqmon_plot.worldmap_plot(params,attributeDics,"Magnitude",plotName)

            plotName = os.path.join(earthquakesDirectory,"worldmap_traveltimes.png")
            pylal.pylal_seismon_eqmon_plot.worldmap_plot(params,attributeDics,"Traveltimes",plotName)

            plotName = os.path.join(earthquakesDirectory,"worldmap_restimates.png")
            pylal.pylal_seismon_eqmon_plot.worldmap_plot(params,attributeDics,"Restimates",plotName)

            plotName = os.path.join(earthquakesDirectory,"restimates.png")
            pylal.pylal_seismon_eqmon_plot.restimates(params,attributeDics,plotName)

        plotName = os.path.join(earthquakesDirectory,"magnitudes.png")
        pylal.pylal_seismon_eqmon_plot.magnitudes(params,attributeDics,plotName)
        plotName = os.path.join(earthquakesDirectory,"magnitudes_latencies.png")
        pylal.pylal_seismon_eqmon_plot.magnitudes_latencies(params,attributeDics,plotName)
        plotName = os.path.join(earthquakesDirectory,"latencies_sent.png")
        pylal.pylal_seismon_eqmon_plot.latencies_sent(params,attributeDics,plotName)
        plotName = os.path.join(earthquakesDirectory,"latencies_written.png")
        pylal.pylal_seismon_eqmon_plot.latencies_written(params,attributeDics,plotName)
        plotName = os.path.join(earthquakesDirectory,"traveltimes%s.png"%params["ifo"])
        pylal.pylal_seismon_eqmon_plot.traveltimes(params,attributeDics,ifo,params["gpsEnd"],plotName)
        plotName = os.path.join(earthquakesDirectory,"worldmap.png")
        pylal.pylal_seismon_eqmon_plot.worldmap_wavefronts(params,attributeDics,params["gpsEnd"],plotName)


    return segmentlist

def get_envelope(start_time,end_time,files):

    time = []
    data = []

    #== loop over frames in cache
    for file in files:

        fileSplit = file.split("/")
        txtFile = fileSplit[-1].replace(".txt","")
        txtFileSplit = txtFile.split("-")
        thisTTStart = int(txtFileSplit[0])
        thisTTEnd = int(txtFileSplit[1])

        if end_time < thisTTStart:
            continue
        if start_time > thisTTEnd:
            continue

        try:
            data_out = np.loadtxt(file)
        except:
            continue
        file_time = data_out[:,0]
        file_data = data_out[:,1]
        data_out[0,1] = data_out[3,1]
        data_out[1,1] = data_out[3,1]
        data_out[2,1] = data_out[3,1]

        for i in range(len(file_data)):
            if file_time[i] <= start_time:  continue
            if file_time[i] >= end_time:  continue
            time.append(file_time[i])
            data.append(file_data[i])

    return time,data

def parse_xml(element):

    subdic = {}

    numChildren = 0
    for subelement in element.iterchildren():
        tag = str(subelement.tag)
        tag = tag.replace("{http://www.usgs.gov/ansseqmsg}","")
        tag = tag.replace("{http://quakeml.org/xmlns/quakeml/1.2}","")
        subdic[tag] = parse_xml(subelement)
        numChildren += 1

    if numChildren == 0:
        value = str(element.text)
        value = value.replace("{http://www.usgs.gov/ansseqmsg}","")
    else:
        value = subdic

    tag = str(element.tag)
    tag = tag.replace("{http://www.usgs.gov/ansseqmsg}","")
    tag = tag.replace("{http://quakeml.org/xmlns/quakeml/1.2}","")
    dic = value

    return dic

def read_eqxml(file,eventName):

    tree = etree.parse(file)
    root = tree.getroot()
    dic = parse_xml(root)

    attributeDic = {}

    if not "Origin" in dic["Event"] or not "Magnitude" in dic["Event"]["Origin"]:
        return attributeDic

    attributeDic["Longitude"] = float(dic["Event"]["Origin"]["Longitude"])
    attributeDic["Latitude"] = float(dic["Event"]["Origin"]["Latitude"])
    attributeDic["Depth"] = float(dic["Event"]["Origin"]["Depth"])
    attributeDic["eventID"] = dic["Event"]["EventID"]
    attributeDic["eventName"] = eventName
    attributeDic["Magnitude"] = float(dic["Event"]["Origin"]["Magnitude"]["Value"])

    if "Region" in dic["Event"]["Origin"]:
        attributeDic["Region"] = dic["Event"]["Origin"]["Region"]
    else:
        attributeDic["Region"] = "N/A"

    attributeDic["Time"] = dic["Event"]["Origin"]["Time"]
    timeString = attributeDic["Time"].replace("T"," ").replace("Z","")
    dt = datetime.strptime(timeString, "%Y-%m-%d %H:%M:%S.%f")
    tm = time.struct_time(dt.timetuple())
    attributeDic['GPS'] = float(XLALUTCToGPS(tm))
    attributeDic['UTC'] = float(dt.strftime("%s"))

    attributeDic["Sent"] = dic["Sent"]
    timeString = attributeDic["Sent"].replace("T"," ").replace("Z","")
    dt = datetime.strptime(timeString, "%Y-%m-%d %H:%M:%S.%f")
    tm = time.struct_time(dt.timetuple())
    attributeDic['SentGPS'] = float(XLALUTCToGPS(tm))
    attributeDic['SentUTC'] = float(dt.strftime("%s"))

    attributeDic["DataSource"] = dic["Source"]
    attributeDic["Version"] = dic["Event"]["Version"]

    if "Type" in dic["Event"]:
        attributeDic["Type"] = dic["Event"]["Type"]
    else:
        attributeDic["Type"] = "N/A"  

    if dic["Event"]["Origin"]["Status"] == "Automatic":
        attributeDic["Review"] = "Automatic"
    else:
        attributeDic["Review"] = "Manual"

    attributeDic = calculate_traveltimes(attributeDic)
    tm = time.struct_time(time.gmtime())
    attributeDic['WrittenGPS'] = float(XLALUTCToGPS(tm))
    attributeDic['WrittenUTC'] = float(time.time())

    return attributeDic

def read_quakeml(file,eventName):

    tree = etree.parse(file)
    root = tree.getroot()
    dic = parse_xml(root)

    attributeDic = {}
    attributeDic["Longitude"] = float(dic["eventParameters"]["event"]["origin"]["longitude"]["value"])
    attributeDic["Latitude"] = float(dic["eventParameters"]["event"]["origin"]["latitude"]["value"])
    attributeDic["Depth"] = float(dic["eventParameters"]["event"]["origin"]["depth"]["value"]) / 1000
    attributeDic["eventID"] = ""
    attributeDic["eventName"] = eventName
    attributeDic["Magnitude"] = float(dic["eventParameters"]["event"]["magnitude"]["mag"]["value"])

    attributeDic["Time"] = dic["eventParameters"]["event"]["origin"]["time"]["value"]
    timeString = attributeDic["Time"].replace("T"," ").replace("Z","")
    dt = datetime.strptime(timeString, "%Y-%m-%d %H:%M:%S.%f")
    tm = time.struct_time(dt.timetuple())
    attributeDic['GPS'] = float(XLALUTCToGPS(tm))
    attributeDic['UTC'] = float(dt.strftime("%s"))

    attributeDic["Sent"] = dic["eventParameters"]["event"]["creationInfo"]["creationTime"]
    timeString = attributeDic["Sent"].replace("T"," ").replace("Z","")
    dt = datetime.strptime(timeString, "%Y-%m-%d %H:%M:%S.%f")
    tm = time.struct_time(dt.timetuple())
    attributeDic['SentGPS'] = float(XLALUTCToGPS(tm))
    attributeDic['SentUTC'] = float(dt.strftime("%s"))

    attributeDic["DataSource"] = dic["eventParameters"]["event"]["creationInfo"]["agencyID"]
    attributeDic["Version"] = float(dic["eventParameters"]["event"]["creationInfo"]["version"])
    attributeDic["Type"] = dic["eventParameters"]["event"]["type"]

    if dic["eventParameters"]["event"]["origin"]["evaluationMode"] == "automatic":
        attributeDic["Review"] = "Automatic"
    else:
        attributeDic["Review"] = "Manual"

    attributeDic = traveltimes(attributeDic)
    tm = time.struct_time(time.gmtime())
    attributeDic['WrittenGPS'] = float(XLALUTCToGPS(tm))
    attributeDic['WrittenUTC'] = float(time.time())

    return attributeDic

def read_eqmon(params,file):

    attributeDic = {}
    tree = etree.parse(file)
    root = tree.getroot()       # get the document root
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
    return attributeDic

def jsonread(event):

    attributeDic = {}
    attributeDic["Longitude"] = event["geometry"]["coordinates"][0]
    attributeDic["Latitude"] = event["geometry"]["coordinates"][1]
    attributeDic["Depth"] = event["geometry"]["coordinates"][2]
    attributeDic["eventID"] = event["properties"]["code"]
    attributeDic["eventName"] = event["properties"]["ids"].replace(",","")
    attributeDic["Magnitude"] = event["properties"]["mag"]
    attributeDic["UTC"] = float(event["properties"]["time"])
    attributeDic["DataSource"] = event["properties"]["sources"].replace(",","")
    attributeDic["Version"] = 1.0
    attributeDic["Type"] = 1.0
    attributeDic['Region'] = event["properties"]["place"]

    if event["properties"]["status"] == "AUTOMATIC":
        attributeDic["Review"] = "Automatic"
    else:
        attributeDic["Review"] = "Manual"

    Time = time.gmtime(attributeDic["UTC"])
    attributeDic['GPS'] = float(XLALUTCToGPS(Time))
    SentTime = time.gmtime()
    attributeDic['SentGPS'] = float(XLALUTCToGPS(SentTime))
    attributeDic['SentUTC'] = time.time()

    attributeDic['Time'] = time.strftime("%Y-%m-%dT%H:%M:%S.000Z", Time)
    attributeDic['Sent'] = time.strftime("%Y-%m-%dT%H:%M:%S.000Z", SentTime)

    attributeDic = calculate_traveltimes(attributeDic)
    tm = time.struct_time(time.gmtime())
    attributeDic['WrittenGPS'] = float(XLALUTCToGPS(tm))
    attributeDic['WrittenUTC'] = float(time.time())

    return attributeDic

def databaseread(event):

    attributeDic = {}
    eventSplit = event.split(",")

    year = int(eventSplit[0])
    month = int(eventSplit[1])
    day = int(eventSplit[2])
    hour = int(eventSplit[3][0:1])
    minute = int(eventSplit[3][0:1])
    second = int(eventSplit[3][0:1])

    timeString = "%d-%02d-%02d %02d:%02d:%02d"%(year,month,day,hour,minute,second)
    dt = datetime.strptime(timeString, "%Y-%m-%d %H:%M:%S")

    eventID = ''.join(random.sample(string.digits,8))
    eventName = ''.join(["db",eventID])

    attributeDic["Longitude"] = float(eventSplit[5])
    attributeDic["Latitude"] = float(eventSplit[4])
    attributeDic["Depth"] = float(eventSplit[7])
    attributeDic["eventID"] = float(eventID)
    attributeDic["eventName"] = eventName
    try:
        attributeDic["Magnitude"] = float(eventSplit[6])
    except:
        attributeDic["Magnitude"] = 0
    tm = time.struct_time(dt.timetuple())
    attributeDic['GPS'] = float(XLALUTCToGPS(tm))
    attributeDic['UTC'] = float(dt.strftime("%s"))
    attributeDic["DataSource"] = "DB"
    attributeDic["Version"] = 1.0
    attributeDic["Type"] = 1.0
    attributeDic['Region'] = "N/A"
    attributeDic["Review"] = "Manual"

    SentTime = time.gmtime()
    attributeDic['SentGPS'] = float(XLALUTCToGPS(SentTime))
    attributeDic['SentUTC'] = time.time()

    attributeDic['Time'] = time.strftime("%Y-%m-%dT%H:%M:%S.000Z", tm)
    attributeDic['Sent'] = time.strftime("%Y-%m-%dT%H:%M:%S.000Z", SentTime)

    attributeDic = calculate_traveltimes(attributeDic)
    tm = time.struct_time(time.gmtime())
    attributeDic['WrittenGPS'] = float(XLALUTCToGPS(tm))
    attributeDic['WrittenUTC'] = float(time.time())

    return attributeDic

def calculate_traveltimes(attributeDic): 

    if not "traveltimes" in attributeDic:
        attributeDic["traveltimes"] = {}

    if not "Latitude" in attributeDic and not "Longitude" in attributeDic:
        return attributeDic

    attributeDic = ifotraveltimes(attributeDic, "LHO", 46.6475, -119.5986)
    attributeDic = ifotraveltimes(attributeDic, "LLO", 30.4986, -90.7483)
    attributeDic = ifotraveltimes(attributeDic, "GEO", 52.246944, 9.808333)
    attributeDic = ifotraveltimes(attributeDic, "VIRGO", 43.631389, 10.505)
    attributeDic = ifotraveltimes(attributeDic, "FortyMeter", 34.1391, -118.1238)
    attributeDic = ifotraveltimes(attributeDic, "Homestake", 44.3465, -103.7574)

    return attributeDic

def do_kdtree(combined_x_y_arrays,points):
    mytree = scipy.spatial.cKDTree(combined_x_y_arrays)
    dist, indexes = mytree.query(points)
    return indexes

def ifotraveltimes(attributeDic,ifo,ifolat,ifolon):

    try:
        from obspy.taup.taup import getTravelTimes
        from obspy.core.util.geodetics import gps2DistAzimuth
    except:
        print "Enable ObsPy if updated earthquake estimates desired...\n"
        return attributeDic

    distance,fwd,back = gps2DistAzimuth(attributeDic["Latitude"],attributeDic["Longitude"],ifolat,ifolon)
    distances = np.linspace(0,distance,1000)
    degrees = (distances/6370000)*(180/np.pi)

    distance_delta = distances[1] - distances[0]

    periods = [25.0,27.0,30.0,32.0,35.0,40.0,45.0,50.0,60.0,75.0,100.0,125.0,150.0,200.0,250.0]
    frequencies = 1 / np.array(periods)
   
    c = 18
    fc = 10**(2.3-(attributeDic["Magnitude"]/2.))
    Q = np.max([500,80/np.sqrt(fc)])

    Rfamp = ((attributeDic["Magnitude"]/fc)*0.0035) * np.exp(-2*math.pi*attributeDic["Depth"]*fc/c) * np.exp(-2*math.pi*(distances[-1]/1000)*(fc/c)*1/Q)/(distances[-1]/1000)
    Pamp = 1e-6
    Samp = 1e-5

    index = np.argmin(np.absolute(frequencies - fc))

    lats = []
    lons = []
    Ptimes = []
    Stimes = []
    Rtimes = []
    Rfamps = []

    velocityFile = '/home/mcoughlin/Seismon/velocity_maps/GR025_1_GDM52.pix'
    velocity_map = np.loadtxt(velocityFile)
    base_velocity = 3.59738 

    for distance, degree in zip(distances, degrees):

        lon, lat, baz = shoot(attributeDic["Longitude"], attributeDic["Latitude"], fwd, distance/1000)
        lats.append(lat)
        lons.append(lon)

    combined_x_y_arrays = np.dstack([velocity_map[:,0],velocity_map[:,1]])[0]
    points_list = np.dstack([lats, lons])

    indexes = do_kdtree(combined_x_y_arrays,points_list)[0]

    time = 0

    for distance, degree, index in zip(distances, degrees,indexes):

        velocity = 1000 * (1 + 0.01*velocity_map[index,3])*base_velocity

        time_delta = distance_delta / velocity
        time = time + time_delta

        #degrees = locations2degrees(lat,lon,attributeDic["Latitude"],attributeDic["Longitude"])
        #distance,fwd,back = gps2DistAzimuth(lat,lon,attributeDic["Latitude"],attributeDic["Longitude"])
        tt = getTravelTimes(delta=degree, depth=attributeDic["Depth"])
        tt.append({'phase_name': 'R', 'dT/dD': 0, 'take-off angle': 0, 'time': time, 'd2T/dD2': 0, 'dT/dh': 0})
        Ptime = -1
        Stime = -1
        Rtime = -1
        for phase in tt:
            if Ptime == -1 and phase["phase_name"][0] == "P":
                Ptime = attributeDic["GPS"]+phase["time"]
            if Stime == -1 and phase["phase_name"][0] == "S":
                Stime = attributeDic["GPS"]+phase["time"]
            if Rtime == -1 and phase["phase_name"][0] == "R":
                Rtime = attributeDic["GPS"]+phase["time"]

        Ptimes.append(Ptime)
        Stimes.append(Stime)
        Rtimes.append(Rtime)

    #if ifo == "LHO":
    #    print time - distance / 3500

    traveltimes = {}
    traveltimes["Latitudes"] = lats
    traveltimes["Longitudes"] = lons
    traveltimes["Distances"] = distances
    traveltimes["Degrees"] = degrees
    traveltimes["Ptimes"] = Ptimes
    traveltimes["Stimes"] = Stimes
    traveltimes["Rtimes"] = Rtimes
    traveltimes["Rfamp"] = [Rfamp] 
    traveltimes["Pamp"] = [Pamp]
    traveltimes["Samp"] = [Samp]

    attributeDic["traveltimes"][ifo] = traveltimes

    return attributeDic

def GPSToUTCDateTime(gps):

    utc = XLALGPSToUTC(LIGOTimeGPS(int(gps)))
    tt = time.strftime("%Y-%jT%H:%M:%S", utc)
    ttUTC = UTCDateTime(tt)

    return ttUTC    

def attribute_array(attributeDics,type):

    array = []
    for attributeDic in attributeDics:
        if len(type) == 1:
            attribute = attributeDic[type[0]]
        elif len(type) == 2:
            attribute = attributeDic[type[0]][type[1]]
        elif len(type) == 3:
            attribute = attributeDic[type[0]][type[1]][type[2]]
        elif len(type) == 4:
            attribute = attributeDic[type[0]][type[1]][type[2]][type[3]]

        array.append(attribute)

    return array

def eventDiff(attributeDics, magnitudeDiff, latitudeDiff, longitudeDiff):

    if len(attributeDics) > 1:
        for i in xrange(len(attributeDics)-1):
            if "Magnitude" in attributeDics[i] and "Magnitude" in attributeDics[i+1] and \
                "Latitude" in attributeDics[i] and "Latitude" in attributeDics[i+1] and\
                "Longitude" in attributeDics[i] and "Longitude" in attributeDics[i+1]:

                magnitudeDiff.append(attributeDics[i]["Magnitude"]-attributeDics[i+1]["Magnitude"])
                latitudeDiff.append(attributeDics[i]["Latitude"]-attributeDics[i+1]["Latitude"])
                longitudeDiff.append(attributeDics[i]["Longitude"]-attributeDics[i+1]["Longitude"])
    return magnitudeDiff, latitudeDiff, longitudeDiff

def great_circle_distance(latlong_a, latlong_b):

    EARTH_CIRCUMFERENCE = 6378.137 # earth circumference in kilometers

    lat1, lon1 = latlong_a
    lat2, lon2 = latlong_b

    dLat = math.radians(lat2 - lat1)
    dLon = math.radians(lon2 - lon1)
    a = (math.sin(dLat / 2) * math.sin(dLat / 2) +
            math.cos(math.radians(lat1)) * math.cos(math.radians(lat2)) *
            math.sin(dLon / 2) * math.sin(dLon / 2))
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    d = EARTH_CIRCUMFERENCE * c
    
    return d

def retrieve_earthquakes(params):

    attributeDics = []
    files = glob.glob(os.path.join(params["eventfilesLocation"],"*.xml"))

#    for numFile in xrange(100):
    for numFile in xrange(len(files)):

        file = files[numFile]

        fileSplit = file.replace(".xml","").split("-")
        gps = float(fileSplit[-1])
        if (gps < params["gpsStart"] - 3600) or (gps > params["gpsEnd"]):
           continue

        attributeDic = read_eqmon(params,file)

        if attributeDic["Magnitude"] >= params["earthquakesMinMag"]:
            attributeDics.append(attributeDic)

    return attributeDics

def retrieve_earthquakes_text(params):
    attributeDics = []
    files = glob.glob(os.path.join(params["eventfilesLocation"],"*.txt"))
    for numFile in xrange(len(files)):

        file = files[numFile]

        fileSplit = file.split("/")
        fileName = fileSplit[-1]

        eventName = fileName.replace(".txt","").split("-")[0]
        gps = int(fileName.replace(".txt","").split("-")[1])

        #if gps + 10000 < params["gps"]:
        #    continue

        attributeDic = {}
        attributeDic["eventName"] = eventName
        attributeDic["traveltimes"] = {}
        with open(file,'r') as f:
            counter = 0
            attributeDic["traveltimes"]["LHO"] = []
            attributeDic["traveltimes"]["LLO"] = []
            for line in f:
                counter += 1
                line_without_return = line.split("\n")
                space_split = line_without_return[0].split(" ")
                colon_split = line_without_return[0].split(";")
                if len(colon_split) > 1:
                    if colon_split[0] == "traveltimes":
                        continue
                    try:
                        attributeDic[colon_split[0]] = float(colon_split[1])
                    except:
                        attributeDic[colon_split[0]] = colon_split[1]
                elif len(space_split) > 1:
                    if counter < 4:
                        attributeDic["traveltimes"]["LHO"].append(float(space_split[-2]))
                    else:
                        attributeDic["traveltimes"]["LLO"].append(float(space_split[-2]))

        magThreshold = 0
        if not "Magnitude" in attributeDic or attributeDic["Magnitude"] < magThreshold:
            continue
        attributeDic["doPlots"] = 0
        for ifoName, traveltimes in attributeDic["traveltimes"].items():
            if params["gps"] <= max(traveltimes):
                attributeDic["doPlots"] = 1
        attributeDics.append(attributeDic)

    return attributeDics

def equi(m, centerlon, centerlat, radius):
    glon1 = centerlon
    glat1 = centerlat
    X = []
    Y = []
    for azimuth in range(0, 360):
        glon2, glat2, baz = shoot(glon1, glat1, azimuth, radius)
        X.append(glon2)
        Y.append(glat2)
    X.append(X[0])
    Y.append(Y[0])

    #m.plot(X,Y,**kwargs) #Should work, but doesn't...
    X,Y = m(X,Y)
    return X,Y

def shoot(lon, lat, azimuth, maxdist=None):
    """Shooter Function
    Original javascript on http://williams.best.vwh.net/gccalc.htm
    Translated to python by Thomas Lecocq
    """
    glat1 = lat * np.pi / 180.
    glon1 = lon * np.pi / 180.
    s = maxdist / 1.852
    faz = azimuth * np.pi / 180.

    EPS= 0.00000000005
    if ((np.abs(np.cos(glat1))<EPS) and not (np.abs(np.sin(faz))<EPS)):
        alert("Only N-S courses are meaningful, starting at a pole!")

    a=6378.13/1.852
    f=1/298.257223563
    r = 1 - f
    tu = r * np.tan(glat1)
    sf = np.sin(faz)
    cf = np.cos(faz)
    if (cf==0):
        b=0.
    else:
        b=2. * np.arctan2 (tu, cf)

    cu = 1. / np.sqrt(1 + tu * tu)
    su = tu * cu
    sa = cu * sf
    c2a = 1 - sa * sa
    x = 1. + np.sqrt(1. + c2a * (1. / (r * r) - 1.))
    x = (x - 2.) / x
    c = 1. - x
    c = (x * x / 4. + 1.) / c
    d = (0.375 * x * x - 1.) * x
    tu = s / (r * a * c)
    y = tu
    c = y + 1
    while (np.abs (y - c) > EPS):

        sy = np.sin(y)
        cy = np.cos(y)
        cz = np.cos(b + y)
        e = 2. * cz * cz - 1.
        c = y
        x = e * cy
        y = e + e - 1.
        y = (((sy * sy * 4. - 3.) * y * cz * d / 6. + x) *
              d / 4. - cz) * sy * d + tu

    b = cu * cy * cf - su * sy
    c = r * np.sqrt(sa * sa + b * b)
    d = su * cy + cu * sy * cf
    glat2 = (np.arctan2(d, c) + np.pi) % (2*np.pi) - np.pi
    c = cu * cy - su * sy * cf
    x = np.arctan2(sy * sf, c)
    c = ((-3. * c2a + 4.) * f + 4.) * c2a * f / 16.
    d = ((e * cy * c + cz) * sy * c + y) * sa
    glon2 = ((glon1 + x - (1. - c) * d * f + np.pi) % (2*np.pi)) - np.pi

    baz = (np.arctan2(sa, b) + np.pi) % (2 * np.pi)

    glon2 *= 180./np.pi
    glat2 *= 180./np.pi
    baz *= 180./np.pi

    return (glon2, glat2, baz)

