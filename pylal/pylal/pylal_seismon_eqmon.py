#!/usr/bin/python

import os, sys, time, glob, numpy, math, matplotlib, random, string
from datetime import datetime
from operator import itemgetter
import xml.dom.minidom
from pylal.xlal.datatypes.ligotimegps import LIGOTimeGPS
from pylal.xlal.date import XLALUTCToGPS, XLALGPSToUTC
from pylal import Fr
from collections import namedtuple
from lxml import etree

import pylal.pylal_seismon_eqmon_plot

def run_earthquakes(params):

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

    attributeDics = retrieve_earthquakes(params)
    attributeDics = sorted(attributeDics, key=itemgetter("Magnitude"), reverse=True)

    f = open(os.path.join(earthquakesDirectory,"earthquakes.txt"),"w+")

    for attributeDic in attributeDics:
        traveltimes = attributeDic["traveltimes"][ifo]

        gpsStart = traveltimes["arrivalMin"] - 200
        gpsEnd = traveltimes["arrivalMax"] + 200

        f.write("%.1f %.1f %.1f %.1f %.1f %.5e\n"%(attributeDic["GPS"],attributeDic["Magnitude"],max(traveltimes["Ptimes"]),max(traveltimes["Stimes"]),max(traveltimes["Rtimes"]),traveltimes["Rfamp"][0]))

    f.close()

    if params["doPlots"]:
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
        pylal.pylal_seismon_eqmon_plot.worldmap_plot(params,attributeDics,params["gpsEnd"],plotName)

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

    attributeDic = traveltimes(attributeDic)
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

    attributeDic = traveltimes(attributeDic)
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

    attributeDic = traveltimes(attributeDic)
    tm = time.struct_time(time.gmtime())
    attributeDic['WrittenGPS'] = float(XLALUTCToGPS(tm))
    attributeDic['WrittenUTC'] = float(time.time())

    return attributeDic

def traveltimes(attributeDic): 

    attributeDic["traveltimes"] = {}

    if not "Latitude" in attributeDic and not "Longitude" in attributeDic:
        return attributeDic

    attributeDic = ifotraveltimes(attributeDic, "LHO", 46.6475, -119.5986)
    attributeDic = ifotraveltimes(attributeDic, "LLO", 30.4986, -90.7483)
    attributeDic = ifotraveltimes(attributeDic, "GEO", 52.246944, 9.808333)
    attributeDic = ifotraveltimes(attributeDic, "VIRGO", 43.631389, 10.505)
    attributeDic = ifotraveltimes(attributeDic, "FortyMeter", 34.1391, -118.1238)
    attributeDic = ifotraveltimes(attributeDic, "Homestake", 44.3465, -103.7574)

    #attributeDic["distanceLHO"] = attributeDic["traveltimes"]["LHO"]["Distances"][-1]
    #attributeDic["distanceLLO"] = attributeDic["traveltimes"]["LLO"]["Distances"][-1]

    return attributeDic

def ifotraveltimes(attributeDic,ifo,ifolat,ifolon):

    distance,fwd,back = gps2DistAzimuth(attributeDic["Latitude"],attributeDic["Longitude"],ifolat,ifolon)
    distances = numpy.linspace(0,distance,100)
    degrees = (distances/6370000)*(180/numpy.pi)

    lats = []
    lons = []
    Ptimes = []
    Stimes = []
    Rtimes = []
    Rfamps = []

    # Pmag = T * 10^(Mb - 5.9 - 0.01*dist)
    # Rmag = T * 10^(Ms - 3.3 - 1.66*log_10(dist))
    T = 20

    for distance, degree in zip(distances, degrees):

        lon, lat, baz = shoot(attributeDic["Longitude"], attributeDic["Latitude"], fwd, distance/1000)
        lats.append(lat)
        lons.append(lon)

        #degrees = locations2degrees(lat,lon,attributeDic["Latitude"],attributeDic["Longitude"])
        #distance,fwd,back = gps2DistAzimuth(lat,lon,attributeDic["Latitude"],attributeDic["Longitude"])
        tt = getTravelTimes(delta=degree, depth=attributeDic["Depth"])
        tt.append({'phase_name': 'R', 'dT/dD': 0, 'take-off angle': 0, 'time': distance/3500, 'd2T/dD2': 0, 'dT/dh': 0})
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

    traveltimes = {}
    traveltimes["Latitudes"] = lats
    traveltimes["Longitudes"] = lons
    traveltimes["Distances"] = distances
    traveltimes["Degrees"] = degrees
    traveltimes["Ptimes"] = Ptimes
    traveltimes["Stimes"] = Stimes
    traveltimes["Rtimes"] = Rtimes

    c = 18
    fc = 10**(2.3-(attributeDic["Magnitude"]/2.))
    Q = numpy.max([500,80/numpy.sqrt(fc)])

    Rfamp = ((attributeDic["Magnitude"]/fc)*0.0035) * numpy.exp(-2*math.pi*attributeDic["Depth"]*fc/c) * numpy.exp(-2*math.pi*(distances[-1]/1000)*(fc/c)*1/Q)/(distances[-1]/1000)

    traveltimes["Rfamp"] = [Rfamp]

    Pamp = 1e-6
    Samp = 1e-5
    
    traveltimes["Pamp"] = [Pamp]
    traveltimes["Samp"] = [Samp]

    attributeDic["traveltimes"][ifo] = traveltimes
    return attributeDic

def event_distance(attributeDic,trace):

    distance = -1
    traveltimes = []

    with open("/home/mcoughlin/git-repo/eqmon/input/eqmon-EARTHWORM-channel_list.txt") as f:

       for line in f:

           line_without_return = line.split("\n")
           line_split = line_without_return[0].split(" ")

           org = line_split[0]
           nw = line_split[1]
           sta = line_split[2]
           lat = float(line_split[3])
           lon = float(line_split[4])

           if trace.stats.network == nw and trace.stats.station == sta:
               latitude = lat
               longitude = lon
               distance,fwd,back = gps2DistAzimuth(lat,lon,attributeDic["Latitude"],attributeDic["Longitude"])

               degrees = locations2degrees(lat,lon,attributeDic["Latitude"],attributeDic["Longitude"])
               traveltimes = getTravelTimes(delta=degrees, depth=attributeDic["Depth"])
               traveltimes.append({'phase_name': 'R', 'dT/dD': 0, 'take-off angle': 0, 'time': distance/3500, 'd2T/dD2': 0, 'dT/dh': 0})

    return distance, traveltimes

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

def html_bgcolor(snr,data):

    data = numpy.append(data,snr)

    # Number of colors in array
    N = 256

    colormap = []
    for i in xrange(N):
        r,g,b,a = matplotlib.pyplot.cm.jet(i)
        r = int(round((r * 255),0))
        g = int(round((g * 255),0))
        b = int(round((b* 255),0))
        colormap.append((r,g,b))

    data = numpy.sort(data)
    itemIndex = numpy.where(data==snr)

    # Determine significance of snr (between 0 and 1)
    snrSig = itemIndex[0][0] / float(len(data)+1)

    # Determine color index of this significance
    index = int(numpy.floor(N * snrSig))

    # Return colors of this index
    thisColor = colormap[index]
    # Return rgb string containing these colors
    bgcolor = "rgb(%d, %d, %d)"%(thisColor[0],thisColor[1],thisColor[2])

    return snrSig, bgcolor

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

def create_cache(ifo,frameType,gpsStart,gpsEnd):

    p = Popen("ligo_data_find -o %s -t %s -s %.0f -e %.0f -u file"%(ifo,frameType,gpsStart,gpsEnd),shell=True,stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    output = p.stdout.read()
    frameList = output.split("\n")
    frames = []
    for frame in frameList:
        thisFrame = frame.replace("file://localhost","")
        if thisFrame is not "":
            frames.append(thisFrame)

    return frames

def read_frames(params,gpsStart,gpsEnd,channel_name,cache):
    time = []
    data = []

    #== loop over frames in cache
    for frame in cache:
        if frame == "No files found!":
            continue
        try:
            frame_data,data_start,_,dt,_,_ = Fr.frgetvect1d(frame,channel_name)
        except:
            continue
        frame_length = float(dt)*len(frame_data)
        frame_time = data_start+dt*numpy.arange(len(frame_data))

        if frame_time[0] >= gpsStart and frame_time[-1] <= gpsEnd:
            for i in range(len(frame_data)):
                time.append(frame_time[i])
                data.append(frame_data[i])
        else:
            for i in range(len(frame_data)):
                if frame_time[i] <= gpsStart:  continue
                if frame_time[i] >= gpsEnd:  continue
                time.append(frame_time[i])
                data.append(frame_data[i])

    return time,data

def channel_struct(channelList):
    # Create channel structure
    structproxy_channel = namedtuple( "structproxy_channel", "name samplef calibration latitude longitude" )

    channel = []

    with open(channelList,'r') as f:

       for line in f:

           line_without_return = line.split("\n")
           line_split = line_without_return[0].split(" ")

           name = line_split[0]
           samplef = float(line_split[1])
           calibration = float(line_split[2])

           if name[0] == "H":
               latitude = 46.6475
               longitude = -119.5986;
           elif name[0] == "L":
               latitude = 30.4986
               longitude = -90.7483
           elif name[0] == "G":
               latitude = 52.246944
               longitude = 9.808333
           elif name[0] == "V":
               latitude = 43.631389
               longitude = 10.505
           elif name[0] == "M":
               latitude = 44.3465
               longitude = -103.7574

           channel.append( structproxy_channel(name,samplef,calibration,latitude,longitude))
    return channel

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

def NLNM(unit):

    PL = [0.1, 0.17, 0.4, 0.8, 1.24, 2.4, 4.3, 5, 6, 10, 12, 15.6, 21.9, 31.6, 45, 70,\
        101, 154, 328, 600, 10000]
    AL = [-162.36, -166.7, -170, -166.4, -168.6, -159.98, -141.1, -71.36, -97.26,\
        -132.18, -205.27, -37.65, -114.37, -160.58, -187.5, -216.47, -185,\
        -168.34, -217.43, -258.28, -346.88]
    BL = [5.64, 0, -8.3, 28.9, 52.48, 29.81, 0, -99.77, -66.49, -31.57, 36.16,\
        -104.33, -47.1, -16.28, 0, 15.7, 0, -7.61, 11.9, 26.6, 48.75]

    PH = [0.1, 0.22, 0.32, 0.8, 3.8, 4.6, 6.3, 7.9, 15.4, 20, 354.8, 10000]
    AH = [-108.73, -150.34, -122.31, -108.48, -116.85, -74.66, 0.66, -93.37, 73.54,\
        -151.52, -206.66, -206.66];
    BH = [-17.23, -80.5, -23.87, 32.51, 18.08, -32.95, -127.18, -22.42, -162.98,\
        10.01, 31.63, 31.63]

    fl = [1/float(e) for e in PL]
    fh = [1/float(e) for e in PH]

    lownoise = []
    highnoise = []
    for i in xrange(len(PL)):
        lownoise.append(10**((AL[i] + BL[i]*math.log10(PL[i]))/20))
    for i in xrange(len(PH)):
        highnoise.append(10**((AH[i] + BH[i]*math.log10(PH[i]))/20))

    for i in xrange(len(PL)):
        if unit == 1:
            lownoise[i] = lownoise[i] * (PL[i]/(2*math.pi))**2
        elif unit==2:
            lownoise[i] = lownoise[i] * (PL[i]/(2*math.pi))
    for i in xrange(len(PH)):
        if unit == 1:
            highnoise[i] = highnoise[i] * (PH[i]/(2*math.pi))**2
        elif unit==2:
            highnoise[i] = highnoise[i] * (PH[i]/(2*math.pi))

    return fl, lownoise, fh, highnoise

def get_KW_triggers(trigger_file,trigger_list):

    with open(trigger_file, 'r') as f:

        for line in f:

            line_without_return = line.split("\n")
            line_split_all = line_without_return[0].split(" ")

            line_split = [d for d in line_split_all \
                             if d != ""]

            gps_time = float(line_split[1])
            SNR = math.sqrt(float(line_split[5]) - float(line_split[6]))
            significance = float(line_split[7])
            trigger_list.append([gps_time,SNR,significance])

    return trigger_list

def run_kleinewelle(params,gpsStart,gpsEnd,channel,frames):

    kleinewelleFolder = os.path.join(params["outputFolder"],"kleinewelle")
    if not os.path.isdir(kleinewelleFolder):
        os.makedirs(kleinewelleFolder)

    kleinewelleParamsFolder = os.path.join(kleinewelleFolder,"params")
    if not os.path.isdir(kleinewelleParamsFolder):
        os.makedirs(kleinewelleParamsFolder)

    kleinewelleParams = """
    stride 128
    basename %s
    segname segments
    significance 20
    threshold 3.0
    decimateFactor -1
    channel %s 1 8
    channel %s 8 64
    channel %s 64 1024
    """ % (channel.replace(":","_"),channel,channel,channel)

    f = open("%s/%s.txt"%(kleinewelleParamsFolder,channel.replace(":","_")),"w")
    f.write(kleinewelleParams)
    f.close()

    kleinewelleCacheFolder = os.path.join(kleinewelleFolder,"cache")
    if not os.path.isdir(kleinewelleCacheFolder):
        os.makedirs(kleinewelleCacheFolder)

    f = open("%s/%s.txt"%(kleinewelleCacheFolder,channel.replace(":","_")),"w")
    for frame in frames:
        f.write(frame + "\n")
    f.close()

    os.chdir(kleinewelleTriggersUnfilteredFolder)
    os.system("kleineWelleM %s/%s.txt -inlist %s/%s.txt"%(kleinewelleParamsFolder,channel.replace(":","_"),kleinewelleCacheFolder,channel.replace(":","_")))
    triggerFolders = glob.glob(os.path.join(kleinewelleTriggersUnfilteredFolder,"%s*"%channel.replace(":","_")))
    trigger_list = []
    for folder in triggerFolders:
        triggerFiles = glob.glob(folder + "/*.trg")
        for file in triggerFiles:
            trigger_list = get_KW_triggers(file,trigger_list)
        os.system("rm -r -f %s"%(folder))

    time = []
    data = []
    for trigger in trigger_list:
        if trigger[0] >= gpsStart and trigger[0] <= gpsEnd:
            time.append(trigger[0])
            data.append(trigger[1])

    return time,data

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
    glat1 = lat * numpy.pi / 180.
    glon1 = lon * numpy.pi / 180.
    s = maxdist / 1.852
    faz = azimuth * numpy.pi / 180.

    EPS= 0.00000000005
    if ((numpy.abs(numpy.cos(glat1))<EPS) and not (numpy.abs(numpy.sin(faz))<EPS)):
        alert("Only N-S courses are meaningful, starting at a pole!")

    a=6378.13/1.852
    f=1/298.257223563
    r = 1 - f
    tu = r * numpy.tan(glat1)
    sf = numpy.sin(faz)
    cf = numpy.cos(faz)
    if (cf==0):
        b=0.
    else:
        b=2. * numpy.arctan2 (tu, cf)

    cu = 1. / numpy.sqrt(1 + tu * tu)
    su = tu * cu
    sa = cu * sf
    c2a = 1 - sa * sa
    x = 1. + numpy.sqrt(1. + c2a * (1. / (r * r) - 1.))
    x = (x - 2.) / x
    c = 1. - x
    c = (x * x / 4. + 1.) / c
    d = (0.375 * x * x - 1.) * x
    tu = s / (r * a * c)
    y = tu
    c = y + 1
    while (numpy.abs (y - c) > EPS):

        sy = numpy.sin(y)
        cy = numpy.cos(y)
        cz = numpy.cos(b + y)
        e = 2. * cz * cz - 1.
        c = y
        x = e * cy
        y = e + e - 1.
        y = (((sy * sy * 4. - 3.) * y * cz * d / 6. + x) *
              d / 4. - cz) * sy * d + tu

    b = cu * cy * cf - su * sy
    c = r * numpy.sqrt(sa * sa + b * b)
    d = su * cy + cu * sy * cf
    glat2 = (numpy.arctan2(d, c) + numpy.pi) % (2*numpy.pi) - numpy.pi
    c = cu * cy - su * sy * cf
    x = numpy.arctan2(sy * sf, c)
    c = ((-3. * c2a + 4.) * f + 4.) * c2a * f / 16.
    d = ((e * cy * c + cz) * sy * c + y) * sa
    glon2 = ((glon1 + x - (1. - c) * d * f + numpy.pi) % (2*numpy.pi)) - numpy.pi

    baz = (numpy.arctan2(sa, b) + numpy.pi) % (2 * numpy.pi)

    glon2 *= 180./numpy.pi
    glat2 *= 180./numpy.pi
    baz *= 180./numpy.pi

    return (glon2, glat2, baz)

