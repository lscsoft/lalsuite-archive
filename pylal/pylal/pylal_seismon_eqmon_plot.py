#!/usr/bin/python

"""
%prog

Michael Coughlin (coughlim@carleton.edu)

This program checks for earthquakes.

"""

import os, time, glob, matplotlib, math
import numpy as np
from mpl_toolkits.basemap import Basemap
matplotlib.use("AGG")
matplotlib.rcParams.update({'font.size': 18})
from subprocess import Popen, PIPE, STDOUT
from pylal import Fr 
from pylab import *
from matplotlib import pyplot as plt

import pylal.pylal_seismon_eqmon

# =============================================================================
#
#                               DEFINITIONS
#
# =============================================================================

def plot_envelope(params,time,data,traveltimes,plotName):

    plotNameSplit = plotName.split("/")
    plotTitle = plotNameSplit[-1].replace(".png","")

    time = np.array(time)
    data = np.array(data)

    startTime = np.min(time)
    endTime = np.max(time)

    time = time - startTime

    Ptime = traveltimes["Ptimes"][-1] - startTime
    Stime = traveltimes["Stimes"][-1] - startTime
    Rtime = traveltimes["Rtimes"][-1] - startTime

    plt.plot(time,data, 'k')
    plt.axvline(x=Ptime,color='r')
    plt.axvline(x=Stime,color='b')
    plt.axvline(x=Rtime,color='g')
    plt.text(Ptime, -0.05, 'P', fontsize=18, ha='center', va='top')
    plt.text(Stime, -0.05, 'S', fontsize=18, ha='center', va='top')
    plt.text(Rtime, -0.05, 'R', fontsize=18, ha='center', va='top')
    #xlim([min(time), max(time)])
    plt.xlim([0, endTime-startTime])
    plt.ylim([0, 1])
    plt.xlabel('Time [s]')
    plt.ylabel("%.0f - %.0f"%(startTime,endTime))
    plt.title(plotTitle)
    plt.show()
    plt.savefig(plotName,dpi=200)
    plt.close('all')

def restimates(params,attributeDics,plotName):

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

    gps = []
    magnitudes = []
    for attributeDic in attributeDics:

        if "Restimate" not in attributeDic["traveltimes"][ifo]:
            continue

        travel_time = attributeDic["traveltimes"][ifo]["Restimate"] - attributeDic["traveltimes"][ifo]["Rtimes"][-1]

        gps.append(travel_time)
        magnitudes.append(attributeDic["Magnitude"])

    if magnitudes == []:
        return

    gps = np.array(gps)
    magnitudes = np.array(magnitudes)

    plt.figure()
    plt.plot(magnitudes,gps, 'k*')
    plt.xlim([min(magnitudes)-0.5, max(magnitudes)+0.5])
    plt.xlabel('Magnitude')
    plt.ylabel("\Delta T")
    plt.show()
    plt.savefig(plotName,dpi=200)
    plt.close('all')

def prediction(data,plotName):

    if len(data["prediction"]["tt"]) == 0:
        return

    timeStart = data["prediction"]["tt"][0]
    t = data["prediction"]["tt"] - timeStart
    t_prediction = t / 86400.0

    amp = np.log10(data["prediction"]["data"])
    indexes = np.isinf(amp)
    amp[indexes] = -250
    amp_prediction = amp

    threshold = -8

    plt.figure()
    plt.plot(t_prediction,amp_prediction,'*',label="Predicted")
    plt.plot(t_prediction,np.ones(t_prediction.shape)*threshold,label='Threshold')

    for key in data["channels"].iterkeys():

        t = data["channels"][key]["tt"] - timeStart
        t = t / 86400.0

        amp = np.log10(data["channels"][key]["data"])
        indexes = np.isinf(amp)
        amp[indexes] = -250
        plt.plot(t,amp,'*',label=key)

        amp_interp = np.interp(t_prediction,t,amp)
        residual = np.absolute(amp_interp - amp_prediction)

        residual_label = "%s residual"%key
        #plt.plot(t_prediction,residual,'*',label=residual_label)

    plt.legend(loc=1,prop={'size':6})
    plt.xlim([np.min(t_prediction),np.max(t_prediction)])
    plt.ylim([-10,-2])
    plt.xlabel('Time [Days] [%d]'%timeStart)
    plt.ylabel('Mean ASD in 0.02-0.1 Hz Band [log10(m/s)]')
    plt.show()
    plt.savefig(plotName,dpi=200)
    plt.close('all')

def residual(data,plotName):

    if len(data["prediction"]["tt"]) == 0:
        return

    timeStart = data["prediction"]["tt"][0]
    t = data["prediction"]["tt"] - timeStart
    t_prediction = t / 86400.0

    amp = np.log10(data["prediction"]["data"])
    indexes = np.isinf(amp)
    amp[indexes] = -250
    amp_prediction = amp

    threshold = -8

    indexes = np.where(amp_prediction >= threshold)
    t_prediction = t_prediction[indexes]
    amp_prediction = amp_prediction[indexes]

    plt.figure()
    #plt.plot(t_prediction,amp_prediction,'*',label="Predicted")
    #plt.plot(t_prediction,np.ones(t_prediction.shape)*threshold,label='Threshold')

    for key in data["channels"].iterkeys():

        t = data["channels"][key]["tt"] - timeStart
        t = t / 86400.0

        amp = np.log10(data["channels"][key]["data"])
        indexes = np.isinf(amp)
        amp[indexes] = -250
        #plt.plot(t,amp,'*',label=key)

        amp_interp = np.interp(t_prediction,t,amp)
        residual = amp_interp - amp_prediction

        residual_label = "%s residual"%key
        plt.plot(t_prediction,residual,'*',label=residual_label)

    plt.legend(loc=1,prop={'size':6})
    plt.xlim([np.min(t_prediction),np.max(t_prediction)])
    plt.ylim([-5,5])
    plt.xlabel('Time [Days] [%d]'%timeStart)
    plt.ylabel('Residual in prediction [log10(m/s)]')
    plt.show()
    plt.savefig(plotName,dpi=200)
    plt.close('all')

def efficiency(data,plotName):

    data = np.array(data)
    t = data[:,0] - data[0,0]

    plt.figure()
    plt.plot(t,np.log10(data[:,2]),marker='*',label="Data")
    plt.plot(t,np.log10(data[:,3]),marker='*',label="Predicted")
    plt.legend(loc=3)
    plt.xlabel('Time [s] [%d]'%data[0,0])
    plt.ylabel('Amplitude')
    plt.show()
    plt.savefig(plotName,dpi=200)
    plt.close('all')

def efficiency_limits(data,plotName):

    data = np.array(data)
    t = data[:,0] - data[0,0]

    plt.figure()
    plt.plot(t,np.log10(data[:,2]),marker='*',label="Data")
    plt.plot(t,np.log10(data[:,3]),marker='*',label="Predicted")
    plt.legend(loc=3)
    plt.ylim([-10,-2])
    plt.xlabel('Time [s] [%d]'%data[0,0])
    plt.ylabel('Amplitude')
    plt.show()
    plt.savefig(plotName,dpi=200)
    plt.close('all')

def plot_timeseries(time,data,traveltimes,gpsStart,gpsEnd,type,plotName):

    plotNameSplit = plotName.split("/")
    plotTitle = plotNameSplit[-1].replace(".png","")

    startTime = gpsStart
    endTime = gpsEnd

    time = np.array(time)
    data = np.array(data)

    minData = min(data)
    maxData = max(data)

    time = time - startTime  
    data = 1/(maxData - minData) * (data - minData)

    Ptime = traveltimes["Ptimes"][-1] - startTime
    Stime = traveltimes["Stimes"][-1] - startTime
    Rtime = traveltimes["Rtimes"][-1] - startTime

    if type == "timeseries":
        plt.plot(time,data, 'k')
    elif type == "kw":
        plot.plot(time,data, 'k*')
    plt.axvline(x=Ptime,color='r')
    plt.axvline(x=Stime,color='b')
    plt.axvline(x=Rtime,color='g')
    plt.text(Ptime, -0.05, 'P', fontsize=18, ha='center', va='top')
    plt.text(Stime, -0.05, 'S', fontsize=18, ha='center', va='top')
    plt.text(Rtime, -0.05, 'R', fontsize=18, ha='center', va='top')
    #xlim([min(time), max(time)])
    plt.xlim([0, endTime-startTime])
    plt.ylim([0, 1])
    plt.xlabel('Time [s]')
    plt.ylabel("%.0f - %.0f"%(startTime,endTime))
    plt.title(plotTitle)
    plt.show()
    plt.savefig(plotName,dpi=200)
    plt.close('all')

def latencies_sent(params,attributeDics,plotName):

    latencies = []
    for attributeDic in attributeDics:
        latencies.append(attributeDic["SentGPS"]-attributeDic["GPS"])

    if latencies == []:
        return

    figure()
    bins=np.logspace(1,5,15)
    hist(latencies, bins=bins, rwidth=1)
    gca().set_xscale("log")
    gca().set_xlim([10**1,10**5])
    xlabel('Latencies [s]')
    #title("Latencies Sent")
    show()
    savefig(plotName,dpi=200)
    close('all')

def latencies_written(params,attributeDics,plotName):

    latencies = []
    for attributeDic in attributeDics:
        latencies.append(attributeDic["WrittenGPS"]-attributeDic["SentGPS"])

    if latencies == []:
        return

    plt.figure()
    bins=np.linspace(0,100,100)
    plt.hist(latencies, bins=bins, rwidth=1)
    plt.gca().set_xscale("linear")
    plt.gca().set_xlim([0,25])
    plt.xlabel('Latencies [s]')
    #title("Latencies Written")
    plt.show()
    plt.savefig(plotName,dpi=200)
    plt.close('all')

def magnitudes(params,attributeDics,plotName):

    gps = []
    magnitudes = []
    for attributeDic in attributeDics:
        gps.append(attributeDic["GPS"])
        magnitudes.append(attributeDic["Magnitude"])

    if magnitudes == []:
        return

    gps = np.array(gps)
    magnitudes = np.array(magnitudes)

    startTime = min(gps)
    #endTime = max(gps)
    endTime = params["gpsEnd"]
    gps = (gps - endTime)/(60)

    plt.figure()
    plt.plot(gps,magnitudes, 'k*')
    plt.ylim([min(magnitudes)-0.5, max(magnitudes)+0.5])
    plt.xlabel('Time [Minutes]')
    plt.ylabel("%.0f - %.0f"%(startTime,endTime))
    plt.xlim([-60, 0])
    plt.title("Magnitude vs. Time")
    plt.show()
    plt.savefig(plotName,dpi=200)
    plt.close('all')

def magnitudes_latencies(params,attributeDics,plotName):

    latencies = []
    magnitudes = []
    for attributeDic in attributeDics:
        latencies.append(attributeDic["SentGPS"]-attributeDic["GPS"])
        magnitudes.append(attributeDic["Magnitude"])

    if latencies == []:
        return

    figure()
    plot(latencies, magnitudes, '*')
    gca().set_xscale("log")
    gca().set_xlim([10**1,10**5])
    xlabel('Latencies [s]')
    #ylabel('Magnitudes')
    #title("Magnitudes vs. Latencies Sent")
    show()
    savefig(plotName,dpi=200)
    close('all')

def variable_magnitudes(attributeDicsDiff,variable,plotName):

    variables = []
    magnitudes = []
    for attributeDicDiff in attributeDicsDiff:
        variables.append(attributeDicDiff[variable])
        magnitudes.append(attributeDicDiff["attributeDic2"]["Magnitude"])

    if variables == []:
        return

    figure()
    plot(magnitudes, variables, '*')
    xlabel('Magnitudes')
    #ylabel(variable)

    xlim([min(magnitudes)-0.5,max(magnitudes)+0.5])
    ylim([min(variables)-0.5,max(variables)+0.5])

    #title("%s vs. Magnitudes"%variable)
    show()
    savefig(plotName,dpi=200)
    close('all')

def traveltimes(params,attributeDics,ifo,currentGPS,plotName):

    traveltimes = []

    for attributeDic in attributeDics:

        if not ifo in attributeDic["traveltimes"]:
            continue

        if traveltimes == []:
            arrivalTimes = [max(attributeDic["traveltimes"][ifo]["Rtimes"]),max(attributeDic["traveltimes"][ifo]["Stimes"]),max(attributeDic["traveltimes"][ifo]["Ptimes"])]
            traveltimes = np.array(arrivalTimes)
        else:
            arrivalTimes = [max(attributeDic["traveltimes"][ifo]["Rtimes"]),max(attributeDic["traveltimes"][ifo]["Stimes"]),max(attributeDic["traveltimes"][ifo]["Ptimes"])]            
            traveltimes = np.vstack([traveltimes,arrivalTimes])

    if traveltimes == []:
        return

    traveltimes = currentGPS-traveltimes

    plotNameSplit = plotName.split("/")
    plotTitle = plotNameSplit[-1].replace(".png","").replace("traveltimes","")

    startTime = np.min(traveltimes)
    endTime = np.max(traveltimes)

    ax = subplot(1,1,1)

    if len(traveltimes.shape) == 1:
        plot(traveltimes[0],1.5 * ones(1), 'r*', label='P', markersize=10.0)
        plot(traveltimes[1],2.0 * ones(1), 'g*', label='S', markersize=10.0)
        plot(traveltimes[2],2.5 * ones(1), 'b*', label='R', markersize=10.0)
        vlines(traveltimes[0], 0, 1.5)
        vlines(traveltimes[1], 0, 2.0)
        vlines(traveltimes[2], 0, 2.5)
    else:
        plot(traveltimes[:,0],1.5 * ones(len(traveltimes[:,0])), 'r*', label='P', markersize=10.0)
        plot(traveltimes[:,1],2.0 * ones(len(traveltimes[:,1])), 'g*', label='S', markersize=10.0)
        plot(traveltimes[:,2],2.5 * ones(len(traveltimes[:,2])), 'b*', label='R', markersize=10.0)
        vlines(traveltimes[:,0], 0, 1.5)
        vlines(traveltimes[:,1], 0, 2.0)
        vlines(traveltimes[:,2], 0, 2.5)
    xlim([startTime-1000, 1000])
    ylim([1, 3])
    xlabel('Countdown [s]')
    title(plotTitle)

    handles, labels = ax.get_legend_handles_labels()
    legend(handles[0:3], labels[0:3])
    show()
    savefig(plotName,dpi=200)
    close('all')

def find_nearest(array,value):
    array = np.array(array)
    index=(np.abs(array-value)).argmin()
    return array[index], index

def worldmap_plot(params,attributeDics,type,plotName):

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

    plt.figure(figsize=(10,5))
    plt.axes([0,0,1,1])

    # lon_0 is central longitude of robinson projection.
    # resolution = 'c' means use crude resolution coastlines.
    m = Basemap(projection='robin',lon_0=0,resolution='c')
    #set a background colour
    m.drawmapboundary(fill_color='#85A6D9')

    # draw coastlines, country boundaries, fill continents.
    m.fillcontinents(color='white',lake_color='#85A6D9')
    m.drawcoastlines(color='#6D5F47', linewidth=.4)
    m.drawcountries(color='#6D5F47', linewidth=.4)

    # draw lat/lon grid lines every 30 degrees.
    m.drawmeridians(np.arange(-180, 180, 30), color='#bbbbbb')
    m.drawparallels(np.arange(-90, 90, 30), color='#bbbbbb')

    if not attributeDics == []:
        traveltimes = attributeDics[0]["traveltimes"][ifo]
        ifolat = traveltimes["Latitudes"][-1]
        ifolng = traveltimes["Longitudes"][-1]
        # compute the native map projection coordinates for cities
        ifox,ifoy = m(ifolng,ifolat)

        m.scatter(
            ifox,
            ifoy,
            s=10, #size
            c='black', #color
            marker='x', #symbol
            alpha=0.5, #transparency
            zorder = 1, #plotting order
            )
        text(
            ifox+50000,
            ifoy+50000,
            ifo,
            color = 'black',
            size='small',
            horizontalalignment='center',
            verticalalignment='center',
            zorder = 2,
            )

    for attributeDic in attributeDics:

        if type == "Restimates" and "Restimate" not in attributeDic["traveltimes"][ifo]:
            continue

        x,y = m(attributeDic["Longitude"], attributeDic["Latitude"])
        if type == "Magnitude":
            color = attributeDic["Magnitude"]
            colorbar_label = "Magnitude"
            vmin = 0
            vmax = 7
        elif type == "Traveltimes":
            travel_time = attributeDic["traveltimes"][ifo]["Rtimes"][-1] - attributeDic["traveltimes"][ifo]["Rtimes"][0]
            travel_time = travel_time / 60
            color = travel_time
            colorbar_label = "Travel times [minutes]"
            vmin = 0
            vmax = 60
        elif type == "Restimates":
            
            travel_time = attributeDic["traveltimes"][ifo]["Restimate"] - attributeDic["traveltimes"][ifo]["Rtimes"][0]
            travel_time = travel_time / 60
            color = travel_time
            colorbar_label = "Travel times [minutes]"
            vmin = 0
            vmax = 60
        m.scatter(
                x,
                y,
                s=10, #size
                marker='o', #symbol
                alpha=0.5, #transparency
                zorder = 3, #plotting order
                c=color, 
                vmin=vmin, 
                vmax=vmax
        )

    try:
       cbar=plt.colorbar()
       cbar.set_label(colorbar_label)
       cbar.set_clim(vmin=vmin,vmax=vmax)
    except:
       pass
    plt.show()
    plt.savefig(plotName,dpi=200)
    plt.close('all')

def worldmap_wavefronts(params,attributeDics,currentGPS,plotName):

    figure(figsize=(10,5))
    axes([0,0,1,1])

    # lon_0 is central longitude of robinson projection.
    # resolution = 'c' means use crude resolution coastlines.
    m = Basemap(projection='robin',lon_0=0,resolution='c')
    #set a background colour
    m.drawmapboundary(fill_color='#85A6D9')

    # draw coastlines, country boundaries, fill continents.
    m.fillcontinents(color='white',lake_color='#85A6D9')
    m.drawcoastlines(color='#6D5F47', linewidth=.4)
    m.drawcountries(color='#6D5F47', linewidth=.4)

    # draw lat/lon grid lines every 30 degrees.
    m.drawmeridians(np.arange(-180, 180, 30), color='#bbbbbb')
    m.drawparallels(np.arange(-90, 90, 30), color='#bbbbbb')

    if not attributeDics == []:
     for ifoName, traveltimes in attributeDics[0]["traveltimes"].items():
        ifolat = traveltimes["Latitudes"][-1]
        ifolng = traveltimes["Longitudes"][-1]
        # compute the native map projection coordinates for cities
        ifox,ifoy = m(ifolng,ifolat)

        m.scatter(
            ifox,
            ifoy,
            s=10, #size
            c='blue', #color
            marker='o', #symbol
            alpha=0.5, #transparency
            zorder = 2, #plotting order
            )
        text(
            ifox+50000,
            ifoy+50000,
            ifoName,
            color = 'black',
            size='small',
            horizontalalignment='center',
            verticalalignment='center',
            zorder = 3,
            )


    for attributeDic in attributeDics:

        if attributeDic["traveltimes"] == {}:
            continue

        ifoNames = []
        ifoDist = []
        for ifoName, traveltimes in attributeDic["traveltimes"].items():
            ifoDist.append(traveltimes["Distances"][-1])
            ifoNames.append(ifoName)
        ifoIndex = np.array(ifoDist).argmax()
        ifo = ifoNames[ifoIndex]

        nearestGPS, Pindex = find_nearest(attributeDic["traveltimes"][ifo]["Ptimes"],currentGPS)
        Pdist = attributeDic["traveltimes"][ifo]["Distances"][Pindex]/1000
        nearestGPS, Sindex = find_nearest(attributeDic["traveltimes"][ifo]["Stimes"],currentGPS)
        Sdist = attributeDic["traveltimes"][ifo]["Distances"][Sindex]/1000
        nearestGPS, Rindex = find_nearest(attributeDic["traveltimes"][ifo]["Ptimes"],currentGPS)
        Rdist = attributeDic["traveltimes"][ifo]["Distances"][Rindex]/1000

        if currentGPS > max([attributeDic["traveltimes"][ifo]["Ptimes"][-1],attributeDic["traveltimes"][ifo]["Stimes"][-1],attributeDic["traveltimes"][ifo]["Rtimes"][-1]]):
            continue

        x,y = m(attributeDic["Longitude"], attributeDic["Latitude"])
        m.scatter(
                x,
                y,
                s=10, #size
                marker='o', #symbol
                alpha=0.5, #transparency
                zorder = 2, #plotting order
        )
        text(
                x,
                y,
                attributeDic["eventName"],
                color = 'black',
                size='small',
                horizontalalignment='center',
                verticalalignment='center',
                zorder = 3,
        )

        X,Y = pylal.pylal_seismon_eqmon.equi(m, attributeDic["Longitude"], attributeDic["Latitude"], Pdist)
        m.plot(
                X,
                Y,
                linewidth = attributeDic["Magnitude"] / 2,
                zorder = 3, #plotting order
                color = 'b'
        )
        X,Y = pylal.pylal_seismon_eqmon.equi(m, attributeDic["Longitude"], attributeDic["Latitude"], Sdist)
        m.plot(
                X,
                Y,
                linewidth = attributeDic["Magnitude"] / 2,
                zorder = 3, #plotting order
                color = 'r'
        )
        X,Y = pylal.pylal_seismon_eqmon.equi(m, attributeDic["Longitude"], attributeDic["Latitude"], Rdist)
        m.plot(
                X,
                Y,
                linewidth = attributeDic["Magnitude"] / 2,
                zorder = 3, #plotting order
                color = 'y'
        )

    show()
    savefig(plotName,dpi=200)
    #savefig(plotNameCounter,dpi=200)
    close('all')

