#!/usr/bin/python

import os, glob, optparse, shutil, warnings, pickle, math, copy, pickle, matplotlib
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal, scipy.stats
from lxml import etree
import pylal.pylal_seismon_NLNM, pylal.pylal_seismon_html
import pylal.pylal_seismon_eqmon, pylal.pylal_seismon_utils

import gwpy.time, gwpy.timeseries
import gwpy.spectrum, gwpy.spectrogram
import gwpy.plotter

__author__ = "Michael Coughlin <michael.coughlin@ligo.org>"
__date__ = "2012/8/26"
__version__ = "0.1"

# =============================================================================
#
#                               DEFINITIONS
#
# =============================================================================

def save_data(params,channel,gpsStart,gpsEnd,data):
    """@saves spectral data in text files.

    @param params
        seismon params dictionary
    @param channel
        seismon channel structure
    @param gpsStart
        start gps
    @param gpsStart
        start gps
    @param gpsEnd
        end gps
    @param data
        spectral data structure 
    """

    psdDirectory = params["dirPath"] + "/Text_Files/PSD/" + channel.station_underscore + "/" + str(params["fftDuration"])
    pylal.pylal_seismon_utils.mkdir(psdDirectory)

    fftDirectory = params["dirPath"] + "/Text_Files/FFT/" + channel.station_underscore + "/" + str(params["fftDuration"])
    pylal.pylal_seismon_utils.mkdir(fftDirectory)

    timeseriesDirectory = params["dirPath"] + "/Text_Files/Timeseries/" + channel.station_underscore + "/" + str(params["fftDuration"])
    pylal.pylal_seismon_utils.mkdir(timeseriesDirectory)

    freq = np.array(data["dataASD"].frequencies)

    psdFile = os.path.join(psdDirectory,"%d-%d.txt"%(gpsStart,gpsEnd))
    f = open(psdFile,"wb")
    for i in xrange(len(freq)):
        f.write("%e %e\n"%(freq[i],data["dataASD"][i]))
    f.close()

    freq = np.array(data["dataFFT"].frequencies)

    fftFile = os.path.join(fftDirectory,"%d-%d.txt"%(gpsStart,gpsEnd))
    f = open(fftFile,"wb")
    for i in xrange(len(freq)):
        f.write("%e %e %e\n"%(freq[i],data["dataFFT"].data[i].real,data["dataFFT"].data[i].imag))
    f.close()

    timeseriesFile = os.path.join(timeseriesDirectory,"%d-%d.txt"%(gpsStart,gpsEnd))
    f = open(timeseriesFile,"wb")
    f.write("%e %e %e\n"%(np.min(data["dataFull"].data),np.median(data["dataFull"].data),np.max(data["dataFull"].data)))
    f.close()

def calculate_spectra(params,channel,dataFull):
    """@calculate spectral data

    @param params
        seismon params dictionary
    @param channel
        seismon channel structure
    @param dataFull
        timeseries data structure
    """

    fs = channel.samplef # 1 ns -> 1 GHz
    cutoff_high = 0.5 # 10 MHz
    cutoff_low = 0.3
    n = 3
    worN = 16384
    B_low, A_low = scipy.signal.butter(n, cutoff_low / (fs / 2.0), btype='lowpass')
    #w_low, h_low = scipy.signal.freqz(B_low,A_low)
    w_low, h_low = scipy.signal.freqz(B_low,A_low,worN=worN)
    B_high, A_high = scipy.signal.butter(n, cutoff_high / (fs / 2.0), btype='highpass') 
    w_high, h_high = scipy.signal.freqz(B_high,A_high,worN=worN)

    w = w_high * (fs / (2.0*np.pi))

    if params["doPlots"]:

        plotDirectory = params["path"] + "/" + channel.station_underscore
        pylal.pylal_seismon_utils.mkdir(plotDirectory)

        pngFile = os.path.join(plotDirectory,"bode.png")
        kwargs = {'logx':True}
        plot = gwpy.plotter.BodePlot(figsize=[14,8],**kwargs)
        plot.add_filter((B_low,A_low),frequencies=w_high,sample_rate=fs,label="lowpass")
        plot.add_filter((B_high,A_high),frequencies=w_high,sample_rate=fs,label="highpass")
        plot.add_legend(loc=1,prop={'size':10})
        plot.save(pngFile,dpi=200)
        plot.close()

    #dataLowpass = dataFull.lowpass(1.0)
    dataLowpass = scipy.signal.lfilter(B_low, A_low, dataFull,
                                       axis=0).view(dataFull.__class__)
    dataLowpass.data[:2*channel.samplef] = dataLowpass.data[2*channel.samplef]
    dataLowpass.data[-2*channel.samplef:] = dataLowpass.data[-2*channel.samplef]
    dataLowpass.sample_rate =  dataFull.sample_rate
    dataLowpass.epoch = dataFull.epoch

    #dataHighpass = dataFull.highpass(1.0)
    dataHighpass = scipy.signal.lfilter(B_high, A_high, dataFull,
                                        axis=0).view(dataFull.__class__)
    dataHighpass.sample_rate =  dataFull.sample_rate
    dataHighpass.epoch = dataFull.epoch

    # calculate spectrum
    NFFT = params["fftDuration"]
    #window = None
    #dataASD = dataFull.asd(NFFT,NFFT,'welch')
    dataPSD = dataFull.psd(NFFT,NFFT,'welch')
    dataASD = dataPSD**(1/2.)
    #dataFFT = np.fft.fft(dataFull.data)
    dataFFT = np.fft.fft(dataFull.data,params["fftDuration"]*channel.samplef)
    #freqFFT = np.fft.fftfreq(dataFull.data.shape[-1],d=1.0/channel.samplef)
    freqFFT = np.fft.fftfreq(int(params["fftDuration"]*channel.samplef),d=1.0/channel.samplef)
    freq = np.array(dataASD.frequencies)

    indexes = np.where((freq >= params["fmin"]) & (freq <= params["fmax"]))[0]
    dataASD = np.array(dataASD.data)

    freq = freq[indexes]
    dataASD = dataASD[indexes]
    dataASD = gwpy.spectrum.Spectrum(dataASD, f0=np.min(freq), df=(freq[1]-freq[0]))

    indexes = np.where((freqFFT >= params["fmin"]) & (freqFFT <= params["fmax"]))[0]
    freqFFT = freqFFT[indexes]
    dataFFT = dataFFT[indexes]
    dataFFT = gwpy.spectrum.Spectrum(dataFFT, f0=np.min(freqFFT), df=(freqFFT[1]-freqFFT[0]))

    # manually set units (units in CIS aren't correct)
    dataASD.unit = 'counts/Hz^(1/2)'
    dataFFT.unit = 'counts/Hz^(1/2)'

    data = {}
    data["dataFull"] = dataFull
    data["dataLowpass"] = dataLowpass
    data["dataHighpass"] = dataHighpass
    data["dataASD"] = dataASD
    data["dataFFT"] = dataFFT

    return data

def apply_calibration(params,channel,data):
    """@applies calibration to necessary channels

    @param params
        seismon params dictionary
    @param channel
        seismon channel structure
    @param data
        spectral data structure
    """  

    if ("L4C" in channel.station) or ("GS13" in channel.station):

        zeros = [0,0,0]
        poles = [0.70711 + 0.70711*1j , 0.70711 - 0.70711*1j]
        gain = 1

        b = [1,0,0,0];
        a = [0,1,-1.414,1];
        w, h = scipy.signal.freqz(b, a)

        f = data["dataASD"].frequencies.data
      
        # Divide by f to get to displacement
        data["dataASD"]/=f
        # Filter spectrum 
        data["dataASD"].filterba(a,b,inplace=True)
        fresp = abs(scipy.signal.freqs(a, b, f)[1])
        # Multiply by f to get to velocity
        data["dataASD"]*=f

        if params["doPlots"]:

            plotDirectory = params["path"] + "/" + channel.station_underscore
            pylal.pylal_seismon_utils.mkdir(plotDirectory)

            pngFile = os.path.join(plotDirectory,"calibration.png")
            plot = gwpy.plotter.Plot(figsize=[14,8])
            plot.add_line(f,fresp)
            plot.xlim = [params["fmin"],params["fmax"]]
            plot.xlabel = "Frequency [Hz]"
            plot.ylabel = "Response"
            plot.title = channel.station.replace("_","\_")
            plot.axes.set_xscale("log")
            plot.axes.set_yscale("log")
            plot.save(pngFile,dpi=200)
            plot.close()

    return data

def retrieve_timeseries(params,channel,segment):
    """@retrieves timeseries for given channel and segment.

    @param params
        seismon params dictionary
    @param channel
        seismon channel structure
    @param segment
        [start,end] gps
    """

    gpsStart = segment[0]
    gpsEnd = segment[1]

    # set the times
    duration = np.ceil(gpsEnd-gpsStart)

    dataFull = []
    if params["ifo"] == "IRIS":
        import obspy.iris
        client = obspy.iris.Client()
        tstart = pylal.pylal_seismon_utils.GPSToUTCDateTime(gpsStart)
        tend = pylal.pylal_seismon_utils.GPSToUTCDateTime(gpsEnd)

        channelSplit = channel.station.split(":")
        st = client.getWaveform(channelSplit[0], channelSplit[1], channelSplit[2], channelSplit[3],\
            tstart, tend)

        data = np.array(st[0].data)
        data = data.astype(float)

        dataFull = gwpy.timeseries.TimeSeries(data, times=None, epoch=gpsStart, channel=channel.station, unit=None,sample_rate=channel.samplef, name=channel.station)

    else:
        # make timeseries
        try:
            dataFull = gwpy.timeseries.TimeSeries.read(params["frame"], channel.station, epoch=gpsStart, duration=duration)
        except:
            print "data read from frames failed... continuing\n"
            return

    return dataFull

def spectra(params, channel, segment):
    """@calculates spectral data for given channel and segment.

    @param params
        seismon params dictionary
    @param channel
        seismon channel structure
    @param segment
        [start,end] gps
    """

    ifo = pylal.pylal_seismon_utils.getIfo(params)

    gpsStart = segment[0]
    gpsEnd = segment[1]
   
    # set the times
    duration = np.ceil(gpsEnd-gpsStart)

    # make timeseries
    dataFull = retrieve_timeseries(params, channel, segment)
    if dataFull == []:
        return 

    dataFull /= channel.calibration
    indexes = np.where(np.isnan(dataFull.data))[0]
    meanSamples = np.mean(np.ma.masked_array(dataFull.data,np.isnan(dataFull.data)))
    for index in indexes:
        dataFull[index] = meanSamples
    dataFull -= np.mean(dataFull.data)

    if np.mean(dataFull) == 0.0:
        print "data only zeroes... continuing\n"
        return

    data = calculate_spectra(params,channel,dataFull)
    data = apply_calibration(params,channel,data)

    save_data(params,channel,gpsStart,gpsEnd,data)

    if params["doEarthquakes"]:
        earthquakesDirectory = os.path.join(params["path"],"earthquakes")
        earthquakesXMLFile = os.path.join(earthquakesDirectory,"earthquakes.xml")
        attributeDics = pylal.pylal_seismon_utils.read_eqmons(earthquakesXMLFile)
    else:
        attributeDics = []

    if params["doPlots"]:

        plotDirectory = params["path"] + "/" + channel.station_underscore
        pylal.pylal_seismon_utils.mkdir(plotDirectory)        

        pngFile = os.path.join(plotDirectory,"timeseries.png")
        plot = gwpy.plotter.TimeSeriesPlot(figsize=[14,8])

        dataHighpass = data["dataHighpass"].resample(16)
        dataFull = data["dataFull"].resample(16)
        dataLowpass = data["dataLowpass"].resample(16)

        #dataHighpass = data["dataHighpass"]
        #dataFull = data["dataFull"]
        #dataLowpass = data["dataLowpass"]

        dataHighpass *= 1e6
        dataFull *= 1e6
        dataLowpass *= 1e6

        plot.add_timeseries(dataHighpass,label="highpass")
        plot.add_timeseries(dataFull,label="data")
        plot.add_timeseries(dataLowpass,label="lowpass")

        xlim = [plot.xlim[0],plot.xlim[1]]
        ylim = [plot.ylim[0],plot.ylim[1]]

        for attributeDic in attributeDics:

            traveltimes = attributeDic["traveltimes"][ifo]

            Ptime = max(traveltimes["Ptimes"])
            Stime = max(traveltimes["Stimes"])
            Rtwotime = max(traveltimes["Rtwotimes"])
            RthreePointFivetime = max(traveltimes["RthreePointFivetimes"])
            Rfivetime = max(traveltimes["Rfivetimes"])

            peak_velocity = traveltimes["Rfamp"][0]
            peak_velocity = peak_velocity * 1e6

            if peak_velocity > ylim[1]:
               ylim[1] = peak_velocity*1.1
            if -peak_velocity < ylim[0]:
               ylim[0] = -peak_velocity*1.1

            kwargs = {"linestyle":"--","color":"r"}
            plot.add_line([Ptime,Ptime],ylim,label="P Est. Arrival",**kwargs)
            kwargs = {"linestyle":"--","color":"g"}
            plot.add_line([Stime,Stime],ylim,label="S Est. Arrival",**kwargs)
            kwargs = {"linestyle":"--","color":"b"}
            plot.add_line([RthreePointFivetime,RthreePointFivetime],ylim,label="Middle R Est. Arrival",**kwargs)            
            kwargs = {"linestyle":"--","color":"b"}
            plot.add_line([Rtwotime,Rtwotime],ylim,label="High R Est. Arrival",**kwargs)
            kwargs = {"linestyle":"--","color":"b"}
            plot.add_line([Rfivetime,Rfivetime],ylim,label="Low R Est. Arrival",**kwargs)
            kwargs = {"linestyle":"--","color":"k"}
            plot.add_line(xlim,[peak_velocity,peak_velocity],label="pred. vel.",**kwargs)
            plot.add_line(xlim,[-peak_velocity,-peak_velocity],**kwargs)

        plot.ylabel = r"Velocity [$\mu$m/s]"
        plot.title = channel.station.replace("_","\_")
        plot.xlim = xlim
        plot.ylim = ylim
        plot.add_legend(loc=1,prop={'size':10})

        plot.save(pngFile)
        plot.close()

        fl, low, fh, high = pylal.pylal_seismon_NLNM.NLNM(2)

        pngFile = os.path.join(plotDirectory,"psd.png")
        label = channel.station.replace("_","\_")

        plot = gwpy.plotter.Plot(figsize=[14,8])
        plot.add_spectrum(data["dataASD"],label=label)
        kwargs = {"linestyle":"-.","color":"k"}
        plot.add_line(fl, low, label="HNM/LNM", **kwargs)
        plot.add_line(fh, high, **kwargs)
        plot.xlim = [params["fmin"],params["fmax"]]
        plot.ylim = [10**-10, 10**-5]
        plot.xlabel = "Frequency [Hz]"
        plot.ylabel = "Amplitude Spectrum [(m/s)/rtHz]"
        plot.title = channel.station.replace("_","\_")
        plot.axes.set_xscale("log")
        plot.axes.set_yscale("log")
        plot.add_legend(loc=1,prop={'size':10})

        plot.save(pngFile,dpi=200)
        plot.close()

def freq_analysis(params,channel,tt,freq,spectra):
    """@frequency analysis of spectral data.

    @param params
        seismon params dictionary
    @param channel
        seismon channel structure
    @param tt
        array of start times
    @param freq
        frequency vector
    @param spectra
        spectrogram structure
    """

    if params["doPlots"]:
        plotDirectory = params["path"] + "/" + channel.station_underscore + "/freq"
        pylal_seismon_utils.mkdir(plotDirectory)
 
    indexes = np.logspace(0,np.log10(len(freq)-1),num=100)
    indexes = list(np.unique(np.ceil(indexes)))
    indexes = range(len(freq))
    #indexes = range(16)

    indexes = np.where(10.0 >= freq)[0]

    deltaT = tt[1] - tt[0]

    n_dist = []
    for j in xrange(1000):
        n_dist.append(scipy.stats.chi2.rvs(2))

    p_chi2_vals = []
    p_ks_vals = []
    ttCoh_vals = []

    for i in indexes:
        vals = spectra[:,i]

        meanPSD = np.mean(vals) 
        stdPSD = np.std(vals)

        vals_norm = 2 * vals / meanPSD

        bins = np.arange(0,10,1)
        (n,bins) = np.histogram(vals_norm,bins=bins)
        n_total = np.sum(n)

        bins = (bins[1:] + bins[:len(bins)-1])/2

        n_expected = []
        for bin in bins:
            expected_val = n_total * scipy.stats.chi2.pdf(bin, 2)
            n_expected.append(expected_val)
        n_expected = np.array(n_expected)

        (stat_chi2,p_chi2) = scipy.stats.mstats.chisquare(n, f_exp=n_expected)
        p_chi2_vals.append(p_chi2)

        (stat_ks,p_ks) = scipy.stats.ks_2samp(vals_norm, n_dist)
        p_ks_vals.append(p_ks)

        acov = np.correlate(vals,vals,"full")
        acov = acov / np.max(acov)

        ttCov = (np.arange(len(acov)) - len(acov)/2) * float(deltaT)

        #ttLimitMin = - 5/freq[i]
        #ttLimitMax = 5 /freq[i]

        ttLimitMin = - float('inf')
        ttLimitMax = float('inf')

        ttIndexes = np.intersect1d(np.where(ttCov >= ttLimitMin)[0],np.where(ttCov <= ttLimitMax)[0])
        #ttCov = ttCov / (60)

        acov_minus_05 = np.absolute(acov[ttIndexes] - 0.66)
        index_min = acov_minus_05.argmin()

        ttCoh = np.absolute(ttCov[ttIndexes[index_min]])
        ttCoh_vals.append(ttCoh)

        print freq[i], ttCoh, len(ttIndexes)

        if len(ttIndexes) == 0:
            continue

        #if freq[i] > 0:
        #    continue

        if params["doPlots"]:
            ax = plt.subplot(111)
            plt.plot(bins,n,label='true')
            plt.plot(bins,n_expected,'k*',label='expected')
            plt.xlabel("2 * data / mean")
            plt.ylabel("Counts")
            plot_title = "p-value: %f"%p_chi2
            plt.title(plot_title)
            plt.legend()
            plt.show()
            plt.savefig(os.path.join(plotDirectory,"%s_dist.png"%str(freq[i])),dpi=200)
            plt.savefig(os.path.join(plotDirectory,"%s_dist.eps"%str(freq[i])),dpi=200)
            plt.close('all')

            ax = plt.subplot(111)
            plt.semilogy(ttCov[ttIndexes],acov[ttIndexes])
            plt.vlines(ttCoh,10**(-3),1,color='r')
            plt.vlines(-ttCoh,10**(-3),1,color='r')
            plt.xlabel("Time [Seconds]")
            plt.ylabel("Correlation")
            plt.show()
            plt.savefig(os.path.join(plotDirectory,"%s_cov.png"%str(freq[i])),dpi=200)
            plt.savefig(os.path.join(plotDirectory,"%s_cov.eps"%str(freq[i])),dpi=200)
            plt.close('all')

    if params["doPlots"]:
        ax = plt.subplot(111)
        plt.loglog(freq[indexes],p_chi2_vals,label='chi2')
        plt.loglog(freq[indexes],p_ks_vals,label='k-s')
        plt.xlabel("Frequency [Hz]")
        plt.ylabel("p-value")
        plt.legend(loc=3)
        plt.show()
        plt.savefig(os.path.join(plotDirectory,"freq_analysis.png"),dpi=200)
        plt.savefig(os.path.join(plotDirectory,"freq_analysis.eps"),dpi=200)
        plt.close('all')      

        ax = plt.subplot(111)
        plt.semilogx(freq[indexes],ttCoh_vals)
        plt.xlabel("Frequency [Hz]")
        plt.ylabel("Coherence Time [s]")
        plt.show()
        plt.savefig(os.path.join(plotDirectory,"ttCohs.png"),dpi=200)
        plt.savefig(os.path.join(plotDirectory,"ttCohs.eps"),dpi=200)
        plt.close('all')

def analysis(params, channel):
    """@analysis of spectral data.

    @param params
        seismon params dictionary
    @param channel
        seismon channel structure
    """

    psdDirectory = params["dirPath"] + "/Text_Files/PSD/" + channel.station_underscore + "/" + str(params["fftDuration"])
    files = glob.glob(os.path.join(psdDirectory,"*.txt"))

    files = sorted(files)

    if not params["doFreqAnalysis"]:
        if len(files) > 1000:
            files = files[-1000:]

    #files = files[-10:]

    tts = []
    spectra = []

    for file in files:

        fileSplit = file.split("/")
        txtFile = fileSplit[-1].replace(".txt","")
        txtFileSplit = txtFile.split("-")
        thisTTStart = int(txtFileSplit[0])
        thisTTEnd = int(txtFileSplit[1])
        tt = thisTTStart
        tts.append(tt)

        spectra_out = gwpy.spectrum.Spectrum.read(file)
        spectra_out.unit = 'counts/Hz^(1/2)'
        spectra.append(spectra_out)

        if tt == params["gpsStart"]:
            spectraNow = spectra_out

    if not 'spectraNow' in locals():
        print "no data at requested time... continuing\n"
        return

    if np.mean(spectraNow.data) == 0.0:
        print "data only zeroes... continuing\n"
        return

    dt = tts[-1] - tts[-2]
    epoch = gwpy.time.Time(tts[0], format='gps')
    specgram = gwpy.spectrogram.Spectrogram.from_spectra(*spectra, dt=dt,epoch=epoch)
    
    freq = np.array(specgram.frequencies)

    bins,specvar = pylal.pylal_seismon_utils.spectral_histogram(specgram)

    if params["doFreqAnalysis"]:
        freq_analysis(params,channel,ttStart,freq,specgram)

    # Calculate percentiles
    spectral_variation_1per = pylal.pylal_seismon_utils.spectral_percentiles(specvar,bins,1)
    spectral_variation_10per = pylal.pylal_seismon_utils.spectral_percentiles(specvar,bins,10)
    spectral_variation_50per = pylal.pylal_seismon_utils.spectral_percentiles(specvar,bins,50)
    spectral_variation_90per = pylal.pylal_seismon_utils.spectral_percentiles(specvar,bins,90)
    spectral_variation_99per = pylal.pylal_seismon_utils.spectral_percentiles(specvar,bins,99)

    textDirectory = params["path"] + "/" + channel.station_underscore
    pylal.pylal_seismon_utils.mkdir(textDirectory)

    f = open(os.path.join(textDirectory,"spectra.txt"),"w")
    for i in xrange(len(freq)):
        f.write("%e %e %e %e %e %e %e\n"%(freq[i],spectral_variation_1per[i],spectral_variation_10per[i],spectral_variation_50per[i],spectral_variation_90per[i],spectral_variation_99per[i],spectraNow[i]))
    f.close()

    sigDict = {}
    # Break up entire frequency band into 6 segments
    ff_ave = [1/float(128), 1/float(64),  0.1, 1, 3, 5, 10]

    f = open(os.path.join(textDirectory,"sig.txt"),"w")
    for i in xrange(len(ff_ave)-1):
        newSpectra = []
        newSpectraNow = []
        newFreq = []

        for j in xrange(len(freq)):
            if ff_ave[i] <= freq[j] and freq[j] <= ff_ave[i+1]:
                newFreq.append(freq[j])
                newSpectraNow.append(spectraNow.data[j])
                if newSpectra == []:
                    newSpectra = specgram.data[:,j]
                else:                 
                    newSpectra = np.vstack([newSpectra,specgram.data[:,j]])

        newSpectra = np.array(newSpectra)
        if len(newSpectra.shape) > 1:
            newSpectra = np.mean(newSpectra, axis = 0)
        sig, bgcolor = pylal.pylal_seismon_utils.html_bgcolor(np.mean(newSpectraNow),newSpectra)

        f.write("%e %e %e %e %s\n"%(ff_ave[i],ff_ave[i+1],np.mean(newSpectraNow),sig,bgcolor))

        key = "%s-%s"%(ff_ave[i],ff_ave[i+1])

        dt = tts[-1] - tts[-2]
        epoch = gwpy.time.Time(tts[0], format='gps')

        timeseries = gwpy.timeseries.TimeSeries(newSpectra, epoch=epoch, sample_rate=1.0/dt)
        
        sigDict[key] = {}
        #timeseries.data = np.log10(timeseries.data) 
        sigDict[key]["data"] = timeseries

    f.close()

    if params["doPlots"]:

        plotDirectory = params["path"] + "/" + channel.station_underscore
        pylal.pylal_seismon_utils.mkdir(plotDirectory)

        fl, low, fh, high = pylal.pylal_seismon_NLNM.NLNM(2)

        pngFile = os.path.join(plotDirectory,"psd.png")

        plot = spectraNow.plot(auto_refresh=True)
        kwargs = {"linestyle":"-","color":"k"}
        plot.add_line(freq, spectral_variation_10per, **kwargs)
        plot.add_line(freq, spectral_variation_50per, **kwargs)
        plot.add_line(freq, spectral_variation_90per, **kwargs)
        kwargs = {"linestyle":"-.","color":"k"}
        plot.add_line(fl, low, **kwargs)
        plot.add_line(fh, high, **kwargs)
        plot.xlim = [params["fmin"],params["fmax"]]
        plot.ylim = [np.min(bins), np.max(bins)]
        plot.xlabel = "Frequency [Hz]"
        plot.ylabel = "Amplitude Spectrum [(m/s)/rtHz]"
        plot.save(pngFile,dpi=200)
        plot.close()

        pngFile = os.path.join(plotDirectory,"disp.png")

        spectraNowDisplacement = spectraNow / freq
        plot = spectraNowDisplacement.plot(auto_refresh=True)
        kwargs = {"linestyle":"-","color":"w"}
        plot.add_line(freq, spectral_variation_10per/freq, **kwargs)
        plot.add_line(freq, spectral_variation_50per/freq, **kwargs)
        plot.add_line(freq, spectral_variation_90per/freq, **kwargs)
        kwargs = {"linestyle":"-.","color":"k"}
        plot.add_line(fl, low/fl, **kwargs)
        plot.add_line(fh, high/fh, **kwargs)
        plot.xlim = [params["fmin"],params["fmax"]]
        plot.ylim = [np.min(bins)/np.max(freq), np.max(bins)/np.min(freq)]
        plot.xlabel = "Frequency [Hz]"
        plot.ylabel = "Displacement Spectrum [m/rtHz]"
        plot.save(pngFile,dpi=200)
        plot.close()

        pngFile = os.path.join(plotDirectory,"tf.png")
        specgramLog = specgram.to_logscale(fmin=np.min(freq),fmax=np.max(freq))
        plot = specgramLog.plot(auto_refresh=True)
        plot.ylim = [params["fmin"],params["fmax"]]
        plot.ylabel = "Frequency [Hz]"
        colorbar_label = "Amplitude Spectrum [(m/s)/rtHz]"
        kwargs = {}
        plot.add_colorbar(location='right', log=True, label=colorbar_label, clim=None, visible=True, **kwargs)
        plot.save(pngFile,dpi=200)
        plot.close()

        pngFile = os.path.join(plotDirectory,"psd-%d-%d.png"%(gpsStart,gpsEnd))
        plot = spectraNow.plot(auto_refresh=True)
        kwargs = {"linestyle":"-","color":"k"}
        plot.add_line(freq, spectral_variation_10per, **kwargs)
        plot.add_line(freq, spectral_variation_50per, **kwargs)
        plot.add_line(freq, spectral_variation_90per, **kwargs)
        kwargs = {"linestyle":"-.","color":"k"}
        plot.add_line(fl, low, **kwargs)
        plot.add_line(fh, high, **kwargs)
        plot.xlim = [params["fmin"],params["fmax"]]
        plot.ylim = [np.min(bins),np.max(bins)]
        plot.xlabel = "Frequency [Hz]"
        plot.ylabel = "Amplitude Spectrum [(m/s)/rtHz]"

        plot.save(pngFile,dpi=200)
        plot.close()

        pngFile = os.path.join(plotDirectory,"specvar.png")
        kwargs = {"linestyle":"-","color":"w"}
        plot = spectraNow.plot(**kwargs)
        kwargs = {"linestyle":"-","color":"k"}
        plot.add_line(freq, spectral_variation_10per, **kwargs)
        plot.add_line(freq, spectral_variation_50per, **kwargs)
        plot.add_line(freq, spectral_variation_90per, **kwargs)
        extent = [np.min(freq), np.max(freq),
                   np.min(bins), np.max(bins)]
        kwargs = {}
        plot.add_image(specvar, extent=extent, **kwargs)
        plot.axes.set_xscale("log")
        plot.axes.set_yscale("log")
        kwargs = {"linestyle":"-.","color":"k"}
        plot.add_line(fl, low, **kwargs)
        plot.add_line(fh, high, **kwargs)
        plot.xlim = [params["fmin"],params["fmax"]]
        plot.ylim = [np.min(bins), np.max(bins)]
        plot.xlabel = "Frequency [Hz]"
        plot.ylabel = "Amplitude Spectrum [(m/s)/rtHz]"

        plot.save(pngFile,dpi=200)
        plot.close()

        pngFile = os.path.join(plotDirectory,"bands.png")
        plot = gwpy.plotter.TimeSeriesPlot(auto_refresh=True)
        for key in sigDict.iterkeys():
            label = key
            plot.add_timeseries(sigDict[key]["data"], label=label)
        plot.axes.set_yscale("log")
        plot.ylabel = "Average Amplitude Spectrum log10[(m/s)/rtHz]"
        plot.add_legend(loc=1,prop={'size':10})
        plot.save(pngFile,dpi=200)
        plot.close()

    htmlPage = pylal.pylal_seismon_html.seismon_page(channel,textDirectory)
    if htmlPage is not None:
        f = open(os.path.join(textDirectory,"psd.html"),"w")
        f.write(htmlPage)
        f.close()

def channel_summary(params, segment):
    """@summary of channels of spectral data.

    @param params
        seismon params dictionary
    @param segment
        [start,end] gps
    """

    gpsStart = segment[0]
    gpsEnd = segment[1]

    data = {}
    for channel in params["channels"]:

        psdDirectory = params["dirPath"] + "/Text_Files/PSD/" + channel.station_underscore + "/" + str(params["fftDuration"])
        file = os.path.join(psdDirectory,"%d-%d.txt"%(gpsStart,gpsEnd))

        if not os.path.isfile(file):
            continue

        spectra_out = gwpy.spectrum.Spectrum.read(file)
        spectra_out.unit = 'counts/Hz^(1/2)'

        if np.sum(spectra_out.data) == 0.0:
            continue

        data[channel.station_underscore] = {}
        data[channel.station_underscore]["data"] = spectra_out

    if data == {}:
        return

    if params["doPlots"]:

        plotDirectory = params["path"] + "/summary"
        pylal.pylal_seismon_utils.mkdir(plotDirectory)

        fl, low, fh, high = pylal.pylal_seismon_NLNM.NLNM(2)

        pngFile = os.path.join(plotDirectory,"psd-%d-%d.png"%(gpsStart,gpsEnd))
        lowBin = np.inf
        highBin = -np.inf
        plot = gwpy.plotter.Plot(figsize=[14,8])
        for key in data.iterkeys():

            label = key.replace("_","\_")

            plot.add_spectrum(data[key]["data"], label=label)
            lowBin = np.min([lowBin,np.min(data[key]["data"])])
            highBin = np.max([highBin,np.max(data[key]["data"])])

        kwargs = {"linestyle":"-.","color":"k"}
        plot.add_line(fl, low, **kwargs)
        plot.add_line(fh, high, **kwargs)
        plot.xlim = [params["fmin"],params["fmax"]]
        plot.ylim = [lowBin, highBin]
        plot.xlabel = "Frequency [Hz]"
        plot.ylabel = "Amplitude Spectrum [(m/s)/rtHz]"
        plot.add_legend(loc=1,prop={'size':10})
        plot.axes.set_xscale("log")
        plot.axes.set_yscale("log")

        plot.save(pngFile,dpi=200)
        plot.close()

        pngFile = os.path.join(plotDirectory,"ratio.png")
        lowBin = np.inf
        highBin = -np.inf
        ref = params["referenceChannel"].replace(":","_")
        plot = gwpy.plotter.Plot(figsize=[14,8])
        for key in data.iterkeys():

            label = key.replace("_","\_")

            plot.add_spectrum(data[key]["data"] / data[ref]["data"], label=label)
            lowBin = np.min([lowBin,np.min(data[key]["data"])])
            highBin = np.max([highBin,np.max(data[key]["data"])])

        kwargs = {"linestyle":"-.","color":"k"}
        #plot.add_line(fl, low, **kwargs)
        #plot.add_line(fh, high, **kwargs)
        plot.xlim = [params["fmin"],params["fmax"]]
        #plot.ylim = [lowBin, highBin]
        plot.xlabel = "Frequency [Hz]"
        label_ref = params["referenceChannel"].replace("_","\_")
        plot.ylabel = "Spectrum / Reference [%s]"%(label_ref)
        plot.add_legend(loc=1,prop={'size':10})
        plot.axes.set_xscale("log")
        plot.axes.set_yscale("log")

        plot.save(pngFile,dpi=200)
        plot.close()


