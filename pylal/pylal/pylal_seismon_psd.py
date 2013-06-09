#!/usr/bin/python

import os, glob, optparse, shutil, warnings, matplotlib, pickle, math, copy, pickle
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal, scipy.stats
from collections import namedtuple
from pylal.xlal.datatypes.ligotimegps import LIGOTimeGPS
from pylal.xlal.date import XLALGPSToUTC
import pylal.seriesutils
from pylal import Fr
from pylal.dq import dqDataUtils
import pylal.pylal_seismon_NLNM, pylal.pylal_seismon_html

__author__ = "Michael Coughlin <michael.coughlin@ligo.org>"
__date__ = "2012/8/26"
__version__ = "0.1"

# =============================================================================
#
#                               DEFINITIONS
#
# =============================================================================

def read_frames(start_time,end_time,channel,cache):
    time = []
    data = []

    #== loop over frames in cache
    for frame in cache:
        frame_data,data_start,_,dt,_,_ = Fr.frgetvect1d(frame.path,channel.station)
        frame_length = float(dt)*len(frame_data)
        frame_time = data_start+dt*np.arange(len(frame_data))

        for i in range(len(frame_data)):
            if frame_time[i] <= start_time:  continue
            if frame_time[i] >= end_time:  continue
            time.append(frame_time[i])
            data.append(frame_data[i])
    data = [e/channel.calibration for e in data]

    return time,data

def mat(params, channel, segment):

    psdLocation = params["dirPath"] + "/Text_Files/PSD/" + channel.station_underscore
    if not os.path.isdir(psdLocation):
        os.makedirs(psdLocation)
    psdLocation = os.path.join(psdLocation,str(params["fftDuration"])) 
    if not os.path.isdir(psdLocation):
        os.makedirs(psdLocation)

    gpsStart = segment[0]
    gpsEnd = segment[1]

    time,data = read_frames(gpsStart,gpsEnd,channel,params["frame"])

    NFFT = params["fftDuration"]*channel.samplef
    spectra, freq = matplotlib.pyplot.psd(data, NFFT=NFFT, Fs=channel.samplef, Fc=0, detrend=matplotlib.mlab.detrend_mean,window=matplotlib.mlab.window_hanning)
    plt.close('all')

    spectra = [math.sqrt(e) for e in spectra]

    newSpectra = []
    newFreq = []

    for i in xrange(len(freq)):
        if freq[i] <= params["fmax"]:
            newFreq.append(freq[i])
            newSpectra.append(spectra[i])

    spectra = newSpectra
    freq = newFreq

    psdFile = os.path.join(psdLocation,"%d-%d.txt"%(gpsStart,gpsEnd))
    f = open(psdFile,"wb")
    for i in xrange(len(freq)):
        f.write("%e %e\n"%(freq[i],spectra[i]))
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

        startTime = np.min(time)
        startTimeUTC = XLALGPSToUTC(LIGOTimeGPS(int(startTime)))
        startTimeUTCString = "%d-%d-%d %d:%d:%d"%(startTimeUTC[0],startTimeUTC[1],startTimeUTC[2],startTimeUTC[3],startTimeUTC[4],startTimeUTC[5])

        time = time - startTime

        norm_pass = 1.0/(channel.samplef/2)
        norm_stop = 1.5*norm_pass
        #(N, Wn) = scipy.signal.buttord(wp=norm_pass, ws=norm_stop, gpass=2, gstop=30, analog=0)
        order = 5
        (b, a) = scipy.signal.butter(order, norm_pass, btype='low', analog=0, output='ba')
        data = np.array(data)
        dataLowpass = scipy.signal.filtfilt(b, a, data) 

        minData = np.min(data)
        maxData = np.max(data)
        data = 1/(maxData - minData) * (data - minData)

        minDataLowpass = np.min(dataLowpass)
        maxDataLowpass = np.max(dataLowpass)
        dataLowpass = 1/(maxDataLowpass - minDataLowpass) * (dataLowpass - minDataLowpass)

        plt.plot(time,data,'k',label='data')
        plt.plot(time,dataLowpass,'b',label='data lowpassed')
        plt.legend(loc=1,prop={'size':10})

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

                plt.text(Ptime, 1.1, 'P', fontsize=18, ha='center', va='top')
                plt.text(Stime, 1.1, 'S', fontsize=18, ha='center', va='top')
                plt.text(Rtime, 1.1, 'R', fontsize=18, ha='center', va='top')

                plt.axvline(x=Ptime,color='r',linewidth=2,zorder = 0,clip_on=False)
                plt.axvline(x=Stime,color='b',linewidth=2,zorder = 0,clip_on=False)
                plt.axvline(x=Rtime,color='g',linewidth=2,zorder = 0,clip_on=False)

        plt.xlabel("Time [s] [%s (%d)]"%(startTimeUTCString,startTime))
        plt.ylabel("Normalized Amplitude")
        plt.xlim([np.min(time),np.max(time)])

        plt.show()
        plt.savefig(os.path.join(plotLocation,"timeseries.png"),dpi=200)
        plt.savefig(os.path.join(plotLocation,"timeseries.eps"),dpi=200)
        plt.close('all')

        fl, low, fh, high = pylal.pylal_seismon_NLNM.NLNM(2)

        try:
            plt.semilogx(freq,spectra, 'k')
            plt.loglog(fl,low,'k-.',fh,high,'k-.')
        except:
            pass
        plt.xlim([params["fmin"],params["fmax"]])
        plt.ylim([10**-10, 10**-5])
        plt.xlabel("Frequency [Hz]")
        plt.ylabel("Seismic Spectrum [(m/s)/rtHz]")
        plt.grid()
        plt.show()
        plt.savefig(os.path.join(plotLocation,"psd.png"),dpi=200)
        plt.savefig(os.path.join(plotLocation,"psd.eps"),dpi=200)
        plt.close('all')

def freq_analysis(params,channel,tt,freq,spectra):

    if params["doPlots"]:
        plotLocation = params["path"] + "/" + channel.station_underscore
        if not os.path.isdir(plotLocation):
            os.makedirs(plotLocation)
        plotLocation = params["path"] + "/" + channel.station_underscore + "/freq" 
        if not os.path.isdir(plotLocation):
            os.makedirs(plotLocation)
 
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
            plt.savefig(os.path.join(plotLocation,"%s_dist.png"%str(freq[i])),dpi=200)
            plt.savefig(os.path.join(plotLocation,"%s_dist.eps"%str(freq[i])),dpi=200)
            plt.close('all')

            ax = plt.subplot(111)
            plt.semilogy(ttCov[ttIndexes],acov[ttIndexes])
            plt.vlines(ttCoh,10**(-3),1,color='r')
            plt.vlines(-ttCoh,10**(-3),1,color='r')
            plt.xlabel("Time [Seconds]")
            plt.ylabel("Correlation")
            plt.show()
            plt.savefig(os.path.join(plotLocation,"%s_cov.png"%str(freq[i])),dpi=200)
            plt.savefig(os.path.join(plotLocation,"%s_cov.eps"%str(freq[i])),dpi=200)
            plt.close('all')

    if params["doPlots"]:
        ax = plt.subplot(111)
        plt.loglog(freq[indexes],p_chi2_vals,label='chi2')
        plt.loglog(freq[indexes],p_ks_vals,label='k-s')
        plt.xlabel("Frequency [Hz]")
        plt.ylabel("p-value")
        plt.legend(loc=3)
        plt.show()
        plt.savefig(os.path.join(plotLocation,"freq_analysis.png"),dpi=200)
        plt.savefig(os.path.join(plotLocation,"freq_analysis.eps"),dpi=200)
        plt.close('all')      

        ax = plt.subplot(111)
        plt.semilogx(freq[indexes],ttCoh_vals)
        plt.xlabel("Frequency [Hz]")
        plt.ylabel("Coherence Time [s]")
        plt.show()
        plt.savefig(os.path.join(plotLocation,"ttCohs.png"),dpi=200)
        plt.savefig(os.path.join(plotLocation,"ttCohs.eps"),dpi=200)
        plt.close('all')

def spectral_histogram(data,bins):
    
    spectral_variation_norm = []
    rows, columns = data.shape

    for i in xrange(columns):
        this_spectral_variation, bin_edges = np.histogram(data[:,i],bins)
        this_spectral_variation = np.array(this_spectral_variation)
        weight = (100/float(sum(this_spectral_variation))) + np.zeros(this_spectral_variation.shape)
        if spectral_variation_norm == []:
            spectral_variation_norm = this_spectral_variation * weight
        else:
            spectral_variation_norm = np.vstack([spectral_variation_norm,this_spectral_variation * weight])
    return spectral_variation_norm

def calculate_percentiles(data,bins,percentile):

    cumsumvals = np.cumsum(data)
    
    abs_cumsumvals_minus_percentile = abs(cumsumvals - percentile)
    minindex = abs_cumsumvals_minus_percentile.argmin()
    val = bins[minindex]

    return val

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

def analysis(params, channel):

    psdLocation = params["dirPath"] + "/Text_Files/PSD/" + channel.station_underscore
    psdLocation = os.path.join(psdLocation,str(params["fftDuration"]))
    files = glob.glob(os.path.join(psdLocation,"*.txt"))

    files = sorted(files)

    if not params["doFreqAnalysis"]:
        if len(files) > 1000:
            files = files[-1000:]

    ttStart = []
    ttEnd = []
    freq = []
    spectra = []

    for file in files:

        fileSplit = file.split("/")
        txtFile = fileSplit[-1].replace(".txt","")
        txtFileSplit = txtFile.split("-")
        thisTTStart = int(txtFileSplit[0])
        thisTTEnd = int(txtFileSplit[1])

        ttStart.append(thisTTStart)
        ttEnd.append(thisTTEnd)

        data = np.loadtxt(file)
        thisSpectra = data[:,1]
        thisFreq = data[:,0]

        if freq == []:
            freq = np.array(thisFreq)
        if spectra == []:
            spectra = np.array(thisSpectra)
        else:
            spectra = np.vstack([spectra,np.array(copy.copy(thisSpectra))])
        if thisTTStart == params["gpsStart"]:
            freqNow = copy.copy(thisFreq)
            spectraNow = copy.copy(thisSpectra)

    # Binning parameters
    nb = 500
    range_binning = np.logspace(-10,-5,num=nb);

    try:
        spectra.max(axis = 1)
    except:
        print "Requires more than one spectra for analysis, continuing..."
        return

    spectra[spectra==float('Inf')] = 0

    # Calculate bin histogram of PSDs
    which_spectra = np.all([spectra.max(axis = 1) <= 10**-3,spectra.max(axis = 1) >= 10**-12],axis=0)

    if len(spectra[which_spectra]) == 0:
        return

    if params["doFreqAnalysis"]:
        freq_analysis(params,channel,ttStart,freq,spectra[which_spectra])

    spectral_variation_norm = spectral_histogram(spectra[which_spectra],range_binning)

    # Initialize arrays
    spectral_variation_norm_1per = []
    spectral_variation_norm_10per = []
    spectral_variation_norm_50per = []
    spectral_variation_norm_90per = []
    spectral_variation_norm_99per = []

    for i in xrange(len(freq)):
        spectral_variation_norm_1per.append(calculate_percentiles(spectral_variation_norm[i,:],range_binning,1))
        spectral_variation_norm_10per.append(calculate_percentiles(spectral_variation_norm[i,:],range_binning,10))
        spectral_variation_norm_50per.append(calculate_percentiles(spectral_variation_norm[i,:],range_binning,50))
        spectral_variation_norm_90per.append(calculate_percentiles(spectral_variation_norm[i,:],range_binning,90))
        spectral_variation_norm_99per.append(calculate_percentiles(spectral_variation_norm[i,:],range_binning,99))

    textLocation = params["path"] + "/" + channel.station_underscore
    if not os.path.isdir(textLocation):
        os.makedirs(textLocation)

    sigDict = []

    f = open(os.path.join(textLocation,"sig.txt"),"w")

    # Break up entire frequency band into 6 segments
    ff_ave = [1/float(128), 1/float(64),  0.1, 1, 3, 5, 10]

    for i in xrange(len(ff_ave)-1):
        newSpectra = []
        newSpectraNow = []
        newFreq = []

        for j in xrange(len(freq)):
            if ff_ave[i] <= freq[j] and freq[j] <= ff_ave[i+1]:
                newFreq.append(freq[j])
                newSpectraNow.append(spectraNow[j])
                if newSpectra == []:
                    newSpectra = spectra[:,j]
                else:                 
                    newSpectra = np.vstack([newSpectra,spectra[:,j]])

        if len(newSpectra.shape) > 1:
            sig, bgcolor = html_bgcolor(np.mean(newSpectraNow),np.mean(newSpectra, axis = 0))
        else:
            sig, bgcolor = html_bgcolor(np.mean(newSpectraNow),newSpectra)

        f.write("%e %e %e %e %s\n"%(ff_ave[i],ff_ave[i+1],np.mean(newSpectraNow),sig,bgcolor))

    f.close()

    if params["doPlots"]:

        plotLocation = params["path"] + "/" + channel.station_underscore
        if not os.path.isdir(plotLocation):
            os.makedirs(plotLocation)

        fl, low, fh, high = pylal.pylal_seismon_NLNM.NLNM(2)

        plt.semilogx(freqNow,spectraNow, 'k', label='Current')
        plt.semilogx(freq,spectral_variation_norm_10per,'b',label='10')
        plt.semilogx(freq,spectral_variation_norm_50per,'r',label='50')
        plt.semilogx(freq,spectral_variation_norm_90per,'g',label='90')
        plt.loglog(fl,low,'k-.')
        plt.loglog(fh,high,'k-.',label='LNM/HNM')
        plt.legend(loc=3,prop={'size':10})
        plt.xlim([params["fmin"],params["fmax"]])
        plt.ylim([10**-10, 10**-5])
        plt.xlabel("Frequency [Hz]")
        plt.ylabel("Seismic Spectrum [(m/s)/rtHz]")
        plt.grid()
        plt.show()
        plt.savefig(os.path.join(plotLocation,"psd.png"),dpi=200)
        plt.savefig(os.path.join(plotLocation,"psd.eps"),dpi=200)
        plt.close('all')

        plt.semilogx(freqNow,[y/x for x,y in zip(freqNow,spectraNow)], 'k', label='Current')
        plt.semilogx(freq,[y/x for x,y in zip(freq,spectral_variation_norm_10per)],'b',label='10')
        plt.semilogx(freq,[y/x for x,y in zip(freq,spectral_variation_norm_50per)],'r',label='50')
        plt.semilogx(freq,[y/x for x,y in zip(freq,spectral_variation_norm_90per)],'g',label='90')
        plt.loglog(fl,[y/x for x,y in zip(fl,low)],'k-.')
        plt.loglog(fh,[y/x for x,y in zip(fl,high)],'k-.',label='LNM/HNM')
        plt.legend(loc=3,prop={'size':10})
        plt.xlim([params["fmin"],params["fmax"]])
        plt.ylim([10**-10, 10**-5])
        plt.xlabel("Frequency [Hz]")
        plt.ylabel("Seismic Spectrum [m/rtHz]")
        plt.grid()
        plt.show()
        plt.savefig(os.path.join(plotLocation,"disp.png"),dpi=200)
        plt.savefig(os.path.join(plotLocation,"disp.eps"),dpi=200)
        plt.close('all')

        indexes = np.unique(np.floor(np.logspace(0, np.log10(len(freq)-1), num=100)))
        indices = [int(x) for x in indexes]

        #X,Y = np.meshgrid(freq, range_binning)
        X,Y = np.meshgrid(freq[indices], range_binning)
        ax = plt.subplot(111)
        #im = plt.pcolor(X,Y,np.transpose(spectral_variation_norm), cmap=plt.cm.jet)
        im = plt.pcolor(X,Y,np.transpose(spectral_variation_norm[indices,:]), cmap=plt.cm.jet)
        ax.set_xscale('log')
        ax.set_yscale('log')
        plt.semilogx(freqNow,spectraNow, 'k', label='Current')
        plt.semilogx(freq,spectral_variation_norm_10per,'w',label='10')
        plt.semilogx(freq,spectral_variation_norm_50per,'w',label='50')
        plt.semilogx(freq,spectral_variation_norm_90per,'w',label='90')
        plt.loglog(fl,low,'k-.')
        plt.loglog(fh,high,'k-.',label='LNM/HNM')  
        plt.xlim([params["fmin"],params["fmax"]])
        plt.ylim([10**-10, 10**-5])
        plt.xlabel("Frequency [Hz]")
        plt.ylabel("Seismic Spectrum [(m/s)/rtHz]")
        plt.grid()
        plt.show()
        plt.savefig(os.path.join(plotLocation,"specvar.png"),dpi=200)
        plt.savefig(os.path.join(plotLocation,"specvar.eps"),dpi=200)
        plt.close('all')

        ttStart = np.array(ttStart)
        indices_ttStart = np.where(ttStart >= params["gpsStart"] - 12*60*60)
        ttStart = ttStart[indices_ttStart]

        spectra = np.squeeze(spectra[indices_ttStart,:])
        ttStartMin = min(ttStart)
        tt = [(c-ttStartMin)/(60*60) for c in ttStart]

        #X,Y = np.meshgrid(freq, tt)
        X,Y = np.meshgrid(freq[indices], tt)
        ax = plt.subplot(111)
        #im = plt.pcolor(X,Y,np.log10(spectra), cmap=plt.cm.jet, vmin=-9, vmax=-5)
        im = plt.pcolor(X,Y,np.log10(spectra[:,indices]), cmap=plt.cm.jet, vmin=-9, vmax=-5)
        ax.set_xscale('log')
        plt.xlim([params["fmin"],params["fmax"]])
        plt.xlabel("Frequency [Hz]")
        plt.ylabel("Time [Hours]")
        cbar=plt.colorbar()
        cbar.set_label('log10(Seismic Spectrum [(m/s)/rtHz])') 
        plt.grid()
        plt.show()
        plt.savefig(os.path.join(plotLocation,"tf.png"),dpi=200)
        plt.savefig(os.path.join(plotLocation,"tf.eps"),dpi=200)
        plt.close('all')

    htmlPage = pylal.pylal_seismon_html.seismon_page(channel,textLocation)
    if htmlPage is not None:
        f = open(os.path.join(textLocation,"psd.html"),"w")
        f.write(htmlPage)
        f.close()


