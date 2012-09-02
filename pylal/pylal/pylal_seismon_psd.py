#!/usr/bin/python

from pylab import *
import os, glob, optparse, shutil, warnings, numpy, matplotlib, pickle, math, copy, pickle
from collections import namedtuple
from subprocess import Popen, PIPE, STDOUT
from pylal.xlal.datatypes.ligotimegps import LIGOTimeGPS
from pylal.xlal.date import XLALGPSToUTC
from pylal import Fr
from pylal.dq import dqDataUtils
import seismon_NLNM, seismon_html

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
        frame_data,data_start,_,dt,_,_ = Fr.frgetvect1d(frame,channel.station)
        frame_length = float(dt)*len(frame_data)
        frame_time = data_start+dt*numpy.arange(len(frame_data))

        for i in range(len(frame_data)):
            time.append(frame_time[i])
            data.append(frame_data[i])
        else:
            for i in range(len(frame_data)):
                if frame_time[i] <= start_time:  continue
                if frame_time[i] >= end_time:  continue
                time.append(frame_time[i])
                data.append(frame_data[i])
        data = [e/channel.calibration for e in data]

    return time,data

def mat(params, channel):

    psdLocation = params["dirPath"] + "/Text_Files/PSD/" + channel.station_underscore
    if not os.path.isdir(psdLocation):
        os.makedirs(psdLocation)

    time,data = read_frames(params["frameGPS"],params["frameGPS"]+params["frameDur"],channel,params["frame"])
    #freq,spectra = dqDataUtils.spectrum(data, channel.samplef, NFFT=len(data), overlap = 0)

    spectra, freq = matplotlib.pyplot.psd(data, NFFT=64*channel.samplef, Fs=channel.samplef, Fc=0, detrend=matplotlib.mlab.detrend_mean,window=matplotlib.mlab.window_hanning)

    spectra = [math.sqrt(e) for e in spectra]

    newSpectra = []
    newFreq = []

    for i in xrange(len(freq)):
        if freq[i] <= params["fmax"]:
            newFreq.append(freq[i])
            newSpectra.append(spectra[i])

    spectra = newSpectra
    freq = newFreq

    psdFile = psdLocation + "/" + str(params["frameGPS"]) + ".pickle"
    itemlist = zip(freq,spectra)

    f = open(psdFile,"wb")
    pickle.dump(itemlist,f)
    f.close()

    if params["doPlots"] == 1:

        plotLocation = params["path"] + "/" + channel.station_underscore
        if not os.path.isdir(plotLocation):
            os.makedirs(plotLocation)        

        fl, low, fh, high = seismon_NLNM.NLNM(2)

        semilogx(freq,spectra, 'k')
        loglog(fl,low,'k-.',fh,high,'k-.')
        xlim([params["fmin"],params["fmax"]])
        ylim([10**-10, 10**-5])
        xlabel("Frequency [Hz]")
        ylabel("Seismic Spectrum [(m/s)/\surd Hz]")
        grid
        show()
        savefig(os.path.join(plotLocation,"psd.png"),dpi=200)
        savefig(os.path.join(plotLocation,"psd.eps"),dpi=200)
        close('all')

def spectral_histogram(data,bins):
    
    spectral_variation_norm = []
    rows, columns = data.shape

    for i in xrange(columns):
        this_spectral_variation, bin_edges = numpy.histogram(data[:,i],bins)
        this_spectral_variation = numpy.array(this_spectral_variation)
        np = (100/float(sum(this_spectral_variation))) + numpy.zeros(this_spectral_variation.shape)
        if spectral_variation_norm == []:
            spectral_variation_norm = this_spectral_variation * np
        else:
            spectral_variation_norm = numpy.vstack([spectral_variation_norm,this_spectral_variation * np])
    return spectral_variation_norm

def calculate_percentiles(data,bins,percentile):

    cumsumvals = numpy.cumsum(data)
    
    abs_cumsumvals_minus_percentile = abs(cumsumvals - percentile)
    minindex = abs_cumsumvals_minus_percentile.argmin()
    val = bins[minindex]

    return val

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

def analysis(params, channel):

    psdLocation = params["dirPath"] + "/Text_Files/PSD/" + channel.station_underscore
    files = glob.glob(os.path.join(psdLocation,"*.pickle"))

    tt = []
    freq = []
    spectra = []

    for file in files:

        fileSplit = file.split("/")
        thisTT = float(fileSplit[-1].replace(".pickle",""))
        tt.append(thisTT)

        f = open(file,"rb")
        itemlist = pickle.load(f)
        f.close()  
        thisFreq, thisSpectra = zip(*itemlist) 

        if freq == []:
            freq = numpy.array(thisFreq)
        if spectra == []:
            spectra = numpy.array(thisSpectra)
        else:
            spectra = numpy.vstack([spectra,numpy.array(copy.copy(thisSpectra))])
        if thisTT == params["frameGPS"]:
            freqNow = copy.copy(thisFreq)
            spectraNow = copy.copy(thisSpectra)

    # Binning parameters
    nb = 500
    range_binning = numpy.logspace(-10,-5,num=nb);

    spectra[spectra==Inf] = 0

    # Calculate bin histogram of PSDs
    which_spectra = numpy.all([spectra.max(axis = 1) <= 10**-4,spectra.max(axis = 1) >= 10**-10],axis=0)

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
                    newSpectra = numpy.vstack([newSpectra,spectra[:,j]])

        if len(newSpectra.shape) > 1:
            sig, bgcolor = html_bgcolor(numpy.mean(newSpectraNow),numpy.mean(newSpectra, axis = 0))
        else:
            sig, bgcolor = html_bgcolor(numpy.mean(newSpectraNow),newSpectra)

        sigDict.append({"flow":ff_ave[i],"fhigh":ff_ave[i+1],"meanPSD":numpy.mean(newSpectraNow),"sig":sig,"bgcolor":bgcolor})  
    f = open(os.path.join(textLocation,"sig.pickle"),"w")
    pickle.dump(sigDict, f)
    f.close()

    if params["doPlots"] == 1:

        plotLocation = params["path"] + "/" + channel.station_underscore
        if not os.path.isdir(plotLocation):
            os.makedirs(plotLocation)

        fl, low, fh, high = seismon_NLNM.NLNM(2)

        semilogx(freqNow,spectraNow, 'k', label='Current')
        semilogx(freq,spectral_variation_norm_10per,'b',label='10')
        semilogx(freq,spectral_variation_norm_50per,'r',label='50')
        semilogx(freq,spectral_variation_norm_90per,'g',label='90')
        loglog(fl,low,'k-.')
        loglog(fh,high,'k-.',label='LNM/HNM')
        legend()
        xlim([params["fmin"],params["fmax"]])
        ylim([10**-10, 10**-5])
        xlabel("Frequency [Hz]")
        ylabel("Seismic Spectrum [(m/s)/\surd Hz]")
        grid
        show()
        savefig(os.path.join(plotLocation,"psd.png"),dpi=200)
        savefig(os.path.join(plotLocation,"psd.eps"),dpi=200)
        close('all')

        X,Y = meshgrid(freq, range_binning)
        ax = subplot(111)
        im = pcolor(X,Y,transpose(spectral_variation_norm), cmap=cm.jet)
        ax.set_xscale('log')
        ax.set_yscale('log')
        #im.set_interpolation('bilinear')
        semilogx(freqNow,spectraNow, 'k', label='Current')
        semilogx(freq,spectral_variation_norm_10per,'w',label='10')
        semilogx(freq,spectral_variation_norm_50per,'w',label='50')
        semilogx(freq,spectral_variation_norm_90per,'w',label='90')
        loglog(fl,low,'k-.')
        loglog(fh,high,'k-.',label='LNM/HNM')  
        xlim([params["fmin"],params["fmax"]])
        ylim([10**-10, 10**-5])
        xlabel("Frequency [Hz]")
        ylabel("Seismic Spectrum [(m/s)/\surd Hz]")        
        grid
        show()
        savefig(os.path.join(plotLocation,"specvar.png"),dpi=200)
        savefig(os.path.join(plotLocation,"specvar.eps"),dpi=200)
        close('all')

        ttStart = min(tt)
        tt = [(c-ttStart)/(60*60) for c in tt]

        X,Y = meshgrid(freq, tt)
        ax = subplot(111)
        im = pcolor(X,Y,numpy.log10(spectra), cmap=cm.jet, vmin=-9, vmax=-5)
        ax.set_xscale('log')
        #im.set_interpolation('bilinear')
        xlim([params["fmin"],params["fmax"]])
        xlabel("Frequency [Hz]")
        ylabel("Time [Hours]")
        cbar=colorbar()
        cbar.set_label('log10(Seismic Spectrum [(m/s)/\surd Hz])') 
        grid
        show()
        savefig(os.path.join(plotLocation,"tf.png"),dpi=200)
        savefig(os.path.join(plotLocation,"tf.eps"),dpi=200)
        close('all')

    htmlPage = seismon_html.seismon_page(channel,textLocation)
    if htmlPage is not None:
        f = open(os.path.join(textLocation,"psd.html"),"w")
        f.write(htmlPage)
        f.close()


