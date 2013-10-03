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

def coherence(params, channel1, channel2, segment):
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
    dataFull1 = pylal.pylal_seismon_utils.retrieve_timeseries(params, channel1, segment)
    if dataFull1 == []:
        return

    dataFull1 = dataFull1 / channel1.calibration
    indexes = np.where(np.isnan(dataFull1.data))[0]
    meanSamples = np.mean(np.ma.masked_array(dataFull1.data,np.isnan(dataFull1.data)))
    for index in indexes:
        dataFull1[index] = meanSamples
    dataFull1 -= np.mean(dataFull1.data)

    if np.mean(dataFull1.data) == 0.0:
        print "data only zeroes... continuing\n"
        return
    if len(dataFull1.data) < 2*channel1.samplef:
        print "timeseries too short for analysis... continuing\n"
        return

    dataFull2 = pylal.pylal_seismon_utils.retrieve_timeseries(params, channel1, segment)
    if dataFull2 == []:
        return

    dataFull2 = dataFull2 / channel2.calibration
    indexes = np.where(np.isnan(dataFull2.data))[0]
    meanSamples = np.mean(np.ma.masked_array(dataFull2.data,np.isnan(dataFull2.data)))
    for index in indexes:
        dataFull2[index] = meanSamples
    dataFull2 -= np.mean(dataFull2.data)

    if np.mean(dataFull2.data) == 0.0:
        print "data only zeroes... continuing\n"
        return
    if len(dataFull2.data) < 2*channel2.samplef:
        print "timeseries too short for analysis... continuing\n"
        return

    gpss = np.arange(gpsStart,gpsEnd,params["fftDuration"])
    fft1 = []
    fft2 = []
    for i in xrange(len(gpss)-1):

        tt1 = np.array(dataFull1.times)
        indexes = np.intersect1d(np.where(tt1 >= gpss[i])[0],np.where(tt1 <= gpss[i+1])[0])
        indexMin = np.min(indexes)
        indexMax = np.max(indexes)
        dataCut1 = dataFull1[indexMin:indexMax]
        dataCut1FFT = dataCut1.fft()

        tt2 = np.array(dataFull2.times)
        indexes = np.intersect1d(np.where(tt2 >= gpss[i])[0],np.where(tt2 <= gpss[i+1])[0])
        indexMin = np.min(indexes)
        indexMax = np.max(indexes)
        dataCut2 = dataFull2[indexMin:indexMax]
        dataCut2FFT = dataCut2.fft()

        print dataCut1FFT
        print penis

        fft1.append(dataCut1FFT)
        fft2.append(dataCut2FFT)

    specgram1 = gwpy.spectrogram.Spectrogram.from_spectra(*fft1, dt=params["fftDuration"],epoch=gpsStart)
    specgram2 = gwpy.spectrogram.Spectrogram.from_spectra(*fft2, dt=params["fftDuration"],epoch=gpsStart)
