
import os, glob, optparse, shutil, warnings, matplotlib, pickle, math, copy, pickle, time
import numpy as np

import pylal.pylal_seismon_utils

import gwpy.time, gwpy.timeseries, gwpy.spectrum, gwpy.plotter
import gwpy.segments

import obspy.core, obspy.signal.array_analysis
import matplotlib.pyplot as plt
import matplotlib.colorbar, matplotlib.cm, matplotlib.colors


__author__ = "Michael Coughlin <michael.coughlin@ligo.org>"
__date__ = "2012/8/26"
__version__ = "0.1"

# =============================================================================
#
#                               DEFINITIONS
#
# =============================================================================

def beamforming(params, segment):
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

    st = obspy.core.Stream()

    tstart = pylal.pylal_seismon_utils.GPSToUTCDateTime(gpsStart)
    tend = pylal.pylal_seismon_utils.GPSToUTCDateTime(gpsEnd)

    for channel in params["channels"]:
        # make timeseries
        dataFull = pylal.pylal_seismon_utils.retrieve_timeseries(params, channel, segment)
        if dataFull == []:
            continue

        dataFull = dataFull / channel.calibration
        indexes = np.where(np.isnan(dataFull.data))[0]
        meanSamples = np.mean(np.ma.masked_array(dataFull.data,np.isnan(dataFull.data)))
        for index in indexes:
            dataFull[index] = meanSamples
        dataFull -= np.mean(dataFull.data)

        trace = obspy.core.Trace()
        trace.stats.network = ""
        trace.stats.station = channel.station
        trace.stats.location = ""
        trace.stats.channel = ""
        trace.stats.sampling_rate = channel.samplef
        trace.stats.delta = 1/channel.samplef
        trace.stats.calib = 1
        trace.data = np.array(dataFull.data)

        trace.stats.starttime = tstart
        trace.stats.npts = len(dataFull.data)

        trace.stats.coordinates = obspy.core.AttribDict({
            'latitude': channel.latitude,
            'longitude': channel.longitude,
            'elevation': 0.0})

        st += trace

    print st

    # Execute sonic
    kwargs = dict(
        # slowness grid: X min, X max, Y min, Y max, Slow Step
        sll_x=-3.0, slm_x=3.0, sll_y=-3.0, slm_y=3.0, sl_s=0.03,
        # sliding window properties
        win_len=1.0, win_frac=0.05,
        # frequency properties
        frqlow=0.0, frqhigh=1.0, prewhiten=0,
        # restrict output
        semb_thres=-1e9, vel_thres=-1e9, timestamp='mlabday',
        #stime=tstart, etime=tstart+300,
        stime=tstart, etime=tend-1,
        coordsys = 'lonlat', method = 0
    )

    out = obspy.signal.array_analysis.array_processing(st, **kwargs)

    if params["doPlots"]:

        plotDirectory = params["path"] + "/beamforming"
        pylal.pylal_seismon_utils.mkdir(plotDirectory)

        pngFile = os.path.join(plotDirectory,"timeseries.png")

        # Plot
        labels = ['rel.power', 'abs.power', 'baz', 'slow']

        fig = plt.figure(figsize=(16, 16))
        for i, lab in enumerate(labels):
            ax = fig.add_subplot(4, 1, i + 1)
            ax.scatter(out[:, 0], out[:, i + 1], c=out[:, 1], alpha=0.6,
                   edgecolors='none')
            ax.set_ylabel(lab)
            ax.set_xlim(out[0, 0], out[-1, 0])
            ax.set_ylim(out[:, i + 1].min(), out[:, i + 1].max())

        fig.autofmt_xdate()
        fig.subplots_adjust(top=0.95, right=0.95, bottom=0.2, hspace=0)
        plt.show()
        plt.savefig(pngFile,dpi=200)
        plt.close('all')

        cmap = matplotlib.cm.hot_r

        # make output human readable, adjust backazimuth to values between 0 and 360
        t, rel_power, abs_power, baz, slow = out.T
        baz[baz < 0.0] += 360

        # choose number of fractions in plot (desirably 360 degree/N is an integer!)
        N = 30
        abins = np.arange(N + 1) * 360. / N
        sbins = np.linspace(0, 3, N + 1)

        # sum rel power in bins given by abins and sbins
        hist, baz_edges, sl_edges = np.histogram2d(baz, slow,
                bins=[abins, sbins], weights=rel_power)

        # transform to gradient
        baz_edges = baz_edges / 180 * np.pi

        pngFile = os.path.join(plotDirectory,"polar.png")

        # add polar and colorbar axes
        fig = plt.figure(figsize=(16, 16))
        cax = fig.add_axes([0.85, 0.2, 0.05, 0.5])
        ax = fig.add_axes([0.10, 0.1, 0.70, 0.7], polar=True)

        dh = abs(sl_edges[1] - sl_edges[0])
        dw = abs(baz_edges[1] - baz_edges[0])

        # circle through backazimuth
        for i, row in enumerate(hist):
            bars = ax.bar(left=(np.pi / 2 - (i + 1) * dw) * np.ones(N),
                          height=dh * np.ones(N),
                          width=dw, bottom=dh * np.arange(N),
                          color=cmap(row / hist.max()))

        ax.set_xticks([np.pi / 2, 0, 3. / 2 * np.pi, np.pi])
        ax.set_xticklabels(['N', 'E', 'S', 'W'])

        # set slowness limits
        ax.set_ylim(0, 3)
        matplotlib.colorbar.ColorbarBase(cax, cmap=cmap,
                     norm=matplotlib.colors.Normalize(vmin=hist.min(), vmax=hist.max()))

        plt.show()
        plt.savefig(pngFile,dpi=200)
        plt.close('all')

