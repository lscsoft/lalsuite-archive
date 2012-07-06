#!/usr/bin/env python

# Copyright (C) 2011 Duncan Macleod
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

"""
This module provides some wrapper functions for plotting TimeSeries and FrequencySeries using classes defined in pylal.plotutils.
"""

# =============================================================================
# Preamble
# =============================================================================

from __future__ import division

import datetime
import re
import itertools
import numpy

import swiglal as lal
import swiglalframe as lalframe

from pylal import plotutils
from pylal import seriesutils
from pylal import git_version

__author__  = "Duncan M. Macleod <duncan.macleod@ligo.org>"
__version__ = git_version.id
__date__    = git_version.date

# =============================================================================
# Plot TimeSeries
# =============================================================================

def plottimeseries(series, outfile, t0=0, zeroline=False, **kwargs):
    """
    Plot a REALXTimeSeries.
    """
    # construct list of series
    if hasattr(series, "__contains__"):
        serieslist = series
    else:
        serieslist = [series]

    # get limits
    xlim = kwargs.pop("xlim", None)
    ylim = kwargs.pop("ylim", None)
    if xlim:
        start,end = xlim
    else:
        start = min(float(s.epoch) for s in serieslist)
        end   = max(float(s.epoch + s.data.length*s.deltaT) for s in serieslist)

    # get axis scales
    logx = kwargs.pop("logx", False)
    logy = kwargs.pop("logy", False)

    # get legend loc
    loc = kwargs.pop("loc", 0)
    alpha = kwargs.pop("alpha", 0.8)

    # get colorbar options
    hidden_colorbar = kwargs.pop("hidden_colorbar", False)

    # get savefig option
    bbox_inches = kwargs.pop("bbox_inches", None)

    # set default params
    if kwargs.has_key("marker"):
        kwargs.setdefault("markersize", 5)
    else:
        kwargs.setdefault("linestyle", "-")
        kwargs.setdefault("linewidth", "1")
    
    #
    # get labels
    #

    xlabel = kwargs.pop("xlabel", None)
    if xlabel:
        unit = 1
    if not xlabel:
        unit, timestr = plotutils.time_axis_unit(end-start)
        if not t0:
            t0 = lal.LIGOTimeGPS(start)
        if int(t0.gpsNanoSeconds) == 0:
            xlabel = datetime.datetime(*lal.XLALGPSToUTC([2000]*6,int(t0))[:6])\
                         .strftime("%B %d %Y, %H:%M:%S %ZUTC")
            xlabel = "Time (%s) since %s (%s)" % (timestr, xlabel, int(t0))
        else:
            xlabel = datetime.datetime(*lal.XLALGPSToUTC([2000]*6,\
                                                         t0.gpsSeconds)[:6])\
                          .strftime("%B %d %Y, %H:%M:%S %ZUTC")
            xlabel = "Time (%s) since %s (%s)"\
                     % (timestr,
                        xlabel.replace(" UTC",".%.3s UTC" % t0.gpsNanoSeconds),\
                        t0)
    ylabel   = kwargs.pop("ylabel", "Amplitude (counts)")
    title    = kwargs.pop("title", "")
    subtitle = kwargs.pop("subtitle", "")

    #
    # make plot
    #

    allnames    = [s.name for s in serieslist]
    namedseries = [s for s in serieslist if not re.search("(min|max)\Z",s.name)]

    plot = plotutils.SimplePlot(xlabel, ylabel, title, subtitle)
    for i,(series,c) in enumerate(itertools.izip(namedseries,\
                                                 plotutils.default_colors())):
        x = numpy.arange(series.data.length) * series.deltaT + series.epoch - t0
        x = x.astype(float)
        x /= unit
        plot.add_content(x, series.data.data, color=c,\
                         label=plotutils.display_name(series.name), **kwargs)
        # find min/max and plot
        for i,name in enumerate(allnames):
            for ext in ["min", "max"]:
                if re.match("%s\w%s" % (name, ext), series.name):
                    series2 = serieslist[i]
                    x2 = numpy.arange(series2.data.length) * series2.deltaT\
                         + series2.epoch - t0
                    x2 /= unit
                    plot.ax.plot(x2, series2.data.data, color=c, linewidth=0.3,\
                                 **kwargs)
                    plot.ax.fill_between(x2, series.data.data,\
                                         series2.data.data, color=c, alpha=0.25)

    # finalize
    plot.finalize(loc=loc, alpha=alpha)
    if hidden_colorbar:
        plotutils.add_colorbar(plot.ax, visible=False)
    
    # add zero line
    axis_lims = plot.ax.get_ylim()
    if zeroline:
        plot.ax.plot([0, 0], [axis_lims[0], axis_lims[1]], 'r--', linewidth=2)
        plot.ax.set_ylim([ axis_lims[0], axis_lims[1] ])

    # set logscale
    if logx:
        plot.ax.set_xscale("log")
    if logy:
        plot.ax.set_yscale("log")

    # format axes
    if xlim:
        xlim = (numpy.asarray(xlim).astype(float)-t0)/unit
        plot.ax.set_xlim(xlim)
    if ylim:
        plot.ax.set_ylim(ylim)

    plotutils.set_time_ticks(plot.ax)
    plotutils.set_minor_ticks(plot.ax, x=False)

    # save and close
    plot.savefig(outfile, bbox_inches=bbox_inches,\
                 bbox_extra_artists=plot.ax.texts)
    plot.close()

# =============================================================================
# Plot FrequencySeries
# =============================================================================

def plotfrequencyseries(series, outfile, **kwargs):
    """
    Plot a (swig)LAL FrequencySeries.
    """
    # construct list of series
    if hasattr(series, "__contains__"):
        serieslist = series
    else:
        serieslist = [series]

    # get limits
    xlim = kwargs.pop("xlim", None)
    ylim = kwargs.pop("ylim", None)
    if not xlim:
        fmin = min(float(s.f0) for s in serieslist)
        fmax = max(float(s.f0 + s.data.length*s.deltaF) for s in serieslist)
        xlim = [fmin, fmax]

    # get axis scales
    logx = kwargs.pop("logx", False)
    logy = kwargs.pop("logy", False)

    # get legend loc
    loc = kwargs.pop("loc", 0)
    alpha = kwargs.pop("alpha", 0.8)

    # get colorbar options
    hidden_colorbar = kwargs.pop("hidden_colorbar", False)

    # get savefig option
    bbox_inches = kwargs.pop("bbox_inches", None)

    # set default params
    if kwargs.has_key("marker"):
        kwargs.setdefault("markersize", 5)
    else:
        kwargs.setdefault("linestyle", "-")
        kwargs.setdefault("linewidth", "1")
    
    #
    # get labels
    #

    xlabel = kwargs.pop("xlabel", "Frequency (Hz)")
    ylabel   = kwargs.pop("ylabel", "Amplitude")
    title    = kwargs.pop("title", "")
    subtitle = kwargs.pop("subtitle", "")

    #
    # make plot
    #

    allnames    = [s.name for s in serieslist]
    namedseries = [s for s in serieslist if not re.search("(min|max)\Z",s.name)]

    plot = plotutils.SimplePlot(xlabel, ylabel, title, subtitle)
    for i,(series,c) in enumerate(itertools.izip(namedseries,\
                                                 plotutils.default_colors())):
        x = numpy.arange(series.data.length) * series.deltaF + series.f0 
        x = x.astype(float)
        plot.add_content(x, series.data.data, color=c,\
                         label=plotutils.display_name(series.name), **kwargs)
        # find min/max and plot
        for i,name in enumerate(allnames):
            for ext in ["min", "max"]:
                if re.match("%s\w%s" % (name, ext), series.name):
                    series2 = serieslist[i]
                    x2 = numpy.arange(series2.data.length) * series2.deltaF\
                         + series2.f0 
                    plot.ax.plot(x2, series2.data.data, color=c, linewidth=0.3,\
                                 **kwargs)
                    plot.ax.fill_between(x2, series.data.data,\
                                         series2.data.data, color=c, alpha=0.25)

    # finalize
    plot.finalize(loc=loc, alpha=alpha)
    if hidden_colorbar:
        plotutils.add_colorbar(plot.ax, visible=False)
    
    # set logscale
    if logx:
        plot.ax.set_xscale("log")
    if logy:
        plot.ax.set_yscale("log")

    # format axes
    if xlim:
        plot.ax.set_xlim(xlim)
    if ylim:
        plot.ax.set_ylim(ylim)

    plotutils.set_minor_ticks(plot.ax)

    # save and close
    plot.savefig(outfile, bbox_inches=bbox_inches,\
                 bbox_extra_artists=plot.ax.texts)
    plot.close()
