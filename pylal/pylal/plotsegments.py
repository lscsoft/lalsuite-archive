# encoding: utf-8
#
# Copyright (C) 2008  Nickolas V Fotopoulos
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
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

__author__ = "Nickolas Fotopoulos <nvf@gravity.phys.uwm.edu>"

import datetime, re

from glue import segments
from pylal import plotutils
from pylal import date
from pylal.datatypes import LIGOTimeGPS
import pylab

##############################################################################
# Plotting class

class PlotSegmentsPlot(object):
  """
  Represent a plotsegments plot.  To use it, one must instantiate a
  PlotSegmentsPlot object, call add_contents(), then call finalize().
  To save, call the savefig() method.  To display to the screen, call
  pylab.show().
  
  Developer note: There is the transformation function _time_transform.
  Follow the convention of applying it before storing internal data so that
  there is no ambiguity about whether data has been transformed or not.
  """
  color_code = {'H1':'r', 'H2':'b', 'L1':'g', 'V1':'m', 'G1':'k'}
  
  def __init__(self, t0=0):
    """
    Create a fresh plot.  Provide t0 to provide a reference time to use as
    zero.
    """
    self.fig = pylab.figure()
    self.ax = self.fig.add_subplot(111)
    self.savefig = self.fig.savefig
    
    self.window = None
    self.ifos = []
    self._time_transform = lambda t: t - t0
    
    self.ax.set_xlabel("time (s)")
    self.ax.set_ylabel("IFO")
  
  def add_contents(self, segdict, ifos=None):
    """
    Add the contents of segdict to the plot.  Provide the list ifos, if you
    wish to specify the order in which they appear or wish to plot a subset of
    the segmentlists in segdict.
    """
    if ifos is None:
      ifos = segdict.keys()
      ifos.sort()
    self.ifos.extend(ifos[::-1])
    for row, ifo in enumerate(ifos[::-1]):
      color = self.color_code[ifo]
      for seg in segdict[ifo]:
        a = self._time_transform(seg[0])
        b = self._time_transform(seg[1])
        self.ax.fill([a, b, b, a, a], [row, row, row+1, row+1, row], color)
  
  def set_window(self, window_seg, padding=0):
    """
    Define a window of interest by setting the x-limits of the plot
    appropriately.  If padding is also present, protract the x-limits by
    that quantity and mark the unpadded window with solid black lines.
    """
    a = self._time_transform(window_seg[0])
    b = self._time_transform(window_seg[1])
    self.window = segments.segment((a - padding, b + padding))
    
    if padding > 0:
      self.ax.axvline(a, color='k', linewidth=2)
      self.ax.axvline(b, color='k', linewidth=2)
  
  def highlight_segment(self, seg):
    """
    Highlight a particular segment with dashed lines.
    """
    self.ax.axvline(self._time_transform(seg[0]), color='k', linestyle='--')
    self.ax.axvline(self._time_transform(seg[1]), color='k', linestyle='--')
  
  def finalize(self):
    """
    Make final changes to the plot with the guarantee that no additional data
    will be added.
    """
    ticks = pylab.arange(len(self.ifos)) + 0.5
    self.ax.set_yticks(ticks)
    self.ax.set_yticklabels(self.ifos)
    
    if self.window is not None:
      self.ax.set_xlim(self.window)
    self.ax.set_ylim((0, len(self.ifos)))

  def close(self):
    pylab.close(self.fig)
  
  def __del__(self):
    self.close()

# =============================================================================
# Plot segments
# =============================================================================

def plotsegmentlistdict(segdict, outfile, keys=None, start=None,\
                        end=None, zero=None, highlight_segments=None,\
                        insetlabels=None, **kwargs):
    """
    Plots a glue.segments.segmentlistdict using the PlotSegmentsPlot class from
    pylal.plotutils.

    Arguments:

        segdict : glue.segments.segmentlistdict
            dictionary of (name, segmentlist) pairs to plot.
        outfile : str
            filepath to write to

    Keyword arguments:

        keys : list
            ordered list of keys to use in this order on the plot
        start : float
            GPS start time for plot
        end : float
            GPS end time for plot
        zero : float
            GPS time around which to zero plot
        highlight_segments : glue.segments.segmentlistdict
            list of segments to highlight with vertical red dashed lines
        insetlabels : [ True | False ]
            write labels inside the plot axis

    Unnamed keyword arguments:

        xlabel : string
            label for x-axis
        ylabel : string
            label for y-axis
        title : string
            title for plot
        subtitle : string
            subtitle for plot
        bbox_inches : str
            use "tight" to get a bounding box tight to the axis.
    """
    # get zeros
    if not start or not end:
        extents = [seg.extent() for seg in segdict.values()]
    if not start:
        start = min(s[0] for s in extents)
    if not end:
        end   = max(s[1] for s in extents)
    if not zero:
        zero = start

    # get unit for plot
    unit, timestr = plotutils.time_axis_unit(end-start)

    # set labels
    zero = LIGOTimeGPS("%.3f" % zero)
    if zero.nanoseconds==0:
        tlabel = datetime.datetime(*date.XLALGPSToUTC(LIGOTimeGPS(zero))[:6])\
                     .strftime("%B %d %Y, %H:%M:%S %ZUTC")
    else:
        tlabel = datetime.datetime(*date.XLALGPSToUTC(LIGOTimeGPS(\
                     zero.seconds))[:6]).strftime("%B %d %Y, %H:%M:%S %ZUTC")
        tlabel = tlabel.replace(" UTC", ".%.3s UTC" % zero.nanoseconds)
    xlabel   = kwargs.pop("xlabel",\
                          "Time (%s) since %s (%s)" % (timestr, tlabel, zero))
    ylabel   = kwargs.pop("ylabel", "")
    title    = kwargs.pop("title", "")
    subtitle = kwargs.pop("subtitle", "")

    # get other parameters
    labels_inset = kwargs.pop("labels_inset", False)
    bbox_inches  = kwargs.pop("bbox_inches", "tight")

    # escape underscores for latex text
    if pylab.rcParams["text.usetex"]:
        if keys:
            keys = [re.sub('(?<!\\\\)_', '\_', key) for key in keys]
        segkeys = segdict.keys()
        newdict = segments.segmentlistdict()
        for key in segkeys:
           newkey = re.sub('(?<!\\\\)_', '\_', key)
           newdict[newkey] = segdict[key]
        segdict = newdict

    # generate plot
    plot = plotutils.PlotSegmentsPlot(xlabel, ylabel, title, subtitle,\
                                       t0=zero, dt=unit)
    plot.add_content(segdict, keys, **kwargs)
    plot.finalize(labels_inset=insetlabels)

    # highligh segments
    if highlight_segments:
        for seg in highlight_segments:
            plot.highlight_segment(seg)
 
    # set axis limits
    xlim = [float(start-zero)/unit, float(end-zero)/unit]
    plot.ax.set_xlim(*xlim)

    # set grid
    plot.ax.grid(True,which='major')
    plot.ax.grid(True,which='majorminor')

    # save
    plot.savefig(outfile, bbox_inches=bbox_inches,\
                 bbox_extra_artists=plot.ax.texts)
    plot.close()

# =============================================================================
# Plot segment histogram
# =============================================================================

def plothistogram(segdict, outfile, keys=None, num_bins=100, **kwargs):
    """
    Plots a histogram of segment duration for each entry in the
    glue.segments.segmentlistdict segdict.

    Arguments:

      segdict : glue.segments.segmentlistdict
          list of segments with which to veto triggers, use dict for multiple
          datasets
      outfile : string 
          string path for output plot

    Keyword arguments:

        keys : [ str | list]
            ordered list of keys to use in this order on the plot
        num_bins : int
            number of bins.

    Unnamed keyword arguments:

        logx : [ True | False ]
            display x-axis in log scale
        logy : [ True | False ]
            display y-axis in log scale
        xlim : tuple
            (xmin,xmax) limits for x-axis 
        ylim : tuple
            (ymin,ymax) limits for y-axis 
        xlabel : string
            label for x-axis
        ylabel : string
            label for y-axis
        title : string
            title for plot
        subtitle : string
            subtitle for plot
        bbox_inches : str
            use "tight" to get a bounding box tight to the axis.
    """
    # get limits
    xlim = kwargs.pop('xlim', None)
    ylim = kwargs.pop('ylim', None)

    # get labels
    xlabel = kwargs.pop('xlabel', 'Length of segment (seconds)')
    ylabel = kwargs.pop('ylabel', 'Number of segments')
    title  = kwargs.pop('title',  'Segment Duration Histogram')
    subtitle = kwargs.pop('subtitle', "")

    # get axis scale
    logx = kwargs.pop('logx', False)
    logy = kwargs.pop('logy', False)

    # get savefig option
    bbox_inches = kwargs.pop('bbox_inches', 'tight')
 
    # escape underscores for latex text
    if pylab.rcParams["text.usetex"]:
        if not keys:
            keys = segdict.keys()
        newdict = segments.segmentlistdict()
        for key in segdict.keys():
           newkey = re.sub('(?<!\\\\)_', '\_', key)
           newdict[newkey] = segdict[key]
        segdict = newdict

    # generate plot object
    plot = plotutils.VerticalBarHistogram(xlabel, ylabel, title, subtitle)

    # add each segmentlist
    for flag in keys:
        plot.add_content([float(abs(seg)) for seg in segdict[flag]],\
                          label=flag, **kwargs)

    # finalize plot with histograms
    plot.finalize(num_bins=num_bins, logx=logx, logy=logy)

    # set limits
    plot.ax.autoscale_view(tight=True, scalex=True, scaley=True)
    if ylim:
      plot.ax.set_ylim(map(float, ylim))
    if xlim:
      plot.ax.set_xlim(map(float, xlim))

    # save figure
    plot.savefig(outfile, bbox_inches=bbox_inches,\
                 bbox_extra_artists=plot.ax.texts)
    plot.close()

# =============================================================================
# Plot segment duration histogram
# =============================================================================

def plotduration(seglist, outfile, flag="_", binlength=3600, start=None,\
                 end=None, zero=None, cumulative=False, normalized=False,\
                 **kwargs):
    """
    Plot the duration of entries in a given glue.segments.segmentlist segs,
    to the given output file, using bins of the given length. Output can be in 
    either format='line' or format='bar', and can be cumulative and normalized
    to show segment number instead of time.
    """
    if not start:
        start = seglist.extent()[0]
    if not end: 
        end   = seglist.extent()[1]
    if not zero:
        zero  = start
    
    # restrict to start and end
    span = segments.segmentlist([segments.segment(start, end)])
    segs = segs & span

    unit,timestr = dqPlotUtils.time_unit(end-start)

    # get labels
    if normalized:
        xlabel  = kwargs.pop("xlabel", "Segment number")
    else:
        zerostr = datetime.datetime(*date.XLALGPSToUTC(\
                                    datatypes.LIGOTimeGPS(zero))[:6])\
                                    .strftime("%B %d %Y, %H:%M:%S %ZUTC")
        xlabel  = kwargs.pop('xlabel', 'Time (%s) since %s (%s)'\
                                       % (timestr, zerostr, zero))
    if cumulative:
        ylabel   = kwargs.pop('ylabel', 'Cumulative segment duration (seconds)')
    else:
        ylabel   = kwargs.pop('ylabel', 'Segment duration (seconds)')
    title    = kwargs.pop('title',  'Segment duration summary')
    subtitle = kwargs.pop("subtitle",\
                          "Running time: %.2f %s. Duty cycle: %.2f\%%"\
                          % ((end-start)/unit, timestr,\
                          livetime/(end-start)*100))

    # get axis scale
    logx = kwargs.pop('logx', False)
    logy = kwargs.pop('logy', False)

    # get savefig option
    bbox_inches = kwargs.pop('bbox_inches', 'tight')

    # get axis limits
    if not normalized:
        xlim = kwargs.pop('xlim', (float(start-zero)/unit,float(end-zero)/unit))
    else:
        xlim = kwargs.pop('xlim', (0.5, len(segs)+0.9))
    ylim = kwargs.pop('ylim', None)


