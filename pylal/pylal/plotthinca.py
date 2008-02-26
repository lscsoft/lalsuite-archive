#!/usr/bin/env python
# encoding: utf-8
# $Id$
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

from __future__ import division

__author__ = "Nickolas Fotopoulos <nvf@gravity.phys.uwm.edu>"
__version__ = "$Revision$"[11:-2]
__date__ = "$Date$"

import bisect
itertools = __import__("itertools")
import sys

from glue import iterutils
from pylal import viz

import numpy
import pylab

epsilon = 1e-16

# helpful warning
try:
  from pylal import rate
except ImportError:
  print >>sys.stderr, "warning: pylal.rate import failed, probably due to "\
        "a scipy import failure.  If you later get a NameError exception "\
        "saying that rate isn't found, it's probably due to scipy not being "\
        "installed."

##############################################################################
# Supporting classes

class BackgroundEstimator(list):
  """
  BackgroundEstimator should hold the SNRs of your background triggers.
  When you've added all the SNRs you care to add, call finalize().
  You can then call your BackgroundEstimator instance with a trigger of
  interest to determine its expected background rate.
  
  It initializes like a normal list, except that it stores bg_fg_ratio
  (a kwarg), which is the ratio of background live-time to foreground
  live-time.  For the inspiral search, this is typically the total number of
  timeslides.
  """
  def __init__(self, *args, **kwargs):
    list.__init__(self, *args)
    if "bg_fg_ratio" in kwargs:
      self.bg_fg_ratio = kwargs["bg_fg_ratio"]
    self.finalized = False
  
  def finalize(self):
    if len(self) == 0:
      raise ValueError, "cannot estimate background with zero triggers"
    self.sort()
    self.finalized = True
  
  def __call__(self, trigger_snr):
    """
    Return p(snr >= trigger_snr | 0), the number of background triggers with
    SNR greater than trigger_snr.  Runs in log2(len(self)) time.
    """
    if not self.finalized:
      raise ValueError, "BackgroundEstimator has not been finalized"
    index = bisect.bisect_left(self, trigger_snr)
    num_as_big = len(self) - index
    effective_num_as_big = num_as_big / self.bg_fg_ratio
    return effective_num_as_big

class FitError(Exception):
    pass

class ExtrapolatingBackgroundEstimator(BackgroundEstimator):
  def finalize(self):
    
    import scipy.optimize
    
    if self.finalized:
      raise ValueError, "ExtrapolatingBackgroundEstimator has already been finalized"
    BackgroundEstimator.finalize(self)
    
    fit_func = lambda p, x: p[0] * numpy.exp(-p[1]*x/2)
    err_func = lambda p, x, y: fit_func(p, x) - y
    
    # initial guess
    p0 = [len(self), 1.]
    
    # fit
    x = numpy.array(self[3*len(self)//4:])
    y = numpy.arange(0, -len(self)//4, -1, dtype=float)
    y += 3*len(self)//4 - 1
    y /= self.bg_fg_ratio
    self.fit_params, err = scipy.optimize.leastsq(err_func, p0[:], args=(x, y))
    if err != 1:
      raise FitError, "background fit did not converge"
    
    # delete SNRs to free memory
    del self[:]
    
    self.finalized = True
  
  def __call__(self, trigger_snr):
    if not self.finalized:
      raise ValueError, "ExtrapolatingBackgroundEstimator has not been finalized"
    return self.fit_params[0] * numpy.exp(-self.fit_params[1]*trigger_snr/2)

##############################################################################
# Plotting classes

class InverseExpectationPlot(object):
  """
  Represent a plotsegments plot.  To use it, one must instantiate a
  InverseExpectationPlot object, call add_contents(), then call
  finalize().  To save, call the savefig() method.  To display to the screen,
  call pylab.show().
  """
  def __init__(self, fig_num=None):
    """
    Create a fresh plot.
    """
    self.fig = pylab.figure(fig_num)
    self.ax = self.fig.add_subplot(111)
    
    self.slide_stats_list = []
    self.bg_estimators = []
    self.label_list = []
    self.expectation_list_list = []
    self.finalized = False
    
    self.ax.set_xscale("log")
    self.ax.set_yscale("log")
    self.ax.set_xlabel(r"$1/\langle N \rangle$")
    self.ax.set_ylabel(r"Cumulative \#")
  
  def add_content(self, trigs, slide_trigs, label=None,
                  zero_lag_playground=False):
    """
    Add the contents to the plot.  The expectation value of each trigger in
    trigs will be computed from the slide_trigs.  At this point, a list of
    these expectation values will be stored here and then plotted only at
    finalize time.
    
    The ability to mix sets of data is a key design feature of this plot.
    """
    if self.finalized:
      raise ValueError, "cannot add data after plot has been finalized"
    
    # determine bg_fg_ratio
    bg_fg_ratio = len(slide_trigs)
    if zero_lag_playground:
      bg_fg_ratio *= 600/6370
    
    # store background stats
    self.slide_stats_list.append([s.getstat() for s in slide_trigs])
    
    # prepare background estimate
    bg_est = BackgroundEstimator(iterutils.flatten(self.slide_stats_list[-1]),
                                 bg_fg_ratio=bg_fg_ratio)
    bg_est.finalize()
    self.bg_estimators.append(bg_est)
    
    # evaluate background rate for each foreground trigger
    self.expectation_list_list.append([bg_est(s) for s in trigs.getstat()])
    self.label_list.append(label)
  
  def finalize(self, min_exp=None, max_exp=None):
    """
    Make final changes to the plot with the guarantee that no additional data
    will be added.
    """
    if self.finalized:
      raise ValueError, "InverseExpectationPlot has already been finalized"
    
    if sum([len(l) for l in self.expectation_list_list]) == 0:
      print >>sys.stderr, "warning: nothing to plot"
      return
    
    # sort
    for expectations in self.expectation_list_list:
      expectations.sort()

    # get the max and min over all data sets
    if (min_exp is None) or (max_exp is None):
      internal_min_exp, internal_max_exp = self._get_extreme_expectations()
      if min_exp is None:
        min_exp = internal_min_exp
      if max_exp is None:
        max_exp = internal_max_exp
    
    colors = itertools.cycle(('b', 'g', 'r', 'c', 'm', 'y', 'k'))
    symbols = itertools.cycle(('^', 'D', 'H', 'o', '1', '+'))
    
    # plot the points
    for expectations, label, color, symbol in \
        zip(self.expectation_list_list, self.label_list, colors, symbols):
      expectations.sort()
      num_zeros = bisect.bisect_right(expectations, 0)
      
      x = 1 / numpy.array(expectations[num_zeros:], dtype=float)
      y = numpy.arange(len(expectations) - num_zeros, dtype=int) + num_zeros + 1
      
      if num_zeros == len(expectations):
        continue
      
      # plot the non-zero rates
      self.ax.plot(x, y, color + symbol, label=label + " (N = %d)" % len(y),
                   markeredgecolor=color, alpha=0.85,
                   scalex=False, scaley=False)
      
      # handle zeros
      if num_zeros > 0:
        self.ax.plot([x[0]], [num_zeros], color + symbol,
                     label=label + " larger than all BG",
                     markerfacecolor='w', markeredgecolor=color,
                     scalex=False, scaley=False)
    
    # plot the background
    x = numpy.logspace(numpy.log10(1 / max_exp), numpy.log10(1 / min_exp), 3)
    y = 1 / x
    self.ax.plot(x, y, 'k--', label="background",
                 scalex=False, scaley=False)
    
    # set plot limits
    self.ax.set_xlim((0.9 / max_exp, 1.1 / min_exp))
    self.ax.set_ylim((0.6,
      1.1 * max([len(l) for l in self.expectation_list_list])))
    
    # add legend
    self.ax.legend()
    
    # mark as finalized
    self.savefig = self.fig.savefig
    self.finalized = True

  def _get_extreme_expectations(self):
    """
    Return min and max expectation numbers (across all data sets) as a tuple.
    Assumes that expectation numbers are sorted!
    """
    internal_min_exp = numpy.inf
    internal_max_exp = -numpy.inf
    for expectations in self.expectation_list_list:
      expectations.sort()
      num_zeros = bisect.bisect_right(expectations, 0)
      if num_zeros == len(expectations):
        continue
      internal_min_exp = min(internal_min_exp, expectations[num_zeros])
      internal_max_exp = max(internal_max_exp, expectations[-1])
    return internal_min_exp, internal_max_exp

  def savefig(self, *args, **kwargs):
    """
    Stub.  Will be replaced by self.fig.savefig after being finalized.
    """
    raise ValueError, "InverseExpectationPlot has not been finalized"

  def close(self):
    """
    Release the memory from the plot, which can be considerable.
    """
    pylab.close(self.fig)

  def __del__(self):
    self.close()


class HistogramInverseExpectationPlot(InverseExpectationPlot):
  """
  A histogrammed version of the InverseExpectationPlot, designed to be more
  like the output of viz.cumhiststat.
  """
  def __init__(self, *args, **kwargs):
    InverseExpectationPlot.__init__(self, *args, **kwargs)
    if "num_bins" in kwargs:
      self.num_bins = kwargs["num_bins"]
    else:
      self.num_bins = 20

  def finalize(self, min_exp=None, max_exp=None, num_bins=20):
    if self.finalized:
      raise ValueError, "HistogramInverseExpectationPlot has already been finalized"

    if sum([len(l) for l in self.expectation_list_list]) == 0:
      print >>sys.stderr, "warning: nothing to plot"
      return

    # sort
    for expectations in self.expectation_list_list:
      expectations.sort()
    
    # get the max and min over all data sets
    if (min_exp is None) or (max_exp is None):
      internal_min_exp, internal_max_exp = self._get_extreme_expectations()
      if min_exp is None:
        min_exp = internal_min_exp
      if max_exp is None:
        max_exp = internal_max_exp

    # plot the points
    histograms = []
    bins = rate.NDBins((rate.LogarithmicBins(1/max_exp, 1/min_exp, num_bins),))
    x = bins.centres()[0]

    colors = itertools.cycle(('b', 'g', 'r', 'c', 'm', 'y', 'k'))
    symbols = itertools.cycle(('^', 'D', 'H', 'o', '1', '+'))
    
    for expectations, label, color, symbol in \
        zip(self.expectation_list_list, self.label_list, colors, symbols):
      num_zeros = bisect.bisect_right(expectations, 0)
      if num_zeros == len(expectations):
        continue
      
      hist = rate.BinnedArray(bins, dtype=float)
      for value in expectations[num_zeros:]:
        hist[1/value,] += 1
      y = hist.array[::-1].cumsum()[::-1]
      y[y == 0] = epsilon
      assert y[0] == len(expectations)
      
      # plot the non-zero rates
      self.ax.plot(x, y, color + symbol, label=label,
                   markeredgecolor=color, alpha=0.65,
                   scalex=False, scaley=False)
      
      # handle zeros
      if num_zeros > 0:
        self.ax.plot([x[-1]], [num_zeros], color + symbol,
                     label=label + " larger than all BG",
                     markerfacecolor='w', markeredgecolor=color, markersize=2,
                     scalex=False, scaley=False)

    # shade background error bars
    # for each data set, histogram each slide separately and take stats
    colors = itertools.cycle(('b', 'g', 'r', 'c', 'm', 'y', 'k'))
    for slide_stats, bg_estimator, color in \
      zip(self.slide_stats_list, self.bg_estimators, colors):
      hist_sum = numpy.zeros(len(x), dtype=float)
      sq_hist_sum = numpy.zeros(len(x), dtype=float)
      for slide in slide_stats:
        hist = rate.BinnedArray(bins, dtype=float)
        for value in slide:
          try:
            hist[1/bg_estimator(value),] += 1
          except IndexError:
            # the histogram bins were chosen based on foreground extrema, so
            # background may exceed limits
            pass
        y = hist.array[::-1].cumsum()[::-1]
        hist_sum += y
        sq_hist_sum += y*y
      
      # get statistics
      N = len(slide_stats)
      means = hist_sum / N
      stds = numpy.sqrt((sq_hist_sum - hist_sum*means) / (N - 1))
      
      # shade in the area
      upper = means + stds
      lower = means - stds
      upper[upper <= epsilon] = epsilon
      lower[lower <= epsilon] = epsilon
      means[means <= epsilon] = epsilon
      tmp_x, tmp_y = viz.makesteps(bins.lower()[0], upper, lower)
      self.ax.fill(tmp_x, tmp_y, facecolor=color, alpha=0.2)
    
    # plot the background
    x = numpy.logspace(numpy.log10(1/max_exp), numpy.log10(1/min_exp), 200)
    y = 1 / x
    self.ax.plot(x, y, 'k--', label="background",
                 scalex=False, scaley=False)

    # set plot limits
    self.ax.set_xlim((0.9 / max_exp, 1.1 / min_exp))
    self.ax.set_ylim((0.6,
      1.15 * max([len(l) for l in self.expectation_list_list])))

    # add legend
    self.ax.legend()

    # mark as finalized
    self.savefig = self.fig.savefig
    self.finalized = True