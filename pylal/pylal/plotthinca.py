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
import operator
import warnings

import numpy
import pylab

import scipy.optimize

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
    
    print self.fit_params
    
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
    
    self.label_list = []
    self.expectation_list_list = []
    self.finalized = False
    
    self.ax.set_xscale("log")
    self.ax.set_yscale("log")
    self.ax.set_xlabel(r"$1/\langle N \rangle$")
    self.ax.set_ylabel(r"Cumulative \#")
  
  def add_content(self, trigs, slide_trigs, label=None):
    """
    Add the contents to the plot.  The expectation value of each trigger in
    trigs will be computed from the slide_trigs.  At this point, a list of
    these expectation values will be stored here and then plotted only at
    finalize time.
    
    The ability to mix sets of data is a key design feature of this plot.
    """
    if self.finalized:
      raise ValueError, "cannot add data after plot has been finalized"
    
    # prepare background estimate
    bg_estimator = BackgroundEstimator(bg_fg_ratio=len(slide_trigs))
    for slides in slide_trigs:
      bg_estimator.extend(slides.getstat())
    bg_estimator.finalize()
    
    # evaluate background rate for each foreground trigger
    self.expectation_list_list.append([bg_estimator(s) for s in trigs.getstat()])
    self.label_list.append(label)
  
  def finalize(self, min_exp=None, max_exp=None):
    """
    Make final changes to the plot with the guarantee that no additional data
    will be added.
    """
    if self.finalized:
      raise ValueError, "InverseExpectationPlot has already been finalized"
    colors = itertools.cycle(['b', 'g', 'r', 'c', 'm', 'y', 'k'])
    symbols = itertools.cycle(['^', 'D', 'H', 'o', '1', '+'])
    
    # plot the points, keeping track of maxes and mins
    internal_min_exp = numpy.inf
    internal_max_exp = -numpy.inf
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

      # contribute to plot limits
      internal_min_exp = min(internal_min_exp, expectations[num_zeros])
      internal_max_exp = max(internal_max_exp, expectations[-1])
      
      # handle zeros
      if num_zeros > 0:
        self.ax.plot([x[0]], [num_zeros], color + symbol,
                     label=label + " larger than all BG",
                     markerfacecolor='w', markeredgecolor=color,
                     scalex=False, scaley=False)
    
    # determine plot limits
    if min_exp is None:
      min_exp = internal_min_exp
    if max_exp is None:
      max_exp = internal_max_exp
    
    # plot the background
    x = numpy.logspace(numpy.log10(1 / max_exp), numpy.log10(1 / min_exp), 200)
    y = 1 / x
    self.ax.plot(x, y, 'k--', label="background",
                 scalex=False, scaley=False)
    
    # set plot limits
    self.ax.set_xlim((0.9 / max_exp, 1.1 / min_exp))
    self.ax.set_ylim((0.9,
      1.1 * max([len(l) for l in self.expectation_list_list])))
    
    # add legend
    self.ax.legend()
    
    # mark as finalized
    self.savefig = self.fig.savefig
    self.finalized = True

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
