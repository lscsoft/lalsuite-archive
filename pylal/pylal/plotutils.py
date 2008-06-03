#!/usr/bin/env python
#
# Copyright (C) 2008  Nickolas Fotopoulos
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
#
"""
This module is intended to store generic, reusable, sub-classable plot classes
to minimize formulaic copying and pasting.
"""

__author__ = "Nickolas Fotopoulos <nvf@gravity.phys.uwm.edu>"
__version__ = "$Revision$"[11:-2]
__date__ = "$Date$"


itertools = __import__("itertools")  # absolute import of system-wide itertools

import numpy
import pylab

from glue import iterutils

from pylal import viz

##############################################################################
# abstract classes

class BasicPlot(object):
    """
    A very default meta-class to almost any plot you might want to make.
    It provides basic initialization, a savefig method, and a close method.
    It is up to developers to subclass BasicPlot and fill in the add_content()
    and finalize() methods.
    """
    def __init__(self, xlabel="", ylabel="", title=""):
        """
        Basic plot initialization.  A subclass can override __init__ and call
        this one (plotutils.BasicPlot.__init__(self, *args, **kwargs)) and
        then initialize variables to hold data to plot and labels.
        """
        self.fig = pylab.figure()
        self.ax = self.fig.add_subplot(111)

        self.ax.set_xlabel(xlabel)
        self.ax.set_ylabel(ylabel)
        self.ax.set_title(title)

        self.ax.grid(True)

    def add_content(self, data, label="_nolabel_"):
        """
        Stub.  Replace with a method that appends values or lists of values
        to self.data_sets and appends labels to self.data_labels.  Feel free
        to accept complicated inputs, but try to store only the raw numbers
        that will enter the plot.
        """
        raise NotImplemented

    def finalize(self):
        """
        Stub.  Replace with a function that creates and makes your plot
        pretty.  Do not do I/O here.
        """
        raise NotImplemented

    def savefig(self, *args, **kwargs):
        self.fig.savefig(*args, **kwargs)

    def close(self):
        """
        Close the plot and release its memory.
        """
        pylab.close(self.fig)

    def add_legend_if_labels_exist(self, *args, **kwargs):
        """
        Create a legend if there are any non-trivial labels.
        """
        for label in self.data_labels:
            if not (label is None or label.startswith("_")):
                self.ax.legend(*args, **kwargs)
                return


##############################################################################
# utility functions

def default_colors():
    return itertools.cycle(('b', 'g', 'r', 'c', 'm', 'y', 'k'))

def default_symbols():
    return itertools.cycle(('x', '^', 'D', 'H', 'o', '1', '+'))

def determine_common_bin_limits(data_sets, default_min=0, default_max=0):
    """
    Given a some nested sequences (e.g. list of lists), determine the largest
    and smallest values over the data sets and determine a common binning.
    """
    max_stat = max(list(iterutils.flatten(data_sets)) + [-numpy.inf])
    min_stat = min(list(iterutils.flatten(data_sets)) + [numpy.inf])
    if numpy.isinf(-max_stat):
        max_stat = default_max
    if numpy.isinf(min_stat):
        min_stat = default_min
    return min_stat, max_stat


##############################################################################
# generic, but usable classes

class VerticalBarPlot(BasicPlot):
    """
    A simple vertical bar plot.  Bars are centered on the x values and have
    height equal to the y values.
    """
    def __init__(self, *args, **kwargs):
        BasicPlot.__init__(self, *args, **kwargs)
        self.x_data_sets = []
        self.y_data_sets = []
        self.data_labels = []

    def add_content(self, x_data, y_data, label="_nolabel_"):
        self.x_data_sets.append(x_data)
        self.y_data_sets.append(y_data)
        self.data_labels.append(label)

    def finalize(self):
        # make plot
        colors = default_colors()

        for x_vals, y_vals, color, label in \
            itertools.izip(self.x_data_sets, self.y_data_sets, colors,
                           self.data_labels):
            self.ax.bar(x_vals, y_vals, color=color, label=label,
                        align="center", linewidth=0)

        # add legend if there are any non-trivial labels
        self.add_legend_if_labels_exist()

        # decrement reference counts
        del self.x_data_sets
        del self.y_data_sets
        del self.data_labels


class VerticalBarHistogram(BasicPlot):
    """
    Histogram data sets with a common binning, then make a vertical bar plot.
    """
    def __init__(self, *args, **kwargs):
        BasicPlot.__init__(self, *args, **kwargs)
        self.data_sets = []
        self.data_labels = []

    def add_content(self, data, label="_nolabel_"):
        self.data_sets.append(data)
        self.data_labels.append(label)

    def finalize(self, num_bins=20):
        # determine binning
        min_stat, max_stat = determine_common_bin_limits(self.data_sets)
        bins = numpy.linspace(min_stat, max_stat, num_bins)

        # determine bar width; gets silly for more than a few data sets
        width = (1 - 0.1 * len(self.data_sets)) * max_stat / num_bins

        # make plot
        colors = default_colors()
        count = itertools.count()
        for i, data_set, color, label in \
            itertools.izip(count, self.data_sets, colors, self.labels):
            # make histogram
            y, x = numpy.histogram(stats, bins=bins)

            # stagger bins for pure aesthetics
            x += 0.1 * i * max_stat / num_bins

            # plot
            self.ax.bar(x, y, color, label=label, width=width, alpha=0.6,
                        align="center")

        # add legend if there are any non-trivial labels
        self.add_legend_if_labels_exist()

        # decrement reference counts
        del self.data_sets
        del self.data_labels

class CumulativeHistogramPlot(BasicPlot):
    """
    Cumulative histogram of foreground that also has a shaded region,
    determined by the mean and standard deviation of the background
    population coincidence statistics.
    """
    def __init__(self, *args, **kwargs):
        BasicPlot.__init__(self, *args, **kwargs)
        self.fg_data_sets = []
        self.fg_labels = []
        self.bg_data_sets = []
        self.bg_label = "_nolegend_"

    def add_content(self, fg_data_set, label="_nolegend_"):
        self.fg_data_sets.append(fg_data_set)
        self.fg_labels.append(label)

    def add_background(self, bg_data_sets, label="_nolegend_"):
        self.bg_data_sets.extend(bg_data_sets)
        self.bg_label = label

    def finalize(self, num_bins=20, normalization=1):
        epsilon = 1e-8

        # determine binning
        min_stat, max_stat = determine_common_bin_limits(\
            self.fg_data_sets + self.bg_data_sets)
        bins = numpy.linspace(min_stat, max_stat, num_bins)

        # plot foreground
        colors = default_colors()
        symbols = default_symbols()
        for data_set, color, symbol, label in \
            itertools.izip(self.fg_data_sets, colors, symbols,
            self.fg_labels):
            # make histogram
            y, x = numpy.histogram(data_set, bins=bins)
            y = y[::-1].cumsum()[::-1]

            # plot
            y = numpy.array(y, dtype=numpy.float32)
            y[y <= epsilon] = epsilon
            self.ax.plot(x, y*normalization, symbol + color, label=label)

        # shade background region
        if len(self.bg_data_sets) > 0:
            # histogram each background instance separately and take stats
            hist_sum = numpy.zeros(len(bins), dtype=float)
            sq_hist_sum = numpy.zeros(len(bins), dtype=float)
            for instance in self.bg_data_sets:
                # make histogram
                y, x = numpy.histogram(instance, bins=bins)
                y = y[::-1].cumsum()[::-1]
                hist_sum += y
                sq_hist_sum += y*y
          
            # get statistics
            N = len(self.bg_data_sets)
            means = hist_sum / N
            stds = numpy.sqrt((sq_hist_sum - hist_sum*means) / (N - 1))
          
            # plot mean
            self.ax.plot(x, means*normalization, 'r+',
                label=r"$\mu_\mathrm{%s}$" % self.bg_label)

            # shade in the area
            upper = means + stds
            lower = means - stds
            upper[upper <= epsilon] = epsilon
            lower[lower <= epsilon] = epsilon
            means[means <= epsilon] = epsilon
            tmp_x, tmp_y = viz.makesteps(bins, upper, lower)
            self.ax.fill(tmp_x, tmp_y*normalization, facecolor='y', alpha=0.3,
                label=r"$\sigma_\mathrm{%s}$" % self.bg_label)

        # make semilogy plot
        self.ax.set_yscale("log")

        # adjust plot range
        self.ax.set_xlim((0.9 * min_stat, 1.1 * max_stat))
        if len(self.bg_data_sets) > 0 and len(self.fg_data_sets) == 0:
            self.ax.set_ylim(ymin=0.6 / N)
        else:
            self.ax.set_ylim(ymin=0.6 * normalization)

        # add legend if there are any non-trivial labels
        self.data_labels = self.fg_labels + [self.bg_label]
        self.add_legend_if_labels_exist()

        # decrement reference counts
        del self.fg_data_sets
        del self.bg_data_sets
        del self.fg_labels
        del self.bg_label
        del self.data_labels

