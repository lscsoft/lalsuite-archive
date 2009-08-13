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

__author__ = "Nickolas Fotopoulos <nvf@gravity.phys.uwm.edu>"
__version__ = "$Revision$"[11:-2]
__date__ = "$Date$"

from glue.segments import segment
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
    
    if t0>0:
      self.ax.set_xlabel("time (s)  offset: %9d"%t0)
    else:
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
    self.window = segment((a - padding, b + padding))
    
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
