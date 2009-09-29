# $Id$
#
# Copyright (C) 2008  Kipp C. Cannon
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
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#


"""
Code to assist in reading and writing LAL time- and frequency series data
encoded in LIGO Light-Weight XML format.
"""


import numpy
import os
import sys


from glue.ligolw import ligolw
from glue.ligolw import array
from glue.ligolw import param
from pylal.xlal.datatypes.ligotimegps import LIGOTimeGPS


__author__ = "Kipp Cannon <kipp@gravity.phys.uwm.edu>"
__version__ = "$Revision$"[11:-2]
__date__ = "$Date$"[7:-2]


#
# =============================================================================
#
#                               REAL8TimeSeries
#
# =============================================================================
#


class REAL8TimeSeries(object):
	# FIXME:  this should be an extension module that exposes LAL's
	# REAL8TimeSeries type, but I can't be bothered yet.
	def __init__(self, name, epoch, f0, deltaT, sampleUnits, n):
		self.name = name
		self.epoch = LIGOTimeGPS(epoch)
		self.f0 = f0
		self.deltaT = deltaT
		self.sampleUnits = sampleUnits
		self.data = numpy.zeros((n,), dtype = "double")


#
# =============================================================================
#
#                                     Body
#
# =============================================================================
#


def seriesiter(xmldoc, laltype = "REAL8TimeSeries"):
	"""
	Iterate over the REAL8TimeSeries objects encoded in a LIGO Light
	Weight XML file.
	"""
	for llw in xmldoc.getElementsByTagName(ligolw.LIGO_LW.tagName):
		try:
			if llw.getAttribute("Name") != laltype:
				continue
		except KeyError:
			# this LIGO_LW doesn't have a Name attribute
			continue

		[t] = llw.getElementsByTagName(ligolw.Time.tagName)
		[a] = llw.getElementsByTagName(ligolw.Array.tagName)
		dims = a.getElementsByTagName(ligolw.Dim.tagName)
		series = REAL8TimeSeries(
			a.getAttribute("Name"),
			LIGOTimeGPS(t.pcdata),
			param.get_pyvalue(llw, "f0"),
			float(dims[0].getAttribute("Scale")),
			a.getAttribute("Unit"),
			0
		)
		series.data = a.array[1]

		yield series
