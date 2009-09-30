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


from glue.ligolw import ligolw
from glue.ligolw import param
from pylal.xlal.datatypes.lalunit import LALUnit
from pylal.xlal.datatypes.ligotimegps import LIGOTimeGPS
from pylal.xlal.datatypes.real8frequencyseries import REAL8FrequencySeries
from pylal.xlal.datatypes.real8timeseries import REAL8TimeSeries


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
__version__ = "$Revision$"[11:-2]
__date__ = "$Date$"[7:-2]


#
# =============================================================================
#
#                                     Body
#
# =============================================================================
#


def parse_REAL8FrequencySeries(elem):
	t, = elem.getElementsByTagName(ligolw.Time.tagName)
	a, = elem.getElementsByTagName(ligolw.Array.tagName)
	dims = a.getElementsByTagName(ligolw.Dim.tagName)
	series = REAL8FrequencySeries(
		name = a.getAttribute("Name"),
		epoch = LIGOTimeGPS(t.pcdata),
		f0 = param.get_pyvalue(elem, "f0"),
		deltaF = float(dims[0].getAttribute("Scale")),
		sampleUnits = LALUnit(a.getAttribute("Unit")),
		data = a.array[1]
	)
	return series


def parse_REAL8TimeSeries(elem):
	t, = elem.getElementsByTagName(ligolw.Time.tagName)
	a, = elem.getElementsByTagName(ligolw.Array.tagName)
	dims = a.getElementsByTagName(ligolw.Dim.tagName)
	series = REAL8TimeSeries(
		name = a.getAttribute("Name"),
		epoch = LIGOTimeGPS(t.pcdata),
		f0 = param.get_pyvalue(elem, "f0"),
		deltaT = float(dims[0].getAttribute("Scale")),
		sampleUnits = LALUnit(a.getAttribute("Unit"))
		data = a.array[1]
	)
	return series


def seriesiter(xmldoc, laltype = "REAL8TimeSeries"):
	"""
	Iterate over the REAL8TimeSeries objects encoded in a LIGO Light
	Weight XML file.
	"""
	for elem in xmldoc.getElementsByTagName(ligolw.LIGO_LW.tagName):
		try:
			if elem.getAttribute("Name") != laltype:
				continue
		except KeyError:
			# this LIGO_LW doesn't have a Name attribute
			continue

		yield parse_REAL8TimeSeries(elem)
