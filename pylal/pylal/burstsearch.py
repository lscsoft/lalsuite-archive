# $Id$
#
# Copyright (C) 2007  Kipp Cannon
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
This module is a wrapper of the xlal.burstsearch module, supplementing the
C code in that module with additional features that are more easily
implemented in Python.  It is recommended that you import this module
rather than importing xlal.burstsearch directly.
"""


from xlal.burstsearch import *
from xlal.burstsearch import XLALEPGetTimingParameters as __XLALEPGetTimingParameters


__author__ = "Kipp Cannon <kipp@gravity.phys.uwm.edu>"
__date__ = "$Date$"[7:-2]
__version__ = "$Revision$"[11:-2]


#
# =============================================================================
#
#                              Timing Parameters
#
# =============================================================================
#


class TimingParameters(object):
	"""
	A place-holder class to hold timing parameter values.
	"""
	pass


def XLALEPGetTimingParameters(window_length, max_tile_length, tile_stride_fraction, psd_length = None):
	"""
	Construct and populate a TimingParameters class instance
	initialized from some input parameter values.
	"""
	#
	# init TimingParameters object
	#

	params = TimingParameters()
	params.window_length = window_length
	params.max_tile_length = max_tile_length
	params.tile_stride_fraction = tile_stride_fraction
	if psd_length is not None:
		params.psd_length = psd_length

	#
	# populate computed parameters from library code
	#

	if psd_length is None:
		params.update(__XLALEPGetTimingParameters(window_length, max_tile_length, tile_stride_fraction))
	else:
		params.update(__XLALEPGetTimingParameters(window_length, max_tile_length, tile_stride_fraction, psd_length))

	#
	# done
	#

	return params
