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


from pylal import git_version
from pylal.xlal.burstsearch import *
from pylal.xlal.burstsearch import XLALEPGetTimingParameters as __XLALEPGetTimingParameters


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date


#
# =============================================================================
#
#                              Timing Parameters
#
# =============================================================================
#


class TimingParameters(object):
	"""
	A class to hold timing parameter values.
	"""
	pass


def XLALEPGetTimingParameters(window_length, max_tile_length, tile_stride_fraction, psd_length = None):
	"""
	Construct and populate a TimingParameters class instance
	initialized from some input parameter values.
	"""
	#
	# Check input
	#

	if not isinstance(window_length, int):
		raise TypeError, window_length

	if not isinstance(max_tile_length, int):
		raise TypeError, max_tile_length

	if psd_length is not None and not isinstance(psd_length, int):
		raise TypeError, psd_length

	#
	# initialize TimingParameters object
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
		params.__dict__.update(__XLALEPGetTimingParameters(window_length, max_tile_length, tile_stride_fraction))
	else:
		params.__dict__.update(__XLALEPGetTimingParameters(window_length, max_tile_length, tile_stride_fraction, psd_length))

	#
	# done
	#

	return params
