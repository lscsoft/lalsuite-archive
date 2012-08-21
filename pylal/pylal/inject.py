# Copyright (C) 2006--2011  Kipp Cannon
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
This module is a wrapper of the xlal.inject module, supplementing the C
code in that module with additional features that are more easily
implemented in Python.  It is recommended that you import this module
rather than importing xlal.inject directly.
"""


import math


from pylal import git_version
from pylal.xlal.tools import *
from pylal.xlal.inject import *
from pylal.xlal.constants import LAL_C_SI


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date


#
# =============================================================================
#
#                                Look-up Tables
#
# =============================================================================
#


name_to_prefix = dict((detector.name, detector.prefix) for detector in cached_detector.itervalues())


prefix_to_name = dict((detector.prefix, detector.name) for detector in cached_detector.itervalues())


#
# =============================================================================
#
#                              Function Wrappers
#
# =============================================================================
#


def light_travel_time(instrument1, instrument2):
	"""
	Compute and return the time required for light to travel through
	free space the distance separating the two instruments.  The inputs
	are two instrument prefixes (e.g., "H1"), and the result is
	returned in seconds.  Note how this differs from LAL's
	XLALLightTravelTime() function, which takes two detector objects as
	input, and returns the time truncated to integer nanoseconds.
	"""
	dx = cached_detector[prefix_to_name[instrument1]].location - cached_detector[prefix_to_name[instrument2]].location
	return math.sqrt((dx * dx).sum()) / LAL_C_SI


def effective_distance_factor(inclination, fp, fc):
	"""
	Returns the ratio of effective distance to physical distance for
	compact binary mergers.  Inclination is the orbital inclination of
	the system in radians, fp and fc are the F+ and Fx antenna factors.
	See XLALComputeDetAMResponse() for a function to compute antenna
	factors.  The effective distance is given by

	Deff = effective_distance_factor * D

	See Equation (4.3) of arXiv:0705.1514.
	"""
	cos2i = math.cos(inclination)**2
	return 1.0 / math.sqrt(fp**2 * (1+cos2i)**2 / 4 + fc**2 * cos2i)
