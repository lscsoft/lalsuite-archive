# Copyright (C) 2006--2011,2013  Kipp Cannon
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
Obsolete module.  Do not use.
"""


from lal import cached_detector_by_prefix, cached_detector_by_name, name_to_prefix, prefix_to_name
from pylal import git_version
from pylal.snglcoinc import light_travel_time


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date


cached_detector = cached_detector_by_name


# FIXME:  this is a hack to allow inject.light_travel_time(), which is used
# internally by snglcoinc.py, to work with the H1H2 coherent and null
# detectors.  the use of light_travel_time() internally by snglcoinc.py
# might be inappropriate and should be re-considered.  in the meantime,
# this will get the sub-solar mass search working.
prefix_to_name["H1H2"] = "LHO_4k"
