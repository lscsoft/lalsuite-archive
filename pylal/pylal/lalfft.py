# Copyright (C) 2009  Kipp Cannon
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
This module is a wrapper of the xlal.fft module, supplementing the C
code in that module with additional features that are more easily
implemented in Python.  It is recommended that you import this module
rather than importing xlal.fft directly.
"""


import git_version
from xlal.datatypes.real8fftplan import REAL8FFTPlan
from xlal.fft import *


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date


#
# =============================================================================
#
#                              Function Wrappers
#
# =============================================================================
#


def XLALCreateForwardREAL8FFTPlan(size, measurelvl = 0):
	return REAL8FFTPlan(size, 1, measurelvl)


def XLALCreateReverseREAL8FFTPlan(size, measurelvl = 0):
	return REAL8FFTPlan(size, 0, measurelvl)
