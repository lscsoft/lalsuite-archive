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


import numpy


import git_version
from xlal.datatypes.complex16frequencyseries import COMPLEX16FrequencySeries
from xlal.datatypes.complex16fftplan import COMPLEX16FFTPlan
from xlal.datatypes.lalunit import lalSecondUnit
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


def XLALCreateForwardCOMPLEX16FFTPlan(size, measurelvl = 0):
	return COMPLEX16FFTPlan(size, 1, measurelvl)


def XLALCreateReverseCOMPLEX16FFTPlan(size, measurelvl = 0):
	return COMPLEX16FFTPlan(size, 0, measurelvl)


def XLALCreateForwardREAL8FFTPlan(size, measurelvl = 0):
	return REAL8FFTPlan(size, 1, measurelvl)


def XLALCreateReverseREAL8FFTPlan(size, measurelvl = 0):
	return REAL8FFTPlan(size, 0, measurelvl)


#
# =============================================================================
#
#                                  Utilities
#
# =============================================================================
#


def prepare_fseries_for_real8tseries(series):
	"""
	Construct a COMPLEX16FrequencySeries object suitable for storing
	the Fourier transform of a REAL8TimeSeries object.
	"""
	n = len(series.data)
	return COMPLEX16FrequencySeries(
		name = series.name,
		epoch = series.epoch,
		f0 = series.f0,	# note: non-zero f0 not supported by LAL
		deltaF = 1.0 / (n * series.deltaT),
		sampleUnits = series.sampleUnits * lalSecondUnit,
		data = numpy.zeros((n / 2 + 1,), dtype = "cdouble")
	)


def prepare_fseries_for_complex16tseries(series):
	"""
	Construct a COMPLEX16FrequencySeries object suitable for storing
	the Fourier transform of a COMPLEX16TimeSeries object.
	"""
	n = len(series.data)
	return COMPLEX16FrequencySeries(
		name = series.name,
		epoch = series.epoch,
		f0 = series.f0,	# note: non-zero f0 not supported by LAL
		deltaF = 1.0 / (n * series.deltaT),
		sampleUnits = series.sampleUnits * lalSecondUnit,
		data = numpy.zeros((n,), dtype = "cdouble")
	)
