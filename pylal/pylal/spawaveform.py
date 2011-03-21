# Copyright (C) 2010  Kipp Cannon,  Drew Keppel
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
This module is a wrapper of the _spawaveform module, supplementing the C
code in that module with additional features that are more easily
implemented in Python.  It is recommended that you import this module
rather than importing _spawaveform directly.
"""


import math


from pylal import git_version
from pylal.xlal import constants as lalconstants
from pylal._spawaveform import *


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date


#
# =============================================================================
#
#                                    Empty
#
# =============================================================================
#


def eta(m1, m2):
	"""
	Compute the symmetric mass ratio, eta.
	"""
	return m1*m2/(m1+m2)**2.


def chirpmass(m1, m2):
	"""
	Compute the chirp mass in seconds.
	"""
	return lalconstants.LAL_MTSUN_SI * (m1+m2) * eta(m1, m2)**.6


def ms2taus(m1, m2, f0 = 40.0):
	"""
	Solve for tau_0 and tau_3 from m1 and m2.
	"""
	tau0 = 5./256./(math.pi*f0)**(8./3.) * chirpmass(m1,m2)**(-5./3.)
	tau3 = math.pi/8./eta(m1,m2)**.6/(math.pi*f0)**(5./3.) * chirpmass(m1,m2)**(-2./3.)
	return tau0, tau3


def taus2ms(tau0, tau3, f0 = 40.0):
	"""
	Solve for m1 and m2 from tau_0 and tau_3.
	"""
	Mc = (5./256./(math.pi*f0)**(8./3.) / tau0)**(3./5.)
	eta = (math.pi/8./(math.pi*f0)**(5./3.) / tau3 / Mc**(2./3.))**(5./3.)

	M = Mc / eta**(3./5.)

	m1 = (1. + abs(1. - 4.*eta)**.5) * M / 2.
	m2 = M - m1

	return m1 / lalconstants.LAL_MTSUN_SI, m2 / lalconstants.LAL_MTSUN_SI
