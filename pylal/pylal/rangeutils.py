# Copyright (C) 2012 Duncan M. Macleod
# 
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
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

"""
This module provides a bunch of user-friendly wrappers to the SWIG-bound REAL8TimeSeries and REAL8FrequencySeries objects and their associated functions.
"""

from __future__ import division

import numpy
from scipy import integrate
import lal
from pylal import git_version

__author__  = "Duncan M. Macleod <duncan.macleod@ligo.org>"
__version__ = git_version.id
__date__    = git_version.date

# =============================================================================
# Inspiral range
# =============================================================================

def inspiralrange(f, S, rho=8, m1=1.4, m2=1.4, fmin=10, fmax=None,\
                  horizon=False):
    """
    Calculate the sensitivity distance to an inspiral with the given masses,
    for the given signal-to-noise ratio.

    See the following reference for information on the integral performed:

    https://dcc.ligo.org/cgi-bin/private/DocDB/ShowDocument?docid=27267
    """
    # compute chirp mass and symmetric mass ratio
    mtot   = m1+m2
    mchirp = (m1*m2)**(0.6) / mtot**(0.4)
  
    # get mchirp in solar masses and compute integral prefactor
    mchirp *= lal.LAL_MSUN_SI
    pre = (5 * lal.LAL_C_SI**(1/3) *\
           (mchirp * lal.LAL_G_SI / lal.LAL_C_SI**2)**(5/3) * 1.77**2) /\
          (96 * lal.LAL_PI ** (4/3) * rho**2)


    # compute ISCO
    fisco = lal.LAL_C_SI**3\
            / (lal.LAL_G_SI * lal.LAL_MSUN_SI * 6**1.5 * lal.LAL_PI * mtot)
    if not fmax:
        fmax = fisco
    elif fmax > fisco:
        warnings.warn("Upper frequency bound greater than %s-%s ISCO "\
                      "frequency of %.2g Hz, using ISCO" % (m1,m2,fisco))
        fmax = fisco

    # integrate
    condition = (f >= fmin) & (f < fmax)
    integrand = f[condition]**(-7/3)/S[condition]
    result = integrate.trapz(integrand, f[condition])

    d = (pre*result) ** (1/2) / (lal.LAL_PC_SI*1e6)
    return d

# =============================================================================
# Burst range
# =============================================================================

def fdependent_burstrange(f, S, rh0=8, E=1e-2):
    """
    Calculate the sensitive distance to a GW burst with the given intrinsic
    energy for the given signal-to-noise ratio rho, as a function of f.
    """
    A = ((lal.LAL_G_SI * E * lal.LAL_MSUN_SI * 0.4)\
         / (lal.LAL_PI**2 * lal.LAL_C_SI))**(1/2) / lal.LAL_PC_SI * 1e6
    return A / (rho * S**(1/2) * frequency)

def burstrange(f, S, rho=8, E=1e-2, fmin=0, fmax=None):
    """
    Calculate the sensitive distance to a GW burst with the given intrinsic
    energy for the given signal-to-noise ratio rho, integrated over frequency.
    """
    if not fmin: 
        fmin = f.min()
    if not fmax:
        fmax = f.max()

    # restrict to band
    condition = (f >= fmin) & (f < fmax)

    # integrate
    integrand = fdependent_burstrange(f[condition], S[condition], rho, E)**3
    result = spectrum.deltaF*integrand.sum()
    result = integrate.trapz(integrand, f[condition])
    
    d = (result / (fmax-fmin))**1/3
    return d
