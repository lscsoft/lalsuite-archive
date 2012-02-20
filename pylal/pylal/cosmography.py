# Copyright (C) 2011 Nickolas Fotopoulos
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
"""
Cosmography based on Hogg 1998 (astro-ph/9905116)
Cosmological parameters taken from WMAP 7 year BAO + H0 ML value (arXiv:1001.4538)
"""

from numpy import pi, sinh, arcsinh
from scipy.integrate import quad
from scipy.optimize import newton

from pylal.xlal.constants import LAL_PC_SI, LAL_C_SI

h0 = 70.3  # km / s / Mpc
H0_SI = h0 * 1000 / (1000000 * LAL_PC_SI)  # Hubble constant in inverse seconds
DH_SI = LAL_C_SI / H0_SI  # Hubble distance in meters
OmegaM = 0.1338 / (h0 * h0)  # fraction of crit. density in matter
OmegaL = 0.729  # fraction of crit. density in dark energy
OmegaR = 1 - OmegaM - OmegaL  # remainder

def Einv(z):
    """
    Integrand used in many cosmography integrals. 1 / E(z) in Hogg's notation.
    """
    return (OmegaM * (1 + z)**3 + OmegaR * (1 + z)**2 + OmegaL)**-0.5

def DM(z):
    """
    Comoving distance (transverse) at z. OmegaR > 0.
    """
    return DH_SI * OmegaR**-0.5 * sinh(OmegaR**0.5 * quad(Einv, 0, z)[0])

def DL(z):
    """
    Luminosity distance
    """
    return DM(z) * (1 + z)

def compute_redshift(DL0):
    """
    Compute the redshift from a luminosity distance.
    Use the Newton-Raphson refinement method on a first-order guess.
    """
    return newton(lambda z: DL(z) - DL0, DL0 / DH_SI)

def dVdz(z):
    """
    Different volume element per unit redshift.
    """
    return DH_SI * DM(z)**2 * Einv(z) * 4 * pi

def V(z):
    """
    Analytical integration of dVdz given in Carroll's astro-ph/0004075.
    Double precision craps out below about 100 kpc (z=2.3e-6).
    """
    DM_DH = OmegaR**-0.5 * sinh(OmegaR**0.5 * quad(Einv, 0, z)[0]) # DM / DH
    return 2 * pi * DH_SI * DH_SI * DH_SI / OmegaR * \
        (DM_DH * (1 + OmegaR * DM_DH * DM_DH)**0.5 \
                  - OmegaR**-0.5 * arcsinh(OmegaR**0.5 * DM_DH))

def V_from_DL_z(D, z=None):
    """
    Analytical integration of dVdz given in Carroll's astro-ph/0004075,
    sped up for the case that you have D rather than (or in addition to) z.
    Double precision craps out below about 100 kpc (z=2.3e-6).
    """
    DM_DH = D / (1 + (z or compute_redshift(D))) / DH_SI
    return 2 * pi * DH_SI * DH_SI * DH_SI / OmegaR * \
        (DM_DH * (1 + OmegaR * DM_DH * DM_DH)**0.5 \
                  - OmegaR**-0.5 * arcsinh(OmegaR**0.5 * DM_DH))
