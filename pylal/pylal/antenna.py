# $Id$
#
# Copyright (C) 2007   Alexander Dietz
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
This module provides functions to calculate antenna factors for a given time, a given sky location and a given detector
"""
import sys
from math import *
from pylal.xlal import date
from pylal import date
from pylal.xlal import inject

__author__ = "Alexander Dietz <Alexander.Dietz@astro.cf.ac.uk>"
__version__ = "$Revision$"[11:-2]
__date__ = "$Date$"[7:-2]

    


def response( gpsTime, ra_deg, dec_deg, iota_deg, psi_deg, det ):

  """
  Calculates the antenna factors for a detector 'det' (e.g. 'H1')
  at a given gps time (as integer) for a given sky location
  (ra_deg, dec_deg) [degree]. This computation also takes into account
  a specific inclination iota_deg and polarization psi_deg, all given
  in degree. 
  The returned values are: (f-plus, f-cross, f-average, q-value). 
  Example: antenna.response( 854378604.780, 11.089, 42.308, 0, 0, 'H1' )
  """

  # convert degree to radians
  ra_rad  = ra_deg * pi / 180.0
  dec_rad = dec_deg * pi / 180.0
  iota_rad= iota_deg * pi / 180.0
  psi_rad = psi_deg * pi / 180.0
  
  # calculate GMST if the GPS time
  gps=date.LIGOTimeGPS( gpsTime )
  gmst_rad = date.XLALGreenwichMeanSiderealTime(gps)

  # create detector-name map
  detMap = {'H1': 'LHO_4k', 'H2': 'LHO_2k', 'L1': 'LLO_4k',
            'G1': 'GEO_600', 'V1': 'VIRGO', 'T1': 'TAMA_300'}
  try:
    detector=detMap[det]
  except KeyError:
    print >>sys.stderr, "ERROR. Key %s is not a valid detector name."\
          % (det)
    sys.exit(0)

  # get detector
  if detector not in inject.cached_detector.keys():
    print >>sys.stderr, "%s is not a cached detector.  "\
          "Cached detectors are: %s" \
          % (det, inject.cached_detector.keys())
    sys.exit(1)

  # get the correct response data
  response = inject.cached_detector[detector].response

  # actual computation of antenna factors
  f_plus, f_cross = inject.XLALComputeDetAMResponse(response, ra_rad, dec_rad,
                                                    psi_rad, gmst_rad)

  f_ave=sqrt( (f_plus*f_plus + f_cross*f_cross)/2.0 );
  ci=cos( iota_rad );
  cc=ci*ci;

  # calculate q-value, e.g. ratio of effective to real distance
  # ref: Duncans PhD, eq. (4.3) on page 57 
  f_q=sqrt( f_plus*f_plus*(1+cc)*(1+cc)/4.0 + f_cross*f_cross*cc ); 

  # output
  return f_plus, f_cross, f_ave, f_q
