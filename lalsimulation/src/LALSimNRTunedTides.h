/*
* Copyright (C) 2017 Tim Dietrich, Sebastiano Bernuzzi, Nathan Johnson-McDaniel, Shasvath J Kapadia, Francesco Pannarale and Sebastian Khan
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

static double NRTunedTidesFDTidalPhase(
    const REAL8 PN_x, /**< PN frequency parameter: PN_x = orb_freq^(2./3.) */
    const REAL8 PN_x_2,
    const REAL8 PN_x_3over2,
    const REAL8 PN_x_5over2,
    const REAL8 Xa, /**< Mass of companion 1 divided by total mass */
    const REAL8 Xb, /**< Mass of companion 2 divided by total mass */
    const REAL8 kappa2T /**< tidal coupling constant. Eq. 2 in arXiv:1706.02969 */
);
