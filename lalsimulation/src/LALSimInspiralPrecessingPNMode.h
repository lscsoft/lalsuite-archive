/*
 * Copyright (C) 2013 Evan Ochsner
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

#ifndef _LALSIMINSPIRALPRECESSINGPNMODE_H
#define _LALSIMINSPIRALPRECESSINGPNMODE_H

#include <lal/TimeSeries.h>

int XLALSimInspiralRadiationFrameToJFrame(
        REAL8TimeSeries **S1x_J,       /**< [returned] J-frame time series */
        REAL8TimeSeries **S1y_J,       /**< [returned] J-frame time series */
        REAL8TimeSeries **S1z_J,       /**< [returned] J-frame time series */
        REAL8TimeSeries **S2x_J,       /**< [returned] J-frame time series */
        REAL8TimeSeries **S2y_J,       /**< [returned] J-frame time series */
        REAL8TimeSeries **S2z_J,       /**< [returned] J-frame time series */
        REAL8TimeSeries **LNhatx_J,    /**< [returned] J-frame time series */
        REAL8TimeSeries **LNhaty_J,    /**< [returned] J-frame time series */
        REAL8TimeSeries **LNhatz_J,    /**< [returned] J-frame time series */
        const REAL8TimeSeries *V,      /**< post-Newtonian parameter */
        const REAL8TimeSeries *S1x,    /**< Spin1 vector x component */
        const REAL8TimeSeries *S1y,    /**< Spin1 vector y component */
        const REAL8TimeSeries *S1z,    /**< Spin1 vector z component */
        const REAL8TimeSeries *S2x,    /**< Spin2 vector x component */
        const REAL8TimeSeries *S2y,    /**< Spin2 vector y component */
        const REAL8TimeSeries *S2z,    /**< Spin2 vector z component */
        const REAL8TimeSeries *LNhatx, /**< unit orbital ang. mom. x comp. */
        const REAL8TimeSeries *LNhaty, /**< unit orbital ang. mom. y comp. */
        const REAL8TimeSeries *LNhatz, /**< unit orbital ang. mom. z comp. */
        const REAL8 m1,                /**< mass of companion 1 (kg) */
        const REAL8 m2,                /**< mass of companion 2 (kg) */
        const REAL8 fref               /**< frame is set by J at this GW freq.*/
        );

COMPLEX16TimeSeries *XLALSimInspiralComputePrecessingPNMode(
        REAL8TimeSeries *V,      /**< post-Newtonian parameter */
        REAL8TimeSeries *Phi,    /**< orbital phase (rad) */
        REAL8TimeSeries *S1x,    /**< spin 1 x-component in J-frame */
        REAL8TimeSeries *S1y,    /**< spin 1 y-component in J-frame */
        REAL8TimeSeries *S1z,    /**< spin 1 z-component in J-frame */
        REAL8TimeSeries *S2x,    /**< spin 2 x-component in J-frame */
        REAL8TimeSeries *S2y,    /**< spin 2 y-component in J-frame */
        REAL8TimeSeries *S2z,    /**< spin 2 z-component in J-frame */
        REAL8TimeSeries *LNhatx, /**< \hat{L}_N x-component in J-frame */
        REAL8TimeSeries *LNhaty, /**< \hat{L}_N y-component in J-frame */
        REAL8TimeSeries *LNhatz, /**< \hat{L}_N z-component in J-frame */
        REAL8 m1,                /**< mass1 (kg) */
        REAL8 m2,                /**< mass2 (kg) */
        REAL8 r,                 /**< distance to binary (m) */
        REAL8 ampO,              /**< twice PN amplitude order (3 == 1.5PN) */
        int l,                   /**< mode number l */
        int m                    /**< mode number m */
        );

COMPLEX16TimeSeries *XLALSimInspiralPrecessingPNMode22(
        const REAL8TimeSeries *V,      /**< post-Newtonian parameter */
        const REAL8TimeSeries *Phi,    /**< orbital phase */
        const REAL8TimeSeries *S1x,    /**< Spin1 vector x component */
        const REAL8TimeSeries *S1y,    /**< Spin1 vector y component */
        const REAL8TimeSeries *S1z,    /**< Spin1 vector z component */
        const REAL8TimeSeries *S2x,    /**< Spin2 vector x component */
        const REAL8TimeSeries *S2y,    /**< Spin2 vector y component */
        const REAL8TimeSeries *S2z,    /**< Spin2 vector z component */
        const REAL8TimeSeries *LNhatx, /**< unit orbital ang. mom. x comp. */
        const REAL8TimeSeries *LNhaty, /**< unit orbital ang. mom. y comp. */
        const REAL8TimeSeries *LNhatz, /**< unit orbital ang. mom. z comp. */
        const REAL8 m1,                /**< mass of companion 1 (kg) */
        const REAL8 m2,                /**< mass of companion 2 (kg) */
        const REAL8 r,                 /**< distance of source (m) */
        const int O                    /**< twice post-Newtonian order */
        );

COMPLEX16TimeSeries *XLALSimInspiralPrecessingPNMode2m2(
        const REAL8TimeSeries *V,      /**< post-Newtonian parameter */
        const REAL8TimeSeries *Phi,    /**< orbital phase */
        const REAL8TimeSeries *S1x,    /**< Spin1 vector x component */
        const REAL8TimeSeries *S1y,    /**< Spin1 vector y component */
        const REAL8TimeSeries *S1z,    /**< Spin1 vector z component */
        const REAL8TimeSeries *S2x,    /**< Spin2 vector x component */
        const REAL8TimeSeries *S2y,    /**< Spin2 vector y component */
        const REAL8TimeSeries *S2z,    /**< Spin2 vector z component */
        const REAL8TimeSeries *LNhatx, /**< unit orbital ang. mom. x comp. */
        const REAL8TimeSeries *LNhaty, /**< unit orbital ang. mom. y comp. */
        const REAL8TimeSeries *LNhatz, /**< unit orbital ang. mom. z comp. */
        const REAL8 m1,                /**< mass of companion 1 (kg) */
        const REAL8 m2,                /**< mass of companion 2 (kg) */
        const REAL8 r,                 /**< distance of source (m) */
        const int O                    /**< twice post-Newtonian order */
        );

COMPLEX16TimeSeries *XLALSimInspiralPrecessingPNMode21(
        const REAL8TimeSeries *V,      /**< post-Newtonian parameter */
        const REAL8TimeSeries *Phi,    /**< orbital phase */
        const REAL8TimeSeries *S1x,    /**< Spin1 vector x component */
        const REAL8TimeSeries *S1y,    /**< Spin1 vector y component */
        const REAL8TimeSeries *S1z,    /**< Spin1 vector z component */
        const REAL8TimeSeries *S2x,    /**< Spin2 vector x component */
        const REAL8TimeSeries *S2y,    /**< Spin2 vector y component */
        const REAL8TimeSeries *S2z,    /**< Spin2 vector z component */
        const REAL8TimeSeries *LNhatx, /**< unit orbital ang. mom. x comp. */
        const REAL8TimeSeries *LNhaty, /**< unit orbital ang. mom. y comp. */
        const REAL8TimeSeries *LNhatz, /**< unit orbital ang. mom. z comp. */
        const REAL8 m1,                /**< mass of companion 1 (kg) */
        const REAL8 m2,                /**< mass of companion 2 (kg) */
        const REAL8 r,                 /**< distance of source (m) */
        const int O                    /**< twice post-Newtonian order */
        );

COMPLEX16TimeSeries *XLALSimInspiralPrecessingPNMode2m1(
        const REAL8TimeSeries *V,      /**< post-Newtonian parameter */
        const REAL8TimeSeries *Phi,    /**< orbital phase */
        const REAL8TimeSeries *S1x,    /**< Spin1 vector x component */
        const REAL8TimeSeries *S1y,    /**< Spin1 vector y component */
        const REAL8TimeSeries *S1z,    /**< Spin1 vector z component */
        const REAL8TimeSeries *S2x,    /**< Spin2 vector x component */
        const REAL8TimeSeries *S2y,    /**< Spin2 vector y component */
        const REAL8TimeSeries *S2z,    /**< Spin2 vector z component */
        const REAL8TimeSeries *LNhatx, /**< unit orbital ang. mom. x comp. */
        const REAL8TimeSeries *LNhaty, /**< unit orbital ang. mom. y comp. */
        const REAL8TimeSeries *LNhatz, /**< unit orbital ang. mom. z comp. */
        const REAL8 m1,                /**< mass of companion 1 (kg) */
        const REAL8 m2,                /**< mass of companion 2 (kg) */
        const REAL8 r,                 /**< distance of source (m) */
        const int O                    /**< twice post-Newtonian order */
        );

COMPLEX16TimeSeries *XLALSimInspiralPrecessingPNMode20(
        const REAL8TimeSeries *V,      /**< post-Newtonian parameter */
        const REAL8TimeSeries *Phi,    /**< orbital phase */
        const REAL8TimeSeries *S1x,    /**< Spin1 vector x component */
        const REAL8TimeSeries *S1y,    /**< Spin1 vector y component */
        const REAL8TimeSeries *S1z,    /**< Spin1 vector z component */
        const REAL8TimeSeries *S2x,    /**< Spin2 vector x component */
        const REAL8TimeSeries *S2y,    /**< Spin2 vector y component */
        const REAL8TimeSeries *S2z,    /**< Spin2 vector z component */
        const REAL8TimeSeries *LNhatx, /**< unit orbital ang. mom. x comp. */
        const REAL8TimeSeries *LNhaty, /**< unit orbital ang. mom. y comp. */
        const REAL8TimeSeries *LNhatz, /**< unit orbital ang. mom. z comp. */
        const REAL8 m1,                /**< mass of companion 1 (kg) */
        const REAL8 m2,                /**< mass of companion 2 (kg) */
        const REAL8 r,                 /**< distance of source (m) */
        const int O                    /**< twice post-Newtonian order */
        );

COMPLEX16TimeSeries *XLALSimInspiralPrecessingPNMode33(
        const REAL8TimeSeries *V,      /**< post-Newtonian parameter */
        const REAL8TimeSeries *Phi,    /**< orbital phase */
        const REAL8TimeSeries *S1x,    /**< Spin1 vector x component */
        const REAL8TimeSeries *S1y,    /**< Spin1 vector y component */
        const REAL8TimeSeries *S1z,    /**< Spin1 vector z component */
        const REAL8TimeSeries *S2x,    /**< Spin2 vector x component */
        const REAL8TimeSeries *S2y,    /**< Spin2 vector y component */
        const REAL8TimeSeries *S2z,    /**< Spin2 vector z component */
        const REAL8TimeSeries *LNhatx, /**< unit orbital ang. mom. x comp. */
        const REAL8TimeSeries *LNhaty, /**< unit orbital ang. mom. y comp. */
        const REAL8TimeSeries *LNhatz, /**< unit orbital ang. mom. z comp. */
        const REAL8 m1,                /**< mass of companion 1 (kg) */
        const REAL8 m2,                /**< mass of companion 2 (kg) */
        const REAL8 r,                 /**< distance of source (m) */
        const int O                    /**< twice post-Newtonian order */
        );

COMPLEX16TimeSeries *XLALSimInspiralPrecessingPNMode3m3(
        const REAL8TimeSeries *V,      /**< post-Newtonian parameter */
        const REAL8TimeSeries *Phi,    /**< orbital phase */
        const REAL8TimeSeries *S1x,    /**< Spin1 vector x component */
        const REAL8TimeSeries *S1y,    /**< Spin1 vector y component */
        const REAL8TimeSeries *S1z,    /**< Spin1 vector z component */
        const REAL8TimeSeries *S2x,    /**< Spin2 vector x component */
        const REAL8TimeSeries *S2y,    /**< Spin2 vector y component */
        const REAL8TimeSeries *S2z,    /**< Spin2 vector z component */
        const REAL8TimeSeries *LNhatx, /**< unit orbital ang. mom. x comp. */
        const REAL8TimeSeries *LNhaty, /**< unit orbital ang. mom. y comp. */
        const REAL8TimeSeries *LNhatz, /**< unit orbital ang. mom. z comp. */
        const REAL8 m1,                /**< mass of companion 1 (kg) */
        const REAL8 m2,                /**< mass of companion 2 (kg) */
        const REAL8 r,                 /**< distance of source (m) */
        const int O                    /**< twice post-Newtonian order */
        );

COMPLEX16TimeSeries *XLALSimInspiralPrecessingPNMode32(
        const REAL8TimeSeries *V,      /**< post-Newtonian parameter */
        const REAL8TimeSeries *Phi,    /**< orbital phase */
        const REAL8TimeSeries *S1x,    /**< Spin1 vector x component */
        const REAL8TimeSeries *S1y,    /**< Spin1 vector y component */
        const REAL8TimeSeries *S1z,    /**< Spin1 vector z component */
        const REAL8TimeSeries *S2x,    /**< Spin2 vector x component */
        const REAL8TimeSeries *S2y,    /**< Spin2 vector y component */
        const REAL8TimeSeries *S2z,    /**< Spin2 vector z component */
        const REAL8TimeSeries *LNhatx, /**< unit orbital ang. mom. x comp. */
        const REAL8TimeSeries *LNhaty, /**< unit orbital ang. mom. y comp. */
        const REAL8TimeSeries *LNhatz, /**< unit orbital ang. mom. z comp. */
        const REAL8 m1,                /**< mass of companion 1 (kg) */
        const REAL8 m2,                /**< mass of companion 2 (kg) */
        const REAL8 r,                 /**< distance of source (m) */
        const int O                    /**< twice post-Newtonian order */
        );

COMPLEX16TimeSeries *XLALSimInspiralPrecessingPNMode3m2(
        const REAL8TimeSeries *V,      /**< post-Newtonian parameter */
        const REAL8TimeSeries *Phi,    /**< orbital phase */
        const REAL8TimeSeries *S1x,    /**< Spin1 vector x component */
        const REAL8TimeSeries *S1y,    /**< Spin1 vector y component */
        const REAL8TimeSeries *S1z,    /**< Spin1 vector z component */
        const REAL8TimeSeries *S2x,    /**< Spin2 vector x component */
        const REAL8TimeSeries *S2y,    /**< Spin2 vector y component */
        const REAL8TimeSeries *S2z,    /**< Spin2 vector z component */
        const REAL8TimeSeries *LNhatx, /**< unit orbital ang. mom. x comp. */
        const REAL8TimeSeries *LNhaty, /**< unit orbital ang. mom. y comp. */
        const REAL8TimeSeries *LNhatz, /**< unit orbital ang. mom. z comp. */
        const REAL8 m1,                /**< mass of companion 1 (kg) */
        const REAL8 m2,                /**< mass of companion 2 (kg) */
        const REAL8 r,                 /**< distance of source (m) */
        const int O                    /**< twice post-Newtonian order */
        );

COMPLEX16TimeSeries *XLALSimInspiralPrecessingPNMode31(
        const REAL8TimeSeries *V,      /**< post-Newtonian parameter */
        const REAL8TimeSeries *Phi,    /**< orbital phase */
        const REAL8TimeSeries *S1x,    /**< Spin1 vector x component */
        const REAL8TimeSeries *S1y,    /**< Spin1 vector y component */
        const REAL8TimeSeries *S1z,    /**< Spin1 vector z component */
        const REAL8TimeSeries *S2x,    /**< Spin2 vector x component */
        const REAL8TimeSeries *S2y,    /**< Spin2 vector y component */
        const REAL8TimeSeries *S2z,    /**< Spin2 vector z component */
        const REAL8TimeSeries *LNhatx, /**< unit orbital ang. mom. x comp. */
        const REAL8TimeSeries *LNhaty, /**< unit orbital ang. mom. y comp. */
        const REAL8TimeSeries *LNhatz, /**< unit orbital ang. mom. z comp. */
        const REAL8 m1,                /**< mass of companion 1 (kg) */
        const REAL8 m2,                /**< mass of companion 2 (kg) */
        const REAL8 r,                 /**< distance of source (m) */
        const int O                    /**< twice post-Newtonian order */
        );

COMPLEX16TimeSeries *XLALSimInspiralPrecessingPNMode3m1(
        const REAL8TimeSeries *V,      /**< post-Newtonian parameter */
        const REAL8TimeSeries *Phi,    /**< orbital phase */
        const REAL8TimeSeries *S1x,    /**< Spin1 vector x component */
        const REAL8TimeSeries *S1y,    /**< Spin1 vector y component */
        const REAL8TimeSeries *S1z,    /**< Spin1 vector z component */
        const REAL8TimeSeries *S2x,    /**< Spin2 vector x component */
        const REAL8TimeSeries *S2y,    /**< Spin2 vector y component */
        const REAL8TimeSeries *S2z,    /**< Spin2 vector z component */
        const REAL8TimeSeries *LNhatx, /**< unit orbital ang. mom. x comp. */
        const REAL8TimeSeries *LNhaty, /**< unit orbital ang. mom. y comp. */
        const REAL8TimeSeries *LNhatz, /**< unit orbital ang. mom. z comp. */
        const REAL8 m1,                /**< mass of companion 1 (kg) */
        const REAL8 m2,                /**< mass of companion 2 (kg) */
        const REAL8 r,                 /**< distance of source (m) */
        const int O                    /**< twice post-Newtonian order */
        );

COMPLEX16TimeSeries *XLALSimInspiralPrecessingPNMode30(
        const REAL8TimeSeries *V,      /**< post-Newtonian parameter */
        const REAL8TimeSeries *Phi,    /**< orbital phase */
        const REAL8TimeSeries *S1x,    /**< Spin1 vector x component */
        const REAL8TimeSeries *S1y,    /**< Spin1 vector y component */
        const REAL8TimeSeries *S1z,    /**< Spin1 vector z component */
        const REAL8TimeSeries *S2x,    /**< Spin2 vector x component */
        const REAL8TimeSeries *S2y,    /**< Spin2 vector y component */
        const REAL8TimeSeries *S2z,    /**< Spin2 vector z component */
        const REAL8TimeSeries *LNhatx, /**< unit orbital ang. mom. x comp. */
        const REAL8TimeSeries *LNhaty, /**< unit orbital ang. mom. y comp. */
        const REAL8TimeSeries *LNhatz, /**< unit orbital ang. mom. z comp. */
        const REAL8 m1,                /**< mass of companion 1 (kg) */
        const REAL8 m2,                /**< mass of companion 2 (kg) */
        const REAL8 r,                 /**< distance of source (m) */
        const int O                    /**< twice post-Newtonian order */
        );

COMPLEX16TimeSeries *XLALSimInspiralPrecessingPNMode44(
        const REAL8TimeSeries *V,      /**< post-Newtonian parameter */
        const REAL8TimeSeries *Phi,    /**< orbital phase */
        const REAL8TimeSeries *S1x,    /**< Spin1 vector x component */
        const REAL8TimeSeries *S1y,    /**< Spin1 vector y component */
        const REAL8TimeSeries *S1z,    /**< Spin1 vector z component */
        const REAL8TimeSeries *S2x,    /**< Spin2 vector x component */
        const REAL8TimeSeries *S2y,    /**< Spin2 vector y component */
        const REAL8TimeSeries *S2z,    /**< Spin2 vector z component */
        const REAL8TimeSeries *LNhatx, /**< unit orbital ang. mom. x comp. */
        const REAL8TimeSeries *LNhaty, /**< unit orbital ang. mom. y comp. */
        const REAL8TimeSeries *LNhatz, /**< unit orbital ang. mom. z comp. */
        const REAL8 m1,                /**< mass of companion 1 (kg) */
        const REAL8 m2,                /**< mass of companion 2 (kg) */
        const REAL8 r,                 /**< distance of source (m) */
        const int O                    /**< twice post-Newtonian order */
        );

COMPLEX16TimeSeries *XLALSimInspiralPrecessingPNMode4m4(
        const REAL8TimeSeries *V,      /**< post-Newtonian parameter */
        const REAL8TimeSeries *Phi,    /**< orbital phase */
        const REAL8TimeSeries *S1x,    /**< Spin1 vector x component */
        const REAL8TimeSeries *S1y,    /**< Spin1 vector y component */
        const REAL8TimeSeries *S1z,    /**< Spin1 vector z component */
        const REAL8TimeSeries *S2x,    /**< Spin2 vector x component */
        const REAL8TimeSeries *S2y,    /**< Spin2 vector y component */
        const REAL8TimeSeries *S2z,    /**< Spin2 vector z component */
        const REAL8TimeSeries *LNhatx, /**< unit orbital ang. mom. x comp. */
        const REAL8TimeSeries *LNhaty, /**< unit orbital ang. mom. y comp. */
        const REAL8TimeSeries *LNhatz, /**< unit orbital ang. mom. z comp. */
        const REAL8 m1,                /**< mass of companion 1 (kg) */
        const REAL8 m2,                /**< mass of companion 2 (kg) */
        const REAL8 r,                 /**< distance of source (m) */
        const int O                    /**< twice post-Newtonian order */
        );

COMPLEX16TimeSeries *XLALSimInspiralPrecessingPNMode43(
        const REAL8TimeSeries *V,      /**< post-Newtonian parameter */
        const REAL8TimeSeries *Phi,    /**< orbital phase */
        const REAL8TimeSeries *S1x,    /**< Spin1 vector x component */
        const REAL8TimeSeries *S1y,    /**< Spin1 vector y component */
        const REAL8TimeSeries *S1z,    /**< Spin1 vector z component */
        const REAL8TimeSeries *S2x,    /**< Spin2 vector x component */
        const REAL8TimeSeries *S2y,    /**< Spin2 vector y component */
        const REAL8TimeSeries *S2z,    /**< Spin2 vector z component */
        const REAL8TimeSeries *LNhatx, /**< unit orbital ang. mom. x comp. */
        const REAL8TimeSeries *LNhaty, /**< unit orbital ang. mom. y comp. */
        const REAL8TimeSeries *LNhatz, /**< unit orbital ang. mom. z comp. */
        const REAL8 m1,                /**< mass of companion 1 (kg) */
        const REAL8 m2,                /**< mass of companion 2 (kg) */
        const REAL8 r,                 /**< distance of source (m) */
        const int O                    /**< twice post-Newtonian order */
        );

COMPLEX16TimeSeries *XLALSimInspiralPrecessingPNMode4m3(
        const REAL8TimeSeries *V,      /**< post-Newtonian parameter */
        const REAL8TimeSeries *Phi,    /**< orbital phase */
        const REAL8TimeSeries *S1x,    /**< Spin1 vector x component */
        const REAL8TimeSeries *S1y,    /**< Spin1 vector y component */
        const REAL8TimeSeries *S1z,    /**< Spin1 vector z component */
        const REAL8TimeSeries *S2x,    /**< Spin2 vector x component */
        const REAL8TimeSeries *S2y,    /**< Spin2 vector y component */
        const REAL8TimeSeries *S2z,    /**< Spin2 vector z component */
        const REAL8TimeSeries *LNhatx, /**< unit orbital ang. mom. x comp. */
        const REAL8TimeSeries *LNhaty, /**< unit orbital ang. mom. y comp. */
        const REAL8TimeSeries *LNhatz, /**< unit orbital ang. mom. z comp. */
        const REAL8 m1,                /**< mass of companion 1 (kg) */
        const REAL8 m2,                /**< mass of companion 2 (kg) */
        const REAL8 r,                 /**< distance of source (m) */
        const int O                    /**< twice post-Newtonian order */
        );

COMPLEX16TimeSeries *XLALSimInspiralPrecessingPNMode42(
        const REAL8TimeSeries *V,      /**< post-Newtonian parameter */
        const REAL8TimeSeries *Phi,    /**< orbital phase */
        const REAL8TimeSeries *S1x,    /**< Spin1 vector x component */
        const REAL8TimeSeries *S1y,    /**< Spin1 vector y component */
        const REAL8TimeSeries *S1z,    /**< Spin1 vector z component */
        const REAL8TimeSeries *S2x,    /**< Spin2 vector x component */
        const REAL8TimeSeries *S2y,    /**< Spin2 vector y component */
        const REAL8TimeSeries *S2z,    /**< Spin2 vector z component */
        const REAL8TimeSeries *LNhatx, /**< unit orbital ang. mom. x comp. */
        const REAL8TimeSeries *LNhaty, /**< unit orbital ang. mom. y comp. */
        const REAL8TimeSeries *LNhatz, /**< unit orbital ang. mom. z comp. */
        const REAL8 m1,                /**< mass of companion 1 (kg) */
        const REAL8 m2,                /**< mass of companion 2 (kg) */
        const REAL8 r,                 /**< distance of source (m) */
        const int O                    /**< twice post-Newtonian order */
        );

COMPLEX16TimeSeries *XLALSimInspiralPrecessingPNMode4m2(
        const REAL8TimeSeries *V,      /**< post-Newtonian parameter */
        const REAL8TimeSeries *Phi,    /**< orbital phase */
        const REAL8TimeSeries *S1x,    /**< Spin1 vector x component */
        const REAL8TimeSeries *S1y,    /**< Spin1 vector y component */
        const REAL8TimeSeries *S1z,    /**< Spin1 vector z component */
        const REAL8TimeSeries *S2x,    /**< Spin2 vector x component */
        const REAL8TimeSeries *S2y,    /**< Spin2 vector y component */
        const REAL8TimeSeries *S2z,    /**< Spin2 vector z component */
        const REAL8TimeSeries *LNhatx, /**< unit orbital ang. mom. x comp. */
        const REAL8TimeSeries *LNhaty, /**< unit orbital ang. mom. y comp. */
        const REAL8TimeSeries *LNhatz, /**< unit orbital ang. mom. z comp. */
        const REAL8 m1,                /**< mass of companion 1 (kg) */
        const REAL8 m2,                /**< mass of companion 2 (kg) */
        const REAL8 r,                 /**< distance of source (m) */
        const int O                    /**< twice post-Newtonian order */
        );

COMPLEX16TimeSeries *XLALSimInspiralPrecessingPNMode41(
        const REAL8TimeSeries *V,      /**< post-Newtonian parameter */
        const REAL8TimeSeries *Phi,    /**< orbital phase */
        const REAL8TimeSeries *S1x,    /**< Spin1 vector x component */
        const REAL8TimeSeries *S1y,    /**< Spin1 vector y component */
        const REAL8TimeSeries *S1z,    /**< Spin1 vector z component */
        const REAL8TimeSeries *S2x,    /**< Spin2 vector x component */
        const REAL8TimeSeries *S2y,    /**< Spin2 vector y component */
        const REAL8TimeSeries *S2z,    /**< Spin2 vector z component */
        const REAL8TimeSeries *LNhatx, /**< unit orbital ang. mom. x comp. */
        const REAL8TimeSeries *LNhaty, /**< unit orbital ang. mom. y comp. */
        const REAL8TimeSeries *LNhatz, /**< unit orbital ang. mom. z comp. */
        const REAL8 m1,                /**< mass of companion 1 (kg) */
        const REAL8 m2,                /**< mass of companion 2 (kg) */
        const REAL8 r,                 /**< distance of source (m) */
        const int O                    /**< twice post-Newtonian order */
        );

COMPLEX16TimeSeries *XLALSimInspiralPrecessingPNMode4m1(
        const REAL8TimeSeries *V,      /**< post-Newtonian parameter */
        const REAL8TimeSeries *Phi,    /**< orbital phase */
        const REAL8TimeSeries *S1x,    /**< Spin1 vector x component */
        const REAL8TimeSeries *S1y,    /**< Spin1 vector y component */
        const REAL8TimeSeries *S1z,    /**< Spin1 vector z component */
        const REAL8TimeSeries *S2x,    /**< Spin2 vector x component */
        const REAL8TimeSeries *S2y,    /**< Spin2 vector y component */
        const REAL8TimeSeries *S2z,    /**< Spin2 vector z component */
        const REAL8TimeSeries *LNhatx, /**< unit orbital ang. mom. x comp. */
        const REAL8TimeSeries *LNhaty, /**< unit orbital ang. mom. y comp. */
        const REAL8TimeSeries *LNhatz, /**< unit orbital ang. mom. z comp. */
        const REAL8 m1,                /**< mass of companion 1 (kg) */
        const REAL8 m2,                /**< mass of companion 2 (kg) */
        const REAL8 r,                 /**< distance of source (m) */
        const int O                    /**< twice post-Newtonian order */
        );

COMPLEX16TimeSeries *XLALSimInspiralPrecessingPNMode40(
        const REAL8TimeSeries *V,      /**< post-Newtonian parameter */
        const REAL8TimeSeries *Phi,    /**< orbital phase */
        const REAL8TimeSeries *S1x,    /**< Spin1 vector x component */
        const REAL8TimeSeries *S1y,    /**< Spin1 vector y component */
        const REAL8TimeSeries *S1z,    /**< Spin1 vector z component */
        const REAL8TimeSeries *S2x,    /**< Spin2 vector x component */
        const REAL8TimeSeries *S2y,    /**< Spin2 vector y component */
        const REAL8TimeSeries *S2z,    /**< Spin2 vector z component */
        const REAL8TimeSeries *LNhatx, /**< unit orbital ang. mom. x comp. */
        const REAL8TimeSeries *LNhaty, /**< unit orbital ang. mom. y comp. */
        const REAL8TimeSeries *LNhatz, /**< unit orbital ang. mom. z comp. */
        const REAL8 m1,                /**< mass of companion 1 (kg) */
        const REAL8 m2,                /**< mass of companion 2 (kg) */
        const REAL8 r,                 /**< distance of source (m) */
        const int O                    /**< twice post-Newtonian order */
        );

#endif /* _LALSIMINSPIRALPRECESSINGPNMODE_H */
