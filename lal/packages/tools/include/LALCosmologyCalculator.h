/* Copyright (C) 2012 Walter Del Pozzo, Tjonnie Li
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
 
#ifndef LALCOSMOLOGYCALCULATOR_H
#define LALCOSMOLOGYCALCULATOR_H

#ifdef  __cplusplus   /* C++ protection. */
extern "C" {
#endif

typedef struct tagLALCosmologicalParameters
{
    double h;
    double om;
    double ol;
    double ok;
    double w0;
    double w1;
    double w2;
}   LALCosmologicalParameters;


double XLALLuminosityDistance(
            LALCosmologicalParameters *omega, 
            double z);

double XLALAngularDistance(
            LALCosmologicalParameters *omega, 
            double z);

double XLALComovingLOSDistance(
            LALCosmologicalParameters *omega, 
            double z);

double XLALComovingTransverseDistance(
            LALCosmologicalParameters *omega, 
            double z);            

double XLALHubbleDistance(
            LALCosmologicalParameters *omega
            );
            
double XLALHubbleParameter(double z,
            void *omega
            );
            
double XLALIntegrateHubbleParameter(
            LALCosmologicalParameters *omega, 
            double z);

double XLALUniformComovingVolumeDensity(
            double z,
            void *omega);

double XLALUniformComovingVolumeDistribution(
            LALCosmologicalParameters *omega, 
            double z,
            double zmin,
            double zmax);

double XLALIntegrateComovingVolumeDensity(LALCosmologicalParameters *omega, double zmin, double zmax);
            
LALCosmologicalParameters *XLALCreateCosmologicalParameters(double h, double om, double ok, double ol, double w0, double w1, double w2);

void XLALDestroyCosmologicalParameters(LALCosmologicalParameters *omega);

#endif
