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

#include <stdio.h>
#include <stdlib.h> 
#include <string.h>
#include <math.h>
#include <lal/LALMalloc.h>
#include <lal/LALConstants.h>
#include <gsl/gsl_integration.h>
#include "LALCosmologyCalculator.h"

/** 
* The set of functions in this module implements the standard cosmological distance measures 
* defined a Friedmann-Robertson-Walker-Lemaitre metric.
* For a detailed reference to the various formulae, see Hogg 1999 ( http://arxiv.org/abs/astro-ph/9905116 ). 
**/ 

/** 
* Computes the luminosity distance as the angular distance divided by (1+z)^2
* Eq. 21 in Hogg 1999 ( http://arxiv.org/abs/astro-ph/9905116 ) 
**/
double XLALLuminosityDistance(
            LALCosmologicalParameters *omega, /** structure holding the cosmological parameters **/
            double z /** redshift **/
            )
{
    double x=1.0+z;
    return XLALAngularDistance(omega,z)*x*x;
}

/** 
* Computes the angular size distance by dividing the comoving transverse distance by 1+z
* Eq. 18 in Hogg 1999 ( http://arxiv.org/abs/astro-ph/9905116 )
**/
 
double XLALAngularDistance(
            LALCosmologicalParameters *omega, /** structure holding the cosmological parameters **/
            double z /** redshift **/
            )
{
    double x=1.0+z;
    return XLALComovingTransverseDistance(omega, z)/x;
}

/** 
* Computes the comoving line-of-sight distance (usually called the proper distance) by integrating the Hubble parameter.  
* Eq. 15 in Hogg 1999 ( http://arxiv.org/abs/astro-ph/9905116 ).
**/
double XLALComovingLOSDistance(
            LALCosmologicalParameters *omega, /** structure holding the cosmological parameters **/
            double z /** redshift **/
            )
{
    double dH = XLALHubbleDistance(omega);
    return dH*XLALIntegrateHubbleParameter(omega,z);
}

/**
* Computes the comoving transverse distance from the comoving line-of-sight distance (proper distance) and
* Eq. 16 in Hogg 1999 ( http://arxiv.org/abs/astro-ph/9905116 )
**/

double XLALComovingTransverseDistance(
            LALCosmologicalParameters *omega, /** structure holding the cosmological parameters **/
            double z /** redshift **/
            )
{
    double ok = omega->ok;
    double dH = XLALHubbleDistance(omega);
    double dC = XLALComovingLOSDistance(omega,z);
    
    if (fabs(ok)<1e-14) 
    {
        return dC;
    }
    else if (ok>0)
    {
        return dH*sinh(sqrt(ok)*dC/dH)/sqrt(ok);
    }
    else if (ok<0) 
    {
        return dH*sin(sqrt(fabs(ok))*dC/dH)/sqrt(fabs(ok));
    }
    else 
    {
        fprintf(stderr,"Something funny happened. Aborting.\n");
        exit(-1);
    }

}

/**
* Computes the Hubble distance, independent of the redshift. This is the distance with sets the untis for all others.
* It returns the distance in Mpc. 
* Eq. 4 in Hogg 1999 ( http://arxiv.org/abs/astro-ph/9905116 )
**/

double XLALHubbleDistance(
            LALCosmologicalParameters *omega /** structure holding the cosmological parameters */
            )
{
    double x = LAL_C_SI/omega->h;

    return 1e-5*x;
}

/**
* Computes the inverse of the Hubble parameter at redshift z 
* Eq. 14 in Hogg 1999 ( http://arxiv.org/abs/astro-ph/9905116 )
**/

double XLALHubbleParameter(double z, /** redshift */
            void *omega /** structure holding the cosmological parameters, voided for integration purposes **/
            )
{
    LALCosmologicalParameters *p = (LALCosmologicalParameters *) omega;

    double E=0.0;
    double x = 1.0+z;
    double om=p->om;
    double ol=p->ol;
    double ok=p->ok;
    double w0=p->w0;
    double w1=p->w1;
    double w2=p->w2;
    w0=-1.0;
    w1=0.0;
    w2=0.0;
    E = sqrt(om*x*x*x+ok*x*x+ol*pow(x,3.*(1.0+w0+w1+w2))*exp(-3.0*((w1+w2)*z/x + w2*z*z/(2.0*x*x))));
    
    return  1.0/E; 
}

/**
* Computes the integral of inverse of the Hubble parameter at redshift z.
* The integral is computed using the built-in function gsl_integration_qag of the gsl library. 
* The integral is performed using a Gauss-Kronrod adaptive method (see http://www.gnu.org/software/gsl/manual/html_node/QAG-adaptive-integration.html ) 
**/

double XLALIntegrateHubbleParameter(
            LALCosmologicalParameters *omega, /** structure holding the cosmological parameters **/
            double z /** redshift **/
            )
{
    double result = 0.0;
    double error;
    double epsabs = 5e-3;
    double epsrel = 1e-3;
    size_t limit = 512;
    int key = 1;
    
    gsl_function F;
    F.function = &XLALHubbleParameter;
    F.params  = omega;
    
    gsl_integration_workspace * w 
    = gsl_integration_workspace_alloc (512);

    gsl_integration_qag (&F, 0.0, z, epsabs, epsrel, 
                    limit, key, w, &result, &error);

    gsl_integration_workspace_free (w);

    return result;
}

/** 
* This function computes the value of a uniform probability distribution over the comoving volume.
* Details of the derivation of these formulae can be found in Coward, Burman 2005 ( http://arxiv.org/abs/astro-ph/0505181 )
**/

double XLALUniformComovingVolumeDistribution(
            LALCosmologicalParameters *omega, /** structure holding the cosmological parameters **/
            double z, /** redshift **/
            double zmax /** UNUSED maximal redshift **/
            )
{
    zmax+=1;
    double norm = XLALIntegrateComovingVolumeDensity(omega, -1.0);
    double unnorm_density = XLALUniformComovingVolumeDensity(z,(void *)omega);

    return unnorm_density/norm;
}

/** 
* This function computes the value of a uniform probability density distribution over the comoving volume.
**/

double XLALUniformComovingVolumeDensity( 
            double z, /** redshift **/
            void *omega /** structure holding the cosmological parameters **/
            )
{
    
    LALCosmologicalParameters *p = (LALCosmologicalParameters *)omega;

    double x = 1.0+z;
    double dc = XLALComovingLOSDistance(p,z);
    double E = XLALHubbleParameter(z,omega);
    double unnorm_density = dc*dc*E/(x);
    
    return unnorm_density;
}

/**
* Function that integrates the uniform in comoving volume density to compute the normalisation factor. 
* Consistently with Coward, Burman 2005 ( http://arxiv.org/abs/astro-ph/0505181 ) the integral is always performed up to infinity.
* The integration is performed using gsl_integration_qagiu from the gsl library (see http://www.gnu.org/software/gsl/manual/html_node/QAGI-adaptive-integration-on-infinite-intervals.html )
**/

double XLALIntegrateComovingVolumeDensity(
            LALCosmologicalParameters *omega, /** structure holding the cosmological parameters **/
            double z /** redshift **/
            )
{
    double result = 0.0;
    double error;
    double epsabs = 5e-3;
    double epsrel = 1e-3;
    size_t limit = 512;
    int key = 1;
    
    gsl_function F;
    F.function = &XLALUniformComovingVolumeDensity;
    F.params  = omega;
    
    gsl_integration_workspace * w 
    = gsl_integration_workspace_alloc (512);

    if (z<0.0) gsl_integration_qagiu (&F, 0.0, epsabs, epsrel, 
                    limit, w, &result, &error);
    
    else gsl_integration_qag (&F, 0.0, z, epsabs, epsrel, 
                    limit, key, w, &result, &error);

    gsl_integration_workspace_free (w);

    return result;
}

/**
* Creates a LALCosmologicalParameters structure from the values of the cosmological parameters.
* Note that the boundary condition \Omega_m + \Omega_k + \Omega_\Lambda = 1 is NOT imposed here, 
* but must be set prior to populating this structure and calling the required functions. 
* Failure to impose the correct buondary condition will result in the Hubble parameter to be ill-defined 
* and thus the code to crash.  
**/

LALCosmologicalParameters *XLALCreateCosmologicalParameters(
                            double h, /** hubble constant in normalised units H0/100 Mpc s */
                            double om, /** matter energy density */
                            double ok, /** curvature ''energy density'' */
                            double ol, /** cosmological constant energy density */
                            double w0, /** zeroth order correction for a (putative) dark enegy equation of state */
                            double w1, /** first order correction for a (putative) dark enegy equation of state */
                            double w2  /** second order correction for a (putative) dark enegy equation of state */
                            ) 
{
    LALCosmologicalParameters *p = (LALCosmologicalParameters *)XLALMalloc(sizeof(LALCosmologicalParameters));
    p->h = h;
    p->om = om;
    p->ok = ok;
    p->ol=ol;
    p->w0=w0;
    p->w1=w1;
    p->w2=w2;
    return p;
}

/**
* Destroys a LALCosmologicalParameters structure.
**/

void XLALDestroyCosmologicalParameters(
                    LALCosmologicalParameters *omega /** structure holding the cosmological parameters */
                    )
{
    XLALFree(omega);
}
