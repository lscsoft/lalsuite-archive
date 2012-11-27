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
 
#include "LALCosmologyCalculator.h"

double XLALLuminosityDistance(
            LALCosmologicalParameters *omega, 
            double z)
{
    double x=1.0+z;
    return XLALAngularDistance(omega,z)*x*x;
}

double XLALAngularDistance(
            LALCosmologicalParameters *omega, 
            double z)
{
    double x=1.0+z;
    return XLALComovingTransverseDistance(omega, z)/x;
}

double XLALComovingLOSDistance(
            LALCosmologicalParameters *omega, 
            double z)
{
    double dH = XLALHubbleDistance(omega);
    return dH*XLALIntegrateHubbleParameter(omega,z);
}
            
double XLALComovingTransverseDistance(
            LALCosmologicalParameters *omega, 
            double z)
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

double XLALHubbleDistance(
            LALCosmologicalParameters *omega
            )
{
    double x = GSL_CONST_MKSA_SPEED_OF_LIGHT/omega->h;

    return 1e-5*x;
}

double XLALHubbleParameter(double z,
            void *omega
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

double XLALIntegrateHubbleParameter(LALCosmologicalParameters *omega, double z)
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

double XLALUniformComovingVolumeDistribution(
            LALCosmologicalParameters *omega, 
            double z,
            double zmax)
{

    double x = 1.0+z;
    
    double norm = XLALIntegrateComovingVolumeDensity(omega, zmax);
    double da = XLALAngularDistance(omega,z);
    double E = XLALHubbleParameter(z,(void *)omega);
          
    double unnorm_density = da*da/(x*E);
    
    return unnorm_density/norm;
}

double XLALUniformComovingVolumeDensity( 
            double z,
            void *omega)
{
    
    LALCosmologicalParameters *p = (LALCosmologicalParameters *)omega;

    double x = 1.0+z;
    double da = XLALAngularDistance(p,z);
    double E = XLALHubbleParameter(z,omega);
    double unnorm_density = da*da/(x*E);
    
    return unnorm_density;
}

double XLALIntegrateComovingVolumeDensity(LALCosmologicalParameters *omega, double z)
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

    gsl_integration_qag (&F, 0.0, z, epsabs, epsrel, 
                    limit, key, w, &result, &error);

    gsl_integration_workspace_free (w);

    return result;
}

LALCosmologicalParameters *XLALFillCosmologicalParameters(double h, double om, double ok, double ol, double w0, double w1, double w2)
{
    LALCosmologicalParameters *p = (LALCosmologicalParameters *)malloc(sizeof(LALCosmologicalParameters));
    p->h = h;
    p->om = om;
    p->ok = ok;
    p->ol=ol;
    p->w0=w0;
    p->w1=w1;
    p->w2=w2;
    return p;
}
