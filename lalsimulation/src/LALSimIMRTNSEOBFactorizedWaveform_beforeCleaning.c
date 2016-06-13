/*
*  Copyright (C) 2010 Craig Robinson, Yi Pan
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


/**
 * \author Swetha Bhagwat 
 *
 */

#include <math.h>
#include <complex.h>
#include "LALSimIMREOBNRv2.h"
/* Include static functions */
#include "LALSimInspiraldEnergyFlux.c"
#include "LALSimIMREOBNewtonianMultipole.c" 
#include "LALSimIMREOBNQCCorrection.c"
#include "LALSimIMRTNSEOB.h"
#ifndef _LALSIMIMRFACTORIZEDWAVEFORM_C
#define _LALSIMIMRFACTORIZEDWAVEFORM_C

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

/**
 * Constant which comes up in some of the EOB models. Its value is
 * (94/3 -41/32*pi*pi)
 */
#define ninty4by3etc 18.687902694437592603


static REAL8
XLALCalculateTNSEOBA_nontidal( const REAL8 r, 
                         const REAL8 eta                     /**<< Orbital separation (in units of total mass M) */
                     //    EOBACoefficients * restrict coeffs /**<< Pre-computed coefficients for the A function */
                       )
{

  REAL8 u,u2,u3,u4,u5,logu;
  REAL8 n1,d1,d2,d3,d4,d5;
  REAL8 a5,a6,a5tot,a6tot,a5tot2,pi2,pi4,eta2;
  REAL8 Num, Den;
  

  /* Note that this function uses pre-computed coefficients,
   * and assumes they have been calculated. Since this is a static function,
   * so only used here, I assume it is okay to neglect error checking
   */

///CHECKME : COMPUTING EOBA WITHOUT ANY PRECOMPUTED COEFF
  u=1./r;
  logu=log(u);
  eta2=eta*eta;
  pi2=LAL_PI*LAL_PI;
  pi4=pi2*pi2;
  
  REAL8 EulerGamma,a5l,a5c0,a5c1;
  EulerGamma=0.5772156649015328606065121;
  a5l=64./5.;
  a5c0=-4237./60.+2275./512.*pi2+256./5.*log(2.)+128./5.*EulerGamma;
  a5c1=-221./6.+41./32.*pi2;
  a5=a5c0+eta*a5c1;
  a5tot=a5+a5l*logu;

//  a6 = -122.+147.*(1.-4.*eta); (Where did u get this from ?)

  a6=(-110.5+129.*(1.-4.*eta))*sqrt(1.-0.000015/((eta-0.26)*(eta-0.26)));
  a6tot  = a6  + (-7004./105. - 144./5.*eta)*logu;
  a5tot2 = a5tot*a5tot;
  u2=u*u;
  u3=u2*u;
  u4=u2*u2;
  u5=u2*u3;

  n1  = (-3.*(-512. - 32.*eta2 + eta*(3520. + 32.*a5tot + 8.*a6tot - 123.*pi2)))/(-768. + eta*(3584. + 24.*a5tot - 123.*pi2));

  d1 = (eta*(-3392. - 48.*a5tot - 24.*a6tot + 96.*eta + 123.*pi2))/(-768. + eta*(3584. + 24.*a5tot - 123.*pi2));
  d2 =  (2.*eta*(-3392. - 48.*a5tot - 24.*a6tot + 96.*eta + 123.*pi2))/(-768. + eta*(3584. + 24.*a5tot - 123.*pi2));
  d3 = (-2.*eta*(6016. + 48.*a6tot + 3392.*eta + 24.*a5tot*(4. + eta) - 246.*pi2 - 123.*eta*pi2))/(-768. + eta*(3584. + 24.*a5tot - 123.*pi2));
  d4 = -(eta*(-4608.*a6tot*(-4. + eta) + a5tot*(36864. + eta*(72192. - 2952.*pi2)) + eta*(2048.*(5582. + 9.*eta) - 834432.*pi2 + 15129.*pi4)))/(96.*(-768. + eta*(3584.+ 24.*a5tot - 123.*pi2)));
  d5 = (eta*(-24.*a6tot*(1536. + eta*(-3776. + 123.*pi2)) + eta*(-2304.*a5tot2 + 96.*a5tot*(-3392. + 123.*pi2) - (-3776. + 123.*pi2)*(-3008. - 96.*eta + 123.*pi2))))/(96.*(-768. + eta*(3584. + 24.*a5tot - 123.*pi2)));


  Num = 1.+ (u * n1);

  Den = 1.+
     + u  * d1
     + u2 * d2
     + u3 * d3
     + u4 * d4
     + u5 * d5;

 // printf("\n%f",Num);
 // printf("Num.....");
//  printf("\n%f",Den);
//  printf("Den...");
  
//  printf("\n%f",u);
//  printf("u...");
  return Num/Den;
}

static REAL8
XLALCalculateTNSEOBA_nontidal_Debugging( const REAL8 r, 
                         const REAL8 eta                     /**<< Orbital separation (in units of total mass M) */
                     //    EOBACoefficients * restrict coeffs /**<< Pre-computed coefficients for the A function */
                       )
{

  REAL8 u,u2,u3,u4,u5,logu;
  REAL8 n1,d1,d2,d3,d4,d5;
  REAL8 a5,a6,a5tot,a6tot,a5tot2,pi2,pi4,eta2;
  REAL8 Num, Den;
  

  /* Note that this function uses pre-computed coefficients,
   * and assumes they have been calculated. Since this is a static function,
   * so only used here, I assume it is okay to neglect error checking
   */

///CHECKME : COMPUTING EOBA WITHOUT ANY PRECOMPUTED COEFF

  u=1./r;
  logu=log(u);
  eta2=eta*eta;
  pi2=LAL_PI*LAL_PI;
  pi4=pi2*pi2;
  
  REAL8 EulerGamma,a5l,a5c0,a5c1;
  EulerGamma=0.5772156649015328606065121;
  a5l=64./5.;
  a5c0=-4237./60.+2275./512.*pi2+256./5.*log(2.)+128./5.*EulerGamma;
  a5c1=-221./6.+41./32.*pi2;
  a5=a5c0+eta*a5c1;
  a5tot=a5+a5l*logu;

//  a6 = -122.+147.*(1.-4.*eta); (Where did u get this from ?)

  a6=(-110.5+129.*(1.-4.*eta))*sqrt(1.-0.000015/((eta-0.26)*(eta-0.26)));
  a6tot  = a6  + (-7004./105. - 144./5.*eta)*logu;
  a5tot2 = a5tot*a5tot;
  u2=u*u;
  u3=u2*u;
  u4=u2*u2;
  u5=u2*u3;


printf("\n eta :%.10f",eta); 
printf("\n u :%.10f",u); 
printf("\n a5tot :%.10f",a5); 
printf("\n a5l:%.10f",a5l); 
printf("\n a6tot :%.10f",a6tot); 

  n1  = (-3.*(-512. - 32.*eta2 + eta*(3520. + 32.*a5tot + 8.*a6tot - 123.*pi2)))/(-768. + eta*(3584. + 24.*a5tot - 123.*pi2));

  d1 = (eta*(-3392. - 48.*a5tot - 24.*a6tot + 96.*eta + 123.*pi2))/(-768. + eta*(3584. + 24.*a5tot - 123.*pi2));
  d2 =  (2.*eta*(-3392. - 48.*a5tot - 24.*a6tot + 96.*eta + 123.*pi2))/(-768. + eta*(3584. + 24.*a5tot - 123.*pi2));
  d3 = (-2.*eta*(6016. + 48.*a6tot + 3392.*eta + 24.*a5tot*(4. + eta) - 246.*pi2 - 123.*eta*pi2))/(-768. + eta*(3584. + 24.*a5tot - 123.*pi2));
  d4 = -(eta*(-4608.*a6tot*(-4. + eta) + a5tot*(36864. + eta*(72192. - 2952.*pi2)) + eta*(2048.*(5582. + 9.*eta) - 834432.*pi2 + 15129.*pi4)))/(96.*(-768. + eta*(3584.+ 24.*a5tot - 123.*pi2)));
  d5 = (eta*(-24.*a6tot*(1536. + eta*(-3776. + 123.*pi2)) + eta*(-2304.*a5tot2 + 96.*a5tot*(-3392. + 123.*pi2) - (-3776. + 123.*pi2)*(-3008. - 96.*eta + 123.*pi2))))/(96.*(-768. + eta*(3584. + 24.*a5tot - 123.*pi2)));


  Num = 1.+ (u * n1);

  Den = 1.+
     + u  * d1
     + u2 * d2
     + u3 * d3
     + u4 * d4
     + u5 * d5;

 
  printf("\n Numerator_A :%.10f",Num);
 // printf("Num.....");
  printf("\n Denominator_A %.10f",Den);
//  printf("Den...");
  
//  printf("\n%f",u);
//  printf("u...");
  return Num/Den;
}

static REAL8
XLALCalculatek2T( const REAL8 Lambda2A,const REAL8 Lambda2B,const REAL8 m1, const REAL8 m2)
{
REAL8 XA,XB,k2A,k2B,k2T;

XA=m1/(m1+m2);
XB=m2/(m1+m2);

k2A=3.*(XB/XA)*Lambda2A*(XA*XA*XA*XA*XA);
k2B=3.*(XA/XB)*Lambda2B*(XB*XB*XB*XB*XB);

k2T=k2A+k2B;

return k2T;                       
}

static REAL8
XLALCalculateLambda3( const REAL8 lambda2)
{REAL8 lambda3;
 REAL8 a3,b3,c3,d3,e3,loglambda2;

a3=-1.15;
b3=1.18;
c3=0.0251;
d3=-0.00131;
e3=0.0000252;

loglambda2=log(lambda2);
lambda3=exp(a3+(b3*loglambda2)+(c3*loglambda2*loglambda2)+(d3*loglambda2*loglambda2*loglambda2)+(e3*loglambda2*loglambda2*loglambda2*loglambda2));

return lambda3;
}

static REAL8
XLALCalculatek3T( const REAL8 lambda2A,const REAL8 lambda2B,const REAL8 m1, const REAL8 m2)
{REAL8 k3T;
REAL8 XA,XB,XA7,XB7,lambda3A,lambda3B,k3A,k3B;
XA=m1/(m1+m2);
XA7=XA*XA*XA*XA*XA*XA*XA;
XB=m2/(m1+m2);
XB7=XB*XB*XB*XB*XB*XB*XB;

lambda3A=XLALCalculateLambda3(lambda2A);
lambda3B=XLALCalculateLambda3(lambda2B);
k3A=15.*(m2/m1)*lambda3A*XA7;
k3B=15.*(m1/m2)*lambda3B*XB7;
k3T=k3A+k3B;

return k3T;
}


static REAL8
XLALCalculateLambda4( const REAL8 lambda2)
{REAL8 lambda4;
 REAL8 a4,b4,c4,d4,e4,loglambda2;

a4=-2.45;
b4=1.43;
c4=0.0395;
d4=-0.00181;
e4=0.0000280;

loglambda2=log(lambda2);
lambda4=exp(a4+(b4*loglambda2)+(c4*loglambda2*loglambda2)+(d4*loglambda2*loglambda2*loglambda2)+(e4*loglambda2*loglambda2*loglambda2*loglambda2));

return lambda4;
}


static REAL8
XLALCalculatek4T( const REAL8 lambda2A,const REAL8 lambda2B,const REAL8 m1, const REAL8 m2)
{REAL8 k4T;
REAL8 XA,XB,XA9,XB9,lambda4A,lambda4B,k4A,k4B;
XA=m1/(m1+m2);
XA9=XA*XA*XA*XA*XA*XA*XA*XA*XA;
XB=m2/(m1+m2);
XB9=XB*XB*XB*XB*XB*XB*XB*XB*XB;

lambda4A=XLALCalculateLambda4(lambda2A);
lambda4B=XLALCalculateLambda4(lambda2B);
k4A=105.*(m2/m1)*lambda4A*XA9;
k4B=105.*(m1/m2)*lambda4B*XB9;
k4T=k4A+k4B;

return k4T;
}


static
REAL8 XLALCalculateTNSEOBA2capTidal( const REAL8 r,
                         const REAL8 m1, const REAL8 m2, const REAL8 Lambda2A, const REAL8 Lambda2B)
{

REAL8 u,u2,XA,XB,k2A,k2B,k2T,A2captidal;
REAL8 star1alpha1_2,star2alpha1_2,star1alpha2_2,star2alpha2_2;
REAL8 alpha1bar_2,alpha2bar_2;


u=1./r;
u2=u*u;

if((Lambda2A==0)&&(Lambda2B==0))
 { A2captidal =0.;
  }
else {
    XA=m1/(m1+m2);
    XB=m2/(m1+m2);
  
    k2A=3.*(XB/XA)*Lambda2A*(XA*XA*XA*XA*XA);
    k2B=3.*(XA/XB)*Lambda2B*(XB*XB*XB*XB*XB);
    k2T=k2A+k2B;
  
    star1alpha1_2=(5./2.)*XA;
    star2alpha1_2=(5./2.)*XB;
  
    star1alpha2_2=((337./28.)*XA*XA)+((1./8.)*XA)+3.;
    star2alpha2_2=((337./28.)*XB*XB)+((1./8.)*XB)+3.;
  
    alpha1bar_2=(k2A*star1alpha1_2+k2B*star2alpha1_2)/k2T;
    alpha2bar_2=(k2A*star1alpha2_2+k2B*star2alpha2_2)/k2T;
  
  
    A2captidal=1.+(alpha1bar_2*u)+(alpha2bar_2 *u2);
  
    }
   return A2captidal;
    }
  
  

static
REAL8 XLALCalculateTNSEOBA3capTidal( const REAL8 r,
                         const REAL8 m1, const REAL8 m2, const REAL8 Lambda2A, const REAL8 Lambda2B)
{

REAL8 u,u2,XA,XB,XA7,XB7,k3A,k3B,k3T,A3captidal;
REAL8 star1alpha1_3,star2alpha1_3,star1alpha2_3,star2alpha2_3; 
REAL8 alpha1bar_3,alpha2bar_3,lambda3A,lambda3B;


u=1./r;
u2=u*u;

if((Lambda2A==0)&&(Lambda2B==0))
 { A3captidal =0.;
  }
else {
    

        XA=m1/(m1+m2);
	XA7=XA*XA*XA*XA*XA*XA*XA;
	XB=m2/(m1+m2);
	XB7=XB*XB*XB*XB*XB*XB*XB;

	lambda3A=XLALCalculateLambda3(Lambda2A);
	lambda3B=XLALCalculateLambda3(Lambda2B);
	k3A=15.*(m2/m1)*lambda3A*XA7;
	k3B=15.*(m1/m2)*lambda3B*XB7;
	k3T=k3A+k3B;
  
  
    	star1alpha1_3=-2.+(15./2.)*XA;
    	star2alpha1_3=-2.+(15./2.)*XB;
  
    	star1alpha2_3=8./3. -((311./24.)*XA)+((110./3.)*XA*XA);
    	star2alpha2_3=8./3. -((311./24.)*XB)+((110./3.)*XB*XB);
  
   	alpha1bar_3=(k3A*star1alpha1_3+k3B*star2alpha1_3)/k3T;
    	alpha2bar_3=(k3A*star1alpha2_3+k3B*star2alpha2_3)/k3T;
  
  
       A3captidal=1.+(alpha1bar_3*u)+(alpha2bar_3 *u2);
  
    }
   return A3captidal;
    }



static
REAL8 XLALCalculatedTNSEOBA3capTidal_du( const REAL8 r,
                         const REAL8 m1, const REAL8 m2, const REAL8 Lambda2A, const REAL8 Lambda2B)
{

REAL8 u,u2,XA,XB,XA7,XB7,k3A,k3B,k3T,dA3captidal_du;
REAL8 star1alpha1_3,star2alpha1_3,star1alpha2_3,star2alpha2_3;
REAL8 alpha1bar_3,alpha2bar_3,lambda3A,lambda3B;


u=1./r;
u2=u*u;

if((Lambda2A==0)&&(Lambda2B==0))
 { dA3captidal_du =0.;
  }
else {


        XA=m1/(m1+m2);
        XA7=XA*XA*XA*XA*XA*XA*XA;
        XB=m2/(m1+m2);
        XB7=XB*XB*XB*XB*XB*XB*XB;

        lambda3A=XLALCalculateLambda3(Lambda2A);
        lambda3B=XLALCalculateLambda3(Lambda2B);
        k3A=15.*(m2/m1)*lambda3A*XA7;
        k3B=15.*(m1/m2)*lambda3B*XB7;
        k3T=k3A+k3B;


        star1alpha1_3=-2.+(15./2.)*XA;
        star2alpha1_3=-2.+(15./2.)*XB;

        star1alpha2_3=8./3. -((311./24.)*XA)+((110./3.)*XA*XA);
        star2alpha2_3=8./3. -((311./24.)*XB)+((110./3.)*XB*XB);

        alpha1bar_3=(k3A*star1alpha1_3+k3B*star2alpha1_3)/k3T;
        alpha2bar_3=(k3A*star1alpha2_3+k3B*star2alpha2_3)/k3T;


       dA3captidal_du=(alpha1bar_3)+(2.*alpha2bar_3);

    }
   return dA3captidal_du;
    }


static
REAL8 XLALCalculated2TNSEOBA3capTidal_d2u( const REAL8 r,
                         const REAL8 m1, const REAL8 m2, const REAL8 Lambda2A, const REAL8 Lambda2B)
{

REAL8 u,u2,XA,XB,XA7,XB7,k3A,k3B,k3T,d2A3captidal_d2u;
REAL8 star1alpha1_3,star2alpha1_3,star1alpha2_3,star2alpha2_3;
REAL8 alpha1bar_3,alpha2bar_3,lambda3A,lambda3B;


u=1./r;
u2=u*u;

if((Lambda2A==0)&&(Lambda2B==0))
 { d2A3captidal_d2u =0.;
  }
else {


        XA=m1/(m1+m2);
        XA7=XA*XA*XA*XA*XA*XA*XA;
        XB=m2/(m1+m2);
        XB7=XB*XB*XB*XB*XB*XB*XB;

        lambda3A=XLALCalculateLambda3(Lambda2A);
        lambda3B=XLALCalculateLambda3(Lambda2B);
        k3A=15.*(m2/m1)*lambda3A*XA7;
        k3B=15.*(m1/m2)*lambda3B*XB7;
        k3T=k3A+k3B;


        star1alpha1_3=-2.+(15./2.)*XA;
        star2alpha1_3=-2.+(15./2.)*XB;

        star1alpha2_3=8./3. -((311./24.)*XA)+((110./3.)*XA*XA);
        star2alpha2_3=8./3. -((311./24.)*XB)+((110./3.)*XB*XB);

        alpha1bar_3=(k3A*star1alpha1_3+k3B*star2alpha1_3)/k3T;
        alpha2bar_3=(k3A*star1alpha2_3+k3B*star2alpha2_3)/k3T;


       d2A3captidal_d2u=(2.*alpha2bar_3);

    }
   return d2A3captidal_d2u;
    }




static
REAL8 XLALCalculatedTNSEOBA2capTidal_du( const REAL8 r,
                         const REAL8 m1, const REAL8 m2, const REAL8 Lambda2A, const REAL8 Lambda2B)
{

REAL8 u,u2,XA,XB,k2A,k2B,k2T,dA2captidal_du;
REAL8 star1alpha1_2,star2alpha1_2,star1alpha2_2,star2alpha2_2;
REAL8 alpha1bar_2,alpha2bar_2;


u=1./r;
u2=u*u;

if((Lambda2A==0)&&(Lambda2B==0))
 { dA2captidal_du =0.;
  }
else {
    XA=m1/(m1+m2);
    XB=m2/(m1+m2);

    k2A=3.*(XB/XA)*Lambda2A*(XA*XA*XA*XA*XA);
    k2B=3.*(XA/XB)*Lambda2B*(XB*XB*XB*XB*XB);
    k2T=k2A+k2B;

    star1alpha1_2=(5./2.)*XA;
    star2alpha1_2=(5./2.)*XB;

    star1alpha2_2=((337./28.)*XA*XA)+((1./8.)*XA)+3.;
    star2alpha2_2=((337./28.)*XB*XB)+((1./8.)*XB)+3.;

    alpha1bar_2=(k2A*star1alpha1_2+k2B*star2alpha1_2)/k2T;
    alpha2bar_2=(k2A*star1alpha2_2+k2B*star2alpha2_2)/k2T;


    dA2captidal_du=(alpha1bar_2)+(2.*alpha2bar_2 *u);

    }
   return dA2captidal_du;
    }





static
REAL8 XLALCalculated2TNSEOBA2capTidal_d2u( const REAL8 r,
                         const REAL8 m1, const REAL8 m2, const REAL8 Lambda2A, const REAL8 Lambda2B)
{

REAL8 u,u2,XA,XB,k2A,k2B,k2T,d2A2captidal_d2u;
REAL8 star1alpha1_2,star2alpha1_2,star1alpha2_2,star2alpha2_2;
REAL8 alpha1bar_2,alpha2bar_2;


u=1./r;
u2=u*u;

if((Lambda2A==0)&&(Lambda2B==0))
 { d2A2captidal_d2u =0.;
  }
else {
    XA=m1/(m1+m2);
    XB=m2/(m1+m2);

    k2A=3.*(XB/XA)*Lambda2A*(XA*XA*XA*XA*XA);
    k2B=3.*(XA/XB)*Lambda2B*(XB*XB*XB*XB*XB);
    k2T=k2A+k2B;

    star1alpha1_2=(5./2.)*XA;
    star2alpha1_2=(5./2.)*XB;

    star1alpha2_2=((337./28.)*XA*XA)+((1./8.)*XA)+3.;
    star2alpha2_2=((337./28.)*XB*XB)+((1./8.)*XB)+3.;

    alpha1bar_2=(k2A*star1alpha1_2+k2B*star2alpha1_2)/k2T;
    alpha2bar_2=(k2A*star1alpha2_2+k2B*star2alpha2_2)/k2T;


    d2A2captidal_d2u=(2.*alpha2bar_2);

    }
   return d2A2captidal_d2u;
    }




static
REAL8 XLALCalculateTNSEOBA_Tidal_nnlo( const REAL8 r,
                         const REAL8 m1, const REAL8 m2, const REAL8 Lambda2A, const REAL8 Lambda2B)
{


REAL8 u,u3,u6,u8,u10,XA,XB,k2T,A2capTidal,k3T,A3capTidal,k4T,A_tidal;
  

if((Lambda2A==0)&&(Lambda2B==0))
 { A_tidal =0.;
  }

else
{
u=1./r;
u3=u*u*u;
u6=u3*u3;
u8=u6*u*u;
u10=u8*u*u;
XA=m1/(m1+m2);
XB=m2/(m1+m2);



k2T=XLALCalculatek2T(Lambda2A,Lambda2B,m1, m2);
k3T=XLALCalculatek3T(Lambda2A,Lambda2B,m1, m2);
k4T=XLALCalculatek4T(Lambda2A,Lambda2B,m1, m2);



A2capTidal=XLALCalculateTNSEOBA2capTidal(r,m1,m2,Lambda2A,Lambda2B);
A3capTidal=XLALCalculateTNSEOBA3capTidal(r,m1,m2,Lambda2A,Lambda2B);

A_tidal=-(k4T*u10)-(k3T*u8*A3capTidal)-(k2T*u6*A2capTidal);
}
return A_tidal;
}


static
REAL8 XLALCalculatedTNSEOBA_Tidal_nnlo_du( const REAL8 r,
                         const REAL8 m1, const REAL8 m2, const REAL8 Lambda2A, const REAL8 Lambda2B)
{

REAL8 u,u3,u6,u8,u10,XA,XB,k2T,dA2capTidal_du,A2capTidal,A3capTidal,k3T,dA3capTidal_du,k4T,dA_tidal_du;

if((Lambda2A==0)&&(Lambda2B==0))
 { dA_tidal_du =0.;
  }


else
{u=1./r;
u3=u*u*u;
u6=u3*u3;
u8=u6*u*u;
u10=u8*u*u;
XA=m1/(m1+m2);
XB=m2/(m1+m2);



k2T=XLALCalculatek2T(Lambda2A,Lambda2B,m1, m2);
k3T=XLALCalculatek3T(Lambda2A,Lambda2B,m1, m2);
k4T=XLALCalculatek4T(Lambda2A,Lambda2B,m1, m2);


A2capTidal=XLALCalculateTNSEOBA2capTidal(r,m1,m2,Lambda2A,Lambda2B);
A3capTidal=XLALCalculateTNSEOBA3capTidal(r,m1,m2,Lambda2A,Lambda2B);
dA2capTidal_du=XLALCalculatedTNSEOBA2capTidal_du(r,m1,m2,Lambda2A,Lambda2B);
dA3capTidal_du=XLALCalculatedTNSEOBA3capTidal_du(r,m1,m2,Lambda2A,Lambda2B);

dA_tidal_du=-(k4T*u8*u*10.)-(k3T*u6*u*8.*A3capTidal)-(k3T*u8*dA3capTidal_du)-(k2T*u3*u*u*6*A2capTidal)-(k2T*u6*dA2capTidal_du);
}
return dA_tidal_du;
}




static
REAL8 XLALCalculated2TNSEOBA_Tidal_nnlo_d2u( const REAL8 r,
                         const REAL8 m1, const REAL8 m2, const REAL8 Lambda2A, const REAL8 Lambda2B)
{

REAL8 u,u3,u6,u8,u10,XA,XB,k2T,dA2capTidal_du,A2capTidal,A3capTidal,k3T,dA3capTidal_du,k4T,d2A_tidal_d2u,d2A2capTidal_d2u,d2A3capTidal_d2u;

if((Lambda2A==0)&&(Lambda2B==0))
 { d2A_tidal_d2u =0.;
  }


else
{
u=1./r;
u3=u*u*u;
u6=u3*u3;
u8=u6*u*u;
u10=u8*u*u;
XA=m1/(m1+m2);
XB=m2/(m1+m2);

k2T=XLALCalculatek2T(Lambda2A,Lambda2B,m1, m2);
k3T=XLALCalculatek3T(Lambda2A,Lambda2B,m1, m2);
k4T=XLALCalculatek4T(Lambda2A,Lambda2B,m1, m2);

A2capTidal=XLALCalculateTNSEOBA2capTidal(r,m1,m2,Lambda2A,Lambda2B);
A3capTidal=XLALCalculateTNSEOBA3capTidal(r,m1,m2,Lambda2A,Lambda2B);
dA2capTidal_du=XLALCalculatedTNSEOBA2capTidal_du(r,m1,m2,Lambda2A,Lambda2B);
dA3capTidal_du=XLALCalculatedTNSEOBA3capTidal_du(r,m1,m2,Lambda2A,Lambda2B);
d2A2capTidal_d2u=XLALCalculated2TNSEOBA2capTidal_d2u(r,m1,m2,Lambda2A,Lambda2B);
d2A3capTidal_d2u=XLALCalculated2TNSEOBA3capTidal_d2u(r,m1,m2,Lambda2A,Lambda2B);

d2A_tidal_d2u=-(k4T*u8*10.*9.)-(k3T*8.*((7.*u6*A3capTidal)+(u6*u*dA3capTidal_du)))-((k3T*8.*u6*u*dA3capTidal_du)+(k3T*u8*d2A3capTidal_d2u))-((k2T*u3*u*6.*5.*A2capTidal)+(k2T*u3*u*u*6.*dA2capTidal_du))-((k2T*u3*u*u*6.*dA2capTidal_du)+(k2T*u6*d2A2capTidal_d2u));
}
return d2A_tidal_d2u;
}



static REAL8
XLALCalculateTNSEOBA( const REAL8 r,
                         const REAL8 m1, const REAL8 m2, const REAL8 Lambda2A, const REAL8 Lambda2B)
{
//REAL8 dAtidaldu;
REAL8 totalMass,eta,A_nonTidal,A_tidal,A;

totalMass = m1 + m2;
eta = m1 * m2 / (totalMass*totalMass);

//printf("\n%f",r);
//printf("at radius....");
A_nonTidal=XLALCalculateTNSEOBA_nontidal(r,eta);
//printf("\n%f",A_nonTidal);
//printf("A_nonTidal....");
A_tidal=XLALCalculateTNSEOBA_Tidal_nnlo( r,m1,m2,Lambda2A,Lambda2B);
//printf("\n%f",A2_tidal);
//printf("A2_tidal....");
A=A_nonTidal+A_tidal;
return A;
}


//CHECK ME : dAdu need to get updated once u add tidal potential. 
static
REAL8 XLALCalculateEOB_Nontidal_dAdu( const REAL8 r, const REAL8 eta                     /**<< Orbital separation (in units of total mass M) */
                          //  EOBACoefficients * restrict coeffs /**<< Pre-computed coefficients for the A function */
                          )
{ REAL8 u,u2, u3, u4, u5,logu;
  REAL8 dn1,dd1,dd2,dd3,dd4,dd5;
  REAL8 n1,d1,d2,d3,d4,d5;
  REAL8 a5,a6,a5tot,a6tot,a5tot2,pi2,eta2,pi4;
  REAL8 A,Num,Den,dNum,dDen,prefactor,dANonTidal_du;

//OLD NR CALIB
  //u=1./r;
  //logu=log(u);
  //pi2=LAL_PI*LAL_PI;
  //pi4=pi2*pi2;
  //eta2=eta*eta;

 // REAL8 EulerGamma,a5l,a5c0,a5c1;
 // EulerGamma=0.5772156649015328606065121;
 // a5l=64./5.;
 // a5c0=-4237./60.+2275./512.*pi2+256./5.*log(2.)+128./5.*EulerGamma;
 // a5c1=-221./6.+41./32.*pi2;
 // a5=a5c0+eta*a5c1;
 // a5tot=a5+a5l*logu;  

 // a6 = -122.+147.*(1.-4.*eta);
  //a6tot  = a6  + (-7004./105. - 144./5.*eta)*logu;
  //a5tot2 = a5tot*a5tot;
 // u2=u*u;
 // u3=u2*u;
 // u4=u2*u2;
//  u5=u2*u3;

  u=1./r;
  logu=log(u);
  eta2=eta*eta;
  pi2=LAL_PI*LAL_PI;
  pi4=pi2*pi2;
  
  REAL8 EulerGamma,a5l,a5c0,a5c1;
  EulerGamma=0.5772156649015328606065121;
  a5l=64./5.;
  a5c0=-4237./60.+2275./512.*pi2+256./5.*log(2.)+128./5.*EulerGamma;
  a5c1=-221./6.+41./32.*pi2;
  a5=a5c0+eta*a5c1;
  a5tot=a5+a5l*logu;

//  a6 = -122.+147.*(1.-4.*eta); (Where did u get this from ?)

  a6=(-110.5+129.*(1.-4.*eta))*sqrt(1.-0.000015/((eta-0.26)*(eta-0.26)));
  a6tot  = a6  + (-7004./105. - 144./5.*eta)*logu;
  a5tot2 = a5tot*a5tot;
  u2=u*u;
  u3=u2*u;
  u4=u2*u2;
  u5=u2*u3;



  n1  = (-3.*(-512. - 32.*eta2 + eta*(3520. + 32.*a5tot + 8.*a6tot - 123.*pi2)))/(-768. + eta*(3584. + 24.*a5tot - 123.*pi2));

  d1 = (eta*(-3392. - 48.*a5tot - 24.*a6tot + 96.*eta + 123.*pi2))/(-768. + eta*(3584. + 24.*a5tot - 123.*pi2));
  d2 =  (2.*eta*(-3392. - 48.*a5tot - 24.*a6tot + 96.*eta + 123.*pi2))/(-768. + eta*(3584. + 24.*a5tot - 123.*pi2));
  d3 = (-2.*eta*(6016. + 48.*a6tot + 3392.*eta + 24.*a5tot*(4. + eta) - 246.*pi2 - 123.*eta*pi2))/(-768. + eta*(3584. + 24.*a5tot - 123.*pi2));
  d4 = -(eta*(-4608.*a6tot*(-4. + eta) + a5tot*(36864. + eta*(72192. - 2952.*pi2)) + eta*(2048.*(5582. + 9.*eta) - 834432.*pi2 + 15129.*pi4)))/(96.*(-768. + eta*(3584.+ 24.*a5tot - 123.*pi2)));
  d5 = (eta*(-24.*a6tot*(1536. + eta*(-3776. + 123.*pi2)) + eta*(-2304.*a5tot2 + 96.*a5tot*(-3392. + 123.*pi2) - (-3776. + 123.*pi2)*(-3008. - 96.*eta + 123.*pi2))))/(96.*(-768. + eta*(3584. + 24.*a5tot - 123.*pi2)));


  Num = 1.+ (u * n1);

  Den = 1.+
     + u  * d1
     + u2 * d2
     + u3 * d3
     + u4 * d4
     + u5 * d5;

  A=Num/Den;
  dn1 = (160.*eta*(-828672. - 32256.*eta2 + 756.*eta*(-768. + eta*(3584. + 24.*a5 - 123.*pi2)) + eta*(5006848. + 42024.*a5 + 8064.*a6 - 174045.*pi2)))/(7.*pow(1536*logu*eta + 5.*(-768. + eta*(3584. + 24.*a5 - 123.*pi2)),2.)*u);

  dd1 = (160.*eta*(-828672. - 32256.*eta2 + 756.*eta*(-768. + eta*(3584. + 24.*a5 - 123.*pi2)) + eta*(5006848. + 42024.*a5 + 8064.*a6 - 174045.*pi2)))/(7.*pow(1536.*logu*eta + 5.*(-768. + eta*(3584. + 24.*a5 - 123.*pi2)),2.)*u);

  dd2 = (320.*eta*(-828672. - 32256.*eta2 + 756.*eta*(-768. + eta*(3584. + 24.*a5 - 123.*pi2)) + eta*(5006848. + 42024.*a5 + 8064.*a6 - 174045.*pi2)))/(7.*pow(1536.*logu*eta + 5.*(-768. + eta*(3584. + 24.*a5 - 123.*pi2)),2.)*u);

  dd3 = (640.*eta*(-828672. - 32256.*eta2 + 756.*eta*(-768. + eta*(3584. + 24.*a5 - 123.*pi2)) + eta*(5006848. + 42024.*a5 + 8064.*a6 - 174045.*pi2)))/(7.*pow(1536.*logu*eta + 5.*(-768. + eta*(3584. + 24.*a5 - 123.*pi2)),2.)*u);

  dd4 = (-320.*(-4. + eta)*eta*(-828672. - 32256.*eta2 + 756.*eta*(-768. + eta*(3584. + 24.*a5 - 123.*pi2)) + eta*(5006848. + 42024.*a5 + 8064.*a6 - 174045.*pi2)))/(7.*pow(1536.*logu*eta + 5.*(-768. + eta*(3584. + 24.*a5 - 123.*pi2)),2.)*u);

  dd5 = (eta*(-8400.*eta*(-24.*(a6 - (4.*logu*(1751. + 756.*eta))/105.)*(1536. + eta*(-3776. + 123.*pi2))+ eta*(-2304.*pow(a5 + (64.*logu)/5.,2.) + 96.*(a5 + (64.*logu)/5.)*(-3392. + 123.*pi2) - (-3776. + 123.*pi2)*(-32.*(94. + 3.*eta) + 123.*pi2)))- (1536.*logu*eta + 5.*(-768. + eta*(3584. + 24.*a5 - 123.*pi2)))*(4128768.*logu*eta + 5.*(-2689536. + eta*(11170624. + 64512.*a5 - 380685.*pi2) - 756.*eta*(1536. + eta*(-3776. + 123.*pi2))))))/(2625.*pow(-768. + eta*(3584. + 24.*(a5 + (64.*logu)/5.) - 123.*pi2),2.)*u);

 //First derivative
  dNum  = dn1*u + n1;
  dDen  = d1 + u*(dd1 + 2*d2) + u2*(dd2 + 3*d3) + u3*(dd3 + 4*d4) + u4*(dd4 + 5*d5) + dd5*u5;

// Derivative of A function with respect to u
  prefactor = A/(Num*Den);
  dANonTidal_du = prefactor*(dNum*Den - dDen*Num);
  
  return dANonTidal_du;
}


static
REAL8 XLALCalculateEOBdAdu( const REAL8 r,
                         const REAL8 m1, const REAL8 m2, const REAL8 Lambda2A, const REAL8 Lambda2B)
{
REAL8 totalMass,eta,dA_tidal_du,dANonTidal_du,dA_du;

totalMass = m1 + m2;
eta = m1 * m2 / (totalMass*totalMass);


dA_tidal_du=XLALCalculatedTNSEOBA_Tidal_nnlo_du(r,m1,m2,Lambda2A,Lambda2B);


dANonTidal_du=XLALCalculateEOB_Nontidal_dAdu(r, eta);
//printf("\n%f",dANonTidal_du);
//printf("dANonTidal_du");


dA_du=dA_tidal_du+dANonTidal_du;
return dA_du;
}
/**
 * Calculated the derivative of the EOB A function with respect to
 * r
 */ 

//THIS IS NOT GETTING SED ANYWHERE. YOU ONLY HAVE NON TIDAL PEICE IN THIS! 
static
REAL8 XLALCalculateEOBdAdr( const REAL8 r, const REAL8 eta                     /**<< Orbital separation (in units of total mass M) */
                          //  EOBACoefficients * restrict coeffs /**<< Pre-computed coefficients for the A function */
                          )
{
  REAL8 u,u2, u3, u4, u5,logu;
  REAL8 dn1,dd1,dd2,dd3,dd4,dd5;
  REAL8 n1,d1,d2,d3,d4,d5;
  REAL8 a5,a6,a5tot,a6tot,a5tot2,pi2,eta2,pi4;
  REAL8 A,Num,Den,dNum,dDen,prefactor,dA_u;


 // u=1./r;
 // logu=log(u);
 // pi2=LAL_PI*LAL_PI;
 // pi4=pi2*pi2;
 // eta2=eta*eta;
  
//  REAL8 EulerGamma,a5l,a5c0,a5c1;
  //EulerGamma=0.5772156649015328606065121;
 // a5l=64./5.;
 // a5c0=-4237./60.+2275./512.*pi2+256./5.*log(2.)+128./5.*EulerGamma;
 // a5c1=-221./6.+41./32.*pi2;
  //a5=a5c0+eta*a5c1;
 // a5tot=a5+a5l*logu;

//  a6 = -122.+147.*(1.-4.*eta);
//  a6tot  = a6  + (-7004./105. - 144./5.*eta)*logu;
//  a5tot2 = a5tot*a5tot;
  
  u=1./r;
  logu=log(u);
  eta2=eta*eta;
  pi2=LAL_PI*LAL_PI;
  pi4=pi2*pi2;
  
  REAL8 EulerGamma,a5l,a5c0,a5c1;
  EulerGamma=0.5772156649015328606065121;
  a5l=64./5.;
  a5c0=-4237./60.+2275./512.*pi2+256./5.*log(2.)+128./5.*EulerGamma;
  a5c1=-221./6.+41./32.*pi2;
  a5=a5c0+eta*a5c1;
  a5tot=a5+a5l*logu;

//  a6 = -122.+147.*(1.-4.*eta); (Where did u get this from ?)

  a6=(-110.5+129.*(1.-4.*eta))*sqrt(1.-0.000015/((eta-0.26)*(eta-0.26)));
  a6tot  = a6  + (-7004./105. - 144./5.*eta)*logu;
  a5tot2 = a5tot*a5tot;
  u2=u*u;
  u3=u2*u;
  u4=u2*u2;
  u5=u2*u3;



  n1  = (-3.*(-512. - 32.*eta2 + eta*(3520. + 32.*a5tot + 8.*a6tot - 123.*pi2)))/(-768. + eta*(3584. + 24.*a5tot - 123.*pi2));

  d1 = (eta*(-3392. - 48.*a5tot - 24.*a6tot + 96.*eta + 123.*pi2))/(-768. + eta*(3584. + 24.*a5tot - 123.*pi2));
  d2 =  (2.*eta*(-3392. - 48.*a5tot - 24.*a6tot + 96.*eta + 123.*pi2))/(-768. + eta*(3584. + 24.*a5tot - 123.*pi2));
  d3 = (-2.*eta*(6016. + 48.*a6tot + 3392.*eta + 24.*a5tot*(4. + eta) - 246.*pi2 - 123.*eta*pi2))/(-768. + eta*(3584. + 24.*a5tot - 123.*pi2));
  d4 = -(eta*(-4608.*a6tot*(-4. + eta) + a5tot*(36864. + eta*(72192. - 2952.*pi2)) + eta*(2048.*(5582. + 9.*eta) - 834432.*pi2 + 15129.*pi4)))/(96.*(-768. + eta*(3584.+ 24.*a5tot - 123.*pi2)));
  d5 = (eta*(-24.*a6tot*(1536. + eta*(-3776. + 123.*pi2)) + eta*(-2304.*a5tot2 + 96.*a5tot*(-3392. + 123.*pi2) - (-3776. + 123.*pi2)*(-3008. - 96.*eta + 123.*pi2))))/(96.*(-768. + eta*(3584. + 24.*a5tot - 123.*pi2)));


  Num = 1.+ (u * n1);

  Den = 1.+
     + u  * d1
     + u2 * d2
     + u3 * d3
     + u4 * d4
     + u5 * d5;

  A=Num/Den;
  dn1 = (160.*eta*(-828672. - 32256.*eta2 + 756.*eta*(-768. + eta*(3584. + 24.*a5 - 123.*pi2)) + eta*(5006848. + 42024.*a5 + 8064.*a6 - 174045.*pi2)))/(7.*pow(1536*logu*eta + 5.*(-768. + eta*(3584. + 24.*a5 - 123.*pi2)),2.)*u);

  dd1 = (160.*eta*(-828672. - 32256.*eta2 + 756.*eta*(-768. + eta*(3584. + 24.*a5 - 123.*pi2)) + eta*(5006848. + 42024.*a5 + 8064.*a6 - 174045.*pi2)))/(7.*pow(1536.*logu*eta + 5.*(-768. + eta*(3584. + 24.*a5 - 123.*pi2)),2.)*u);

  dd2 = (320.*eta*(-828672. - 32256.*eta2 + 756.*eta*(-768. + eta*(3584. + 24.*a5 - 123.*pi2)) + eta*(5006848. + 42024.*a5 + 8064.*a6 - 174045.*pi2)))/(7.*pow(1536.*logu*eta + 5.*(-768. + eta*(3584. + 24.*a5 - 123.*pi2)),2.)*u);

  dd3 = (640.*eta*(-828672. - 32256.*eta2 + 756.*eta*(-768. + eta*(3584. + 24.*a5 - 123.*pi2)) + eta*(5006848. + 42024.*a5 + 8064.*a6 - 174045.*pi2)))/(7.*pow(1536.*logu*eta + 5.*(-768. + eta*(3584. + 24.*a5 - 123.*pi2)),2.)*u);

  dd4 = (-320.*(-4. + eta)*eta*(-828672. - 32256.*eta2 + 756.*eta*(-768. + eta*(3584. + 24.*a5 - 123.*pi2)) + eta*(5006848. + 42024.*a5 + 8064.*a6 - 174045.*pi2)))/(7.*pow(1536.*logu*eta + 5.*(-768. + eta*(3584. + 24.*a5 - 123.*pi2)),2.)*u);

  dd5 = (eta*(-8400.*eta*(-24.*(a6 - (4.*logu*(1751. + 756.*eta))/105.)*(1536. + eta*(-3776. + 123.*pi2))+ eta*(-2304.*pow(a5 + (64.*logu)/5.,2.) + 96.*(a5 + (64.*logu)/5.)*(-3392. + 123.*pi2) - (-3776. + 123.*pi2)*(-32.*(94. + 3.*eta) + 123.*pi2)))- (1536.*logu*eta + 5.*(-768. + eta*(3584. + 24.*a5 - 123.*pi2)))*(4128768.*logu*eta + 5.*(-2689536. + eta*(11170624. + 64512.*a5 - 380685.*pi2) - 756.*eta*(1536. + eta*(-3776. + 123.*pi2))))))/(2625.*pow(-768. + eta*(3584. + 24.*(a5 + (64.*logu)/5.) - 123.*pi2),2.)*u);

 //First derivative
  dNum  = dn1*u + n1;
  dDen  = d1 + u*(dd1 + 2*d2) + u2*(dd2 + 3*d3) + u3*(dd3 + 4*d4) + u4*(dd4 + 5*d5) + dd5*u5;

// Derivative of A function with respect to u
  prefactor = A/(Num*Den);
  dA_u = prefactor*(dNum*Den - dDen*Num);

// Derivative of A w



  return -u2*dA_u;
}



static
REAL8 XLALCalculateEOBd2Adu2_NonTidal( const REAL8 r, const REAL8 eta                     /**<< Orbital separation (in units of total mass M) */
                          //  EOBACoefficients * restrict coeffs /**<< Pre-computed coefficients for the A function */
                          )
{
  REAL8 u,u2, u3, u4, u5,logu;
  REAL8 dn1,dd1,dd2,dd3,dd4,dd5;
  REAL8 ddn1,ddd1,ddd2,ddd3,ddd4,ddd5;
  REAL8 n1,d1,d2,d3,d4,d5;
  REAL8 a5,a6,a5tot,a6tot,a5tot2,pi2,eta2,pi4;
  REAL8 A,Num,Den,dNum,dDen,ddNum, ddDen,prefactor,dA_u,ddA_u;


  //u=1./r;
  //logu=log(u);
  //pi2=LAL_PI*LAL_PI;
  //pi4=pi2*pi2;
  //eta2=eta*eta;
  
  //REAL8 EulerGamma,a5l,a5c0,a5c1;
  //EulerGamma=0.5772156649015328606065121;
  //a5l=64./5.;
  //a5c0=-4237./60.+2275./512.*pi2+256./5.*log(2.)+128./5.*EulerGamma;
  //a5c1=-221./6.+41./32.*pi2;
  //a5=a5c0+eta*a5c1;
  //a5tot=a5+a5l*logu;

  //a6 = -122.+147.*(1.-4.*eta);
  //a6tot  = a6  + (-7004./105. - 144./5.*eta)*logu;
  //a5tot2 = a5tot*a5tot;
  //u2=u*u;
  //u3=u2*u;
  //u4=u2*u2;
  //u5=u2*u3;

//
  u=1./r;
  logu=log(u);
  eta2=eta*eta;
  pi2=LAL_PI*LAL_PI;
  pi4=pi2*pi2;
  
  REAL8 EulerGamma,a5l,a5c0,a5c1;
  EulerGamma=0.5772156649015328606065121;
  a5l=64./5.;
  a5c0=-4237./60.+2275./512.*pi2+256./5.*log(2.)+128./5.*EulerGamma;
  a5c1=-221./6.+41./32.*pi2;
  a5=a5c0+eta*a5c1;
  a5tot=a5+a5l*logu;

//  a6 = -122.+147.*(1.-4.*eta); (Where did u get this from ?)

  a6=(-110.5+129.*(1.-4.*eta))*sqrt(1.-0.000015/((eta-0.26)*(eta-0.26)));
  a6tot  = a6  + (-7004./105. - 144./5.*eta)*logu;
  a5tot2 = a5tot*a5tot;
  u2=u*u;
  u3=u2*u;
  u4=u2*u2;
  u5=u2*u3;


  n1  = (-3.*(-512. - 32.*eta2 + eta*(3520. + 32.*a5tot + 8.*a6tot - 123.*pi2)))/(-768. + eta*(3584. + 24.*a5tot - 123.*pi2));

  d1 = (eta*(-3392. - 48.*a5tot - 24.*a6tot + 96.*eta + 123.*pi2))/(-768. + eta*(3584. + 24.*a5tot - 123.*pi2));
  d2 =  (2.*eta*(-3392. - 48.*a5tot - 24.*a6tot + 96.*eta + 123.*pi2))/(-768. + eta*(3584. + 24.*a5tot - 123.*pi2));
  d3 = (-2.*eta*(6016. + 48.*a6tot + 3392.*eta + 24.*a5tot*(4. + eta) - 246.*pi2 - 123.*eta*pi2))/(-768. + eta*(3584. + 24.*a5tot - 123.*pi2));
  d4 = -(eta*(-4608.*a6tot*(-4. + eta) + a5tot*(36864. + eta*(72192. - 2952.*pi2)) + eta*(2048.*(5582. + 9.*eta) - 834432.*pi2 + 15129.*pi4)))/(96.*(-768. + eta*(3584.+ 24.*a5tot - 123.*pi2)));
  d5 = (eta*(-24.*a6tot*(1536. + eta*(-3776. + 123.*pi2)) + eta*(-2304.*a5tot2 + 96.*a5tot*(-3392. + 123.*pi2) - (-3776. + 123.*pi2)*(-3008. - 96.*eta + 123.*pi2))))/(96.*(-768. + eta*(3584. + 24.*a5tot - 123.*pi2)));


  Num = 1.+ (u * n1);

  Den = 1.+
     + u  * d1
     + u2 * d2
     + u3 * d3
     + u4 * d4
     + u5 * d5;

  A=Num/Den;
  dn1 = (160.*eta*(-828672. - 32256.*eta2 + 756.*eta*(-768. + eta*(3584. + 24.*a5 - 123.*pi2)) + eta*(5006848. + 42024.*a5 + 8064.*a6 - 174045.*pi2)))/(7.*pow(1536*logu*eta + 5.*(-768. + eta*(3584. + 24.*a5 - 123.*pi2)),2.)*u);

  dd1 = (160.*eta*(-828672. - 32256.*eta2 + 756.*eta*(-768. + eta*(3584. + 24.*a5 - 123.*pi2)) + eta*(5006848. + 42024.*a5 + 8064.*a6 - 174045.*pi2)))/(7.*pow(1536.*logu*eta + 5.*(-768. + eta*(3584. + 24.*a5 - 123.*pi2)),2.)*u);

  dd2 = (320.*eta*(-828672. - 32256.*eta2 + 756.*eta*(-768. + eta*(3584. + 24.*a5 - 123.*pi2)) + eta*(5006848. + 42024.*a5 + 8064.*a6 - 174045.*pi2)))/(7.*pow(1536.*logu*eta + 5.*(-768. + eta*(3584. + 24.*a5 - 123.*pi2)),2.)*u);

  dd3 = (640.*eta*(-828672. - 32256.*eta2 + 756.*eta*(-768. + eta*(3584. + 24.*a5 - 123.*pi2)) + eta*(5006848. + 42024.*a5 + 8064.*a6 - 174045.*pi2)))/(7.*pow(1536.*logu*eta + 5.*(-768. + eta*(3584. + 24.*a5 - 123.*pi2)),2.)*u);

  dd4 = (-320.*(-4. + eta)*eta*(-828672. - 32256.*eta2 + 756.*eta*(-768. + eta*(3584. + 24.*a5 - 123.*pi2)) + eta*(5006848. + 42024.*a5 + 8064.*a6 - 174045.*pi2)))/(7.*pow(1536.*logu*eta + 5.*(-768. + eta*(3584. + 24.*a5 - 123.*pi2)),2.)*u);

  dd5 = (eta*(-8400.*eta*(-24.*(a6 - (4.*logu*(1751. + 756.*eta))/105.)*(1536. + eta*(-3776. + 123.*pi2))+ eta*(-2304.*pow(a5 + (64.*logu)/5.,2.) + 96.*(a5 + (64.*logu)/5.)*(-3392. + 123.*pi2) - (-3776. + 123.*pi2)*(-32.*(94. + 3.*eta) + 123.*pi2)))- (1536.*logu*eta + 5.*(-768. + eta*(3584. + 24.*a5 - 123.*pi2)))*(4128768.*logu*eta + 5.*(-2689536. + eta*(11170624. + 64512.*a5 - 380685.*pi2) - 756.*eta*(1536. + eta*(-3776. + 123.*pi2))))))/(2625.*pow(-768. + eta*(3584. + 24.*(a5 + (64.*logu)/5.) - 123.*pi2),2.)*u);

 //First derivative
  dNum  = dn1*u + n1;
  dDen  = d1 + u*(dd1 + 2.*d2) + u2*(dd2 + 3.*d3) + u3*(dd3 + 4.*d4) + u4*(dd4 + 5.*d5) + dd5*u5;

// Derivative of A function with respect to u
  prefactor = A/(Num*Den);
  dA_u = prefactor*(dNum*Den - dDen*Num);

//2nd derivative

  ddn1 = (160.*eta*(-3840. + 1536.*logu*eta + eta*(20992. + 120.*a5 - 615.*pi2))*(828672. + eta*(-42024.*a5 - 8064.*a6 + 3584.*(-1397. + 9.*eta) + 174045. *pi2) + 756. *eta*(768. + eta*(-3584. - 24.*a5 + 123.*pi2))))/(7.*pow (1536.*logu*eta + 5.*(-768. + eta*(3584. + 24.*a5 - 123.*pi2)),3.)*u2);

  ddd1 = (160.*eta*(-3840. + 1536.*logu*eta + eta*(20992. + 120.*a5 - 615.*pi2))*(828672. + eta*(-42024.*a5 - 8064.*a6 + 3584.*(-1397. + 9.*eta)+ 174045.*pi2) + 756.*eta*(768. + eta*(-3584. - 24.*a5 + 123.*pi2))))/(7.*pow (1536.*logu*eta + 5.*(-768. + eta*(3584. + 24.*a5 - 123.*pi2)),3)*u2);

  ddd2 = (320.*eta*(-3840. + 1536.*logu*eta + eta*(20992. + 120.*a5 - 615.*pi2))*(828672. + eta*(-42024.*a5 - 8064.*a6 + 3584.*(-1397. + 9.*eta) + 174045.*pi2) + 756.*eta*(768. + eta*(-3584. - 24.*a5 + 123.*pi2))))/(7.*pow (1536.*logu*eta + 5.*(-768. + eta*(3584. + 24.*a5 - 123.*pi2)),3.)*u2);
  ddd3 = (640.*eta*(-3840. + 1536.*logu*eta + eta*(20992. + 120.*a5 - 615.*pi2))*(828672. + eta*(-42024.*a5 - 8064.*a6 + 3584.*(-1397. + 9.*eta)+ 174045.*pi2) + 756.*eta*(768. + eta*(-3584. - 24.*a5 + 123.*pi2))))/(7.*pow (1536.*logu*eta + 5.*(-768. + eta*(3584. + 24.*a5 - 123.*pi2)),3.)*u2);
  ddd4 = (320.*(-4. + eta)*eta*(-828672. + 756.*eta*(-768. + eta*(3584. + 24.*a5 - 123.*pi2))+ eta*(5006848. + 42024.*a5 + 8064.*a6 - 32256.*eta - 174045.*pi2))*(-3840. + 1536.*logu*eta + eta*(20992. + 120.*a5 - 615.*pi2)))/(7.*pow (1536.*logu*eta + 5.*(-768. + eta*(3584. + 24.*a5 - 123.*pi2)),3)*u2);
  ddd5 = (eta*(pow (1536.*logu*eta + 5.*(-768. + eta*(3584. + 24.*a5 - 123.*pi2)),2.)*(4128768.*logu*eta - 7680.*(1751. + 756.*eta) + eta*(64.*(808193.+ 5040.*a5 + 223020.*eta) - 615.*(3095. + 756.*eta)*pi2)) + 3072.*eta*(1536.*logu*eta + 5.*(-768. + eta*(3584. + 24.*a5 - 123.*pi2)))*(4128768.*logu*eta - 7680.*(1751. + 756.*eta) + 5.*eta*(64.*(174541. + 1008.*a5 + 44604.*eta) - 123.*(3095. + 756.*eta)*pi2))+ 25804800.*eta2*(-24.*(a6 - (4.*logu*(1751. + 756.*eta))/105.)*(1536. + eta*(-3776. + 123.*pi2)) + eta*(-2304.*pow(a5 + (64.*logu)/5.,2) + 96.*(a5 + (64.*logu)/5.)*(-3392. + 123.*pi2) - (-3776. + 123.*pi2)*(-32.*(94. + 3.*eta) + 123.*pi2)))+ 42000.*eta*(-768. + eta*(3584. + 24.*(a5 + (64.*logu)/5.) - 123.*pi2))*(-24.*(a6 - (4.*logu*(1751. + 756.*eta))/105.)*(1536. + eta*(-3776. + 123.*pi2))+ eta*(-2304.*pow(a5 + (64.*logu)/5.,2) + 96.*(a5 + (64.*logu)/5.)*(-3392. + 123.*pi2) - (-3776. + 123.*pi2)*(-32.*(94. + 3.*eta) + 123.*pi2)))))/(13125.*pow (-768 + eta*(3584 + 24*(a5 + (64*logu)/5.) - 123*pi2),3)*u2);

 ddNum = 2.*dn1 + ddn1*u;
 ddDen = 2.*(d2 + dd1) + u*(6.*d3 + 4.*dd2 + ddd1) + u2*(12.*d4 + 6.*dd3 + ddd2)+ u3 *(20.*d5 + 8.*dd4 + ddd3) + u4*(10.*dd5 + ddd4) + u5*ddd5;
 
  ddA_u    = prefactor*(2.*pow(dDen,2)*A - 2.*dNum*dDen + Den*ddNum - ddDen*Num);
            
  return ddA_u;
}



UNUSED static
REAL8 XLALCalculateEOBd2Adu2( const REAL8 r, const REAL8 eta,const REAL8 m1, const REAL8 m2, const REAL8 Lambda2A, const REAL8 Lambda2B)
{
REAL8 d2Adu2_NonTidal,d2Adu2_Tidal,d2Adu2;

d2Adu2_NonTidal= XLALCalculateEOBd2Adu2_NonTidal(r,eta);
d2Adu2_Tidal=XLALCalculated2TNSEOBA_Tidal_nnlo_d2u(r,m1,m2,Lambda2A,Lambda2B);
d2Adu2=d2Adu2_NonTidal+d2Adu2_Tidal;

return d2Adu2;
}

/**
 * Calculated the derivative of the EOB A function with respect to
 * r, using the pre-computed A coefficients
 */
//static//
//REAL8 XLALCalculateEOBdAdr( const REAL8 r,                     /**<< Orbital separation (in units of total mass M) */
//                            EOBACoefficients * restrict coeffs /**<< Pre-computed coefficients for the A function */
//                          )
//{
//  REAL8 r2, r3, r4, r5;

//  REAL8 NA, DA, dNA, dDA, dA;

//  r2 = r*r;
 // r3 = r2 * r;
//  r4 = r2*r2;
//  r5 = r4*r;

//  NA = r4 * coeffs->n4
//     + r5 * coeffs->n5;
//
//  DA = coeffs->d0
//     + r  * coeffs->d1
//     + r2 * coeffs->d2
//     + r3 * coeffs->d3
//     + r4 * coeffs->d4
//     + r5 * coeffs->d5;

//  dNA = 4. * coeffs->n4 * r3
//      + 5. * coeffs->n5 * r4;

//  dDA = coeffs->d1
//      + 2. * coeffs->d2 * r
//      + 3. * coeffs->d3 * r2
//      + 4. * coeffs->d4 * r3
//     + 5. * coeffs->d5 * r4;

//  dA = dNA * DA - dDA * NA;

//  return dA / (DA*DA);
//}



//static
//REAL8 XLALCalculateEOBdAdr( const REAL8 r, const REAL8 eta                     /**<< Orbital separation (in units of total mass M) */
                          //  EOBACoefficients * restrict coeffs /**<< Pre-computed coefficients for the A function */
//                          )
//{
//  REAL8 u,u2, u3, u4, u5,logu;
//  REAL8 dn1,dd1,dd2,dd3,dd4,dd5;
//  REAL8 n1,d1,d2,d3,d4,d5;
//  REAL8 a5,a6,a5tot,a6tot,a5tot2,pi2,eta2,pi4;
//  REAL8 A,Num,Den,dNum,dDen,prefactor,dA_u;


//  u=1./r;
//  logu=log(u);
//  pi2=LAL_PI*LAL_PI;
//  pi4=pi2*pi2;
//  eta2=eta*eta;
//  a5 =  23.5;
//  a6 = -122.+147.*(1.-4.*eta);
//  a5tot  = a5  + 64./5.*logu;
//  a6tot  = a6  + (-7004./105. - 144./5.*eta)*logu;
//  a5tot2 = a5tot*a5tot;
//  u2=u*u;
//  u3=u2*u;
//  u4=u2*u2;
//  u5=u2*u3;

//  n1  = (-3.*(-512. - 32.*eta2 + eta*(3520. + 32.*a5tot + 8.*a6tot - 123.*pi2)))/(-768. + eta*(3584. + 24.*a5tot - 123.*pi2));

//  d1 = (eta*(-3392. - 48.*a5tot - 24.*a6tot + 96.*eta + 123.*pi2))/(-768. + eta*(3584. + 24.*a5tot - 123.*pi2));
//  d2 =  (2.*eta*(-3392. - 48.*a5tot - 24.*a6tot + 96.*eta + 123.*pi2))/(-768. + eta*(3584. + 24.*a5tot - 123.*pi2));
//  d3 = (-2.*eta*(6016. + 48.*a6tot + 3392.*eta + 24.*a5tot*(4. + eta) - 246.*pi2 - 123.*eta*pi2))/(-768. + eta*(3584. + 24.*a5tot - 123.*pi2));
//  d4 = -(eta*(-4608.*a6tot*(-4. + eta) + a5tot*(36864. + eta*(72192. - 2952.*pi2)) + eta*(2048.*(5582. + 9.*eta) - 834432.*pi2 + 15129.*pi4)))/(96.*(-768. + eta*(3584.+ 24.*a5tot - 123.*pi2)));
//  d5 = (eta*(-24.*a6tot*(1536. + eta*(-3776. + 123.*pi2)) + eta*(-2304.*a5tot2 + 96.*a5tot*(-3392. + 123.*pi2) - (-3776. + 123.*pi2)*(-3008. - 96.*eta + 123.*pi2))))/(96.*(-768. + eta*(3584. + 24.*a5tot - 123.*pi2)));


 // Num = 1.+ (u * n1);

 // Den = 1.+
   //  + u  * d1
   //  + u2 * d2
    // + u3 * d3
    // + u4 * d4
     //+ u5 * d5;

//  A=Num/Den;
//  dn1 = (160.*eta*(-828672. - 32256.*eta2 + 756.*eta*(-768. + eta*(3584. + 24.*a5 - 123.*pi2)) + eta*(5006848. + 42024.*a5 + 8064.*a6 - 174045.*pi2)))/(7.*pow(1536*logu*eta + 5.*(-768. + eta*(3584. + 24.*a5 - 123.*pi2)),2.)*u);

//  dd1 = (160.*eta*(-828672. - 32256.*eta2 + 756.*eta*(-768. + eta*(3584. + 24.*a5 - 123.*pi2)) + eta*(5006848. + 42024.*a5 + 8064.*a6 - 174045.*pi2)))/(7.*pow(1536.*logu*eta + 5.*(-768. + eta*(3584. + 24.*a5 - 123.*pi2)),2.)*u);

//  dd2 = (320.*eta*(-828672. - 32256.*eta2 + 756.*eta*(-768. + eta*(3584. + 24.*a5 - 123.*pi2)) + eta*(5006848. + 42024.*a5 + 8064.*a6 - 174045.*pi2)))/(7.*pow(1536.*logu*eta + 5.*(-768. + eta*(3584. + 24.*a5 - 123.*pi2)),2.)*u);

//  dd3 = (640.*eta*(-828672. - 32256.*eta2 + 756.*eta*(-768. + eta*(3584. + 24.*a5 - 123.*pi2)) + eta*(5006848. + 42024.*a5 + 8064.*a6 - 174045.*pi2)))/(7.*pow(1536.*logu*eta + 5.*(-768. + eta*(3584. + 24.*a5 - 123.*pi2)),2.)*u);

//  dd4 = (-320.*(-4. + eta)*eta*(-828672. - 32256.*eta2 + 756.*eta*(-768. + eta*(3584. + 24.*a5 - 123.*pi2)) + eta*(5006848. + 42024.*a5 + 8064.*a6 - 174045.*pi2)))/(7.*pow(1536.*logu*eta + 5.*(-768. + eta*(3584. + 24.*a5 - 123.*pi2)),2.)*u);

//  dd5 = (eta*(-8400.*eta*(-24.*(a6 - (4.*logu*(1751. + 756.*eta))/105.)*(1536. + eta*(-3776. + 123.*pi2))+ eta*(-2304.*pow(a5 + (64.*logu)/5.,2.) + 96.*(a5 + (64.*logu)/5.)*(-3392. + 123.*pi2) - (-3776. + 123.*pi2)*(-32.*(94. + 3.*eta) + 123.*pi2)))- (1536.*logu*eta + 5.*(-768. + eta*(3584. + 24.*a5 - 123.*pi2)))*(4128768.*logu*eta + 5.*(-2689536. + eta*(11170624. + 64512.*a5 - 380685.*pi2) - 756.*eta*(1536. + eta*(-3776. + 123.*pi2))))))/(2625.*pow(-768. + eta*(3584. + 24.*(a5 + (64.*logu)/5.) - 123.*pi2),2.)*u);

 //First derivative
//  dNum  = dn1*u + n1;
//  dDen  = d1 + u*(dd1 + 2*d2) + u2*(dd2 + 3*d3) + u3*(dd3 + 4*d4) + u4*(dd4 + 5*d5) + dd5*u5;

// Derivative of A function with respect to u
//  prefactor = A/(Num*Den);
//  dA_u = prefactor*(dNum*Den - dDen*Num);

// Derivative of A w



//  return -u2*dA_u;
//}



/**
 * Calculate the EOB D function.
 */
static REAL8 XLALCalculateEOBD( REAL8   r, /**<< Orbital separation (in units of total mass M) */
                         REAL8 eta  /**<< Symmetric mass ratio */
                       )
{
	REAL8  u, u2, u3;

	u = 1./r;
	u2 = u*u;
	u3 = u2*u;

	return 1./(1.+6.*eta*u2+2.*eta*(26.-3.*eta)*u3);
}


/**
 * Function to calculate the EOB effective Hamiltonian for the
 * given values of the dynamical variables. The coefficients in the
 * A potential function should already have been computed.
 * Note that the pr used here is the tortoise co-ordinate.
 */
static
REAL8 XLALEffectiveHamiltonian( const REAL8 eta,
				const REAL8 mass1, 
                                const REAL8 mass2,
                                const REAL8 lambda1,
                                const REAL8 lambda2,            /**<< Symmetric mass ratio */
                                const REAL8 r,            /**<< Orbital separation */
                                const REAL8 pr,           /**<< Tortoise co-ordinate */
                                const REAL8 pp           /**<< Momentum pphi */
                              /**<< Pre-computed coefficients in A function */
                              )
{

        /* The pr used in here is the tortoise co-ordinate */
        REAL8 r2, pr2, pp2, z3, eoba;

        r2   = r * r;
        pr2  = pr * pr;
        pp2  = pp * pp;

        eoba = XLALCalculateTNSEOBA( r, mass1, mass2, lambda1, lambda2 );
        z3   = 2. * ( 4. - 3. * eta ) * eta;
        return sqrt( pr2 + eoba * ( 1.  + pp2/r2 + z3*pr2*pr2/r2 ) );
}


/**
 * Function which calculates the various coefficients used in the generation
 * of the factorized waveform. These coefficients depend only on the symmetric
 * mass ratio eta. It should be noted that this function calculates the
 * coefficients used in calculating the flux. For generating the waveforms
 * themselves, the coefficients have additional terms added which are calculated
 * using XLALModifyFacWaveformCoefficients(). THe non-spinning parts of these
 * coefficients can be found in Pan et al, arXiv:1106.1021v1 [gr-qc].
 */
UNUSED static int XLALSimIMRTNSEOBCalcFacWaveformCoefficients(
          TNSFacWaveformCoeffs * const coeffs, /**<< Structure containing coefficients (populated in function) */
          const REAL8               eta     /**<< Symmetric mass ratio */
          )
{

  REAL8 eta2 = eta*eta;
  REAL8 eta3 = eta2 * eta;

  REAL8 dM, dM2; //dM3;

  REAL8 a = 0.;
  REAL8 a2 = 0.;
  REAL8 a3 = 0.;
  REAL8 chiS = 0.;
  REAL8 chiA = 0.;

  /* Combination which appears a lot */
  REAL8 m1Plus3eta, m1Plus3eta2, m1Plus3eta3;

  dM2 = 1. - 4.*eta;
  
  /* Check that deltaM has a reasonable value */
  if ( dM2 < 0 )
  {
    XLALPrintError( "eta seems to be < 0.25 - this isn't allowed!\n" );
    XLAL_ERROR( XLAL_EINVAL );
  }

  dM  = sqrt( dM2 );
  //dM3 = dM2 * dM;

  m1Plus3eta  = - 1. + 3.*eta;
  m1Plus3eta2 = m1Plus3eta * m1Plus3eta;
  m1Plus3eta3 = m1Plus3eta * m1Plus3eta2;

  /* Initialize all coefficients to zero */
  /* This is important, as we will not set some if dM is zero */
  memset( coeffs, 0, sizeof( TNSFacWaveformCoeffs ) );


  /* l = 2 */

  coeffs->delta22vh3 = 7./3.;
  coeffs->delta22vh6 = (-4.*a)/3. + (428.*LAL_PI)/105.;
  coeffs->delta22v8 = (20.*a)/63.;
  coeffs->delta22vh9 = -2203./81. + (1712.*LAL_PI*LAL_PI)/315.;
  coeffs->delta22v5  = - 24.*eta;

  coeffs->rho22v2   = -43./42. + (55.*eta)/84.;
  coeffs->rho22v3   = (-2.*(chiS + chiA*dM - chiS*eta))/3.;
  coeffs->rho22v4   = -20555./10584. + (chiS*chiS + 2.*chiA*chiS*dM + chiA*chiA*dM2)/2.
       - (33025.*eta)/21168. + (19583.*eta2)/42336.;
  coeffs->rho22v5   = (-34.*a)/21.;
  coeffs->rho22v6   = 1556919113./122245200. + (89.*a2)/252. - (48993925.*eta)/9779616. 
       - (6292061.*eta2)/3259872. + (10620745.*eta3)/39118464.
       + (41.*eta*LAL_PI*LAL_PI)/192.;
  coeffs->rho22v6l  = - 428./105.;
  coeffs->rho22v7   = (18733.*a)/15876. + a*a2/3.;
   coeffs->rho22v8   = -387216563023./160190110080. + (18353.*a2)/21168. - a2*a2/8.;
  coeffs->rho22v8l  =  9202./2205.;
  coeffs->rho22v10  = -16094530514677./533967033600.;
  coeffs->rho22v10l =  439877./55566.;


//  coeffs->rho22v2   = -1.0238095238095237 + (0.6547619047619048* eta);
//  coeffs->rho22v4   = -1.94208238851096 - (1.5601379440665155*eta) + (0.4625614134542706*eta2);
//  coeffs->rho22v6   = 12.736034731834051 - (2.902228713904598*eta) - (1.9301558466099282*eta2) + (0.2715020968103451*eta3);
//  coeffs->rho22v6l  = - 4.076190476190476;
//  coeffs->rho22v8   = -2.4172313935587004;
//  coeffs->rho22v8l  =  4.173242630385488;
//  coeffs->rho22v10  = -30.14143102836864;
//  coeffs->rho22v10l =  7.916297736025627;




  if ( dM2 )
  {
    coeffs->delta21vh3 = 2./3.;
    coeffs->delta21vh6 = (-17.*a)/35. + (107.*LAL_PI)/105.;
    coeffs->delta21vh7 = (3.*a2)/140.;
    coeffs->delta21vh9 = -272./81. + (214.*LAL_PI*LAL_PI)/315.;
    coeffs->delta21v5  = - 493. * eta /42.;

    coeffs->rho21v1   = (-3.*(chiS+chiA/dM))/(4.);
    //coeffs->rho21v2   = -59./56 - (9.*chiAPlusChiSdM*chiAPlusChiSdM)/(32.*dM2) + (23.*eta)/84.;
    /*coeffs->rho21v3   = (-567.*chiA*chiA*chiA - 1701.*chiA*chiA*chiS*dM
                        + chiA*(-4708. + 1701.*chiS*chiS - 2648.*eta)*(-1. + 4.*eta)
                        + chiS* dM3 *(4708. - 567.*chiS*chiS
                        + 1816.*eta))/(2688.*dM3);*/
    coeffs->rho21v2   = -59./56. + (23.*eta)/84. - 9./32.*a2;
    coeffs->rho21v3   = 1177./672.*a - 27./128.*a3;
    coeffs->rho21v4   = -47009./56448.- (865.*a2)/1792. - (405.*a2*a2)/2048. - (10993.*eta)/14112.
                        + (617.*eta2)/4704.;
    coeffs->rho21v5   = (-98635.*a)/75264. + (2031.*a*a2)/7168. - (1701.*a2*a3)/8192.;
    coeffs->rho21v6   = 7613184941./2607897600.+ (9032393.*a2)/1806336. + (3897.*a2*a2)/16384.
                        - (15309.*a3*a3)/65536.; 
    coeffs->rho21v6l  = - 107./105.;
    coeffs->rho21v7   = (-3859374457.*a)/1159065600. - (55169.*a3)/16384.
                        + (18603.*a2*a3)/65536. - (72171.*a2*a2*a3)/262144.;
    coeffs->rho21v7l  =  107.*a/140.;
    coeffs->rho21v8   = -1168617463883./911303737344.;
    coeffs->rho21v8l  = 6313./5880.;
    coeffs->rho21v10  = -63735873771463./16569158860800.; 
    coeffs->rho21v10l = 5029963./5927040.;
  }

  /* l = 3 */
  if ( dM2 )
  {
    coeffs->delta33vh3 = 13./10.;
    coeffs->delta33vh6 = (-81.*a)/20. + (39.*LAL_PI)/7.;
    coeffs->delta33vh9 = -227827./3000. + (78.*LAL_PI*LAL_PI)/7.;
    coeffs->delta33v5  = - 80897.*eta / 2430.;

    coeffs->rho33v2 = -7./6. + (2.*eta)/3.;
    coeffs->rho33v3 = (chiS*dM*(-4. + 5.*eta) + chiA*(-4. + 19.*eta))/(6.*dM);
    coeffs->rho33v4 = -6719./3960. + a2/2. - (1861.*eta)/990. + (149.*eta2)/330.;
    coeffs->rho33v5 = (-4.*a)/3.;
    coeffs->rho33v6 = 3203101567./227026800. + (5.*a2)/36.;
    coeffs->rho33v6l = - 26./7.;
    coeffs->rho33v7 = (5297.*a)/2970. + a*a2/3.;
    coeffs->rho33v8 = -57566572157./8562153600.;
    coeffs->rho33v8l = 13./3.;

    coeffs->rho33v10 = 903823148417327./30566888352000.;
    coeffs->rho33v10l = 87347./13860.;




  }

  coeffs->delta32vh3 = (10. + 33.*eta)/(-15.*m1Plus3eta);
  coeffs->delta32vh4 = 4.*a;
  coeffs->delta32vh6 = (-136.*a)/45. + (52.*LAL_PI)/21.;
  coeffs->delta32vh9 = -9112./405. + (208.*LAL_PI*LAL_PI)/63.;

  coeffs->rho32v   = (4.*chiS*eta)/(-3.*m1Plus3eta);
  coeffs->rho32v2  = (-4.*a2*eta2)/(9.*m1Plus3eta2) + (328. - 1115.*eta
                        + 320.*eta2)/(270.*m1Plus3eta);
  coeffs->rho32v3  = (2.*(45.*a*m1Plus3eta3
                        - a*eta*(328. - 2099.*eta + 5.*(733. + 20.*a2)*eta2
                        - 960.*eta3)))/(405.*m1Plus3eta3);
  coeffs->rho32v4  = a2/3. + (-1444528.
                        + 8050045.*eta - 4725605.*eta2 - 20338960.*eta3
                        + 3085640.*eta2*eta2)/(1603800.*m1Plus3eta2);
  coeffs->rho32v5  = (-2788.*a)/1215.;
  coeffs->rho32v6  = 5849948554./940355325. + (488.*a2)/405.;
  coeffs->rho32v6l =  - 104./63.;
  coeffs->rho32v8  = -10607269449358./3072140846775.;
  coeffs->rho32v8l = 17056./8505.;

  if ( dM2 )
  {
    coeffs->delta31vh3 = 13./30.;
    coeffs->delta31vh6 = (61.*a)/20. + (13.*LAL_PI)/21.;
    coeffs->delta31vh7 = (-24.*a2)/5.;
    coeffs->delta31vh9 = -227827./81000. + (26.*LAL_PI*LAL_PI)/63.;
    coeffs->delta31v5  = - 17.*eta/10.; 
 
    coeffs->rho31v2  = -13./18. - (2.*eta)/9.;
    coeffs->rho31v3  = (chiA*(-4. + 11.*eta) + chiS*dM*(-4. + 13.*eta))/(6.*dM);
    coeffs->rho31v4  = 101./7128.
                        - (5.*a2)/6. - (1685.*eta)/1782. - (829.*eta2)/1782.;
    coeffs->rho31v5  = (4.*a)/9.;
    coeffs->rho31v6  = 11706720301./6129723600. - (49.*a2)/108.;
    coeffs->rho31v6l =  - 26./63.;
    coeffs->rho31v7  = (-2579.*a)/5346. + a*a2/9.;
    coeffs->rho31v8  = 2606097992581./4854741091200.;
    coeffs->rho31v8l = 169./567.;


   coeffs->rho31v10 = 430750057673539./297110154781440.;
   coeffs->rho31v10l = -1313./224532.;
  }

  /* l = 4 */
  
  coeffs->delta44vh3 = (112. + 219.*eta)/(-120.*m1Plus3eta);
  coeffs->delta44vh6 = (-464.*a)/75. + (25136.*LAL_PI)/3465.;

  coeffs->rho44v2 = (1614. - 5870.*eta + 2625.*eta2)/(1320.*m1Plus3eta);
  coeffs->rho44v3 = (chiA*(10. - 39.*eta)*dM + chiS*(10. - 41.*eta
                        + 42.*eta2))/(15.*m1Plus3eta);
  coeffs->rho44v4 = a2/2. + (-511573572.
                        + 2338945704.*eta - 313857376.*eta2 - 6733146000.*eta3
                        + 1252563795.*eta2*eta2)/(317116800.*m1Plus3eta2);
  coeffs->rho44v5 = (-69.*a)/55.;
  coeffs->rho44v6 = 16600939332793./1098809712000. + (217.*a2)/3960.;
  coeffs->rho44v6l = - 12568./3465.;

  coeffs->rho44v8 = 845198./190575.;
  coeffs->rho44v8l= -172066910136202271./19426955708160000.;

  if ( dM2 )
  {
    coeffs->delta43vh3 = (486. + 4961.*eta)/(810.*(1. - 2.*eta));
    coeffs->delta43vh4 = (11.*a)/4.;
    coeffs->delta43vh6 = 1571.*LAL_PI/385.;

    coeffs->rho43v   = (5.*(chiA - chiS*dM)*eta)/(8.*dM*(-1. + 2.*eta));
    coeffs->rho43v2  = (222. - 547.*eta + 160.*eta2)/(176.*(-1. + 2.*eta));
    coeffs->rho43v4  = -6894273./7047040. + (3.*a2)/8.;
    coeffs->rho43v5  = (-12113.*a)/6160.;
    coeffs->rho43v6  = 1664224207351./195343948800.;
    coeffs->rho43v6l = - 1571./770.;

    coeffs->rho43v8 = -2465107182496333./460490801971200.;
    coeffs ->rho43v8l = -174381./67760;

  }

  coeffs->delta42vh3 = (7.*(1. + 6.*eta))/(-15.*m1Plus3eta);
  coeffs->delta42vh6 = (212.*a)/75. + (6284.*LAL_PI)/3465.;

  coeffs->rho42v2  = (1146. - 3530.*eta + 285.*eta2)/(1320.*m1Plus3eta);
  coeffs->rho42v3  = (chiA*(10. - 21.*eta)*dM + chiS*(10. - 59.*eta
                        + 78.*eta2))/(15.*m1Plus3eta);
  coeffs->rho42v4  = a2/2. + (-114859044. + 295834536.*eta + 1204388696.*eta2 - 3047981160.*eta3
                        - 379526805.*eta2*eta2)/(317116800.*m1Plus3eta2);
  coeffs->rho42v5  = (-7.*a)/110.;
  coeffs->rho42v6  = 848238724511./219761942400. + (2323.*a2)/3960.;
  coeffs->rho42v6l = - 3142./3465.;

  coeffs->rho42v8 = -12864377174485679./19426955708160000.;
  coeffs -> rho42v8l =300061./381150.;

  if ( dM2 )
  {
    coeffs->delta41vh3 = (2. + 507.*eta)/(10.*(1. - 2.*eta));
    coeffs->delta41vh4 = (11.*a)/12.;
    coeffs->delta41vh6 = 1571.*LAL_PI/3465.;

    coeffs->rho41v   = (5.*(chiA - chiS*dM)*eta)/(8.*dM*(-1. + 2.*eta));
    coeffs->rho41v2  = (602. - 1385.*eta + 288.*eta2)/(528.*(-1. + 2.*eta));
    coeffs->rho41v4  = -7775491./21141120. + (3.*a2)/8.;
    coeffs->rho41v5  = (-20033.*a)/55440. - (5*a*a2)/6.;
    coeffs->rho41v6  = 1227423222031./1758095539200.;
    coeffs->rho41v6l = - 1571./6930.;


   coeffs->rho41v8 = -9584392078751453./37299754959667200.;
   coeffs->rho41v8l= 67553./261360.;
  }

  /* l = 5 */
  if ( dM2 )
  {
    coeffs->delta55vh3 = (96875. + 857528.*eta)/(131250.*(1 - 2*eta));

    coeffs->rho55v2 = (487. - 1298.*eta + 512.*eta2)/(390.*(-1. + 2.*eta));
    coeffs->rho55v3 = (-2.*a)/3.;
    coeffs->rho55v4 = -3353747./2129400. + a2/2.;
    coeffs->rho55v5 = - 241. * a / 195.;

    coeffs->rho55v6= 190606537999247./11957879934000.;
    coeffs -> rho55v6l = -1546./429.;
    coeffs ->rho55v8 = -1213641959949291437./118143853747920000.;
    coeffs ->rho55v8l = 376451./83655.;
  }

  coeffs->delta54vh3 = 8./15.;
  coeffs->delta54vh4 = 12.*a/5.;

  coeffs->rho54v2 = (-17448. + 96019.*eta - 127610.*eta2
                        + 33320.*eta3)/(13650.*(1. - 5.*eta + 5.*eta2));
  coeffs->rho54v3 = (-2.*a)/15.;
  coeffs->rho54v4 = -16213384./15526875. + (2.*a2)/5.;

 
  coeffs->rho54v6 = 6704294638171892./653946558890625.;
  coeffs ->rho54v6l = -24736./10725.;
   
  if ( dM2 )
  {
    coeffs->delta53vh3 = 31./70.;

    coeffs->rho53v2 = (375. - 850.*eta + 176.*eta2)/(390.*(-1. + 2.*eta));
    coeffs->rho53v3 = (-2.*a)/3.;
    coeffs->rho53v4 = -410833./709800. + a2/2.;
    coeffs->rho53v5 = - 103.*a/325.;


    coeffs->rho53v6 =7618462680967./1328653326000.;
    coeffs->rho53v6l=-4638./3575.;
    coeffs->rho53v8= -77082121019870543./39381284582640000.;
    coeffs->rho53v8l=2319./1859.;
  }

  coeffs->delta52vh3 = 4./15.;
  coeffs->delta52vh4 = 6.*a/5.;

  coeffs->rho52v2 = (-15828. + 84679.*eta - 104930.*eta2
                        + 21980.*eta3)/(13650.*(1. - 5.*eta + 5.*eta2));
  coeffs->rho52v3 = (-2.*a)/15.;
  coeffs->rho52v4 = -7187914./15526875. + (2.*a2)/5.;

  coeffs->rho52v6 = 1539689950126502./653946558890625.;
  coeffs->rho52v6l = -6184./10725.; 

  if ( dM2 )
  {
    coeffs->delta51vh3 = 31./210.;

    coeffs->rho51v2 = (319. - 626.*eta + 8.*eta2)/(390.*(-1. + 2.*eta));
    coeffs->rho51v3 = (-2.*a)/3.;
    coeffs->rho51v4 = -31877./304200. + a2/2.;
    coeffs->rho51v5 = 139.*a/975.;

    coeffs ->rho51v6 = 7685351978519./11957879934000.;
    coeffs-> rho51v6l = -1546./10725.;
    coeffs->rho51v8 = -821807362819271./10740350340720000.;
    coeffs->rho51v8l= 22417./190125.;

  }

  /* l = 6 */

  coeffs->delta66vh3 = 43./70.;
  
  coeffs->rho66v2 = (-106. + 602.*eta - 861.*eta2
                        + 273.*eta3)/(84.*(1. - 5.*eta + 5.*eta2));
  coeffs->rho66v3 = (-2.*a)/3.;
  coeffs->rho66v4 = -1025435./659736. + a2/2.;

  coeffs -> rho66v6 =  610931247213169./36701493028200.;
  coeffs ->rho66v6l = -3604./1001.;

  if ( dM2 )
  {
    coeffs->delta65vh3 = 10./21.;
    
    coeffs->rho65v2 = (-185. + 838.*eta - 910.*eta2
                        + 220.*eta3)/(144.*(dM2 + 3.*eta2));
    coeffs->rho65v3 = - 2.*a/9.;



    coeffs->rho65v4 =-5954065./54286848.;
    coeffs ->rho65v6 = 67397117912549267./5798416452820992.;
    coeffs -> rho65v6l = -22525./9009.;    
  }

  coeffs->delta64vh3 = 43./105.;
  
  coeffs->rho64v2 = (-86. + 462.*eta - 581.*eta2
                        + 133.*eta3)/(84.*(1. - 5.*eta + 5.*eta2));
  coeffs->rho64v3 = (-2.*a)/3.;
  coeffs->rho64v4 = -476887./659736. + a2/2.;
  


  coeffs-> rho64v6 = 180067034480351./24467662018800.;
  coeffs-> rho64v6l = -14416./9009.;

  if ( dM2 )
  {
    coeffs->delta63vh3 = 2./7.;

    coeffs->rho63v2 = (-169. + 742.*eta - 750.*eta2
                        + 156.*eta3)/(144.*(dM2 + 3.*eta2));
    coeffs->rho63v3 = - 2.*a/9.;


    coeffs->rho63v4 = -152153941./271434240.;
    coeffs->rho63v6 = 116042497264681103./28992082264104960.;
    coeffs-> rho63v6l = -901./1001.;

  }

  coeffs->delta62vh3 = 43./210.;

  coeffs->rho62v2 = (-74. + 378.*eta - 413.*eta2
                        + 49.*eta3)/(84.*(1. - 5.*eta + 5.*eta2));
  coeffs->rho62v3 = (-2.*a)/3.;
  coeffs->rho62v4 = -817991./3298680. + a2/2.;


  coeffs->rho62v6 = 812992177581./453104852200.;
  coeffs->rho62v6l = -3604./9009.;


  if ( dM2 )
  {
    coeffs->delta61vh3 = 2./21.;

    coeffs->rho61v2 = (-161. + 694.*eta - 670.*eta2
                        + 124.*eta3)/(144.*(dM2 + 3.*eta2));
    coeffs->rho61v3 = - 2. * a / 9.;



    coeffs -> rho61v4 = -79192261./271434240.;
    coeffs -> rho61v6 = 6277796663889319./28992082264104960.;
    coeffs -> rho61v6l = -901./9009.;
  }

  /* l = 7 */
  if ( dM2 )
  {
    coeffs->delta77vh3 = 19./36.;

    coeffs->rho77v2 = (-906. + 4246.*eta - 4963.*eta2
                        + 1380.*eta3)/(714.*(dM2 + 3.*eta2));
    coeffs->rho77v3 = - 2.*a/3.;

    
    coeffs -> rho77v4 = -32358125./20986602.;
    coeffs ->rho77v6 = 66555794049401803./3856993267327200;
    coeffs ->rho77v6l = -11948./3315.;
  }

  coeffs->rho76v2 = (2144. - 16185.*eta + 37828.*eta2 - 29351.*eta3
                        + 6104.*eta2*eta2) / (1666.*(-1 + 7*eta - 14*eta2
                        + 7*eta3));




  coeffs -> rho76v4 = -195441224./171390583.;


  if ( dM2 )
  {
    coeffs->delta75vh3 = 95./252.;

    coeffs->rho75v2 = (-762. + 3382.*eta - 3523.*eta2
                        + 804.*eta3)/(714.*(dM2 + 3.*eta2));
    coeffs->rho75v3 = - 2.*a/3.;



    coeffs ->rho75v4 = -17354227./20986602.;
    coeffs ->rho75v6 = 192862646381533./22039961527584.;
    coeffs ->rho75v6l = - 59740./32487.;
  }

  coeffs->rho74v2 = (17756. - 131805.*eta + 298872.*eta2 - 217959.*eta3
                        + 41076.*eta2*eta2) / (14994.*(-1. + 7.*eta - 14.*eta2
                        + 7.*eta3));


 coeffs->rho74v4 = -2995755988./4627545741.;


  if ( dM2 )
  {
    coeffs->delta73vh3 = 19./84.;

    coeffs->rho73v2 = (-666. + 2806.*eta - 2563.*eta2
                        + 420.*eta3)/(714.*(dM2 + 3.*eta2));
    coeffs->rho73v3 = - 2.*a/3.;


    coeffs->rho73v4 = -7804375./20986602.;
    coeffs->rho73v6 = 1321461327981547./428554807480800.;
    coeffs -> rho73v6l = -35844./54145.;
  }

  coeffs->rho72v2 = (16832. - 123489.*eta + 273924.*eta2 - 190239.*eta3
                        + 32760.*eta2*eta2) /(14994.*(-1. + 7.*eta - 14.*eta2
                        + 7.*eta3));


  coeffs->rho72v4 = -1625746984./4627545741.;

  if ( dM2 )
  {
    coeffs->delta71vh3 = 19./252.;

    coeffs->rho71v2 = (-618. + 2518.*eta - 2083.*eta2
                        + 228.*eta3)/(714.*(dM2 + 3.*eta2));
    coeffs->rho71v3 = - 2.*a/3.;
    
  

    coeffs ->rho71v4 = - 1055091./6995534.;
    coeffs ->rho71v6 = 142228318411021./550999038189600.;
    coeffs ->rho71v6l = -11948./162435.;

  }

  /* l = 8 */
  
  coeffs->rho88v2 = (3482. - 26778.*eta + 64659.*eta2 - 53445.*eta3
                        + 12243.*eta2*eta2) / (2736.*(-1. + 7.*eta - 14.*eta2
                        + 7.*eta3));


  coeffs->rho88v4 = -1.5337092502821381;

  if ( dM2 )
  {
    coeffs->rho87v2 = (23478. - 154099.*eta + 309498.*eta2 - 207550.*eta3
                        + 38920*eta2*eta2) / (18240.*(-1 + 6*eta - 10*eta2
                        + 4*eta3));

    coeffs->rho87v4 = -1.175404252991305;
  }

  coeffs->rho86v2 = (1002. - 7498.*eta + 17269.*eta2 - 13055.*eta3
                        + 2653.*eta2*eta2) / (912.*(-1. + 7.*eta - 14.*eta2
                        + 7.*eta3));

  coeffs -> rho86v4 = - 0.9061610303170207;

  if ( dM2 )
  {
    coeffs->rho85v2 = (4350. - 28055.*eta + 54642.*eta2 - 34598.*eta3
                        + 6056.*eta2*eta2) / (3648.*(-1. + 6.*eta - 10.*eta2
                        + 4.*eta3));

    coeffs -> rho85v4 = -0.7220789990670207 ;
  }

  coeffs->rho84v2 = (2666. - 19434.*eta + 42627.*eta2 - 28965.*eta3
                        + 4899.*eta2*eta2) / (2736.*(-1. + 7.*eta - 14.*eta2
                        + 7.*eta3));

  coeffs -> rho84v4 = -0.47652059150068155;

  if ( dM2 )
  {
    coeffs->rho83v2 = (20598. - 131059.*eta + 249018.*eta2 - 149950.*eta3
                        + 24520.*eta2*eta2) / (18240.*(-1. + 6.*eta - 10.*eta2
                        + 4.*eta3));

   coeffs -> rho83v4 = -0.4196774909106648;
  }

  coeffs->rho82v2 = (2462. - 17598.*eta + 37119.*eta2 - 22845.*eta3
                        + 3063.*eta2*eta2) / (2736.*(-1. + 7.*eta - 14.*eta2
                        + 7.*eta3));

  coeffs ->rho82v4 = -0.2261796441029474 ;

  if ( dM2 )
  {
    coeffs->rho81v2 = (20022. - 126451.*eta + 236922.*eta2 - 138430.*eta3
                        + 21640.*eta2*eta2) / (18240.*(-1. + 6.*eta - 10.*eta2
                        + 4.*eta3));

    coeffs->rho81v4 = -0.26842133517043704 ;
  }

  /* All relevant coefficients should be set, so we return */

  return XLAL_SUCCESS;
}


/**
 * Function which adds the additional terms required for waveform generation
 * to the factorized waveform coefficients. Note that this function only calculates
 * additional terms not present in the flux, so the factorized waveform coefficients
 * SHOULD ALREADY HAVE BEEN CALCULATED using XLALCalcFacWaveformCoefficients() prior
 * to calling this function.
 */
//UNUSED static int XLALSimIMREOBModifyFacWaveformCoefficients( 
//                                       TNSFacWaveformCoeffs * const coeffs, /**<< Structure containing coefficients */
//                                       const REAL8 eta                   /**<< Symmetric mass ratio */
//                                     )
//{
//
//  if ( !coeffs )
//  {
//    XLAL_ERROR( XLAL_EINVAL );
//  }

  /* Tweak the relevant coefficients for the generation of the waveform */
//  coeffs->rho21v6 += -5. * eta;
//  coeffs->rho33v6 += -20. * eta;
//  coeffs->rho44v6 += -15. * eta;
//  coeffs->rho55v6 += 4. * eta;

//  coeffs->delta21v7 += 30. * eta;
//  coeffs->delta33v7 += -10. * eta;
//  coeffs->delta44v5 += -70. * eta;
//  coeffs->delta55v5 += 40. * eta;

//  return XLAL_SUCCESS;
//}

/**
 * Computes the non-Keplerian correction to the velocity as determined from the
 * frequency obtained assuming a circular orbit. In the early stages of the evolution,
 * this should be a number close to 1.
 */
static REAL8
nonKeplerianCoefficient(
                   REAL8Vector * restrict values, /**<< Dynamics r, phi, pr, pphi */
                   const REAL8       eta,         /**<< Symmetric mass ratio */
                  const REAL8       mass1,
                  const REAL8       mass2,
                  const REAL8       lambda1,
                  const REAL8       lambda2                                         
 // EOBACoefficients *coeffs       /**<< Pre-computed A coefficients */
                   )
{

  REAL8 r    = values->data[0];
  REAL8 pphi = values->data[3];

  REAL8 A  = XLALCalculateTNSEOBA( r, mass1,mass2,lambda1, lambda2 );
  REAL8 dA = XLALCalculateEOBdAdr( r, eta );

  return 2. * (1. + 2. * eta * ( -1. + sqrt( (1. + pphi*pphi/(r*r)) * A ) ) )
          / ( r*r * dA );
}

UNUSED static REAL8
nonKeplerianCoefficient_Debugging(
                   REAL8Vector * restrict values, /**<< Dynamics r, phi, pr, pphi */
                   const REAL8       eta,         /**<< Symmetric mass ratio */
                  const REAL8       mass1,
                  const REAL8       mass2,
                  const REAL8       lambda1,
                  const REAL8       lambda2                                         
 // EOBACoefficients *coeffs       /**<< Pre-computed A coefficients */
                   )
{

  REAL8 r    = values->data[0];
  REAL8 pphi = values->data[3];

  REAL8 A  = XLALCalculateTNSEOBA( r, mass1,mass2,lambda1, lambda2 );
  REAL8 A_nontidal = XLALCalculateTNSEOBA_nontidal_Debugging(r,eta ); 
   printf("\n A_nontidal :%.10f",A_nontidal);
  REAL8 dA = XLALCalculateEOBdAdr( r, eta );
  printf("\n A :%.10f",A);
  printf("\n dA : %.10f", dA);
  REAL8 W=( -1. + sqrt( (1. + pphi*pphi/(r*r)) * A ) );
  printf("\n W : %.10f",W);
  return 2. * (1. + 2. * eta * ( -1. + sqrt( (1. + pphi*pphi/(r*r)) * A ) ) )
          / ( r*r * dA );
}

/**
 * Computes the factorized waveform according to the prescription
 * given in Pan et al, arXiv:1106.1021v1 [gr-qc], for a given
 * mode l,m, for the given values of the dynamics at that point.
 * The function returns XLAL_SUCCESS if everything works out properly,
 * otherwise XLAL_FAILURE will be returned.

*/

static int
XLALSimIMRTNSEOBCalculateTidal_hlm(
                 COMPLEX16 *hlm_tidal, /**<< OUTPUT, Newtonian multipole */
                 REAL8 v,              /**<< Dimensionless parameter \f$\equiv v^2\f$ */
                 UNUSED REAL8 r,       /**<< Orbital separation (units of total mass M) */
                 UNUSED  REAL8 phi,            /**<< Orbital phase (in radians) */
                 UINT4  l,             /**<< Mode l */
                 INT4  m,              /**<< Mode m */
                 TNSEOBParams *params,
                 REAL8Vector * restrict values     /**<< Pre-computed coefficients, parameters, etc. */
                 )
{

//==========IMPORTANT NOTE==================
//Please check what v,vPhi etc whould be passed to this. Right now in v this function gets vPhi. Once the NONtidal part works, you have to CHECK the what vel goes to the tidal parts
//===========================================

REAL8 v2,v3,v10,X1,X2,m1,m2,eta,lambda1,lambda2,Omega;
REAL8 star1k2,star2k2,betaA_22,betaB_22,prefactorA_22,prefactorB_22,totalMass,x1,x2;
INT4 epsilon;
INT4 sign;
REAL8 vPhi;
REAL8 c,tidal_peice;
INT4 status;
COMPLEX16 hNewton;
  /* Pre-computed coefficients */
//TNSFacWaveformCoeffs *hCoeffs = params->hCoeffs;


eta = params->eta;
m1 = params->mass1;
m2 = params->mass2;
lambda1 = params -> lambda1;
lambda2 = params -> lambda2;

//printf("Calculating the tidal contribution to h_lm.....");
v2=v*v;
v3=v2*v;
v10=v3*v3*v2*v2;

X1=m1/(m1+m2);
X2=m2/(m1+m2);

if((lambda2==0.)&&(lambda1==0.))
 { 
  //printf("Both the lambdas are zero. hlm_tidal to zero");
  *hlm_tidal=(COMPLEX16)0.;
 }
else
{
star1k2=3.*(X2/X1)*lambda1*(X1*X1*X1*X1*X1);
star2k2=3.*(X1/X2)*lambda2*(X2*X2*X2*X2*X2);

betaA_22 = (-202. + 560.*X1 - 340.*X1*X1 + 45*X1*X1*X1)/(42.*(3.-2.*X1));
betaB_22 = (-202. + 560.*X2 - 340.*X2*X2 + 45*X2*X2*X2)/(42.*(3.-2.*X2));

prefactorA_22=star1k2*(X1/X2 + 3.)*v10;
prefactorB_22=star2k2*(X2/X1 + 3.)*v10;


//Calculating C_(l+eps).... I have written this structure this way so we can add higher mode correction if needed. But right now this WHOLE function gets called only for l=2... so, basically we could hardcode l=2 if we need it to be faster.  
   
  totalMass = m1 + m2;

   epsilon = ( l + m )  % 2;

   x1 = m1 / totalMass;
   x2 = m2 / totalMass;

   eta = m1*m2/(totalMass*totalMass);

   if  ( abs( m % 2 ) == 0 )
   {
     sign = 1;
   }
   else
   {
     sign = -1;
   }
   if  ( m1 != m2 || sign == 1 )
   {
     c = pow( x2, l + epsilon - 1 ) + sign * pow(x1, l + epsilon - 1 );
   }
   else
   {
     switch( l )
     {
       case 2:
         c = -1.0;
         break;
       case 3:
         c = -1.0;
         break;
       case 4:
         c = -0.5;
         break;
       default:
         c = 0.0;
         break;
     }
   }

  Omega = v2 * v;

//Calculating Newtonian Prefactor PSI
 vPhi = nonKeplerianCoefficient( values, eta,m1,m2,lambda1,lambda2 );
 
  /* Assign rOmega value temporarily to vPhi */
  
 vPhi  = r * cbrt(vPhi);
  /* Assign rOmega * Omega to vPhi */
  vPhi *= Omega;
status = XLALSimIMRTNSEOBCalculateNewtonianMultipole( &hNewton, vPhi * vPhi, vPhi/Omega,
            values->data[1], (UINT4)l, m, params );
  if ( status == XLAL_FAILURE )
  {
    XLAL_ERROR( XLAL_EFUNC );
  }

hNewton=hNewton/(COMPLEX16)c;

tidal_peice=((prefactorA_22*(1.+betaA_22*v2))+(prefactorB_22*(1.+betaB_22*v2)));

*hlm_tidal =hNewton*(COMPLEX16)tidal_peice;


}
  return XLAL_SUCCESS;
}

UNUSED static int  XLALSimIMRTNSEOBGetFactorizedWaveform( 
                                COMPLEX16   * restrict hlm,    /**<< The value of hlm (populated by the function) */
                                REAL8Vector * restrict values, /**<< Vector containing dynamics r, phi, pr, pphi for a given point */
                                const REAL8 v,                 /**<< Velocity (in geometric units) */
                                const INT4  l,                 /**<< Mode l */
                                const INT4  m,                 /**<< Mode m */
                                TNSEOBParams   * restrict params 
                                )
{

//Notice :  I am feeding in v=v_phi (the non keplerian velocity for all the parts of h_lm unlike in EOBNR)


  /* Status of function calls */
  INT4 status;
  INT4 i;

  REAL8 eta, mass1, mass2, lambda1, lambda2, pi,pi2,eta2;
  REAL8 r, pr, pp, Omega, v2, vh, vh3,v3, k, hathatk, eulerlogxabs,v4,v5;
  REAL8 Hreal, Heff, Slm, deltalm, rholm, rholmPwrl;
  REAL8 delta22_leading,delta21_leading,delta33_leading,delta31_leading,delta22_Num,delta21_Num,delta33_Num,delta31_Num,delta22_Den, delta21_Den, delta33_Den,delta31_Den;
  COMPLEX16 Tlm;
  COMPLEX16 hNewton,hlm_tidal;
  gsl_sf_result lnr1, arg1, z2;

  /* Non-Keplerian velocity */
  REAL8 psi,rOmega;
  /* Pre-computed coefficients */
  TNSFacWaveformCoeffs *hCoeffs = params->hCoeffs;

  if ( abs(m) > (INT4) l )
  {
    XLAL_ERROR( XLAL_EINVAL );
  }


  eta = params->eta;
  mass1 = params->mass1;
  mass2 = params->mass2;
  lambda1 = params -> lambda1;
  lambda2 = params -> lambda2;

  /* Check our eta was sensible */
  if ( eta > 0.25 )
  {
    XLALPrintError("Eta seems to be > 0.25 - this isn't allowed!\n" );
    XLAL_ERROR( XLAL_EINVAL );
  }
  else if ( eta == 0.25 && m % 2 )
  {
    /* If m is odd and dM = 0, hLM will be zero */
    memset( hlm, 0, sizeof( COMPLEX16 ) );
    return XLAL_SUCCESS;
  }

  r  = values->data[0];
  pr = values->data[2];
  pp = values->data[3];

  //printf("accesing the data r,pr,phi");
//  printf("\n%f",r);
  Heff  = XLALEffectiveHamiltonian( eta, mass1, mass2, lambda1,lambda2, r, pr, pp);
  Hreal = sqrt( 1.0 + 2.0 * eta * ( Heff - 1.0) );
  v2    = v * v;
  v3  = v*v2;
  v4 = v2*v2;
  v5 = v3*v2;
  Omega = v2 * v;
  vh3   = Hreal * Omega;
  vh    = cbrt(vh3);
  eulerlogxabs = LAL_GAMMA + log( 2.0 * (REAL8)m * v );
  pi=LAL_PI;
  pi2=pi*pi;
  eta2=eta*eta;

  /* Calculate the non-Keplerian velocity */
  /* given by Eq. (18) of Pan et al, PRD84, 124052(2011) */
  /* psi given by Eq. (19) of Pan et al, PRD84, 124052(2011) */
  /* Assign temporarily to vPhi */
  //nonKeplerianCoefficient gives psi
  psi = nonKeplerianCoefficient( values, eta,mass1,mass2,lambda1,lambda2 );
  /* Assign rOmega value temporarily to vPhi */
  rOmega  = r * cbrt(psi);


  /* Calculate the newtonian multipole */
  status = XLALSimIMRTNSEOBCalculateNewtonianMultipole( &hNewton, v * v, r,
            values->data[1], (UINT4)l, m, params );
  if ( status == XLAL_FAILURE )
  {
    XLAL_ERROR( XLAL_EFUNC );
  }
  
  //printf("I calculated newtonian fl");

 /*Calculate the tidal part of h_lm*/

  status = XLALSimIMRTNSEOBCalculateTidal_hlm( &hlm_tidal, v , v/Omega,
            values->data[1], (UINT4)l, m, params, values );
  if ( status == XLAL_FAILURE )
  {
    XLAL_ERROR( XLAL_EFUNC );
  }

  /* Calculate the source term */
  if ( ( (l+m)%2 ) == 0)
  {
    Slm = Heff;
  }
  else
  {

  //In gets changed in TNSEOB! 
  // We have  Slm = v * pp in EOBNR but in TNSEOB. 
  Slm = pp/(rOmega*v);
  }

  /* Calculate the Tail term */
  k  = m * Omega;
  hathatk = Hreal * k;
  XLAL_CALLGSL( status = gsl_sf_lngamma_complex_e( l+1.0, -2.0*hathatk, &lnr1, &arg1 ) );
  if (status != GSL_SUCCESS)
  {
    XLALPrintError("Error in GSL function\n" );
    XLAL_ERROR( XLAL_EFUNC );
  }
  XLAL_CALLGSL( status = gsl_sf_fact_e( l, &z2 ) );
  if ( status != GSL_SUCCESS)
  {
    XLALPrintError("Error in GSL function\n" );
    XLAL_ERROR( XLAL_EFUNC );
  }
  Tlm = cexp( ( lnr1.val + LAL_PI * hathatk ) + I * (
        arg1.val + 2.0 * hathatk * log(4.0*k/sqrt(LAL_E)) ) );
  Tlm /= z2.val;

  /* Calculate the residue phase and amplitude terms */
if (eta==0)
   {
	   deltalm=0.0;
	   rholm=0.0;
   }
	//This is the point test limit 
	//
	//

else
 {	//dealing with "non-test-particle"
  switch( l )
  {
    case 2:
      switch( abs(m) )
      {
        case 2:
		  // CHECKME : Changing to Pade resummed coefficients 

          // deltalm = vh3*(hCoeffs->delta22vh3 + vh3*(hCoeffs->delta22vh6
          //  + vh*vh*(hCoeffs->delta22vh9*vh)))
          //  + hCoeffs->delta22v5 *v*v2*v2 + hCoeffs->delta22v8 *v2*v2*v2*v2;
          delta22_leading=(7./3.)*v3;
		  delta22_Num=(808920*eta*pi*v + 137388.*pi2*v2 + 35.*eta2*(136080. + (154975. - 1359276.*eta)*v2));
		  delta22_Den=(808920.*eta*pi*v + 137388.*pi2*v2 + 35.*eta2*(136080. + (154975. + 40404.*eta)*v2));
		  deltalm= delta22_leading *(delta22_Den/delta22_Den);
         // rholm  = 1. + v2*(hCoeffs->rho22v2 + v*(hCoeffs->rho22v3
         //   + v*(hCoeffs->rho22v4
         //   + v*(hCoeffs->rho22v5 + v*(hCoeffs->rho22v6
         //   + hCoeffs->rho22v6l*eulerlogxabs + v*(hCoeffs->rho22v7
         //   + v*(hCoeffs->rho22v8 + hCoeffs->rho22v8l*eulerlogxabs
         //   + (hCoeffs->rho22v10 + hCoeffs->rho22v10l * eulerlogxabs)*v2)))))));
         
 	rholm=1.+(v2*hCoeffs->rho22v2)+(hCoeffs->rho22v4*v4)+(hCoeffs->rho22v6*v2*v2*v2 + hCoeffs->rho22v6l*eulerlogxabs*v2*v2*v2)+(hCoeffs->rho22v8*v4*v4 + hCoeffs->rho22v8l*eulerlogxabs*v4*v4)+(hCoeffs->rho22v10*v4*v4*v2 + hCoeffs->rho22v10l * eulerlogxabs*v4*v4*v2);
          break;
        case 1:
         // deltalm = vh3*(hCoeffs->delta21vh3 + vh3*(hCoeffs->delta21vh6
         //   + vh*(hCoeffs->delta21vh7 + (hCoeffs->delta21vh9)*vh*vh)))
         //   + hCoeffs->delta21v5*v*v2*v2 + hCoeffs->delta21v7*v2*v2*v2*v;
		  delta21_leading=(2./3.)*v3;
		  delta21_Num=69020.*eta + 5992.*pi*v;
		  delta21_Den=5992.*pi*v + 2456.*eta*(28.+493.*eta*v2);
		  deltalm=delta21_leading*(delta21_Num/delta21_Den);

          rholm  = 1. + v*(hCoeffs->rho21v1
            + v*( hCoeffs->rho21v2 + v*(hCoeffs->rho21v3 + v*(hCoeffs->rho21v4
            + v*(hCoeffs->rho21v5 + v*(hCoeffs->rho21v6 + hCoeffs->rho21v6l*eulerlogxabs
            + v*(hCoeffs->rho21v7 + hCoeffs->rho21v7l * eulerlogxabs
            + v*(hCoeffs->rho21v8 + hCoeffs->rho21v8l * eulerlogxabs
            + (hCoeffs->rho21v10 + hCoeffs->rho21v10l * eulerlogxabs)*v2))))))));
          break;
        default:
          XLAL_ERROR( XLAL_EINVAL );
          break;
      }
      break;
    case 3:
      switch (m)
      {
        case 3:
          //deltalm = vh3*(hCoeffs->delta33vh3 + vh3*(hCoeffs->delta33vh6 + hCoeffs->delta3//3vh9*vh3))
		  delta33_leading=(13./10.)*v3;
		  delta33_Num = 1. + 94770.*pi*v/(566279.*eta);
		  delta33_Den = delta33_Num + (80897.*eta*v2/3159.);
		  deltalm=delta33_leading*(delta33_Num/delta33_Den);

            //+ hCoeffs->delta33v5*v*v2*v2 + hCoeffs->delta33v7*v2*v2*v2*v;
          rholm  = 1. + v2*(hCoeffs->rho33v2 + v*(hCoeffs->rho33v3 + v*(hCoeffs->rho33v4
            + v*(hCoeffs->rho33v5 + v*(hCoeffs->rho33v6 + hCoeffs->rho33v6l*eulerlogxabs
            + v*(hCoeffs->rho33v7 + (hCoeffs->rho33v8 + hCoeffs->rho33v8l*eulerlogxabs)*v))))));
          rholm=rholm+ ((hCoeffs->rho33v10 + hCoeffs->rho33v10l*eulerlogxabs)*v5*v5);
          break;
        case 2:
         // deltalm = vh3*(hCoeffs->delta32vh3 + vh*(hCoeffs->delta32vh4 + vh*vh*(hCoeffs->delta32vh6
           // + hCoeffs->delta32vh9*vh3)));
		   
		  deltalm = (v3*hCoeffs->delta32vh3)+(v3*v3*hCoeffs->delta32vh6);
          rholm  = 1. + v*(hCoeffs->rho32v
            + v*(hCoeffs->rho32v2 + v*(hCoeffs->rho32v3 + v*(hCoeffs->rho32v4 + v*(hCoeffs->rho32v5
            + v*(hCoeffs->rho32v6 + hCoeffs->rho32v6l*eulerlogxabs
            + (hCoeffs->rho32v8 + hCoeffs->rho32v8l*eulerlogxabs)*v2))))));
          break;
        case 1:
        //  deltalm = vh3*(hCoeffs->delta31vh3 + vh3*(hCoeffs->delta31vh6
        //    + vh*(hCoeffs->delta31vh7 + hCoeffs->delta31vh9*vh*vh)))
       //     + hCoeffs->delta31v5*v*v2*v2;
	      delta31_leading = (13./30.)*v3;
		  delta31_Num = 4641.*eta + 1690.*pi*v;
		  delta31_Den = delta31_Num + 18207.*eta2*v2;
		  deltalm=delta31_leading*(delta31_Num/delta31_Den);

          rholm  = 1. + v2*(hCoeffs->rho31v2 + v*(hCoeffs->rho31v3 + v*(hCoeffs->rho31v4
            + v*(hCoeffs->rho31v5 + v*(hCoeffs->rho31v6 + hCoeffs->rho31v6l*eulerlogxabs
            + v*(hCoeffs->rho31v7 + (hCoeffs->rho31v8 + hCoeffs->rho31v8l*eulerlogxabs)*v))))));
          rholm=rholm+ ((hCoeffs->rho31v10 + hCoeffs->rho31v10l*eulerlogxabs)*v5*v5);
          break;
        default:
          XLAL_ERROR( XLAL_EINVAL );
          break;
      }
      break;
    case 4:
      switch (m)
      {
        case 4:
          deltalm = v3*(hCoeffs->delta44vh3 + hCoeffs->delta44vh6 *v3);
		 //deltalm = vh3*(hCoeffs->delta44vh3 + hCoeffs->delta44vh6 *vh3)
		 //1034             + hCoeffs->delta44v5*v2*v2*v;
		 //
          rholm  = 1. + v2*(hCoeffs->rho44v2
            + v*( hCoeffs->rho44v3 + v*(hCoeffs->rho44v4
            + v*(hCoeffs->rho44v5 + (hCoeffs->rho44v6
            + hCoeffs->rho44v6l*eulerlogxabs)*v))));
          rholm=rholm + ((hCoeffs->rho44v8 + hCoeffs->rho44v8l*eulerlogxabs)*v4*v4);
          break;
        case 3:
          //deltalm = v3*(hCoeffs->delta43vh3 + vh*(hCoeffs->delta43vh4
          //  + hCoeffs->delta43vh6*vh*vh));
		  deltalm = (v3*hCoeffs->delta43vh3) + (v3*v3*hCoeffs->delta43vh6);

          rholm  = 1. + v*(hCoeffs->rho43v
            + v*(hCoeffs->rho43v2
            + v2*(hCoeffs->rho43v4 + v*(hCoeffs->rho43v5
            + (hCoeffs->rho43v6 + hCoeffs->rho43v6l*eulerlogxabs)*v))));
          rholm=rholm + ((hCoeffs->rho43v8 + hCoeffs->rho43v8l*eulerlogxabs)*v4*v4);
          break;
        case 2:
          //deltalm = vh3*(hCoeffs->delta42vh3 + hCoeffs->delta42vh6*vh3);
          deltalm = (v3*hCoeffs->delta42vh3) + (v3*v3*hCoeffs->delta42vh6);

		  rholm  = 1. + v2*(hCoeffs->rho42v2
            + v*(hCoeffs->rho42v3 + v*(hCoeffs->rho42v4 + v*(hCoeffs->rho42v5
            + (hCoeffs->rho42v6 + hCoeffs->rho42v6l*eulerlogxabs)*v))));
            rholm=rholm + ((hCoeffs->rho42v8 + hCoeffs->rho42v8l*eulerlogxabs)*v4*v4);
          break;
        case 1:
       //   deltalm = vh3*(hCoeffs->delta41vh3 + vh*(hCoeffs->delta41vh4
       //     + hCoeffs->delta41vh6*vh*vh));
	      deltalm = (v3*hCoeffs->delta41vh3)+(v3*v3*+ hCoeffs->delta41vh6);

          rholm  = 1. + v*(hCoeffs->rho41v
            + v*(hCoeffs->rho41v2
            + v2*(hCoeffs->rho41v4 + v*(hCoeffs->rho41v5
            + (hCoeffs->rho41v6 +  hCoeffs->rho41v6l*eulerlogxabs)*v))));
          rholm=rholm + ((hCoeffs->rho41v8 + hCoeffs->rho41v8l*eulerlogxabs)*v4*v4);
          break;
        default:
          XLAL_ERROR( XLAL_EINVAL );
          break;
      }
      break;
    case 5:
      switch (m)
      {
        case 5:
          deltalm = hCoeffs->delta55vh3*v3;
          rholm  = 1. + v2*( hCoeffs->rho55v2
            + v*(hCoeffs->rho55v3 + v*(hCoeffs->rho55v4
            + v*(hCoeffs->rho55v5 + hCoeffs->rho55v6*v))));
          rholm = rholm + (hCoeffs->rho55v6l*eulerlogxabs*v3*v3)+ ((hCoeffs->rho55v8 + hCoeffs->rho55v8l*eulerlogxabs)*v4*v4);
          break;
        case 4:
          deltalm =0.;
			 // vh3*(hCoeffs->delta54vh3 + hCoeffs->delta54vh4*vh);
          rholm  = 1. + v2*(hCoeffs->rho54v2 + v*(hCoeffs->rho54v3
            + hCoeffs->rho54v4*v));
          rholm=rholm + ((hCoeffs->rho54v6 + hCoeffs->rho54v6l*eulerlogxabs)*v3*v3);
          break;
        case 3:
          deltalm =0.;
			 // hCoeffs->delta53vh3 * vh3;
          rholm  = 1. + v2*(hCoeffs->rho53v2
            + v*(hCoeffs->rho53v3 + v*(hCoeffs->rho53v4 + hCoeffs->rho53v5*v)));
          rholm=rholm +  ((hCoeffs->rho53v6 + hCoeffs->rho53v6l*eulerlogxabs)*v3*v3) + ((hCoeffs->rho53v8 + hCoeffs->rho53v8l*eulerlogxabs)*v4*v4);
          break;
        case 2:
          deltalm = 0.;
			  //vh3*(hCoeffs->delta52vh3 + hCoeffs->delta52vh4*vh);
          rholm  = 1. + v2*(hCoeffs->rho52v2 + v*(hCoeffs->rho52v3
            + hCoeffs->rho52v4*v));
          rholm=rholm + ((hCoeffs->rho52v6 + hCoeffs->rho52v6l*eulerlogxabs)*v3*v3);
          break;
        case 1:
          deltalm = 0.;
			  //hCoeffs->delta51vh3 * vh3;
          rholm  = 1. + v2*(hCoeffs->rho51v2
            + v*(hCoeffs->rho51v3 + v*(hCoeffs->rho51v4 + hCoeffs->rho51v5*v)));
           rholm = rholm + ((hCoeffs->rho51v6 +hCoeffs->rho51v6l*eulerlogxabs)*v3*v3)+ ((hCoeffs->rho51v8 + hCoeffs->rho51v8l*eulerlogxabs)*v4*v4);
          break;
        default:
          XLAL_ERROR( XLAL_EINVAL );
          break;
      }
      break;
    case 6:
      switch (m)
      {
        case 6:
          deltalm = 0.;
			  //hCoeffs->delta66vh3*vh3;
          rholm  = 1. + v2*(hCoeffs->rho66v2 + v*(hCoeffs->rho66v3
            + hCoeffs->rho66v4*v));
          rholm=rholm + ((hCoeffs->rho66v6 + hCoeffs->rho66v6l*eulerlogxabs)*v3*v3);
          break;
        case 5:
          deltalm = 0.;
		 // hCoeffs->delta65vh3*vh3;
          rholm  = 1. + v2*(hCoeffs->rho65v2 + hCoeffs->rho65v3*v);

          rholm=rholm + ((hCoeffs->rho65v6 + hCoeffs->rho65v6l*eulerlogxabs)*v3*v3) + (hCoeffs->rho65v4 *v2 *v2);
          break;
        case 4:
          deltalm = 0.;
			 // hCoeffs->delta64vh3 * vh3;
          rholm  = 1. + v2*(hCoeffs->rho64v2 + v*(hCoeffs->rho64v3
            + hCoeffs->rho64v4*v));
          rholm=rholm + ((hCoeffs->rho64v6 + hCoeffs->rho64v6l*eulerlogxabs)*v3*v3);
          break;
        case 3:
          deltalm = 0.;
			  //hCoeffs->delta63vh3 * vh3;
          rholm  = 1. + v2*(hCoeffs->rho63v2 + hCoeffs->rho63v3*v);
          rholm=rholm + ((hCoeffs->rho63v6 + hCoeffs->rho63v6l*eulerlogxabs)*v3*v3) + (hCoeffs->rho63v4 *v2 *v2);
          break;
        case 2:
          deltalm = 0.;
			 // hCoeffs->delta62vh3 * vh3;
          rholm  = 1. + v2*(hCoeffs->rho62v2 + v*(hCoeffs->rho62v3
            + hCoeffs->rho62v4 * v));
          rholm=rholm + ((hCoeffs->rho62v6 + hCoeffs->rho62v6l*eulerlogxabs)*v3*v3);
          break;
        case 1:
          deltalm = 0.;
		 // hCoeffs->delta61vh3 * vh3;
          rholm  = 1. + v2*(hCoeffs->rho61v2 + hCoeffs->rho61v3*v);
          rholm=rholm + ((hCoeffs->rho61v6 + hCoeffs->rho61v6l*eulerlogxabs)*v3*v3) + (hCoeffs->rho61v4 *v2 *v2);
          break;
        default:
          XLAL_ERROR( XLAL_EINVAL );
          break;
      }
      break;
    case 7:
      switch (m)
      {
        case 7:
          deltalm = 0.;
			 // hCoeffs->delta77vh3 * vh3;
          rholm   = 1. + v2*(hCoeffs->rho77v2 + hCoeffs->rho77v3 * v);
          rholm=rholm + ((hCoeffs->rho77v6 + hCoeffs->rho77v6l*eulerlogxabs)*v3*v3) + (hCoeffs->rho77v4 *v2 *v2);
          break;
        case 6:
          deltalm = 0.0;
          rholm   = 1. + hCoeffs->rho76v2 * v2;
          rholm=rholm + (hCoeffs->rho76v4 *v2 *v2);
          break;
        case 5:
          deltalm = 0.0;
		  //hCoeffs->delta75vh3 * vh3;
          rholm   = 1. + v2*(hCoeffs->rho75v2 + hCoeffs->rho75v3*v);
          rholm=rholm + ((hCoeffs->rho75v6 + hCoeffs->rho75v6l*eulerlogxabs)*v3*v3) + (hCoeffs->rho75v4 *v2 *v2);
          break;
        case 4:
          deltalm = 0.0;
          rholm   = 1. + hCoeffs->rho74v2 * v2;
          rholm=rholm + (hCoeffs->rho74v4 *v2 *v2);
          break;
        case 3:
          deltalm = 0.;
		 // hCoeffs->delta73vh3 *vh3;
          rholm   = 1. + v2*(hCoeffs->rho73v2 + hCoeffs->rho73v3 * v);
          rholm=rholm + ((hCoeffs->rho73v6 + hCoeffs->rho73v6l*eulerlogxabs)*v3*v3) + (hCoeffs->rho73v4 *v2 *v2); 
          break;
        case 2:
          deltalm = 0.0;
          rholm   = 1. + hCoeffs->rho72v2 * v2;
          rholm=rholm + (hCoeffs->rho72v4 *v2 *v2);
          break;
        case 1:
          deltalm = 0.0;
		 // hCoeffs->delta71vh3 * vh3;
          rholm   = 1. + v2*(hCoeffs->rho71v2 +hCoeffs->rho71v3 * v);
          rholm=rholm + ((hCoeffs->rho71v6 + hCoeffs->rho71v6l*eulerlogxabs)*v3*v3) + (hCoeffs->rho71v4 *v2 *v2);
          break;
        default:
          XLAL_ERROR( XLAL_EINVAL );
          break;
      }
      break;
    case 8:
      switch (m)
      {
        case 8:
          deltalm = 0.0;
          rholm   = 1. + hCoeffs->rho88v2 * v2;
          rholm=rholm + (hCoeffs->rho88v4 *v2 *v2);
          break;
        case 7:
          deltalm = 0.0;
          rholm   = 1. + hCoeffs->rho87v2 * v2;
          rholm=rholm + (hCoeffs->rho87v4 *v2 *v2);
          break;
        case 6:
          deltalm = 0.0;
          rholm   = 1. + hCoeffs->rho86v2 * v2;
          rholm=rholm + (hCoeffs->rho86v4 *v2 *v2);
          break;
        case 5:
          deltalm = 0.0;
          rholm   = 1. + hCoeffs->rho85v2 * v2;
          rholm=rholm + (hCoeffs->rho85v4 *v2 *v2);
          break;
        case 4:
          deltalm = 0.0;
          rholm  = 1. + hCoeffs->rho84v2 * v2;
          rholm=rholm + (hCoeffs->rho84v4 *v2 *v2);
          break;
        case 3:
          deltalm = 0.0;
          rholm  = 1. + hCoeffs->rho83v2 * v2;
          rholm=rholm + (hCoeffs->rho83v4 *v2 *v2);
          break;
        case 2:
          deltalm = 0.0;
          rholm  = 1. + hCoeffs->rho82v2 * v2;
          rholm=rholm + (hCoeffs->rho82v4 *v2 *v2);
          break;
        case 1:
          deltalm = 0.0;
          rholm  = 1. + hCoeffs->rho81v2 * v2;
          rholm=rholm + (hCoeffs->rho81v4 *v2 *v2);
          break;
        default:
          XLAL_ERROR( XLAL_EINVAL );
          break;
      }
  
      break;
    default:
      XLAL_ERROR( XLAL_EINVAL );
      break;
  }}
  /* Raise rholm to the lth pow */
  rholmPwrl = 1.0;
  i = l;
  while ( i-- )
  {
    rholmPwrl *= rholm;
  }
//CHANGE ME!!!!!!!!!!!!!!!!!!!!!!!!WILL INSPRIAL VERY FAST OTHERWISE
  *hlm = Tlm * cexp(I * deltalm) * Slm * rholmPwrl;
  *hlm *= hNewton;
  
// *hlm = Tlm * cexp(I * deltalm) * Slm * rholmPwrl;
// *hlm = 0;

//  *hlm = hNewton;

//Add tidal corrections to l=2 and m=2 modes
if  ( l == 2 )
   {
     if( m ==2)
          {
            *hlm  = *hlm+ hlm_tidal;
          }
    }
//printf("done with making hlm......");
  return XLAL_SUCCESS;
} 


UNUSED static int  XLALSimIMRTNSEOBGetFactorizedWaveform_withoutNewtonianPart( 
                                COMPLEX16   * restrict hlm,    /**<< The value of hlm (populated by the function) */
                                REAL8Vector * restrict values, /**<< Vector containing dynamics r, phi, pr, pphi for a given point */
                                const REAL8 v,                 /**<< Velocity (in geometric units) */
                                const INT4  l,                 /**<< Mode l */
                                const INT4  m,                 /**<< Mode m */
                                TNSEOBParams   * restrict params 
                                )
{

//Notice :  I am feeding in v=v_phi (the non keplerian velocity for all the parts of h_lm unlike in EOBNR)


  /* Status of function calls */
  INT4 status;
  INT4 i;

  REAL8 eta, mass1, mass2, lambda1, lambda2, pi,pi2,eta2;
  REAL8 r, pr, pp, Omega, v2, vh, vh3,v3, k, hathatk, eulerlogxabs,v4,v5;
  REAL8 Hreal, Heff, Slm, deltalm, rholm, rholmPwrl;
  REAL8 delta22_leading,delta21_leading,delta33_leading,delta31_leading,delta22_Num,delta21_Num,delta33_Num,delta31_Num,delta22_Den, delta21_Den, delta33_Den,delta31_Den;
  COMPLEX16 Tlm;
  COMPLEX16 hlm_tidal;
  gsl_sf_result lnr1, arg1, z2;

  /* Non-Keplerian velocity */
  REAL8 psi,rOmega;
  /* Pre-computed coefficients */
  TNSFacWaveformCoeffs *hCoeffs = params->hCoeffs;

  if ( abs(m) > (INT4) l )
  {
    XLAL_ERROR( XLAL_EINVAL );
  }


  eta = params->eta;
  mass1 = params->mass1;
  mass2 = params->mass2;
  lambda1 = params -> lambda1;
  lambda2 = params -> lambda2;

  /* Check our eta was sensible */
  if ( eta > 0.25 )
  {
    XLALPrintError("Eta seems to be > 0.25 - this isn't allowed!\n" );
    XLAL_ERROR( XLAL_EINVAL );
  }
  else if ( eta == 0.25 && m % 2 )
  {
    /* If m is odd and dM = 0, hLM will be zero */
    memset( hlm, 0, sizeof( COMPLEX16 ) );
    return XLAL_SUCCESS;
  }

  r  = values->data[0];
  pr = values->data[2];
  pp = values->data[3];

  //printf("accesing the data r,pr,phi");
//  printf("\n%f",r);
  Heff  = XLALEffectiveHamiltonian( eta, mass1, mass2, lambda1,lambda2, r, pr, pp);
  Hreal = sqrt( 1.0 + 2.0 * eta * ( Heff - 1.0) );
  v2    = v * v;
  v3  = v*v2;
  v4 = v2*v2;
  v5 = v3*v2;
  eulerlogxabs = LAL_GAMMA + log( 2.0 * (REAL8)m * v );
  pi=LAL_PI;
  pi2=pi*pi;
  eta2=eta*eta;

  /* Calculate the non-Keplerian velocity */
  /* given by Eq. (18) of Pan et al, PRD84, 124052(2011) */
  /* psi given by Eq. (19) of Pan et al, PRD84, 124052(2011) */
  /* Assign temporarily to vPhi */
  //nonKeplerianCoefficient gives psi
  psi = nonKeplerianCoefficient( values, eta,mass1,mass2,lambda1,lambda2 );
  /* Assign rOmega value temporarily to vPhi */
  rOmega  = r * cbrt(psi);
 Omega = v/rOmega;
   vh3   = Hreal * Omega;
     vh    = cbrt(vh3);



//REAL8 xNewt=cbrt(Omega*Omega);

printf("Omega in Newtonian Flux : %g\n", (double) Omega);
  /* Calculate the newtonian multipole */
  
  
  //printf("I calculated newtonian fl");

 /*Calculate the tidal part of h_lm*/

  status = XLALSimIMRTNSEOBCalculateTidal_hlm( &hlm_tidal, v , v/Omega,
            values->data[1], (UINT4)l, m, params, values );
  if ( status == XLAL_FAILURE )
  {
    XLAL_ERROR( XLAL_EFUNC );
  }

  /* Calculate the source term */
  if ( ( (l+m)%2 ) == 0)
  {
    Slm = Heff;
  }
  else
  {

  //In gets changed in TNSEOB! 
  // We have  Slm = v * pp in EOBNR but in TNSEOB. 
  Slm = pp/(rOmega*v);
  }
    printf("Lets check why Slm is slightly off for odd l+m");
    printf("pp : %g\n", (double) pp);
  printf("rOmega : %g\n", (double) rOmega);
  printf("v : %g\n", (double) v);

  /* Calculate the Tail term */
  k  = m * Omega;
  hathatk = Hreal * k;
  XLAL_CALLGSL( status = gsl_sf_lngamma_complex_e( l+1.0, -2.0*hathatk, &lnr1, &arg1 ) );
  if (status != GSL_SUCCESS)
  {
    XLALPrintError("Error in GSL function\n" );
    XLAL_ERROR( XLAL_EFUNC );
  }
  XLAL_CALLGSL( status = gsl_sf_fact_e( l, &z2 ) );
  if ( status != GSL_SUCCESS)
  {
    XLALPrintError("Error in GSL function\n" );
    XLAL_ERROR( XLAL_EFUNC );
  }
  Tlm = cexp( ( lnr1.val + LAL_PI * hathatk ) + I * ( 2.0 * hathatk * log(4.0*k/sqrt(LAL_E)) ) );

  
 // Tlm = cexp( ( lnr1.val + LAL_PI * hathatk ) + I * (
   //     arg1.val + 2.0 * hathatk * log(4.0*k/sqrt(LAL_E)) ) );
  Tlm /= z2.val;

REAL8 omega_check=v/rOmega;

  printf("Construction of Tlm... hathatk : %g\n", (double) hathatk);
  printf("Omega : %g\n", (double) Omega);
  printf("omega_check : %g\n", (double) omega_check);
printf("Hreal  : %g\n", (double) Hreal );
   printf("lnr1.val  : %g\n", (double) lnr1.val );
   printf("arg1.val  : %g\n", (double) arg1.val );
     printf(" z2.val  : %g\n", (double)  z2.val );




REAL8 modTlm=sqrt(creal(Tlm)*creal(Tlm) + cimag(Tlm)*cimag(Tlm));


  /* Calculate the residue phase and amplitude terms */
if (eta==0)
   {
     deltalm=0.0;
     rholm=0.0;
   }
  //This is the point test limit 
  //
  //

else
 {  //dealing with "non-test-particle"
  switch( l )
  {
    case 2:
      switch( abs(m) )
      {
        case 2:
      // CHECKME : Changing to Pade resummed coefficients 

          // deltalm = vh3*(hCoeffs->delta22vh3 + vh3*(hCoeffs->delta22vh6
          //  + vh*vh*(hCoeffs->delta22vh9*vh)))
          //  + hCoeffs->delta22v5 *v*v2*v2 + hCoeffs->delta22v8 *v2*v2*v2*v2;
          delta22_leading=(7./3.)*v3;
      delta22_Num=(808920*eta*pi*v + 137388.*pi2*v2 + 35.*eta2*(136080. + (154975. - 1359276.*eta)*v2));
      delta22_Den=(808920.*eta*pi*v + 137388.*pi2*v2 + 35.*eta2*(136080. + (154975. + 40404.*eta)*v2));
      deltalm= delta22_leading *(delta22_Den/delta22_Den);
         // rholm  = 1. + v2*(hCoeffs->rho22v2 + v*(hCoeffs->rho22v3
         //   + v*(hCoeffs->rho22v4
         //   + v*(hCoeffs->rho22v5 + v*(hCoeffs->rho22v6
         //   + hCoeffs->rho22v6l*eulerlogxabs + v*(hCoeffs->rho22v7
         //   + v*(hCoeffs->rho22v8 + hCoeffs->rho22v8l*eulerlogxabs
         //   + (hCoeffs->rho22v10 + hCoeffs->rho22v10l * eulerlogxabs)*v2)))))));
         
  rholm=1.+(v2*hCoeffs->rho22v2)+(hCoeffs->rho22v4*v4)+(hCoeffs->rho22v6*v2*v2*v2 + hCoeffs->rho22v6l*eulerlogxabs*v2*v2*v2)+(hCoeffs->rho22v8*v4*v4 + hCoeffs->rho22v8l*eulerlogxabs*v4*v4)+(hCoeffs->rho22v10*v4*v4*v2 + hCoeffs->rho22v10l * eulerlogxabs*v4*v4*v2);
         printf(" rholm=1.+(v2*hCoeffs->rho22v2)+everything else; : %g\n", (double) rholm);
     printf("hCoeffs->rho22v2 : %g\n", (double) hCoeffs->rho22v2);
     printf("v2 : %g\n", (double) v2);
     printf("eta used inside the Waveform generation : %g\n", (double) eta);
     break;
        case 1:
         // deltalm = vh3*(hCoeffs->delta21vh3 + vh3*(hCoeffs->delta21vh6
         //   + vh*(hCoeffs->delta21vh7 + (hCoeffs->delta21vh9)*vh*vh)))
         //   + hCoeffs->delta21v5*v*v2*v2 + hCoeffs->delta21v7*v2*v2*v2*v;
      delta21_leading=(2./3.)*v3;
      delta21_Num=69020.*eta + 5992.*pi*v;
      delta21_Den=5992.*pi*v + 2456.*eta*(28.+493.*eta*v2);
      deltalm=delta21_leading*(delta21_Num/delta21_Den);

          rholm  = 1. + v*(hCoeffs->rho21v1
            + v*( hCoeffs->rho21v2 + v*(hCoeffs->rho21v3 + v*(hCoeffs->rho21v4
            + v*(hCoeffs->rho21v5 + v*(hCoeffs->rho21v6 + hCoeffs->rho21v6l*eulerlogxabs
            + v*(hCoeffs->rho21v7 + hCoeffs->rho21v7l * eulerlogxabs
            + v*(hCoeffs->rho21v8 + hCoeffs->rho21v8l * eulerlogxabs
            + (hCoeffs->rho21v10 + hCoeffs->rho21v10l * eulerlogxabs)*v2))))))));
          break;
        default:
          XLAL_ERROR( XLAL_EINVAL );
          break;
      }
      break;
    case 3:
      switch (m)
      {
        case 3:
          //deltalm = vh3*(hCoeffs->delta33vh3 + vh3*(hCoeffs->delta33vh6 + hCoeffs->delta3//3vh9*vh3))
      delta33_leading=(13./10.)*v3;
      delta33_Num = 1. + 94770.*pi*v/(566279.*eta);
      delta33_Den = delta33_Num + (80897.*eta*v2/3159.);
      deltalm=delta33_leading*(delta33_Num/delta33_Den);

            //+ hCoeffs->delta33v5*v*v2*v2 + hCoeffs->delta33v7*v2*v2*v2*v;
          rholm  = 1. + v2*(hCoeffs->rho33v2 + v*(hCoeffs->rho33v3 + v*(hCoeffs->rho33v4
            + v*(hCoeffs->rho33v5 + v*(hCoeffs->rho33v6 + hCoeffs->rho33v6l*eulerlogxabs
            + v*(hCoeffs->rho33v7 + (hCoeffs->rho33v8 + hCoeffs->rho33v8l*eulerlogxabs)*v))))));
          rholm=rholm+ ((hCoeffs->rho33v10 + hCoeffs->rho33v10l*eulerlogxabs)*v5*v5);
          break;
        case 2:
         // deltalm = vh3*(hCoeffs->delta32vh3 + vh*(hCoeffs->delta32vh4 + vh*vh*(hCoeffs->delta32vh6
           // + hCoeffs->delta32vh9*vh3)));
       
      deltalm = (v3*hCoeffs->delta32vh3)+(v3*v3*hCoeffs->delta32vh6);
          rholm  = 1. + v*(hCoeffs->rho32v
            + v*(hCoeffs->rho32v2 + v*(hCoeffs->rho32v3 + v*(hCoeffs->rho32v4 + v*(hCoeffs->rho32v5
            + v*(hCoeffs->rho32v6 + hCoeffs->rho32v6l*eulerlogxabs
            + (hCoeffs->rho32v8 + hCoeffs->rho32v8l*eulerlogxabs)*v2))))));
          break;
        case 1:
        //  deltalm = vh3*(hCoeffs->delta31vh3 + vh3*(hCoeffs->delta31vh6
        //    + vh*(hCoeffs->delta31vh7 + hCoeffs->delta31vh9*vh*vh)))
       //     + hCoeffs->delta31v5*v*v2*v2;
        delta31_leading = (13./30.)*v3;
      delta31_Num = 4641.*eta + 1690.*pi*v;
      delta31_Den = delta31_Num + 18207.*eta2*v2;
      deltalm=delta31_leading*(delta31_Num/delta31_Den);

          rholm  = 1. + v2*(hCoeffs->rho31v2 + v*(hCoeffs->rho31v3 + v*(hCoeffs->rho31v4
            + v*(hCoeffs->rho31v5 + v*(hCoeffs->rho31v6 + hCoeffs->rho31v6l*eulerlogxabs
            + v*(hCoeffs->rho31v7 + (hCoeffs->rho31v8 + hCoeffs->rho31v8l*eulerlogxabs)*v))))));
          rholm=rholm+ ((hCoeffs->rho31v10 + hCoeffs->rho31v10l*eulerlogxabs)*v5*v5);
          break;
        default:
          XLAL_ERROR( XLAL_EINVAL );
          break;
      }
      break;
    case 4:
      switch (m)
      {
        case 4:
          deltalm = v3*(hCoeffs->delta44vh3 + hCoeffs->delta44vh6 *v3);
     //deltalm = vh3*(hCoeffs->delta44vh3 + hCoeffs->delta44vh6 *vh3)
     //1034             + hCoeffs->delta44v5*v2*v2*v;
     //
          rholm  = 1. + v2*(hCoeffs->rho44v2
            + v*( hCoeffs->rho44v3 + v*(hCoeffs->rho44v4
            + v*(hCoeffs->rho44v5 + (hCoeffs->rho44v6
            + hCoeffs->rho44v6l*eulerlogxabs)*v))));
          rholm=rholm + ((hCoeffs->rho44v8 + hCoeffs->rho44v8l*eulerlogxabs)*v4*v4);
          break;
        case 3:
          //deltalm = v3*(hCoeffs->delta43vh3 + vh*(hCoeffs->delta43vh4
          //  + hCoeffs->delta43vh6*vh*vh));
      deltalm = (v3*hCoeffs->delta43vh3) + (v3*v3*hCoeffs->delta43vh6);

          rholm  = 1. + v*(hCoeffs->rho43v
            + v*(hCoeffs->rho43v2
            + v2*(hCoeffs->rho43v4 + v*(hCoeffs->rho43v5
            + (hCoeffs->rho43v6 + hCoeffs->rho43v6l*eulerlogxabs)*v))));
          rholm=rholm + ((hCoeffs->rho43v8 + hCoeffs->rho43v8l*eulerlogxabs)*v4*v4);
          break;
        case 2:
          //deltalm = vh3*(hCoeffs->delta42vh3 + hCoeffs->delta42vh6*vh3);
          deltalm = (v3*hCoeffs->delta42vh3) + (v3*v3*hCoeffs->delta42vh6);

      rholm  = 1. + v2*(hCoeffs->rho42v2
            + v*(hCoeffs->rho42v3 + v*(hCoeffs->rho42v4 + v*(hCoeffs->rho42v5
            + (hCoeffs->rho42v6 + hCoeffs->rho42v6l*eulerlogxabs)*v))));
            rholm=rholm + ((hCoeffs->rho42v8 + hCoeffs->rho42v8l*eulerlogxabs)*v4*v4);
          break;
        case 1:
       //   deltalm = vh3*(hCoeffs->delta41vh3 + vh*(hCoeffs->delta41vh4
       //     + hCoeffs->delta41vh6*vh*vh));
        deltalm = (v3*hCoeffs->delta41vh3)+(v3*v3*+ hCoeffs->delta41vh6);

          rholm  = 1. + v*(hCoeffs->rho41v
            + v*(hCoeffs->rho41v2
            + v2*(hCoeffs->rho41v4 + v*(hCoeffs->rho41v5
            + (hCoeffs->rho41v6 +  hCoeffs->rho41v6l*eulerlogxabs)*v))));
          rholm=rholm + ((hCoeffs->rho41v8 + hCoeffs->rho41v8l*eulerlogxabs)*v4*v4);
          break;
        default:
          XLAL_ERROR( XLAL_EINVAL );
          break;
      }
      break;
    case 5:
      switch (m)
      {
        case 5:
          deltalm = hCoeffs->delta55vh3*v3;
          rholm  = 1. + v2*( hCoeffs->rho55v2
            + v*(hCoeffs->rho55v3 + v*(hCoeffs->rho55v4
            + v*(hCoeffs->rho55v5 + hCoeffs->rho55v6*v))));
          rholm = rholm + (hCoeffs->rho55v6l*eulerlogxabs*v3*v3)+ ((hCoeffs->rho55v8 + hCoeffs->rho55v8l*eulerlogxabs)*v4*v4);
          break;
        case 4:
          deltalm =0.;
       // vh3*(hCoeffs->delta54vh3 + hCoeffs->delta54vh4*vh);
          rholm  = 1. + v2*(hCoeffs->rho54v2 + v*(hCoeffs->rho54v3
            + hCoeffs->rho54v4*v));
          rholm=rholm + ((hCoeffs->rho54v6 + hCoeffs->rho54v6l*eulerlogxabs)*v3*v3);
          break;
        case 3:
          deltalm =0.;
       // hCoeffs->delta53vh3 * vh3;
          rholm  = 1. + v2*(hCoeffs->rho53v2
            + v*(hCoeffs->rho53v3 + v*(hCoeffs->rho53v4 + hCoeffs->rho53v5*v)));
          rholm=rholm +  ((hCoeffs->rho53v6 + hCoeffs->rho53v6l*eulerlogxabs)*v3*v3) + ((hCoeffs->rho53v8 + hCoeffs->rho53v8l*eulerlogxabs)*v4*v4);
          break;
        case 2:
          deltalm = 0.;
        //vh3*(hCoeffs->delta52vh3 + hCoeffs->delta52vh4*vh);
          rholm  = 1. + v2*(hCoeffs->rho52v2 + v*(hCoeffs->rho52v3
            + hCoeffs->rho52v4*v));
          rholm=rholm + ((hCoeffs->rho52v6 + hCoeffs->rho52v6l*eulerlogxabs)*v3*v3);
          break;
        case 1:
          deltalm = 0.;
        //hCoeffs->delta51vh3 * vh3;
          rholm  = 1. + v2*(hCoeffs->rho51v2
            + v*(hCoeffs->rho51v3 + v*(hCoeffs->rho51v4 + hCoeffs->rho51v5*v)));
           rholm = rholm + ((hCoeffs->rho51v6 +hCoeffs->rho51v6l*eulerlogxabs)*v3*v3)+ ((hCoeffs->rho51v8 + hCoeffs->rho51v8l*eulerlogxabs)*v4*v4);
          break;
        default:
          XLAL_ERROR( XLAL_EINVAL );
          break;
      }
      break;
    case 6:
      switch (m)
      {
        case 6:
          deltalm = 0.;
        //hCoeffs->delta66vh3*vh3;
          rholm  = 1. + v2*(hCoeffs->rho66v2 + v*(hCoeffs->rho66v3
            + hCoeffs->rho66v4*v));
          rholm=rholm + ((hCoeffs->rho66v6 + hCoeffs->rho66v6l*eulerlogxabs)*v3*v3);
          break;
        case 5:
          deltalm = 0.;
     // hCoeffs->delta65vh3*vh3;
          rholm  = 1. + v2*(hCoeffs->rho65v2 + hCoeffs->rho65v3*v);

          rholm=rholm + ((hCoeffs->rho65v6 + hCoeffs->rho65v6l*eulerlogxabs)*v3*v3) + (hCoeffs->rho65v4 *v2 *v2);
          break;
        case 4:
          deltalm = 0.;
       // hCoeffs->delta64vh3 * vh3;
          rholm  = 1. + v2*(hCoeffs->rho64v2 + v*(hCoeffs->rho64v3
            + hCoeffs->rho64v4*v));
          rholm=rholm + ((hCoeffs->rho64v6 + hCoeffs->rho64v6l*eulerlogxabs)*v3*v3);
          break;
        case 3:
          deltalm = 0.;
        //hCoeffs->delta63vh3 * vh3;
          rholm  = 1. + v2*(hCoeffs->rho63v2 + hCoeffs->rho63v3*v);
          rholm=rholm + ((hCoeffs->rho63v6 + hCoeffs->rho63v6l*eulerlogxabs)*v3*v3) + (hCoeffs->rho63v4 *v2 *v2);
          break;
        case 2:
          deltalm = 0.;
       // hCoeffs->delta62vh3 * vh3;
          rholm  = 1. + v2*(hCoeffs->rho62v2 + v*(hCoeffs->rho62v3
            + hCoeffs->rho62v4 * v));
          rholm=rholm + ((hCoeffs->rho62v6 + hCoeffs->rho62v6l*eulerlogxabs)*v3*v3);
          break;
        case 1:
          deltalm = 0.;
     // hCoeffs->delta61vh3 * vh3;
          rholm  = 1. + v2*(hCoeffs->rho61v2 + hCoeffs->rho61v3*v);
          rholm=rholm + ((hCoeffs->rho61v6 + hCoeffs->rho61v6l*eulerlogxabs)*v3*v3) + (hCoeffs->rho61v4 *v2 *v2);
          break;
        default:
          XLAL_ERROR( XLAL_EINVAL );
          break;
      }
      break;
    case 7:
      switch (m)
      {
        case 7:
          deltalm = 0.;
       // hCoeffs->delta77vh3 * vh3;
          rholm   = 1. + v2*(hCoeffs->rho77v2 + hCoeffs->rho77v3 * v);
          rholm=rholm + ((hCoeffs->rho77v6 + hCoeffs->rho77v6l*eulerlogxabs)*v3*v3) + (hCoeffs->rho77v4 *v2 *v2);
          break;
        case 6:
          deltalm = 0.0;
          rholm   = 1. + hCoeffs->rho76v2 * v2;
          rholm=rholm + (hCoeffs->rho76v4 *v2 *v2);
          break;
        case 5:
          deltalm = 0.0;
      //hCoeffs->delta75vh3 * vh3;
          rholm   = 1. + v2*(hCoeffs->rho75v2 + hCoeffs->rho75v3*v);
          rholm=rholm + ((hCoeffs->rho75v6 + hCoeffs->rho75v6l*eulerlogxabs)*v3*v3) + (hCoeffs->rho75v4 *v2 *v2);
          break;
        case 4:
          deltalm = 0.0;
          rholm   = 1. + hCoeffs->rho74v2 * v2;
          rholm=rholm + (hCoeffs->rho74v4 *v2 *v2);
          break;
        case 3:
          deltalm = 0.;
     // hCoeffs->delta73vh3 *vh3;
          rholm   = 1. + v2*(hCoeffs->rho73v2 + hCoeffs->rho73v3 * v);
          rholm=rholm + ((hCoeffs->rho73v6 + hCoeffs->rho73v6l*eulerlogxabs)*v3*v3) + (hCoeffs->rho73v4 *v2 *v2); 
          break;
        case 2:
          deltalm = 0.0;
          rholm   = 1. + hCoeffs->rho72v2 * v2;
          rholm=rholm + (hCoeffs->rho72v4 *v2 *v2);
          break;
        case 1:
          deltalm = 0.0;
     // hCoeffs->delta71vh3 * vh3;
          rholm   = 1. + v2*(hCoeffs->rho71v2 +hCoeffs->rho71v3 * v);
          rholm=rholm + ((hCoeffs->rho71v6 + hCoeffs->rho71v6l*eulerlogxabs)*v3*v3) + (hCoeffs->rho71v4 *v2 *v2);
          break;
        default:
          XLAL_ERROR( XLAL_EINVAL );
          break;
      }
      break;
    case 8:
      switch (m)
      {
        case 8:
          deltalm = 0.0;
          rholm   = 1. + hCoeffs->rho88v2 * v2;
          rholm=rholm + (hCoeffs->rho88v4 *v2 *v2);
          break;
        case 7:
          deltalm = 0.0;
          rholm   = 1. + hCoeffs->rho87v2 * v2;
          rholm=rholm + (hCoeffs->rho87v4 *v2 *v2);
          break;
        case 6:
          deltalm = 0.0;
          rholm   = 1. + hCoeffs->rho86v2 * v2;
          rholm=rholm + (hCoeffs->rho86v4 *v2 *v2);
          break;
        case 5:
          deltalm = 0.0;
          rholm   = 1. + hCoeffs->rho85v2 * v2;
          rholm=rholm + (hCoeffs->rho85v4 *v2 *v2);
          break;
        case 4:
          deltalm = 0.0;
          rholm  = 1. + hCoeffs->rho84v2 * v2;
          rholm=rholm + (hCoeffs->rho84v4 *v2 *v2);
          break;
        case 3:
          deltalm = 0.0;
          rholm  = 1. + hCoeffs->rho83v2 * v2;
          rholm=rholm + (hCoeffs->rho83v4 *v2 *v2);
          break;
        case 2:
          deltalm = 0.0;
          rholm  = 1. + hCoeffs->rho82v2 * v2;
          rholm=rholm + (hCoeffs->rho82v4 *v2 *v2);
          break;
        case 1:
          deltalm = 0.0;
          rholm  = 1. + hCoeffs->rho81v2 * v2;
          rholm=rholm + (hCoeffs->rho81v4 *v2 *v2);
          break;
        default:
          XLAL_ERROR( XLAL_EINVAL );
          break;
      }
  
      break;
    default:
      XLAL_ERROR( XLAL_EINVAL );
      break;
  }}
  /* Raise rholm to the lth pow */
  rholmPwrl = 1.0;
  i = l;
  while ( i-- )
  {
    rholmPwrl *= rholm;
  }
  printf("We are dealing with l : %g\n", (double) l);
  printf("and m : %g\n", (double) m);
  printf("psi : %g\n", (double) psi);
  printf("modTlm : %g\n", (double) modTlm);
  printf("arg1.val : %g\n", (double) arg1.val);
  printf("Slm : %g\n", (double) Slm);
  printf("rholmPwrl : %g\n", (double) rholmPwrl);
//CHANGE ME!!!!!!!!!!!!!!!!!!!!!!!!WILL INSPRIAL VERY FAST OTHERWISE
  *hlm = Tlm * cexp(I * deltalm) * Slm * rholmPwrl;
  
  
// *hlm = Tlm * cexp(I * deltalm) * Slm * rholmPwrl;
// *hlm = 0;

//  *hlm = hNewton;

//Add tidal corrections to l=2 and m=2 modes
if  ( l == 2 )
   {
     if( m ==2)
          {
            *hlm  = *hlm+ hlm_tidal;
          }
    }
//printf("done with making hlm......");
  return XLAL_SUCCESS;
   
} 


UNUSED static int  XLALSimIMRTNSEOBGetFactorizedWaveform_Debugging(
  COMPLEX16   * restrict hlm,    /**<< The value of hlm (populated by the function) */
                                REAL8Vector * restrict values, /**<< Vector containing dynamics r, phi, pr, pphi for a given point */
                                const REAL8 v,                 /**<< Velocity (in geometric units) */
                                const INT4  l,                 /**<< Mode l */
                                const INT4  m,                 /**<< Mode m */
                                TNSEOBParams   * restrict params 
                                )
{

//Notice :  I am feeding in v=v_phi (the non keplerian velocity for all the parts of h_lm unlike in EOBNR)


  /* Status of function calls */
  INT4 status;
  INT4 i;

  REAL8 eta, mass1, mass2, lambda1, lambda2, pi,pi2,eta2;
  REAL8 r, pr, pp, Omega, v2, vh, vh3,v3, k, hathatk, eulerlogxabs,v4,v5;
  REAL8 Hreal, Heff, Slm, deltalm, rholm, rholmPwrl;
  REAL8 delta22_leading,delta21_leading,delta33_leading,delta31_leading,delta22_Num,delta21_Num,delta33_Num,delta31_Num,delta22_Den, delta21_Den, delta33_Den,delta31_Den;
  COMPLEX16 Tlm;
  COMPLEX16 hNewton,hlm_tidal;
  gsl_sf_result lnr1, arg1, z2;

  /* Non-Keplerian velocity */
  REAL8 psi,rOmega;
  /* Pre-computed coefficients */
  TNSFacWaveformCoeffs *hCoeffs = params->hCoeffs;

  if ( abs(m) > (INT4) l )
  {
    XLAL_ERROR( XLAL_EINVAL );
  }


  eta = params->eta;
  mass1 = params->mass1;
  mass2 = params->mass2;
  lambda1 = params -> lambda1;
  lambda2 = params -> lambda2;

  /* Check our eta was sensible */
  if ( eta > 0.25 )
  {
    XLALPrintError("Eta seems to be > 0.25 - this isn't allowed!\n" );
    XLAL_ERROR( XLAL_EINVAL );
  }
  else if ( eta == 0.25 && m % 2 )
  {
    /* If m is odd and dM = 0, hLM will be zero */
    memset( hlm, 0, sizeof( COMPLEX16 ) );
    return XLAL_SUCCESS;
  }

  r  = values->data[0];
  pr = values->data[2];
  pp = values->data[3];

  //printf("accesing the data r,pr,phi");
//  printf("\n%f",r);
  Heff  = XLALEffectiveHamiltonian( eta, mass1, mass2, lambda1,lambda2, r, pr, pp);
  Hreal = sqrt( 1.0 + 2.0 * eta * ( Heff - 1.0) );
  v2    = v * v;
  v3  = v*v2;
  v4 = v2*v2;
  v5 = v3*v2;
  eulerlogxabs = LAL_GAMMA + log( 2.0 * (REAL8)m * v );
  pi=LAL_PI;
  pi2=pi*pi;
  eta2=eta*eta;

  /* Calculate the non-Keplerian velocity */
  /* given by Eq. (18) of Pan et al, PRD84, 124052(2011) */
  /* psi given by Eq. (19) of Pan et al, PRD84, 124052(2011) */
  /* Assign temporarily to vPhi */
  //nonKeplerianCoefficient gives psi
  psi = nonKeplerianCoefficient( values, eta,mass1,mass2,lambda1,lambda2 );
  /* Assign rOmega value temporarily to vPhi */
  rOmega  = r * cbrt(psi);
 Omega = v/rOmega;
   vh3   = Hreal * Omega;
     vh    = cbrt(vh3);



//REAL8 xNewt=cbrt(Omega*Omega);

printf("Omega in Newtonian Flux : %g\n", (double) Omega);
  /* Calculate the newtonian multipole */
  status =XLALSimIMRTNSEOBCalculateNewtonianMultipole( &hNewton, v*v, rOmega,
            values->data[1], (UINT4)l, m, params );
  if ( status == XLAL_FAILURE )
  {
    XLAL_ERROR( XLAL_EFUNC );
  }
  
  //printf("I calculated newtonian fl");

 /*Calculate the tidal part of h_lm*/

  status = XLALSimIMRTNSEOBCalculateTidal_hlm( &hlm_tidal, v , v/Omega,
            values->data[1], (UINT4)l, m, params, values );
  if ( status == XLAL_FAILURE )
  {
    XLAL_ERROR( XLAL_EFUNC );
  }

  /* Calculate the source term */
  if ( ( (l+m)%2 ) == 0)
  {
    Slm = Heff;
  }
  else
  {

  //In gets changed in TNSEOB! 
  // We have  Slm = v * pp in EOBNR but in TNSEOB. 
  Slm = pp/(rOmega*v);
  }
    printf("Lets check why Slm is slightly off for odd l+m");
    printf("pp : %g\n", (double) pp);
	printf("rOmega : %g\n", (double) rOmega);
	printf("v : %g\n", (double) v);

  /* Calculate the Tail term */
  k  = m * Omega;
  hathatk = Hreal * k;
  XLAL_CALLGSL( status = gsl_sf_lngamma_complex_e( l+1.0, -2.0*hathatk, &lnr1, &arg1 ) );
  if (status != GSL_SUCCESS)
  {
    XLALPrintError("Error in GSL function\n" );
    XLAL_ERROR( XLAL_EFUNC );
  }
  XLAL_CALLGSL( status = gsl_sf_fact_e( l, &z2 ) );
  if ( status != GSL_SUCCESS)
  {
    XLALPrintError("Error in GSL function\n" );
    XLAL_ERROR( XLAL_EFUNC );
  }
  Tlm = cexp( ( lnr1.val + LAL_PI * hathatk ) + I * ( 2.0 * hathatk * log(4.0*k/sqrt(LAL_E)) ) );

  
 // Tlm = cexp( ( lnr1.val + LAL_PI * hathatk ) + I * (
   //     arg1.val + 2.0 * hathatk * log(4.0*k/sqrt(LAL_E)) ) );
  Tlm /= z2.val;

REAL8 omega_check=v/rOmega;

  printf("Construction of Tlm... hathatk : %g\n", (double) hathatk);
  printf("Omega : %g\n", (double) Omega);
  printf("omega_check : %g\n", (double) omega_check);
printf("Hreal  : %g\n", (double) Hreal );
   printf("lnr1.val  : %g\n", (double) lnr1.val );
   printf("arg1.val  : %g\n", (double) arg1.val );
     printf(" z2.val  : %g\n", (double)  z2.val );




REAL8 modTlm=sqrt(creal(Tlm)*creal(Tlm) + cimag(Tlm)*cimag(Tlm));


  /* Calculate the residue phase and amplitude terms */
if (eta==0)
   {
     deltalm=0.0;
     rholm=0.0;
   }
  //This is the point test limit 
  //
  //

else
 {  //dealing with "non-test-particle"
  switch( l )
  {
    case 2:
      switch( abs(m) )
      {
        case 2:
      // CHECKME : Changing to Pade resummed coefficients 

          // deltalm = vh3*(hCoeffs->delta22vh3 + vh3*(hCoeffs->delta22vh6
          //  + vh*vh*(hCoeffs->delta22vh9*vh)))
          //  + hCoeffs->delta22v5 *v*v2*v2 + hCoeffs->delta22v8 *v2*v2*v2*v2;
          delta22_leading=(7./3.)*v3;
      delta22_Num=(808920*eta*pi*v + 137388.*pi2*v2 + 35.*eta2*(136080. + (154975. - 1359276.*eta)*v2));
      delta22_Den=(808920.*eta*pi*v + 137388.*pi2*v2 + 35.*eta2*(136080. + (154975. + 40404.*eta)*v2));
      deltalm= delta22_leading *(delta22_Den/delta22_Den);
         // rholm  = 1. + v2*(hCoeffs->rho22v2 + v*(hCoeffs->rho22v3
         //   + v*(hCoeffs->rho22v4
         //   + v*(hCoeffs->rho22v5 + v*(hCoeffs->rho22v6
         //   + hCoeffs->rho22v6l*eulerlogxabs + v*(hCoeffs->rho22v7
         //   + v*(hCoeffs->rho22v8 + hCoeffs->rho22v8l*eulerlogxabs
         //   + (hCoeffs->rho22v10 + hCoeffs->rho22v10l * eulerlogxabs)*v2)))))));
         
  rholm=1.+(v2*hCoeffs->rho22v2)+(hCoeffs->rho22v4*v4)+(hCoeffs->rho22v6*v2*v2*v2 + hCoeffs->rho22v6l*eulerlogxabs*v2*v2*v2)+(hCoeffs->rho22v8*v4*v4 + hCoeffs->rho22v8l*eulerlogxabs*v4*v4)+(hCoeffs->rho22v10*v4*v4*v2 + hCoeffs->rho22v10l * eulerlogxabs*v4*v4*v2);
         printf(" rholm=1.+(v2*hCoeffs->rho22v2)+everything else; : %g\n", (double) rholm);
		 printf("hCoeffs->rho22v2 : %g\n", (double) hCoeffs->rho22v2);
		 printf("v2 : %g\n", (double) v2);
		 printf("eta used inside the Waveform generation : %g\n", (double) eta);
		 break;
        case 1:
         // deltalm = vh3*(hCoeffs->delta21vh3 + vh3*(hCoeffs->delta21vh6
         //   + vh*(hCoeffs->delta21vh7 + (hCoeffs->delta21vh9)*vh*vh)))
         //   + hCoeffs->delta21v5*v*v2*v2 + hCoeffs->delta21v7*v2*v2*v2*v;
      delta21_leading=(2./3.)*v3;
      delta21_Num=69020.*eta + 5992.*pi*v;
      delta21_Den=5992.*pi*v + 2456.*eta*(28.+493.*eta*v2);
      deltalm=delta21_leading*(delta21_Num/delta21_Den);

          rholm  = 1. + v*(hCoeffs->rho21v1
            + v*( hCoeffs->rho21v2 + v*(hCoeffs->rho21v3 + v*(hCoeffs->rho21v4
            + v*(hCoeffs->rho21v5 + v*(hCoeffs->rho21v6 + hCoeffs->rho21v6l*eulerlogxabs
            + v*(hCoeffs->rho21v7 + hCoeffs->rho21v7l * eulerlogxabs
            + v*(hCoeffs->rho21v8 + hCoeffs->rho21v8l * eulerlogxabs
            + (hCoeffs->rho21v10 + hCoeffs->rho21v10l * eulerlogxabs)*v2))))))));
          break;
        default:
          XLAL_ERROR( XLAL_EINVAL );
          break;
      }
      break;
    case 3:
      switch (m)
      {
        case 3:
          //deltalm = vh3*(hCoeffs->delta33vh3 + vh3*(hCoeffs->delta33vh6 + hCoeffs->delta3//3vh9*vh3))
      delta33_leading=(13./10.)*v3;
      delta33_Num = 1. + 94770.*pi*v/(566279.*eta);
      delta33_Den = delta33_Num + (80897.*eta*v2/3159.);
      deltalm=delta33_leading*(delta33_Num/delta33_Den);

            //+ hCoeffs->delta33v5*v*v2*v2 + hCoeffs->delta33v7*v2*v2*v2*v;
          rholm  = 1. + v2*(hCoeffs->rho33v2 + v*(hCoeffs->rho33v3 + v*(hCoeffs->rho33v4
            + v*(hCoeffs->rho33v5 + v*(hCoeffs->rho33v6 + hCoeffs->rho33v6l*eulerlogxabs
            + v*(hCoeffs->rho33v7 + (hCoeffs->rho33v8 + hCoeffs->rho33v8l*eulerlogxabs)*v))))));
          rholm=rholm+ ((hCoeffs->rho33v10 + hCoeffs->rho33v10l*eulerlogxabs)*v5*v5);
          break;
        case 2:
         // deltalm = vh3*(hCoeffs->delta32vh3 + vh*(hCoeffs->delta32vh4 + vh*vh*(hCoeffs->delta32vh6
           // + hCoeffs->delta32vh9*vh3)));
       
      deltalm = (v3*hCoeffs->delta32vh3)+(v3*v3*hCoeffs->delta32vh6);
          rholm  = 1. + v*(hCoeffs->rho32v
            + v*(hCoeffs->rho32v2 + v*(hCoeffs->rho32v3 + v*(hCoeffs->rho32v4 + v*(hCoeffs->rho32v5
            + v*(hCoeffs->rho32v6 + hCoeffs->rho32v6l*eulerlogxabs
            + (hCoeffs->rho32v8 + hCoeffs->rho32v8l*eulerlogxabs)*v2))))));
          break;
        case 1:
        //  deltalm = vh3*(hCoeffs->delta31vh3 + vh3*(hCoeffs->delta31vh6
        //    + vh*(hCoeffs->delta31vh7 + hCoeffs->delta31vh9*vh*vh)))
       //     + hCoeffs->delta31v5*v*v2*v2;
        delta31_leading = (13./30.)*v3;
      delta31_Num = 4641.*eta + 1690.*pi*v;
      delta31_Den = delta31_Num + 18207.*eta2*v2;
      deltalm=delta31_leading*(delta31_Num/delta31_Den);

          rholm  = 1. + v2*(hCoeffs->rho31v2 + v*(hCoeffs->rho31v3 + v*(hCoeffs->rho31v4
            + v*(hCoeffs->rho31v5 + v*(hCoeffs->rho31v6 + hCoeffs->rho31v6l*eulerlogxabs
            + v*(hCoeffs->rho31v7 + (hCoeffs->rho31v8 + hCoeffs->rho31v8l*eulerlogxabs)*v))))));
          rholm=rholm+ ((hCoeffs->rho31v10 + hCoeffs->rho31v10l*eulerlogxabs)*v5*v5);
          break;
        default:
          XLAL_ERROR( XLAL_EINVAL );
          break;
      }
      break;
    case 4:
      switch (m)
      {
        case 4:
          deltalm = v3*(hCoeffs->delta44vh3 + hCoeffs->delta44vh6 *v3);
     //deltalm = vh3*(hCoeffs->delta44vh3 + hCoeffs->delta44vh6 *vh3)
     //1034             + hCoeffs->delta44v5*v2*v2*v;
     //
          rholm  = 1. + v2*(hCoeffs->rho44v2
            + v*( hCoeffs->rho44v3 + v*(hCoeffs->rho44v4
            + v*(hCoeffs->rho44v5 + (hCoeffs->rho44v6
            + hCoeffs->rho44v6l*eulerlogxabs)*v))));
          rholm=rholm + ((hCoeffs->rho44v8 + hCoeffs->rho44v8l*eulerlogxabs)*v4*v4);
          break;
        case 3:
          //deltalm = v3*(hCoeffs->delta43vh3 + vh*(hCoeffs->delta43vh4
          //  + hCoeffs->delta43vh6*vh*vh));
      deltalm = (v3*hCoeffs->delta43vh3) + (v3*v3*hCoeffs->delta43vh6);

          rholm  = 1. + v*(hCoeffs->rho43v
            + v*(hCoeffs->rho43v2
            + v2*(hCoeffs->rho43v4 + v*(hCoeffs->rho43v5
            + (hCoeffs->rho43v6 + hCoeffs->rho43v6l*eulerlogxabs)*v))));
          rholm=rholm + ((hCoeffs->rho43v8 + hCoeffs->rho43v8l*eulerlogxabs)*v4*v4);
          break;
        case 2:
          //deltalm = vh3*(hCoeffs->delta42vh3 + hCoeffs->delta42vh6*vh3);
          deltalm = (v3*hCoeffs->delta42vh3) + (v3*v3*hCoeffs->delta42vh6);

      rholm  = 1. + v2*(hCoeffs->rho42v2
            + v*(hCoeffs->rho42v3 + v*(hCoeffs->rho42v4 + v*(hCoeffs->rho42v5
            + (hCoeffs->rho42v6 + hCoeffs->rho42v6l*eulerlogxabs)*v))));
            rholm=rholm + ((hCoeffs->rho42v8 + hCoeffs->rho42v8l*eulerlogxabs)*v4*v4);
          break;
        case 1:
       //   deltalm = vh3*(hCoeffs->delta41vh3 + vh*(hCoeffs->delta41vh4
       //     + hCoeffs->delta41vh6*vh*vh));
        deltalm = (v3*hCoeffs->delta41vh3)+(v3*v3*+ hCoeffs->delta41vh6);

          rholm  = 1. + v*(hCoeffs->rho41v
            + v*(hCoeffs->rho41v2
            + v2*(hCoeffs->rho41v4 + v*(hCoeffs->rho41v5
            + (hCoeffs->rho41v6 +  hCoeffs->rho41v6l*eulerlogxabs)*v))));
          rholm=rholm + ((hCoeffs->rho41v8 + hCoeffs->rho41v8l*eulerlogxabs)*v4*v4);
          break;
        default:
          XLAL_ERROR( XLAL_EINVAL );
          break;
      }
      break;
    case 5:
      switch (m)
      {
        case 5:
          deltalm = hCoeffs->delta55vh3*v3;
          rholm  = 1. + v2*( hCoeffs->rho55v2
            + v*(hCoeffs->rho55v3 + v*(hCoeffs->rho55v4
            + v*(hCoeffs->rho55v5 + hCoeffs->rho55v6*v))));
          rholm = rholm + (hCoeffs->rho55v6l*eulerlogxabs*v3*v3)+ ((hCoeffs->rho55v8 + hCoeffs->rho55v8l*eulerlogxabs)*v4*v4);
          break;
        case 4:
          deltalm =0.;
       // vh3*(hCoeffs->delta54vh3 + hCoeffs->delta54vh4*vh);
          rholm  = 1. + v2*(hCoeffs->rho54v2 + v*(hCoeffs->rho54v3
            + hCoeffs->rho54v4*v));
          rholm=rholm + ((hCoeffs->rho54v6 + hCoeffs->rho54v6l*eulerlogxabs)*v3*v3);
          break;
        case 3:
          deltalm =0.;
       // hCoeffs->delta53vh3 * vh3;
          rholm  = 1. + v2*(hCoeffs->rho53v2
            + v*(hCoeffs->rho53v3 + v*(hCoeffs->rho53v4 + hCoeffs->rho53v5*v)));
          rholm=rholm +  ((hCoeffs->rho53v6 + hCoeffs->rho53v6l*eulerlogxabs)*v3*v3) + ((hCoeffs->rho53v8 + hCoeffs->rho53v8l*eulerlogxabs)*v4*v4);
          break;
        case 2:
          deltalm = 0.;
        //vh3*(hCoeffs->delta52vh3 + hCoeffs->delta52vh4*vh);
          rholm  = 1. + v2*(hCoeffs->rho52v2 + v*(hCoeffs->rho52v3
            + hCoeffs->rho52v4*v));
          rholm=rholm + ((hCoeffs->rho52v6 + hCoeffs->rho52v6l*eulerlogxabs)*v3*v3);
          break;
        case 1:
          deltalm = 0.;
        //hCoeffs->delta51vh3 * vh3;
          rholm  = 1. + v2*(hCoeffs->rho51v2
            + v*(hCoeffs->rho51v3 + v*(hCoeffs->rho51v4 + hCoeffs->rho51v5*v)));
           rholm = rholm + ((hCoeffs->rho51v6 +hCoeffs->rho51v6l*eulerlogxabs)*v3*v3)+ ((hCoeffs->rho51v8 + hCoeffs->rho51v8l*eulerlogxabs)*v4*v4);
          break;
        default:
          XLAL_ERROR( XLAL_EINVAL );
          break;
      }
      break;
    case 6:
      switch (m)
      {
        case 6:
          deltalm = 0.;
        //hCoeffs->delta66vh3*vh3;
          rholm  = 1. + v2*(hCoeffs->rho66v2 + v*(hCoeffs->rho66v3
            + hCoeffs->rho66v4*v));
          rholm=rholm + ((hCoeffs->rho66v6 + hCoeffs->rho66v6l*eulerlogxabs)*v3*v3);
          break;
        case 5:
          deltalm = 0.;
     // hCoeffs->delta65vh3*vh3;
          rholm  = 1. + v2*(hCoeffs->rho65v2 + hCoeffs->rho65v3*v);

          rholm=rholm + ((hCoeffs->rho65v6 + hCoeffs->rho65v6l*eulerlogxabs)*v3*v3) + (hCoeffs->rho65v4 *v2 *v2);
          break;
        case 4:
          deltalm = 0.;
       // hCoeffs->delta64vh3 * vh3;
          rholm  = 1. + v2*(hCoeffs->rho64v2 + v*(hCoeffs->rho64v3
            + hCoeffs->rho64v4*v));
          rholm=rholm + ((hCoeffs->rho64v6 + hCoeffs->rho64v6l*eulerlogxabs)*v3*v3);
          break;
        case 3:
          deltalm = 0.;
        //hCoeffs->delta63vh3 * vh3;
          rholm  = 1. + v2*(hCoeffs->rho63v2 + hCoeffs->rho63v3*v);
          rholm=rholm + ((hCoeffs->rho63v6 + hCoeffs->rho63v6l*eulerlogxabs)*v3*v3) + (hCoeffs->rho63v4 *v2 *v2);
          break;
        case 2:
          deltalm = 0.;
       // hCoeffs->delta62vh3 * vh3;
          rholm  = 1. + v2*(hCoeffs->rho62v2 + v*(hCoeffs->rho62v3
            + hCoeffs->rho62v4 * v));
          rholm=rholm + ((hCoeffs->rho62v6 + hCoeffs->rho62v6l*eulerlogxabs)*v3*v3);
          break;
        case 1:
          deltalm = 0.;
     // hCoeffs->delta61vh3 * vh3;
          rholm  = 1. + v2*(hCoeffs->rho61v2 + hCoeffs->rho61v3*v);
          rholm=rholm + ((hCoeffs->rho61v6 + hCoeffs->rho61v6l*eulerlogxabs)*v3*v3) + (hCoeffs->rho61v4 *v2 *v2);
          break;
        default:
          XLAL_ERROR( XLAL_EINVAL );
          break;
      }
      break;
    case 7:
      switch (m)
      {
        case 7:
          deltalm = 0.;
       // hCoeffs->delta77vh3 * vh3;
          rholm   = 1. + v2*(hCoeffs->rho77v2 + hCoeffs->rho77v3 * v);
          rholm=rholm + ((hCoeffs->rho77v6 + hCoeffs->rho77v6l*eulerlogxabs)*v3*v3) + (hCoeffs->rho77v4 *v2 *v2);
          break;
        case 6:
          deltalm = 0.0;
          rholm   = 1. + hCoeffs->rho76v2 * v2;
          rholm=rholm + (hCoeffs->rho76v4 *v2 *v2);
          break;
        case 5:
          deltalm = 0.0;
      //hCoeffs->delta75vh3 * vh3;
          rholm   = 1. + v2*(hCoeffs->rho75v2 + hCoeffs->rho75v3*v);
          rholm=rholm + ((hCoeffs->rho75v6 + hCoeffs->rho75v6l*eulerlogxabs)*v3*v3) + (hCoeffs->rho75v4 *v2 *v2);
          break;
        case 4:
          deltalm = 0.0;
          rholm   = 1. + hCoeffs->rho74v2 * v2;
          rholm=rholm + (hCoeffs->rho74v4 *v2 *v2);
          break;
        case 3:
          deltalm = 0.;
     // hCoeffs->delta73vh3 *vh3;
          rholm   = 1. + v2*(hCoeffs->rho73v2 + hCoeffs->rho73v3 * v);
          rholm=rholm + ((hCoeffs->rho73v6 + hCoeffs->rho73v6l*eulerlogxabs)*v3*v3) + (hCoeffs->rho73v4 *v2 *v2); 
          break;
        case 2:
          deltalm = 0.0;
          rholm   = 1. + hCoeffs->rho72v2 * v2;
          rholm=rholm + (hCoeffs->rho72v4 *v2 *v2);
          break;
        case 1:
          deltalm = 0.0;
     // hCoeffs->delta71vh3 * vh3;
          rholm   = 1. + v2*(hCoeffs->rho71v2 +hCoeffs->rho71v3 * v);
          rholm=rholm + ((hCoeffs->rho71v6 + hCoeffs->rho71v6l*eulerlogxabs)*v3*v3) + (hCoeffs->rho71v4 *v2 *v2);
          break;
        default:
          XLAL_ERROR( XLAL_EINVAL );
          break;
      }
      break;
    case 8:
      switch (m)
      {
        case 8:
          deltalm = 0.0;
          rholm   = 1. + hCoeffs->rho88v2 * v2;
          rholm=rholm + (hCoeffs->rho88v4 *v2 *v2);
          break;
        case 7:
          deltalm = 0.0;
          rholm   = 1. + hCoeffs->rho87v2 * v2;
          rholm=rholm + (hCoeffs->rho87v4 *v2 *v2);
          break;
        case 6:
          deltalm = 0.0;
          rholm   = 1. + hCoeffs->rho86v2 * v2;
          rholm=rholm + (hCoeffs->rho86v4 *v2 *v2);
          break;
        case 5:
          deltalm = 0.0;
          rholm   = 1. + hCoeffs->rho85v2 * v2;
          rholm=rholm + (hCoeffs->rho85v4 *v2 *v2);
          break;
        case 4:
          deltalm = 0.0;
          rholm  = 1. + hCoeffs->rho84v2 * v2;
          rholm=rholm + (hCoeffs->rho84v4 *v2 *v2);
          break;
        case 3:
          deltalm = 0.0;
          rholm  = 1. + hCoeffs->rho83v2 * v2;
          rholm=rholm + (hCoeffs->rho83v4 *v2 *v2);
          break;
        case 2:
          deltalm = 0.0;
          rholm  = 1. + hCoeffs->rho82v2 * v2;
          rholm=rholm + (hCoeffs->rho82v4 *v2 *v2);
          break;
        case 1:
          deltalm = 0.0;
          rholm  = 1. + hCoeffs->rho81v2 * v2;
          rholm=rholm + (hCoeffs->rho81v4 *v2 *v2);
          break;
        default:
          XLAL_ERROR( XLAL_EINVAL );
          break;
      }
  
      break;
    default:
      XLAL_ERROR( XLAL_EINVAL );
      break;
  }}
  /* Raise rholm to the lth pow */
  rholmPwrl = 1.0;
  i = l;
  while ( i-- )
  {
    rholmPwrl *= rholm;
  }
  printf("We are dealing with l : %g\n", (double) l);
  printf("and m : %g\n", (double) m);
  printf("psi : %g\n", (double) psi);
  printf("modTlm : %g\n", (double) modTlm);
  printf("arg1.val : %g\n", (double) arg1.val);
  printf("Slm : %g\n", (double) Slm);
  printf("rholmPwrl : %g\n", (double) rholmPwrl);
  printf("hNewton : %g\n", (double) hNewton);
//CHANGE ME!!!!!!!!!!!!!!!!!!!!!!!!WILL INSPRIAL VERY FAST OTHERWISE
  *hlm = Tlm * cexp(I * deltalm) * Slm * rholmPwrl;
  *hlm *= hNewton;
  
// *hlm = Tlm * cexp(I * deltalm) * Slm * rholmPwrl;
// *hlm = 0;

//  *hlm = hNewton;

//Add tidal corrections to l=2 and m=2 modes
if  ( l == 2 )
   {
     if( m ==2)
          {
            *hlm  = *hlm+ hlm_tidal;
          }
    }
//printf("done with making hlm......");
  return XLAL_SUCCESS;
} 

#endif /*_LALSIMIMRFACTORIZEDWAVEFORM_C*/
