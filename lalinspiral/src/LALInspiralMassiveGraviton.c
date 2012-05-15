    /*
*  Copyright (C) 2007 Jolien Creighton, B.S. Sathyaprakash, Thomas Cokelaer, Tjonnie G.F. Li, Walter Del Pozzo
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
#include "LALInspiral.h"

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

// NB: PROTOTYPE RESIDES IN LALINSPIRAL.H

/*  <lalVerbatim file="LALInspiralMassiveGraviton"> */
void LALInspiralMassiveGravitonPhasing(
                          InspiralTemplate *params,
                          REAL8 f,
                          REAL8 *psif);
int
XLALInspiralMassiveGraviton (REAL4Vector      *signalvec,
                             InspiralTemplate *params)
{ /* </lalVerbatim>  */
   REAL8 Oneby3, UNUSED h1, UNUSED h2, pimmc, f, v, df, shft, phi, amp0, amp, psif, psi;
   INT4 n, nby2, i, f0, fn;
   expnCoeffs ak;
   expnFunc func;

   /* Perform some initial checks */
   if (signalvec == NULL)
     XLAL_ERROR(XLAL_EFAULT);
   if (signalvec->data == NULL)
     XLAL_ERROR(XLAL_EFAULT);
   if (params == NULL)
     XLAL_ERROR(XLAL_EFAULT);  
   if (signalvec->length<=2)
     XLAL_ERROR(XLAL_EBADLEN);
     
   n = signalvec->length;
   nby2 = n/2;
   memset( &ak, 0, sizeof( ak ) );
   if ( XLALInspiralSetup(&ak, params) == XLAL_FAILURE )
     XLAL_ERROR(XLAL_EFUNC);
   if ( XLALInspiralChooseModel(&func, &ak, params) == XLAL_FAILURE)
     XLAL_ERROR(XLAL_EFUNC);

   
   Oneby3 = 1.L/3.L;
   df = params->tSampling/signalvec->length;
   pimmc = LAL_PI * params->totalMass * LAL_MTSUN_SI;
   f0 = params->fLower;
   fn = (params->fCutoff < ak.fn) ? params->fCutoff : ak.fn;
   v = cbrt(pimmc*f0);

   /* If we want to pad with zeroes in the beginning then the instant of
    * coalescence will be the chirp time + the duration for which padding
    * is needed. Thus, in the equation below nStartPad occurs with a +ve sign.
    * This code doesn't support non-zero start-time. i.e. params->startTime
    * should be necessarily zero.
    */
   shft = 2.L*LAL_PI * (ak.tn + params->nStartPad/params->tSampling + params->startTime);
   phi =  params->startPhase + LAL_PI/4.L;
   amp0 = params->signalAmplitude * ak.totalmass * sqrt(LAL_PI/12.L) * df;
/*
   Compute the standard stationary phase approximation.
*/
   h1 = signalvec->data[0] = 0.L;
   h2 = signalvec->data[nby2] = 0.L;
    

//    for (int k=0;k<10;k++) fprintf(stderr,"dphi%i = %e\n",k,dphis[k]);
//	FILE* model_output;
//	model_output=fopen("output_TF2T.dat","w");

//	fprintf(model_output,"Sampling frequency: %lf\n",params->tSampling);

//	fprintf(stderr,"Mass 1: %lf\n",params->mass1);
//	fprintf(stderr,"Mass 2: %lf\n",params->mass2);

   for (i=1; i<nby2; i++) {
      f = i * df;
      if (f < f0 || f > fn)
      {
	      /*
	       * All frequency components below f0 and above fn are set to zero
	       */
	      signalvec->data[i] = 0.;
	      signalvec->data[n-i] = 0.;
      }
      else
      {
	      v = pow(pimmc * f, Oneby3);
	      LALInspiralMassiveGravitonPhasing(params, f, &psif);
	      psi = shft * f + phi + psif;
	      /*
		 dEnergybyFlux computes 1/(dv/dt) while what we need is 1/(dF/dt):
		 dF/dt=(dF/dv)(dv/dt)=-dEnergybyFlux/(dF/dv)=-dEnergybyFlux 3v^2/(pi m^2)
		 Note that our energy is defined as E=-eta v^2/2 and NOT -eta m v^2/2.
		 This is the reason why there is an extra m in the last equation above
		 amp = amp0 * pow(-dEnergybyFlux(v)/v^2, 0.5) * v^2;
		     = amp0 * pow(-dEnergybyFlux(v), 0.5) * v;
	      */
	      amp = amp0 * pow(-func.dEnergy(v,&ak)/func.flux(v,&ak),0.5L) * v;
	      signalvec->data[i] = (REAL4) (amp * cos(psi));
	      signalvec->data[n-i] = (REAL4) (-amp * sin(psi));
//          fprintf(model_output,"%e\t %e\t %e\t %e\n",i*df,signalvec->data[i],signalvec->data[n-i],psif);  
      }
      	
      /*
	 printf ("%e %e \n", v, psif);
	 printf ("%e %e %e %e %e\n", f, pow(h1,2.)+pow(h2,2.), h2, psi, psif);
	 printf ("&\n");
       */

   }
//   fclose(model_output);
//   exit(0);
   params->fFinal = fn;
   DETATCHSTATUSPTR(status);
   RETURN(status);
}

void LALInspiralMassiveGravitonPhasing(
                          InspiralTemplate *params,
                          REAL8 f,
                          REAL8 *psif)
{
	REAL8 phaseParams[10]={0.0};
	REAL8 mtot = params->totalMass;
	REAL8 lambdaG;
	REAL8 distance;
    REAL8 eta = params->eta;
    REAL8 pimtot = LAL_PI*mtot*LAL_MTSUN_SI;
    REAL8 comprefac = 3.0/(128.0*eta);
    // x is an alias for f^1/3
    REAL8 x=cbrt(f);
    REAL8 psif_loc=0.0;
    REAL8 pimtot1by3=cbrt(pimtot);
    lambdaG=pow(10.0,params->loglambdaG);
    distance=params->distance*LAL_PC_SI*1e6;
    
    // SEE arXiv:1005.0304
    
    phaseParams[0] = comprefac*(1.0/(pimtot1by3*pimtot1by3*pimtot1by3*pimtot1by3*pimtot1by3)); //phi0
    phaseParams[1] = comprefac*(1.0/(pimtot1by3*pimtot1by3*pimtot1by3*pimtot1by3))* 0.0; //phi1
    phaseParams[2] = comprefac*(1.0/pimtot)* (3715.0/756.0 + 55.0/9.0*eta); //phi2
    phaseParams[3] = comprefac*(1.0/(pimtot1by3*pimtot1by3))* -16.0*LAL_PI; //phi3
    phaseParams[4] = comprefac*(1.0/pimtot1by3)* (15293365.0/508032.0 + 27145.0/504.0*eta + 3085.0/72.0*eta*eta); // phi4
    phaseParams[5] = comprefac*LAL_PI*((38645.0/756.0 - 65.0/9.0*eta)+((38645.0/756.0 - 65.0/9.0*eta)*log(pimtot*pow(6.0, 1.5)))); //phi5
    phaseParams[6] = comprefac*LAL_PI*(38645.0/756.0 - 65.0/9.0*eta); //phi5l
    phaseParams[7] = comprefac*pimtot1by3* ((11583231236531.0/4694215680.0 - 640.0/3.0*(LAL_PI*LAL_PI) - 6848.0/21.0*LAL_GAMMA) + eta*(-15335597827.0/3048192.0 + 2255.0/12.0*(LAL_PI*LAL_PI) + 47324.0/63.0-7948.0/9.0) + 76055.0/1728.0*eta*eta - 127825.0/1296.0*eta*eta*eta + -6848.0/21.0*log(4.0*pimtot1by3)); //phi6
    phaseParams[8] = comprefac*pimtot1by3* -6848.0/63.0; //phi6l
    phaseParams[9] = comprefac*pimtot1by3*pimtot1by3* LAL_PI*(77096675.0/254016.0 + 378515.0/1512.0*eta - 74045.0/756.0*eta*eta); //phi7
    // FILL PHASE PARAMETERS
    switch (params->order)
    {
            case LAL_PNORDER_THREE_POINT_FIVE:
                psif_loc +=phaseParams[9]*x*x;
            case LAL_PNORDER_THREE:
                psif_loc +=phaseParams[8]*x*log(f)
                         + phaseParams[7]*x;
            case LAL_PNORDER_TWO_POINT_FIVE:
                psif_loc +=phaseParams[6]*log(f)
                         + phaseParams[5];
            case LAL_PNORDER_TWO:
                psif_loc += phaseParams[4]/x;
            case LAL_PNORDER_ONE_POINT_FIVE:
                psif_loc += phaseParams[3]/(x*x);
            case LAL_PNORDER_ONE:
                psif_loc += (phaseParams[2]-LAL_PI*distance/(lambdaG*lambdaG))/f; // here goes the massive graviton term
            case LAL_PNORDER_HALF:
                psif_loc += phaseParams[1]/(f*x);
            case LAL_PNORDER_NEWTONIAN:
                psif_loc += phaseParams[0]/(f*x*x);
                break;
            default:
                printf("INVALID PN ORDER!");
                exit(-1);
    }
    *psif=psif_loc;    
    return;
}

