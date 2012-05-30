/*
*  Copyright (C) 2007 Jolien Creighton, B.S. Sathyaprakash, Thomas Cokelaer, Tjonnie G.F. Li, Walter Del Pozzo, Michalis Agathos
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

/**** <lalVerbatim file="LALInspiralStationaryPhaseApprox2CV">
 * Author: B.S. Sathyaprakash
 **** </lalVerbatim> */

/**** <lalLaTeX>
 *
 *
 * \subsection{Module \texttt{LALInspiralStationaryPhaseApprox2.c}}
 * %% A one-line description of the function(s) defined in this module.
 * This module computes the usual stationary phase approximation to the
 * Fourier transform of a chirp waveform
 * given by Eq.~(\ref{eq:InspiralFourierPhase:f2}).
 *
 * \subsubsection*{Prototypes}
 * \input{LALInspiralStationaryPhaseApprox2CP}
 * \idx{LALInspiralStationaryPhaseApprox2()}
 * \begin{itemize}
 * \item {\tt signalvec:} Output containing the inspiral waveform.
 * \item {\tt params:} Input containing binary chirp parameters.
 * \end{itemize}
 *
 * \subsubsection*{Description}
 *
 * %% A description of the data analysis task performed by this function;
 * %% this is the main place to document the module.
 * Computes the Fourier transform of the chirp signal in the stationary
 * phase approximation and returns the real and imagninary parts of the
 * Fourier domain signal in the convention of fftw. For a signal vector
 * of length {\tt n=signalvec->length} ({\tt n} even):
 * \begin{itemize}
 * \item {\tt signalvec->data[0]} is the {\it real} 0th frequency component of the Fourier transform.
 * \item {\tt signalvec->data[n/2]} is the {\it real} Nyquist frequency component of the Fourier transform.
 * \item {\tt signalvec->data[k]} and {\tt signalvec->data[n-k],} for {\tt k=1,\ldots, n/2-1,} are
 * the real and imaginary parts of the Fourier transform at a frequency $k\Delta f=k/T,$ $T$ being
 * the duration of the signal and $\Delta f=1/T$ is the frequency resolution.
 * \end{itemize}
 *
 * \subsubsection*{Algorithm}
 *
 * %% A description of the method used to perform the calculation.
 *
 * The standard SPA is given by Eq.~(\ref{eq:InspiralFourierPhase:f2}).
 * We define a variable function pointer {\tt LALInspiralTaylorF2Phasing} and point
 * it to one of the {\texttt static} functions defined within this function
 * that explicitly calculates the Fourier phase at the PN order chosen by the user.
 * The reference points are chosen so that on inverse Fourier transforming
 * the time-domain waveform will
 * \begin{itemize}
 * \item be padded with zeroes in the first {\tt params->nStartPad} bins,
 * \item begin with a phase shift of {\tt params->nStartPhase} radians,
 * \item have an amplitude of ${\tt n} v^2.$
 * \end{itemize}
 * \subsubsection*{Uses}
 * \begin{verbatim}
   LALInspiralSetup
   LALInspiralChooseModel
   LALInspiralTaylorF2Phasing[0234567]PN
 * \end{verbatim}
 *
 * %% List of any external functions called by this function.
 * \begin{verbatim}
 * None
 * \end{verbatim}
 * \subsubsection*{Notes}
 *
 * %% Any relevant notes.
 *
 * If it is required to compare the output of this module with a time domain
 * signal one should use an inverse Fourier transform routine that packs data
 * in the same way as fftw. Moreover, one should divide the resulting inverse
 * Fourier transform by a factor ${\tt n}/2$ to be consistent with the
 * amplitude used in time-domain signal models.
 *
 * \vfill{\footnotesize\input{LALInspiralStationaryPhaseApprox2CV}}
 *
 **** </lalLaTeX> */

#include "LALInspiral.h"
#include "LALInspiralStationaryPhaseApprox2Test.h"

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

// NB: PROTOTYPE RESIDES IN LALINSPIRAL.H

NRCSID (LALINSPIRALSTATIONARYPHASEAPPROX2TESTC, "$Id$");

/*  <lalVerbatim file="LALInspiralStationaryPhaseApprox2CP"> */
void
LALInspiralStationaryPhaseApprox2Test (
                                       LALStatus        *status,
                                       REAL4Vector      *signalvec,
                                       InspiralTemplate *params,
                                       REAL8 *dphis,
                                       REAL4 cutoff)
{ /* </lalVerbatim>  */
   REAL8 Oneby3, UNUSED h1, UNUSED h2, pimmc, f, v, df, shft, phi, amp0, amp, psif, psi;
   INT4 n, nby2, i, f0, fn;
   expnCoeffs ak;
   expnFunc func;
   //void (*LALInspiralTaylorF2PhasingTest)(InspiralTemplate *, REAL8 , REAL8 *) = NULL;

   INITSTATUS (status, "LALInspiralStationaryPhaseApprox2Test", LALINSPIRALSTATIONARYPHASEAPPROX2TESTC);
   ATTATCHSTATUSPTR(status);
   ASSERT (signalvec,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (signalvec->data,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (signalvec->length>2,  status, LALINSPIRALH_ECHOICE, LALINSPIRALH_MSGECHOICE);
   n = signalvec->length;
   nby2 = n/2;
   memset( &ak, 0, sizeof( ak ) );
   LALInspiralSetup(status->statusPtr, &ak, params);
   CHECKSTATUSPTR(status);
   LALInspiralChooseModel(status->statusPtr, &func, &ak, params);
   CHECKSTATUSPTR(status);

   
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
    
    /* FILL PHASE COEFFICIENTS */
    REAL8 phaseParams[10] = {0.0};
    
    TaylorF2fillPhaseParams(params, phaseParams, dphis);
//    for (int k=0;k<10;k++) fprintf(stderr,"dphi%i = %e\n",k,dphis[k]);
//	FILE* model_output;
//	model_output=fopen("output_TF2T.dat","w");

//	fprintf(model_output,"Sampling frequency: %lf\n",params->tSampling);

//	fprintf(stderr,"Mass 1: %lf\n",params->mass1);
//	fprintf(stderr,"Mass 2: %lf\n",params->mass2);

   for (i=1; i<nby2; i++) {
      f = i * df;
      if (f < f0 || f > fn || ( cutoff > f0 && f > cutoff ))
      {
	      /*
	       * All frequency components below f0 and above min{fn,cutoff} are set to zero
	       */
	      signalvec->data[i] = 0.;
	      signalvec->data[n-i] = 0.;
      }
      else
      {
	      v = pow(pimmc * f, Oneby3);
	      LALInspiralTaylorF2PhasingTest(params, f, phaseParams, &psif);
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

void LALInspiralTaylorF2PhasingTest(
                          InspiralTemplate *params,
                          REAL8 f,
                          REAL8 *phaseParams,
                          REAL8 *psif)
{
    // x is an alias for f^1/3
    REAL8 x=cbrt(f);
    REAL8 psif_loc=0.0;
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
                psif_loc += phaseParams[2]/f;
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

void TaylorF2fillPhaseParams(
                                         InspiralTemplate *params,
                                         REAL8 *phaseParams,
                                         REAL8 *dphis)
{
    
    // SYSTEM DEPENDENT PARAMETER - DUMMIES, NEED TO GET FROM PARAMETER STRUCTURES
    REAL8 mtot = params->totalMass;
    REAL8 eta = params->eta;
    // We only need spin magnitudes and signs and we use the 1st component of each spin vector for that.
    REAL8 spin1 = params->spin1[0];
    REAL8 spin2 = params->spin2[0];
    
    UINT4 i;
    REAL8 pimtot = LAL_PI*mtot*LAL_MTSUN_SI;
    REAL8 comprefac = 3.0/(128.0*eta);
    //REAL8 comprefac = 3.0/(256.0*eta);
    
    // POPULATE INDIVIDUAL PHASE PARAMETERS
    // SEE arXiv:gr-qc/0411146
    /*
    phaseParams[0] = 3.0/(128.0*eta)*pow(LAL_PI*LAL_MTSUN_SI*mtot,-5.0/3.0)* 1.0; //phi0
    phaseParams[1] = 3.0/(128.0*eta)*pow(LAL_PI*LAL_MTSUN_SI*mtot,-4.0/3.0)* 0.0; //phi1
    phaseParams[2] = 3.0/(128.0*eta)*pow(LAL_PI*LAL_MTSUN_SI*mtot,-3.0/3.0)* 20.0/9.0*(743.0/336.0 + 11.0/4.0*pow(eta,2.0)); //phi2
    phaseParams[3] = 3.0/(128.0*eta)*pow(LAL_PI*LAL_MTSUN_SI*mtot,-2.0/3.0)* -16.0*LAL_PI; //phi3
    phaseParams[4] = 3.0/(128.0*eta)*pow(LAL_PI*LAL_MTSUN_SI*mtot,-1.0/3.0)* 10.0*(3058673.0/1016064.0 + 5429.0/1008.0*eta + 617.0/144.0*pow(eta,2.0)); // phi4
    phaseParams[5] = 3.0/(128.0*eta)*pow(LAL_PI*LAL_MTSUN_SI*mtot,-0.0/3.0)* LAL_PI*(38645.0/756.0 - 65.0/9.0*eta); //phi5
    phaseParams[6] = 3.0/(128.0*eta)*pow(LAL_PI*LAL_MTSUN_SI*mtot,-0.0/3.0)* 3.0*LAL_PI*(38645.0/756.0 - 65.0/9.0*eta); //phi5l
    phaseParams[7] = 3.0/(128.0*eta)*pow(LAL_PI*LAL_MTSUN_SI*mtot,1.0/3.0)* ((11583231236531.0/4694215680.0 - 640.0/3.0*pow(LAL_PI, 2.0) - 6848.0/21.0*LAL_GAMMA) + eta*(-15335597827.0/3048192.0 + 2255.0/12.0*pow(LAL_PI, 2.0) - 1760.0/3.0*LAL_THETA + 12320.0/9.0*LAL_LAMBDA) + 76055.0/1728.0*pow(eta, 2.0) - 127825.0/1296.0*pow(eta, 3.0) + -6848.0/21.0*log(4*pow(LAL_PI*mtot, 1.0/3.0))); //phi6 + some factors originally coming from phi6l
    phaseParams[8] = 3.0/(128.0*eta)*pow(LAL_PI*LAL_MTSUN_SI*mtot,1.0/3.0)* -6848.0/21.0; //phi6l
    phaseParams[9] = 3.0/(128.0*eta)*pow(LAL_PI*LAL_MTSUN_SI*mtot,2.0/3.0)* LAL_PI*(77096675.0/254016.0 + 378515.0/1512.0*eta - 74045.0/756.0*pow(eta, 2.0)); //phi7
     */
     // pimtot1by3 is an alias for (pi*m)^(1/3)
    REAL8 pimtot1by3=cbrt(pimtot);

    // Spin related parameters 
    REAL8 spinbeta, spinsigma, spingamma, massdelta, spin_a, spin_s;

    spinbeta = 1./12.*((113.*pow((params->mass1/(params->mass1+params->mass2)),2.) + 
		  75.*params->eta)*spin1 + (113.*pow((params->mass2/(params->mass1+params->mass2)),2.)+
		  75.*params->eta)*spin2);
    spinsigma = params->eta/48.*474.*spin1*spin2;
    spin_a = 0.5*(spin1-spin2);
    spin_s = 0.5*(spin1+spin2);
    massdelta = (params->mass1 - params->mass2)/(params->mass1+params->mass2);
    spingamma = (732985./2268. - 24260./81.*params->eta - 340./9.*params->eta*params->eta)*spin_s + 
          (732985./2268. + 140./9.*params->eta)*massdelta*spin_a;

    
    // SEE arXiv:0810.5336  NOTE: arXiv version (v3) includes errata!
    phaseParams[0] = comprefac*(1.0/(pimtot1by3*pimtot1by3*pimtot1by3*pimtot1by3*pimtot1by3)); //phi0
    /* this is needed otherwise phi1 is identically 0 regardless of the dphi1 */
    phaseParams[1] = comprefac*(1.0/(pimtot1by3*pimtot1by3*pimtot1by3*pimtot1by3))* dphis[1]; //phi1
    phaseParams[2] = comprefac*(1.0/pimtot)* (3715.0/756.0 + 55.0/9.0*eta); //phi2
    phaseParams[3] = comprefac*(1.0/(pimtot1by3*pimtot1by3))* -16.0*LAL_PI + 4.*spinbeta; //phi3
    phaseParams[4] = comprefac*(1.0/pimtot1by3)* (15293365.0/508032.0 + 27145.0/504.0*eta + 3085.0/72.0*eta*eta - 10.0*spinsigma); // phi4
    phaseParams[5] = comprefac*LAL_PI*((38645.0/756.0 - 65.0/9.0*eta - spingamma/LAL_PI)+((38645.0/756.0 - 65.0/9.0*eta - spingamma/LAL_PI)*log(pimtot*pow(6.0, 1.5)))); //phi5
    phaseParams[6] = comprefac*LAL_PI*(38645.0/756.0 - 65.0/9.0*eta - spingamma/LAL_PI); //phi5l
    phaseParams[7] = comprefac*pimtot1by3* ((11583231236531.0/4694215680.0 - 640.0/3.0*(LAL_PI*LAL_PI) - 6848.0/21.0*LAL_GAMMA) + eta*(-15335597827.0/3048192.0 + 2255.0/12.0*(LAL_PI*LAL_PI) + 47324.0/63.0-7948.0/9.0) + 76055.0/1728.0*eta*eta - 127825.0/1296.0*eta*eta*eta + -6848.0/21.0*log(4.0*pimtot1by3)); //phi6
    phaseParams[8] = comprefac*pimtot1by3* (-6848.0)/63.0; //phi6l
    phaseParams[9] = comprefac*pimtot1by3*pimtot1by3* LAL_PI*(77096675.0/254016.0 + 378515.0/1512.0*eta - 74045.0/756.0*eta*eta); //phi7
    
    for(i=0;i<10;i++) {if (i!=1) {phaseParams[i]*=(1.0+dphis[i]);}}
    
//    fprintf(stdout,"%e %e %e %e %e %e %e %e\n", params->mass1, params->mass2, spin1, spin2, phaseParams[3], phaseParams[4], phaseParams[5], phaseParams[6]);
    
    return;
}
