/*
*  Copyright (C) 2013 Tjonnie G.F. Li, Michalis Agathos
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

/* <lalVerbatim file="LALInspiralSpinTaylorT4TestCV ">
Author: Li, T.G.F.
$Id: LALInspiralSpinTaylorT4Test.c ,v 1.0 2013/04/01 20:06:23 
</lalVerbatim>*/

/* <lalLaTeX>

\subsection{Test program \texttt{LALInspiralSpinTaylorT4Test.c }}
\label{ss:LALInspiralSpinTaylorT4Test.c }

Module to test the SpinTaylorT4 waveform in
./lalinspiral/src/LALInspiralSpinTaylorT4Test.c 

\subsection*{Usage}

\texttt{LALInspiralSpinTaylorT4Test [switch]}

</lalLaTeX> */

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>
#include <lal/SeqFactories.h>

NRCSID(LALINSPIRALSPINTAYLORT4TESTC , "$Id: LALInspiralSpinTaylorT4Test.c ,v 1.0 2013/04/01 20:06:23 ");

int main(void) {

	/* CREATING NECESSARY STRUCTURES TO CALL LALINSPIRALINTERFACESPINTAYLORT4 */
	static LALStatus mystatus;
	CoherentGW waveform;
	InspiralTemplate params;
	PPNParamStruc ppnParams;
	REAL8 dxis[10]={0.0};
	REAL8 cutoff = 40.;
	LALInspiralEOS eos = PP;
	
	REAL8 m1, m2, mtot, eta, dist, incl, phi_c, psi, lat, lon;
	REAL8 spin1x, spin1y, spin1z, spin2x, spin2y, spin2z;
    char fname[100];

    dist = 100.0; // in Mpc
    phi_c = 0.;
    psi = 0.;
	lat = 0.;
	lon = 0.;

    int i,j;
    UINT4 k;
    
    for (i=0;i<10;i++){
	  m1 = 1.2 + 0.15*i;
	  m2 = 1.0 + 0.1 *i;
	  for (j=0;j<10;j++){
		spin1x = .05 + 0.03*j;  
		spin1y = .1 * (1. - j/10.);  
		spin1z = .01*j;  
		spin2x = .05 + 0.03*i;   
		spin2y = .1 * (1. - i/10.);  
		spin2z = .01*i; 
		incl = LAL_PI*i*0.05;  
	mtot = m1 + m2;
	eta = m1*m2/(mtot*mtot);

	/* ENSURING THAT THE STRUCTURES ARE EMPTY */
	memset( &mystatus, 0, sizeof(LALStatus) );
	memset( &waveform, 0, sizeof(CoherentGW) );
	memset( &params, 0, sizeof(InspiralTemplate) );
	memset( &ppnParams, 0, sizeof(PPNParamStruc) );

    /* Filling in the PPNparam structure */
    ppnParams.deltaT = 1./4096.;
    ppnParams.lengthIn = 0;
    ppnParams.ppn = NULL;
    ppnParams.mTot = m1+m2;
    ppnParams.eta = eta;
    ppnParams.d = dist * 1.0e6 * LAL_PC_SI; // from Mpc to m
    ppnParams.inc = incl;
    ppnParams.phi = phi_c;
    ppnParams.ampOrder = 0;
    ppnParams.fStartIn = 20.;
    ppnParams.fStopIn = -1.0/(6.0*sqrt(6.0) * LAL_PI * mtot * LAL_MTSUN_SI);
    ppnParams.position.longitude = lon;
    ppnParams.position.latitude = lat;
    ppnParams.position.system = COORDINATESYSTEM_EQUATORIAL;
    ppnParams.psi = psi;
    ppnParams.epoch.gpsSeconds = 0;
    ppnParams.epoch.gpsNanoSeconds = 0;

	/* FILLING THE PARAMETER STRUCTURE */
	params.approximant = SpinTaylorT4;
	params.order = LAL_PNORDER_THREE_POINT_FIVE;
	params.ampOrder = LAL_PNORDER_NEWTONIAN;
	params.mass1 = m1;
	params.mass2 = m2;
	params.fCutoff = 1./(ppnParams.deltaT)/1. - 1.;
	params.fLower = 20.0;
	params.tSampling = 4096.;
	params.distance = dist * 1.e6 * LAL_PC_SI; // from Mpc to m (not s!)
	params.signalAmplitude = 1.0;
	params.startPhase = phi_c;
	params.startTime = 0.0;
	params.ieta = 1;
	params.inclination = incl;
	params.orbitTheta0 = 0.0;
	params.orbitPhi0 = 0.0;
	params.spin1[0] = spin1x;
	params.spin1[1] = spin1y;
	params.spin1[2] = spin1z;
	params.spin2[0] = spin2x;
	params.spin2[1] = spin2y;
	params.spin2[2] = spin2z;
	params.sourceTheta = LAL_PI/2. - lat;
	params.sourcePhi = lon;
	params.polarisationAngle = psi;
	params.qmParameter[0] = 0.0;
	params.qmParameter[1] = 0.0;
	params.spinInteraction = LAL_AllInter;

	LALInspiralInterfaceSpinTaylorT4(&mystatus, &waveform, &params, &ppnParams, dxis, cutoff, eos);

	/* 	PRINT THE WAVEFORM */

    sprintf(fname, "SpinTaylorT4_test_%d.dat", 10*i + j);
    FILE *outInj4=fopen(fname ,"w");
    for (k=0; k<waveform.h->data->length; k++) {
            fprintf(outInj4, "%lf %e %e\n", k*ppnParams.deltaT, waveform.h->data->data[2*k], waveform.h->data->data[2*k + 1]);
    }
    fclose(outInj4);
	printf("duration is %f \n", ppnParams.tc);
    // printf("final frequency is %e \n", );
	  }
	}

  /* Destroy data structures */
  LALSDestroyVectorSequence( &mystatus, &(waveform.h->data));
  //CHECKSTATUSPTR(&mystatus);
  /* LALFree(waveform.h);
  LALFree(wavefrom);
  LALFree(params);
  LALFree(ppnParams); */
  LALCheckMemoryLeaks();

  /* EXIT THE FUNCTION */

  return 0;
}
