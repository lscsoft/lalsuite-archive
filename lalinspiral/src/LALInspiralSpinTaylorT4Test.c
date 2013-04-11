/*
*  Copyright (C) 2013 Tjonnie G.F. Li
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

NRCSID(LALInspiralSpinTaylorT4TestC , "$Id: LALInspiralSpinTaylorT4Test.c ,v 1.0 2013/04/01 20:06:23 ");

int main(int argc,char UNUSED **argv) {

	/* CREATING NECESSARY STRUCTURES TO CALL LALINSPIRALINTERFACESPINTAYLORT4 */
	static LALStatus mystatus;
	CoherentGW waveform;
	InspiralTemplate params;
	PPNParamStruc ppnParams;
	REAL8 dxis[10]={0.0};
	REAL8 cutoff = 40.;
	LALInspiralEOS eos = PP;

	/* ENSURING THAT THE STRUCTURES ARE EMPTY */
	memset( &mystatus, 0, sizeof(LALStatus) );
	memset( &waveform, 0, sizeof(CoherentGW) );
	memset( &params, 0, sizeof(InspiralTemplate) );
	memset( &ppnParams, 0, sizeof(PPNParamStruc) );

	/* FILLING THE PARAMETER STRUCTURE */
	params.approximant = SpinTaylorT4;
	params.order = LAL_PNORDER_THREE_POINT_FIVE;
	params.ampOrder = LAL_PNORDER_NEWTONIAN;
	params.mass1 = 1.4;
	params.mass2 = 1.4;
	params.fCutoff = 0.0;
	params.fLower = 20.0;
	params.tSampling = 0.0;
	params.distance = 100.0;
	params.signalAmplitude = 0.0;
	params.startPhase = 0.0;
	params.startTime = 0.0;
	params.ieta = 0;
	params.inclination = 0.0;
	params.orbitTheta0 = 0.0;
	params.orbitPhi0 = 0.0;
	params.spin1[0] = 0.0;
	params.spin1[1] = 0.0;
	params.spin1[2] = 0.0;
	params.spin2[0] = 0.0;
	params.spin2[1] = 0.0;
	params.spin2[2] = 0.0;
	params.sourceTheta = 0.0;
	params.sourcePhi = 0.0;
	params.polarisationAngle = 0.0;
	params.qmParameter[0] = 0.0;
	params.qmParameter[1] = 0.0;
	params.spinInteraction = LAL_AllInter;

	/* WHY DO WE NEED TWO PARAMETER STRUCTURES? */
	/* FILL THE PPNPARAM STRUCTURE - TO BE DONE*/

	LALInspiralInterfaceSpinTaylorT4(&mystatus, &waveform, &params, &ppnParams, dxis, cutoff, eos);

	/* 	PRINT THE WAVEFORM */

	/* EXIT THE FUNCTION */

	return 0;
}
