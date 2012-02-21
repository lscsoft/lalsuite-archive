/**	@file   LALSimInspiralSpinQuadTaylor.c
 *	@author László Veréb
 *	@date   21 Feb 2012
 *	@brief  
 */

#include <lal/XLALError.h>
#include "LALSimInspiralSpinQuadTaylor.h"

#define UNUSED(expr) do { (void)(expr); } while (0)

int XLALSimInspiralSpinQuadTaylorEvolveWaveform(REAL8TimeSeries **hp, REAL8TimeSeries **hc,
		REAL8 mass1, REAL8 mass2, REAL8 qm1, REAL8 qm2, REAL8 chi1x, REAL8 chi1y, REAL8 chi1z,
		REAL8 chi2x, REAL8 chi2y, REAL8 chi2z, REAL8 lnhatx, REAL8 lnhaty, REAL8 lnhatz, REAL8 e1x,
		REAL8 e1y, REAL8 e1z, REAL8 endPhase, REAL8 initialFrequency, REAL8 samplingTime,
		INT4 orderOfPhase, INT4 orderOfAmplitude, LALSimInspiralInteraction interactionFlags) {
	UNUSED(hp);
	UNUSED(hc);
	UNUSED(mass1);
	UNUSED(mass2);
	UNUSED(qm1);
	UNUSED(qm2);
	UNUSED(chi1x);
	UNUSED(chi1y);
	UNUSED(chi1z);
	UNUSED(chi2x);
	UNUSED(chi2y);
	UNUSED(chi2z);
	UNUSED(lnhatx);
	UNUSED(lnhaty);
	UNUSED(lnhatz);
	UNUSED(e1x);
	UNUSED(e1y);
	UNUSED(e1z);
	UNUSED(endPhase);
	UNUSED(initialFrequency);
	UNUSED(samplingTime);
	UNUSED(orderOfPhase);
	UNUSED(orderOfAmplitude);
	UNUSED(interactionFlags);
	return XLAL_SUCCESS;
}

int XLALSimInspiralSpinQuadTaylorEvolveOrbit(REAL8TimeSeries **V, REAL8TimeSeries **Phi,
		REAL8TimeSeries **S1x, REAL8TimeSeries **S1y, REAL8TimeSeries **S1z, REAL8TimeSeries **S2x,
		REAL8TimeSeries **S2y, REAL8TimeSeries **S2z, REAL8TimeSeries **LNhatx,
		REAL8TimeSeries **LNhaty, REAL8TimeSeries **LNhatz, REAL8TimeSeries **E1x,
		REAL8TimeSeries **E1y, REAL8TimeSeries **E1z, REAL8 mass1, REAL8 mass2, REAL8 qm1,
		REAL8 qm2, REAL8 chi1x, REAL8 chi1y, REAL8 chi1z, REAL8 chi2x, REAL8 chi2y, REAL8 chi2z,
		REAL8 lnhatx, REAL8 lnhaty, REAL8 lnhatz, REAL8 e1x, REAL8 e1y, REAL8 e1z, REAL8 endPhase,
		REAL8 initialFrequency, REAL8 samplingTime, INT4 orderOfPhase, INT4 orderOfAmplitude,
		LALSimInspiralInteraction interactionFlags) {
	UNUSED(V);
	UNUSED(Phi);
	UNUSED(S1x);
	UNUSED(S1y);
	UNUSED(S1z);
	UNUSED(S2x);
	UNUSED(S2y);
	UNUSED(S2z);
	UNUSED(LNhatx);
	UNUSED(LNhaty);
	UNUSED(LNhatz);
	UNUSED(E1x);
	UNUSED(E1y);
	UNUSED(E1z);
	UNUSED(mass1);
	UNUSED(mass2);
	UNUSED(qm1);
	UNUSED(qm2);
	UNUSED(chi1x);
	UNUSED(chi1y);
	UNUSED(chi1z);
	UNUSED(chi2x);
	UNUSED(chi2y);
	UNUSED(chi2z);
	UNUSED(lnhatx);
	UNUSED(lnhaty);
	UNUSED(lnhatz);
	UNUSED(e1x);
	UNUSED(e1y);
	UNUSED(e1z);
	UNUSED(endPhase);
	UNUSED(initialFrequency);
	UNUSED(samplingTime);
	UNUSED(orderOfPhase);
	UNUSED(orderOfAmplitude);
	UNUSED(interactionFlags);
	return XLAL_SUCCESS;
}

int XLALSimInspiralSpinQuadTaylorEvolveAll(REAL8TimeSeries **hp, REAL8TimeSeries **hc,
		REAL8TimeSeries **V, REAL8TimeSeries **Phi, REAL8TimeSeries **S1x, REAL8TimeSeries **S1y,
		REAL8TimeSeries **S1z, REAL8TimeSeries **S2x, REAL8TimeSeries **S2y, REAL8TimeSeries **S2z,
		REAL8TimeSeries **LNhatx, REAL8TimeSeries **LNhaty, REAL8TimeSeries **LNhatz,
		REAL8TimeSeries **E1x, REAL8TimeSeries **E1y, REAL8TimeSeries **E1z, REAL8 mass1,
		REAL8 mass2, REAL8 qm1, REAL8 qm2, REAL8 chi1x, REAL8 chi1y, REAL8 chi1z, REAL8 chi2x,
		REAL8 chi2y, REAL8 chi2z, REAL8 lnhatx, REAL8 lnhaty, REAL8 lnhatz, REAL8 e1x, REAL8 e1y,
		REAL8 e1z, REAL8 endPhase, REAL8 initialFrequency, REAL8 samplingTime, INT4 orderOfPhase,
		INT4 orderOfAmplitude, LALSimInspiralInteraction interactionFlags) {
	UNUSED(hp);
	UNUSED(hc);
	UNUSED(V);
	UNUSED(Phi);
	UNUSED(S1x);
	UNUSED(S1y);
	UNUSED(S1z);
	UNUSED(S2x);
	UNUSED(S2y);
	UNUSED(S2z);
	UNUSED(LNhatx);
	UNUSED(LNhaty);
	UNUSED(LNhatz);
	UNUSED(E1x);
	UNUSED(E1y);
	UNUSED(E1z);
	UNUSED(mass1);
	UNUSED(mass2);
	UNUSED(qm1);
	UNUSED(qm2);
	UNUSED(chi1x);
	UNUSED(chi1y);
	UNUSED(chi1z);
	UNUSED(chi2x);
	UNUSED(chi2y);
	UNUSED(chi2z);
	UNUSED(lnhatx);
	UNUSED(lnhaty);
	UNUSED(lnhatz);
	UNUSED(e1x);
	UNUSED(e1y);
	UNUSED(e1z);
	UNUSED(endPhase);
	UNUSED(initialFrequency);
	UNUSED(samplingTime);
	UNUSED(orderOfPhase);
	UNUSED(orderOfAmplitude);
	UNUSED(interactionFlags);
	return XLAL_SUCCESS;
}
