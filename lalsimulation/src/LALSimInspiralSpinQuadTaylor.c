/**	@file   LALSimInspiralSpinQuadTaylor.c
 *	@author László Veréb
 *	@date   21 Feb 2012
 *	@brief  
 */

#include <math.h>
#include <string.h>

#include <lal/AVFactories.h>
#include <lal/LALConstants.h>
#include <lal/XLALError.h>

#include "LALSimInspiralSpinQuadTaylor.h"

#define UNUSED(expr) do { (void)(expr); } while (0)

enum {
	HP = 0, HC, WAVEFORM, X = 0, Y, Z, DIMENSIONS,
};

enum {
	PNALL = -1, PN0_0, PN0_5, PN1_0, PN1_5, PN2_0, PN2_5, PN3_0, PN3_5, PN4_0, PNORDER,
};

/**	@brief Structure containing the coefficients used in evolution equations
 */
typedef struct {
	REAL8 omega[PNORDER]; ///< for orbital frequency
	REAL8 omegaGlobal; ///< global coefficient for orbital frequency
	REAL8 MECO[PNORDER]; ///< for MECO
} Coefficient;

/**	@brief Structure containing initial and calculated parameters.
 * The masses are in (s)
 */
typedef struct {
	REAL8 mass[2]; ///< masses of the system
	REAL8 totalMass; ///< total mass of the system
	REAL8 eta; ///< symmetric mass ratio
	REAL8 chirpMass; ///< chirp mass
	REAL8 qm[2];
	REAL8 chih[2][DIMENSIONS]; ///< unit spin parameters
	REAL8 chiAmp[2]; ///< spin parameter magnitude
	REAL8 lnhat[3]; ///< unit orbital angular momentum
	REAL8 e1[3]; ///< orbital plane basis vector
	INT4 orderOfPhase; ///< twice phase post-Newtonian order
	LALSimInspiralInteraction interactionFlags; ///< flag to control spin effects
	Coefficient coeff; ///< calculated coefficients used in evolution equations
} Parameters;

/**	@brief Structure to produce less error prone code.
 */
typedef struct {
	REAL8TimeSeries *waveform[WAVEFORM]; ///< generated waveform
	REAL8TimeSeries *pnParameter; ///< evolved PN-parameter
	REAL8TimeSeries *phase; ///< evolved phase
	REAL8TimeSeries *chi1[DIMENSIONS]; ///< evolved first spin parameter
	REAL8TimeSeries *chi2[DIMENSIONS]; ///< evolved second spin parameter
	REAL8TimeSeries *e1[DIMENSIONS]; ///< evolved orbital plane basis vector
	REAL8TimeSeries *e3[DIMENSIONS]; ///< evolved unit orbital angular momentum
} Output;

/**	Calculates constants values from current initial values.
 *
 * @param[out] parameter		: the calculated and latter going to be used parameters
 * @param[in]  mass1			: mass of companion 1 (kg)
 * @param[in]  mass2			: mass of companion 2 (kg)
 * @param[in]  chi1x			: initial value of x comp. of the first spin parameter
 * @param[in]  chi1y			: initial value of y comp. of the first spin parameter
 * @param[in]  chi1z			: initial value of z comp. of the first spin parameter
 * @param[in]  chi2x			: initial value of x comp. of the second spin parameter
 * @param[in]  chi2y			: initial value of y comp. of the second spin parameter
 * @param[in]  chi2z			: initial value of z comp. of the second spin parameter
 * @param[in]  lnhatx			: initial value of x comp. of unit orbital angular momentum
 * @param[in]  lnhaty			: initial value of y comp. of unit orbital angular momentum
 * @param[in]  lnhatz			: initial value of z comp. of unit orbital angular momentum
 * @param[in]  e1x				: initial value of x comp. of orbital plane basis vector
 * @param[in]  e1y				: initial value of y comp. of orbital plane basis vector
 * @param[in]  e1z				: initial value of z comp. of orbital plane basis vector
 * @param[in]  endPhase			: orbital phase at last sample
 * @param[in]  initialFrequency	: initial frequency (1/s)
 * @param[in]  samplingTime		: sampling interval (s)
 * @param[in]  orderOfPhase		: twice phase post-Newtonian order
 * @param[in]  interactionFlags : flag to control spin effects
 * @return
 */
static int XLALCalculateConstantParameters(Parameters *param, REAL8 mass1, REAL8 mass2, REAL8 qm1,
		REAL8 qm2, REAL8 chi1x, REAL8 chi1y, REAL8 chi1z, REAL8 chi2x, REAL8 chi2y, REAL8 chi2z,
		REAL8 lnhatx, REAL8 lnhaty, REAL8 lnhatz, REAL8 e1x, REAL8 e1y, REAL8 e1z,
		INT4 orderOfPhase, LALSimInspiralInteraction interactionFlags);

/**	Evolves the parameters of the orbit.
 *
 * @param[out]    out			: evolved value
 * @param[in,out] length		: length of the evolution
 * @param[in]     initial		: initial values
 * @param[in]     samplingTime	: sampling time
 * @return
 */
static int XLALEvolveParameters(REAL8Array **out, REAL8 *length, Parameters *initial,
		REAL8 samplingTime);

/**	Allocates memoory for output structures of orbital parameters.
 *
 * @param[out] param		: structure containing pointers to the allocated memory
 * @param[in]  length		: length of allocated data
 * @param[in]  samplingTime	: sampling time of the data
 * @return
 */
static int XLALCreateOrbitOutput(Output *param, UINT4 length, REAL8 samplingTime);

/**	Allocates memoory for output structures of waveforms.
 *
 * @param[out] param		: structure containing pointers to the allocated memory
 * @param[in]  length		: length of allocated data
 * @param[in]  samplingTime	: sampling time of the data
 * @return
 */
static int XLALCreateWaveformOutput(Output *param, UINT4 length, REAL8 samplingTime);

/**	Calculates the current value of the waveforms.
 *
 * @param[out] hp				: \f$h_+\f$ polarised waveform
 * @param[out] hc				: \f$h_\times\f$ polarised waveform
 * @param[in]  e1				: orbital plane basis vector
 * @param[in]  e3				: unit orbital angular momentum
 * @param[in]  phase			: phase of the waveform
 * @param[in]  V				: PN-parameter
 * @param[in]  params			: parameters
 * @param[in]  orderOfAmplitude	: twice amplitude post-Newtonian order
 * @return
 */
static int XLALCalculateWaveform(REAL8 *hp, REAL8 *hc, REAL8 e1[], REAL8 e3[], REAL8 phase, REAL8 V,
		Parameters *params, INT4 orderOfAmplitude);

int XLALSimInspiralSpinQuadTaylorEvolveWaveform(REAL8TimeSeries **hp, REAL8TimeSeries **hc,
		REAL8 mass1, REAL8 mass2, REAL8 qm1, REAL8 qm2, REAL8 chi1x, REAL8 chi1y, REAL8 chi1z,
		REAL8 chi2x, REAL8 chi2y, REAL8 chi2z, REAL8 lnhatx, REAL8 lnhaty, REAL8 lnhatz, REAL8 e1x,
		REAL8 e1y, REAL8 e1z, REAL8 endPhase, REAL8 initialFrequency, REAL8 samplingTime,
		INT4 orderOfPhase, INT4 orderOfAmplitude, LALSimInspiralInteraction interactionFlags) {
	Parameters parameter;
	memset(&parameter, 0, sizeof(Parameters));
	XLALCalculateConstantParameters(&parameter, mass1, mass2, qm1, qm2, chi1x, chi1y, chi1z, chi2x,
			chi2y, chi2z, lnhatx, lnhaty, lnhatz, e1x, e1y, e1z, orderOfPhase, interactionFlags);
	REAL8 length = (5.0 / 256.0) * pow(LAL_PI, -8.0 / 3.0)
			* pow(parameter.chirpMass * initialFrequency, -5.0 / 3.0) / initialFrequency;
	REAL8Array *out = NULL;
	XLALEvolveParameters(&out, &length, &parameter, samplingTime);
	Output evolved;
	memset(&evolved, 0, sizeof(Output));
	UINT4 size = length;
	XLALCreateWaveformOutput(&evolved, size, samplingTime);
	*hp = evolved.waveform[HP];
	*hc = evolved.waveform[HC];
	REAL8 phase, v;
	REAL8 e1[DIMENSIONS];
	REAL8 e3[DIMENSIONS];
	REAL8 phiShift = endPhase - out->data[2 * size - 1];
	for (UINT4 i = 0; i < size; i++) {
		phase = out->data[size + i] + phiShift;
		v = cbrt(out->data[2 * size + i]);
		e3[X] = out->data[3 * size + i];
		e3[Y] = out->data[4 * size + i];
		e3[Z] = out->data[5 * size + i];
		e1[X] = out->data[12 * size + i];
		e1[Y] = out->data[13 * size + i];
		e1[Z] = out->data[14 * size + i];
		XLALCalculateWaveform(&((*hp)->data->data[i]), &((*hc)->data->data[i]), e1, e3, phase, v,
				&parameter, orderOfAmplitude);
	}
	XLALDestroyREAL8Array(out);
	return XLAL_SUCCESS;
}

int XLALSimInspiralSpinQuadTaylorEvolveOrbit(REAL8TimeSeries **V, REAL8TimeSeries **Phi,
		REAL8TimeSeries **S1x, REAL8TimeSeries **S1y, REAL8TimeSeries **S1z, REAL8TimeSeries **S2x,
		REAL8TimeSeries **S2y, REAL8TimeSeries **S2z, REAL8TimeSeries **LNhatx,
		REAL8TimeSeries **LNhaty, REAL8TimeSeries **LNhatz, REAL8TimeSeries **E1x,
		REAL8TimeSeries **E1y, REAL8TimeSeries **E1z, REAL8 mass1, REAL8 mass2, REAL8 qm1,
		REAL8 qm2, REAL8 chi1x, REAL8 chi1y, REAL8 chi1z, REAL8 chi2x, REAL8 chi2y, REAL8 chi2z,
		REAL8 lnhatx, REAL8 lnhaty, REAL8 lnhatz, REAL8 e1x, REAL8 e1y, REAL8 e1z, REAL8 endPhase,
		REAL8 initialFrequency, REAL8 samplingTime, INT4 orderOfPhase,
		LALSimInspiralInteraction interactionFlags) {
	Parameters parameter;
	memset(&parameter, 0, sizeof(Parameters));
	XLALCalculateConstantParameters(&parameter, mass1, mass2, qm1, qm2, chi1x, chi1y, chi1z, chi2x,
			chi2y, chi2z, lnhatx, lnhaty, lnhatz, e1x, e1y, e1z, orderOfPhase, interactionFlags);
	REAL8 length = (5.0 / 256.0) * pow(LAL_PI, -8.0 / 3.0)
			* pow(parameter.chirpMass * initialFrequency, -5.0 / 3.0) / initialFrequency;
	REAL8Array *out = NULL;
	XLALEvolveParameters(&out, &length, &parameter, samplingTime);
	Output evolved;
	memset(&evolved, 0, sizeof(Output));
	UINT4 size = length;
	XLALCreateWaveformOutput(&evolved, size, samplingTime);
	XLALCreateOrbitOutput(&evolved, size, samplingTime);
	*V = evolved.pnParameter;
	*Phi = evolved.phase;
	*S1x = evolved.chi1[X];
	*S1y = evolved.chi1[Y];
	*S1z = evolved.chi1[Z];
	*S2x = evolved.chi2[X];
	*S2y = evolved.chi2[Y];
	*S2z = evolved.chi2[Z];
	*LNhatx = evolved.e3[X];
	*LNhaty = evolved.e3[Y];
	*LNhatz = evolved.e3[Z];
	*E1x = evolved.e1[X];
	*E1y = evolved.e1[Y];
	*E1z = evolved.e1[Z];
	REAL8 phase, v;
	REAL8 e1[DIMENSIONS];
	REAL8 e3[DIMENSIONS];
	REAL8 phiShift = endPhase - out->data[2 * size - 1];
	for (UINT4 i = 0; i < size; i++) {
		(*Phi)->data->data[i] = phase = out->data[size + i] + phiShift;
		(*V)->data->data[i] = v = cbrt(out->data[2 * size + i]);
		(*LNhatx)->data->data[i] = e3[X] = out->data[3 * size + i];
		(*LNhaty)->data->data[i] = e3[Y] = out->data[4 * size + i];
		(*LNhatz)->data->data[i] = e3[Z] = out->data[5 * size + i];
		(*S1x)->data->data[i] = out->data[6 * size + i];
		(*S1y)->data->data[i] = out->data[7 * size + i];
		(*S1z)->data->data[i] = out->data[8 * size + i];
		(*S2x)->data->data[i] = out->data[9 * size + i];
		(*S2y)->data->data[i] = out->data[10 * size + i];
		(*S2z)->data->data[i] = out->data[11 * size + i];
		(*E1x)->data->data[i] = e1[X] = out->data[12 * size + i];
		(*E1y)->data->data[i] = e1[Y] = out->data[13 * size + i];
		(*E1z)->data->data[i] = e1[Z] = out->data[14 * size + i];
	}
	XLALDestroyREAL8Array(out);
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
	Parameters parameter;
	memset(&parameter, 0, sizeof(Parameters));
	XLALCalculateConstantParameters(&parameter, mass1, mass2, qm1, qm2, chi1x, chi1y, chi1z, chi2x,
			chi2y, chi2z, lnhatx, lnhaty, lnhatz, e1x, e1y, e1z, orderOfPhase, interactionFlags);
	REAL8 length = (5.0 / 256.0) * pow(LAL_PI, -8.0 / 3.0)
			* pow(parameter.chirpMass * initialFrequency, -5.0 / 3.0) / initialFrequency;
	REAL8Array *out = NULL;
	XLALEvolveParameters(&out, &length, &parameter, samplingTime);
	Output evolved;
	memset(&evolved, 0, sizeof(Output));
	UINT4 size = length;
	XLALCreateWaveformOutput(&evolved, size, samplingTime);
	XLALCreateOrbitOutput(&evolved, size, samplingTime);
	*hp = evolved.waveform[HP];
	*hc = evolved.waveform[HC];
	*V = evolved.pnParameter;
	*Phi = evolved.phase;
	*S1x = evolved.chi1[X];
	*S1y = evolved.chi1[Y];
	*S1z = evolved.chi1[Z];
	*S2x = evolved.chi2[X];
	*S2y = evolved.chi2[Y];
	*S2z = evolved.chi2[Z];
	*LNhatx = evolved.e3[X];
	*LNhaty = evolved.e3[Y];
	*LNhatz = evolved.e3[Z];
	*E1x = evolved.e1[X];
	*E1y = evolved.e1[Y];
	*E1z = evolved.e1[Z];
	REAL8 phase, v;
	REAL8 e1[DIMENSIONS];
	REAL8 e3[DIMENSIONS];
	REAL8 phiShift = endPhase - out->data[2 * size - 1];
	for (UINT4 i = 0; i < size; i++) {
		(*Phi)->data->data[i] = phase = out->data[size + i] + phiShift;
		(*V)->data->data[i] = v = cbrt(out->data[2 * size + i]);
		(*LNhatx)->data->data[i] = e3[X] = out->data[3 * size + i];
		(*LNhaty)->data->data[i] = e3[Y] = out->data[4 * size + i];
		(*LNhatz)->data->data[i] = e3[Z] = out->data[5 * size + i];
		(*S1x)->data->data[i] = out->data[6 * size + i];
		(*S1y)->data->data[i] = out->data[7 * size + i];
		(*S1z)->data->data[i] = out->data[8 * size + i];
		(*S2x)->data->data[i] = out->data[9 * size + i];
		(*S2y)->data->data[i] = out->data[10 * size + i];
		(*S2z)->data->data[i] = out->data[11 * size + i];
		(*E1x)->data->data[i] = e1[X] = out->data[12 * size + i];
		(*E1y)->data->data[i] = e1[Y] = out->data[13 * size + i];
		(*E1z)->data->data[i] = e1[Z] = out->data[14 * size + i];
		XLALCalculateWaveform(&((*hp)->data->data[i]), &((*hc)->data->data[i]), e1, e3, phase, v,
				&parameter, orderOfAmplitude);
	}
	XLALDestroyREAL8Array(out);
	return XLAL_SUCCESS;
}

static int XLALCalculateConstantParameters(Parameters *param, REAL8 mass1, REAL8 mass2, REAL8 qm1,
		REAL8 qm2, REAL8 chi1x, REAL8 chi1y, REAL8 chi1z, REAL8 chi2x, REAL8 chi2y, REAL8 chi2z,
		REAL8 lnhatx, REAL8 lnhaty, REAL8 lnhatz, REAL8 e1x, REAL8 e1y, REAL8 e1z,
		INT4 orderOfPhase, LALSimInspiralInteraction interactionFlags) {
	UNUSED(param);
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
	UNUSED(orderOfPhase);
	UNUSED(interactionFlags);
	return XLAL_SUCCESS;
}

static int XLALEvolveParameters(REAL8Array **out, REAL8 *length, Parameters *initial,
		REAL8 samplingTime) {
	UNUSED(out);
	UNUSED(initial);
	UNUSED(length);
	UNUSED(samplingTime);
	// initialize integrator
	// run integrator
	// clean integrator
	return XLAL_SUCCESS;
}

static int XLALCreateOrbitOutput(Output *param, UINT4 length, REAL8 samplingTime) {
	UNUSED(param);
	UNUSED(length);
	UNUSED(samplingTime);
	return XLAL_SUCCESS;
}

static int XLALCreateWaveformOutput(Output *param, UINT4 length, REAL8 samplingTime) {
	UNUSED(param);
	UNUSED(length);
	UNUSED(samplingTime);
	return XLAL_SUCCESS;
}

static int XLALCalculateWaveform(REAL8 *hp, REAL8 *hc, REAL8 e1[], REAL8 e3[], REAL8 phase, REAL8 V,
		Parameters *params, INT4 orderOfAmplitude) {
	UNUSED(hp);
	UNUSED(hc);
	UNUSED(e1);
	UNUSED(e3);
	UNUSED(phase);
	UNUSED(V);
	UNUSED(params);
	UNUSED(orderOfAmplitude);
	return XLAL_SUCCESS;
}
