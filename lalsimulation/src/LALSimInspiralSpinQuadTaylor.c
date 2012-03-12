/**	@file   LALSimInspiralSpinQuadTaylor.c
 *	@author László Veréb
 *	@date   21 Feb 2012
 *	@brief  
 */

#include <math.h>

#include <lal/LALAdaptiveRungeKutta4.h>
#include <lal/LALConstants.h>
#include <lal/LALSimInspiral.h>
#include <lal/TimeSeries.h>
#include <lal/Units.h>

#include "check_series_macros.h"

#define UNUSED(expr) do { (void)(expr); } while (0)

#define ABSOLUTE_TOLERANCE 1.e-12
#define RELATIVE_TOLERANCE 1.e-12

/**	Error codes for the integrator.
 */
enum {
	TEST_ENERGY = 1025, TEST_OMEGADOT, TEST_COORDINATE, TEST_OMEGANAN, TEST_NYQUIST,
};

/**	Enumeration of the evolving quantities.
 */
enum {
	PHASE, ///< PHASE
	OMEGA, ///< OMEGA
	LNHATX, ///< LNHATX
	LNHATY, ///< LNHATY
	LNHATZ, ///< LNHATZ
	CHI1X, ///< CHI1X
	CHI1Y, ///< CHI1Y
	CHI1Z, ///< CHI1Z
	CHI2X, ///< CHI2X
	CHI2Y, ///< CHI2Y
	CHI2Z, ///< CHI2Z
	E1X, ///< E1X
	E1Y, ///< E1Y
	E1Z, ///< E1Z
	EVOLVING_VARIABLES,
///< EVOLVING_VARIABLES
};

/**	Various constants.
 */
enum {
	HP = 0, HC, WAVEFORM, X = 0, Y, Z, DIMENSIONS,
};

/**	Enumeration of the post-newtonian orders.
 */
enum {
	PNALL = -1, PN0_0, PN0_5, PN1_0, PN1_5, PN2_0, PN2_5, PN3_0, PN3_5, PN4_0, PNORDER,
};

/**	@brief Structure containing the coefficients used in evolution equations
 */
typedef struct {
	REAL8 omega[PNORDER]; ///< for orbital frequency
	REAL8 omegaSO[2]; ///< orbital frequency spin-orbit contribution
	REAL8 omegaSS[2]; ///< orbital frequency spin-spin contribution
	REAL8 omegaQM[2]; ///< orbital frequency quad-mono contribution
	REAL8 omegaSELF[2]; ///< orbital frequency spin-self contribution
	REAL8 omegaGlobal; ///< global coefficient for orbital frequency
	REAL8 chihSO[2]; ///< spin parameter spin-orbit contribution
	REAL8 chihSS[2]; ///< spin parameter spin-spin contribution
	REAL8 chihQM[2]; ///< spin parameter quad-mono contribution
	REAL8 lnhat[2]; ///< unit orbital angular momentum
	REAL8 meco[PNORDER]; ///< for MECO
	REAL8 mecoSO[2]; ///< MECO spin-orbit contribution
	REAL8 mecoSS; ///< MECO spin-spin contribution
	REAL8 mecoQM; ///< MECO quad-mono contribution
	REAL8 omegaPowi_3[PNORDER]; ///< current powers of the PN-parameter \f$v^{(i/3)}\f$
	REAL8 lnhchih[2];
	REAL8 chih1chih2;
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
	REAL8 samplingTime; ///< sampling time
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
 * @param[in]  distance			: distance of source (m)
 * @param[in]  samplingTime		: sampling interval (s)
 * @param[in]  orderOfPhase		: twice phase post-Newtonian order
 * @param[in]  interactionFlags : flag to control spin effects
 * @return
 */
static int XLALCalculateConstantParameters(Parameters *param, REAL8 mass1, REAL8 mass2, REAL8 qm1,
		REAL8 qm2, REAL8 chi1x, REAL8 chi1y, REAL8 chi1z, REAL8 chi2x, REAL8 chi2y, REAL8 chi2z,
		REAL8 lnhatx, REAL8 lnhaty, REAL8 lnhatz, REAL8 e1x, REAL8 e1y, REAL8 e1z,
		REAL8 samplingTime, INT4 orderOfPhase, LALSimInspiralInteraction interactionFlags);

/**	Evolves the parameters of the orbit.
 *
 * @param[out]    out				: evolved value
 * @param[in,out] length			: length of the evolution
 * @param[in]     param				: initial values
 * @param[in]     initialFrequency	: initial frequency
 * @return
 */
static int XLALEvolveParameters(REAL8Array **out, REAL8 *length, Parameters *param,
		REAL8 initialFrequency);

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
		Parameters *params, REAL8 distance, INT4 orderOfAmplitude);

/**	Computes the derivatives of the dynamical variables at a given time.
 *
 * @param[in]  t		: time at which the derivatives are calculated
 * @param[in]  values	:
 * @param[out] dvalues	: computed derivatives
 * @param[in]  mparams	: used parameters
 * @return
 */
static int XLALDerivator(double t, const double values[], double dvalues[], void *mparams);

static int XLALStop(double t, const double values[], double dvalues[], void *mparams);

int XLALSimInspiralSpinQuadTaylorEvolveWaveform(REAL8TimeSeries **hp, REAL8TimeSeries **hc,
		REAL8 mass1, REAL8 mass2, REAL8 qm1, REAL8 qm2, REAL8 chi1x, REAL8 chi1y, REAL8 chi1z,
		REAL8 chi2x, REAL8 chi2y, REAL8 chi2z, REAL8 lnhatx, REAL8 lnhaty, REAL8 lnhatz, REAL8 e1x,
		REAL8 e1y, REAL8 e1z, REAL8 distance, REAL8 endPhase, REAL8 initialFrequency,
		REAL8 samplingTime, INT4 orderOfPhase, INT4 orderOfAmplitude,
		LALSimInspiralInteraction interactionFlags) {
	Parameters parameter;
	memset(&parameter, 0, sizeof(Parameters));
	int error;
	error = XLALCalculateConstantParameters(&parameter, mass1, mass2, qm1, qm2, chi1x, chi1y, chi1z,
			chi2x, chi2y, chi2z, lnhatx, lnhaty, lnhatz, e1x, e1y, e1z, samplingTime, orderOfPhase,
			interactionFlags);
	if (error == XLAL_FAILURE) {
		XLAL_ERROR(XLAL_EFUNC);
	}
	REAL8 length = (5.0 / 256.0) * pow(LAL_PI, -8.0 / 3.0)
			* pow(parameter.chirpMass * initialFrequency, -5.0 / 3.0) / initialFrequency;
	REAL8Array *out = NULL;
	error = XLALEvolveParameters(&out, &length, &parameter, initialFrequency);
	if (error == XLAL_FAILURE) {
		XLAL_ERROR(XLAL_EFUNC);
	}
	Output evolved;
	memset(&evolved, 0, sizeof(Output));
	UINT4 size = length;
	error = XLALCreateWaveformOutput(&evolved, size, samplingTime);
	if (error == XLAL_FAILURE) {
		XLAL_ERROR(XLAL_EFUNC);
	}
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
		error = XLALCalculateWaveform(&((*hp)->data->data[i]), &((*hc)->data->data[i]), e1, e3,
				phase, v, &parameter, distance, orderOfAmplitude);
		if (error == XLAL_FAILURE) {
			XLAL_ERROR(XLAL_EFUNC);
		}
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
	int error;
	error = XLALCalculateConstantParameters(&parameter, mass1, mass2, qm1, qm2, chi1x, chi1y, chi1z,
			chi2x, chi2y, chi2z, lnhatx, lnhaty, lnhatz, e1x, e1y, e1z, samplingTime, orderOfPhase,
			interactionFlags);
	if (error == XLAL_FAILURE) {
		XLAL_ERROR(XLAL_EFUNC);
	}
	REAL8 length = (5.0 / 256.0) * pow(LAL_PI, -8.0 / 3.0)
			* pow(parameter.chirpMass * initialFrequency, -5.0 / 3.0) / initialFrequency;
	REAL8Array *out = NULL;
	error = XLALEvolveParameters(&out, &length, &parameter, initialFrequency);
	if (error == XLAL_FAILURE) {
		XLAL_ERROR(XLAL_EFUNC);
	}
	Output evolved;
	memset(&evolved, 0, sizeof(Output));
	UINT4 size = length;
	error = XLALCreateOrbitOutput(&evolved, size, samplingTime);
	if (error == XLAL_FAILURE) {
		XLAL_ERROR(XLAL_EFUNC);
	}
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
		REAL8 e1z, REAL8 distance, REAL8 endPhase, REAL8 initialFrequency, REAL8 samplingTime,
		INT4 orderOfPhase, INT4 orderOfAmplitude, LALSimInspiralInteraction interactionFlags) {
	Parameters parameter;
	memset(&parameter, 0, sizeof(Parameters));
	int error;
	error = XLALCalculateConstantParameters(&parameter, mass1, mass2, qm1, qm2, chi1x, chi1y, chi1z,
			chi2x, chi2y, chi2z, lnhatx, lnhaty, lnhatz, e1x, e1y, e1z, samplingTime, orderOfPhase,
			interactionFlags);
	if (error == XLAL_FAILURE) {
		XLAL_ERROR(XLAL_EFUNC);
	}
	REAL8 length = (5.0 / 256.0) * pow(LAL_PI, -8.0 / 3.0)
			* pow(parameter.chirpMass * initialFrequency, -5.0 / 3.0) / initialFrequency;
	REAL8Array *out = NULL;
	error = XLALEvolveParameters(&out, &length, &parameter, initialFrequency);
	if (error == XLAL_FAILURE) {
		XLAL_ERROR(XLAL_EFUNC);
	}
	Output evolved;
	memset(&evolved, 0, sizeof(Output));
	UINT4 size = length;
	error = XLALCreateWaveformOutput(&evolved, size, samplingTime);
	if (error == XLAL_FAILURE) {
		XLAL_ERROR(XLAL_EFUNC);
	}
	error = XLALCreateOrbitOutput(&evolved, size, samplingTime);
	if (error == XLAL_FAILURE) {
		XLAL_ERROR(XLAL_EFUNC);
	}
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
		error = XLALCalculateWaveform(&((*hp)->data->data[i]), &((*hc)->data->data[i]), e1, e3,
				phase, v, &parameter, distance, orderOfAmplitude);
		if (error == XLAL_FAILURE) {
			XLAL_ERROR(XLAL_EFUNC);
		}
	}
	XLALDestroyREAL8Array(out);
	return XLAL_SUCCESS;
}

inline static REAL8 sq(REAL8 value) {
	return (value * value);
}

inline static int isSet(LALSimInspiralInteraction flag, LALSimInspiralInteraction bit) {
	int set = 0;
	if ((flag & bit) == bit) {
		set = 1;
	}
	return (set);
}

static int XLALCalculateConstantParameters(Parameters *param, REAL8 mass1, REAL8 mass2, REAL8 qm1,
		REAL8 qm2, REAL8 chi1x, REAL8 chi1y, REAL8 chi1z, REAL8 chi2x, REAL8 chi2y, REAL8 chi2z,
		REAL8 lnhatx, REAL8 lnhaty, REAL8 lnhatz, REAL8 e1x, REAL8 e1y, REAL8 e1z,
		REAL8 samplingTime, INT4 orderOfPhase, LALSimInspiralInteraction interactionFlags) {
	memset(param, 0, sizeof(Parameters));
	param->mass[0] = mass1;
	param->mass[1] = mass2;
	param->totalMass = param->mass[0] + param->mass[1];
	param->eta = param->mass[0] * param->mass[1] / sq(param->totalMass);
	param->chirpMass = param->totalMass * pow(param->eta, 3.0 / 5.0);
	param->qm[0] = qm1;
	param->qm[1] = qm2;
	REAL8 chi[2][DIMENSIONS] = { { chi1x, chi1y, chi1z }, { chi2x, chi2y, chi2z } };
	for (UINT2 current = 0; current < 2; current++) {
		for (UINT2 dim = X; dim < DIMENSIONS; dim++) {
			param->chiAmp[current] += sq(chi[current][dim]);
		}
		if (param->chiAmp[current]) {
			param->chiAmp[current] = sqrt(param->chiAmp[current]);
			for (UINT2 dim = X; dim < DIMENSIONS; dim++) {
				param->chih[current][dim] = chi[current][dim] / param->chiAmp[current];
			}
		}
	}
	param->lnhat[X] = lnhatx;
	param->lnhat[Y] = lnhaty;
	param->lnhat[Z] = lnhatz;
	param->e1[X] = e1x;
	param->e1[Y] = e1y;
	param->e1[Z] = e1z;
	param->samplingTime = samplingTime;
	param->orderOfPhase = orderOfPhase;
	param->interactionFlags = interactionFlags;
	REAL8 piPow2 = sq(LAL_PI);
	REAL8 etaPow2 = sq(param->eta);
	REAL8 etaPow3 = param->eta * etaPow2;
	REAL8 mj_mi[2] = { param->mass[1] / param->mass[0], param->mass[0] / param->mass[1] };
	REAL8 spin_MPow2[2];
	for (UINT2 i = 0; i < 2; i++) {
		spin_MPow2[i] = param->chiAmp[i] * sq(param->mass[i]) / sq(param->totalMass);
	}
	param->coeff.omegaGlobal = param->eta * 96.0 / 5.0;
	switch (param->orderOfPhase) {
	case PNALL:
	case PN4_0:
	case PN3_5:
		param->coeff.omega[PN3_5] = (-13245.0 + 717350.0 * param->eta + 731960.0 * etaPow2) * LAL_PI
				/ 12096.0;
	case PN3_0:
		param->coeff.omega[PN3_0] = 16447322263.0 / 139708800.0 - LAL_GAMMA * 1712.0 / 105.0
				+ piPow2 * 16.0 / 3.0 - log(16.0) * 856.0 / 105.0
				+ (-56198689.0 / 217728.0 + piPow2 * 451.0 / 48.0) * param->eta
				+ etaPow2 * 541.0 / 896.0 - etaPow3 * 5605.0 / 2592.0;
		param->coeff.meco[PN3_0] = (54675.0 * param->eta + (11070 * LAL_PI - 310005) * etaPow2
				+ 8370 * etaPow3 + 35.0 * etaPow3 * param->eta) / 3888.0;
	case PN2_5:
		param->coeff.omega[PN2_5] = -LAL_PI * (4159.0 + 15876.0 * param->eta) / 672.0;
	case PN2_0:
		param->coeff.omega[PN2_0] = (34103.0 + 122949.0 * param->eta + 59472.0 * etaPow2) / 18144.0;
		param->coeff.meco[PN2_0] = (81.0 * param->eta - 57.0 * etaPow2 + etaPow3) / 24.0;
		if (isSet(param->interactionFlags, LAL_SIM_INSPIRAL_INTERACTION_SPIN_SPIN_2PN)) {
			param->coeff.chihSS[0] = spin_MPow2[1] / 2.0;
			param->coeff.chihSS[1] = spin_MPow2[0] / 2.0;
			param->coeff.omegaSS[0] = 721.0 * param->eta * param->chiAmp[0] * param->chiAmp[1]
					/ 48.0;
			param->coeff.omegaSS[1] = -247.0 * param->coeff.omegaSS[0] / 721.0;
			param->coeff.mecoSS = -spin_MPow2[0] * spin_MPow2[1];
		}
		if (isSet(param->interactionFlags, LAL_SIM_INSPIRAL_INTERACTION_QUAD_MONO_2PN)) {
			for (UINT2 i = 0; i < 2; i++) {
				param->coeff.omegaQM[i] = spin_MPow2[i] * param->chiAmp[i] * param->qm[i] * 2.5;
				param->coeff.chihQM[i] = -param->qm[i] * param->eta * param->chiAmp[i] * 1.5;
			}
			param->coeff.mecoQM = param->eta / 5.0;
		}
		if (isSet(param->interactionFlags, LAL_SIM_INSPIRAL_INTERACTION_SPIN_SPIN_SELF_2PN)) {
			for (UINT2 i = 0; i < 2; i++) {
				param->coeff.omegaSELF[i] = -spin_MPow2[i] * param->chiAmp[i] / 96.0;
			}
		}
	case PN1_5:
		param->coeff.omega[PN1_5] = 4.0 * LAL_PI;
		if (isSet(param->interactionFlags, LAL_SIM_INSPIRAL_INTERACTION_SPIN_ORBIT_15PN)) {
			for (UINT2 i = 0; i < 2; i++) {
				param->coeff.omegaSO[i] = -spin_MPow2[i] * (113.0 + 75.0 * mj_mi[i]) / 12.0;
				param->coeff.mecoSO[i] = -spin_MPow2[i] * 5.0 * param->eta * (4.0 + 3.0 * mj_mi[i])
						/ 9.0;
				param->coeff.chihSO[i] = (4.0 + 3.0 * mj_mi[i]) * param->eta / 2.0;
			}
		}
		if (param->interactionFlags != 0) {
			for (UINT2 i = 0; i < 2; i++) {
				param->coeff.lnhat[i] = -spin_MPow2[i] / param->eta;
			}
		}
	case PN1_0:
		param->coeff.omega[PN1_0] = -(743.0 + 924.0 * param->eta) / 336.0;
		param->coeff.meco[PN1_0] = (9.0 * param->eta + etaPow2) / 18.0;
	case PN0_5:
	case PN0_0:
		param->coeff.omega[PN0_0] = 1.0;
		param->coeff.meco[PN0_0] = -param->eta / 3.0;
		break;
	default:
		XLALPrintError("XLAL Error - %s: Invalid phase. PN order %s\n", __func__,
				param->orderOfPhase);
		XLAL_ERROR(XLAL_EINVAL);
		break;
	}
	return XLAL_SUCCESS;
}

static int XLALEvolveParameters(REAL8Array **out, REAL8 *length, Parameters *param,
		REAL8 initialFrequency) {
	REAL8 input[EVOLVING_VARIABLES];
	input[PHASE] = 0.0;
	input[OMEGA] = LAL_PI * param->totalMass * initialFrequency;
	for (UINT2 dimension = X; dimension < DIMENSIONS; dimension++) {
		input[LNHATX + dimension] = param->lnhat[dimension];
		input[CHI1X + dimension] = param->chih[0][dimension];
		input[CHI2X + dimension] = param->chih[1][dimension];
		input[E1X + dimension] = param->e1[dimension];
	}
	ark4GSLIntegrator *integrator = XLALAdaptiveRungeKutta4Init(EVOLVING_VARIABLES, XLALDerivator,
			XLALStop, ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE);
	if (!integrator) {
		XLALPrintError("XLAL Error - %s: Cannot allocate integrator\n", __func__);
		XLAL_ERROR(XLAL_EFUNC);
	}
	integrator->stopontestonly = 1;
	// run the integration; note: time is measured in \hat{t} = t / M
	*length = XLALAdaptiveRungeKutta4(integrator, (void *) param, input, 0.0,
			*length / param->totalMass, param->samplingTime / param->totalMass, out);
	INT4 intreturn = integrator->returncode;
	XLALAdaptiveRungeKutta4Free(integrator);
	if (!*length) {
		XLALPrintError("XLAL Error - %s: integration failed with errorcode %d.\n", __func__,
				intreturn);
		XLAL_ERROR(XLAL_EFUNC);
	}
	if (intreturn != 0 && intreturn != TEST_ENERGY && intreturn != TEST_OMEGADOT) {
		XLALPrintWarning(
				"XLAL Warning - %s: integration terminated with code %d.\n Waveform parameters were m1 = %e, m2 = %e, s1 = (%e,%e,%e), s2 = (%e,%e,%e), inc = %e.\n",
				__func__, intreturn, param->mass[0], param->mass[1],
				param->chih[0][X] * param->chiAmp[0], param->chih[0][Y] * param->chiAmp[0],
				param->chih[0][Z] * param->chiAmp[0], param->chih[1][X] * param->chiAmp[1],
				param->chih[1][Y] * param->chiAmp[1], param->chih[1][Z] * param->chiAmp[1],
				acos(param->lnhat[Z]));
	}
	return XLAL_SUCCESS;
}

static int XLALCreateOrbitOutput(Output *param, UINT4 length, REAL8 samplingTime) {
	LIGOTimeGPS tStart = LIGOTIMEGPSZERO;
	XLALGPSAdd(&tStart, -1.0 * (length - 1) * samplingTime);
	param->pnParameter = XLALCreateREAL8TimeSeries("PN_EXPANSION_PARAMETER", &tStart, 0.0,
			samplingTime, &lalDimensionlessUnit, length);
	param->phase = XLALCreateREAL8TimeSeries("ORBITAL_PHASE", &tStart, 0.0, samplingTime,
			&lalDimensionlessUnit, length);
	if (!param->pnParameter || !param->phase) {
		XLAL_ERROR(XLAL_EFUNC);
	} else {
		const char *chi1[] = { "CHI1_X_COMPONENT", "CHI1_Y_COMPONENT", "CHI1_Z_COMPONENT" };
		const char *chi2[] = { "CHI2_X_COMPONENT", "CHI2_Y_COMPONENT", "CHI2_Z_COMPONENT" };
		const char *lnhat[] = { "E1_BASIS_X_COMPONENT", "E1_BASIS_Y_COMPONENT",
				"E1_BASIS_Z_COMPONENT" };
		const char *e1[] = { "CHI1_X_COMPONENT", "CHI1_Y_COMPONENT", "CHI1_Z_COMPONENT" };
		for (UINT2 dim = X; dim < DIMENSIONS; dim++) {
			param->chi1[dim] = XLALCreateREAL8TimeSeries(chi1[dim], &tStart, 0.0, samplingTime,
					&lalDimensionlessUnit, length);
			param->chi2[dim] = XLALCreateREAL8TimeSeries(chi2[dim], &tStart, 0.0, samplingTime,
					&lalDimensionlessUnit, length);
			param->e3[dim] = XLALCreateREAL8TimeSeries(lnhat[dim], &tStart, 0.0, samplingTime,
					&lalDimensionlessUnit, length);
			param->e1[dim] = XLALCreateREAL8TimeSeries(e1[dim], &tStart, 0.0, samplingTime,
					&lalDimensionlessUnit, length);
			if (!param->chi1[dim] || !param->chi2[dim] || !param->e3[dim] || !param->e1[dim]) {
				XLAL_ERROR(XLAL_EFUNC);
			}
		}
	}
	return XLAL_SUCCESS;
}

static int XLALCreateWaveformOutput(Output *param, UINT4 length, REAL8 samplingTime) {
	LIGOTimeGPS tStart = LIGOTIMEGPSZERO;
	XLALGPSAdd(&tStart, -1.0 * (length - 1) * samplingTime);
	;
	param->waveform[HP] = XLALCreateREAL8TimeSeries("H_PLUS", &tStart, 0.0, samplingTime,
			&lalDimensionlessUnit, length);
	;
	param->waveform[HC] = XLALCreateREAL8TimeSeries("H_CROSS", &tStart, 0.0, samplingTime,
			&lalDimensionlessUnit, length);
	;
	if (!param->waveform[HP] || !param->waveform[HC]) {
		XLAL_ERROR(XLAL_EFUNC);
	}
	memset(param->waveform[HP]->data->data, 0,
			param->waveform[HP]->data->length * sizeof(*param->waveform[HP]->data->data));
	memset(param->waveform[HC]->data->data, 0,
			param->waveform[HC]->data->length * sizeof(*param->waveform[HC]->data->data));
	return XLAL_SUCCESS;
}

inline static void crossProduct(REAL8 out[], const REAL8 left[], const REAL8 right[]) {
	out[0] = left[1] * right[2] - left[2] * right[1];
	out[1] = left[2] * right[0] - left[0] * right[2];
	out[2] = left[0] * right[1] - left[1] * right[0];
}

static int XLALCalculateWaveform(REAL8 *hp, REAL8 *hc, REAL8 e1[], REAL8 e3[], REAL8 phase, REAL8 V,
		Parameters *params, REAL8 distance, INT4 orderOfAmplitude) {
	REAL8 h[2];
	hp = &h[HP];
	hc = &h[HC];
	enum {
		TEMP = 4,
	};
	REAL8 vPowi_3[TEMP] = { 1.0, V, sq(V), vPowi_3[2] * V };
	REAL8 siniPhi[TEMP];
	REAL8 cosiPhi[TEMP];
	for (INT2 order = 1; order < orderOfAmplitude + 1; order++) {
		siniPhi[order] = sin((REAL8) order * phase);
		cosiPhi[order] = cos((REAL8) order * phase);
	}
	REAL8 e2[DIMENSIONS];
	crossProduct(e2, e3, e1);
	REAL8 sine[PNORDER][WAVEFORM], cosine[PNORDER][WAVEFORM];
	switch (orderOfAmplitude) {
	case PNALL:
	case PN0_0:
		cosine[PN0_0][HP] = 0.5 * (sq(e1[X]) - sq(e1[Y]) - sq(e2[X]) + sq(e2[Y]));
		cosine[PN0_0][HC] = e1[X] * e1[Y] - e2[X] * e2[Y];
		sine[PN0_0][HP] = e1[X] * e2[X] - e1[Y] * e2[Y];
		sine[PN0_0][HC] = e1[Y] * e2[X] + e1[X] * e2[Y];
		for (UINT2 wave = HP; wave < WAVEFORM; wave++) {
			h[wave] += vPowi_3[2] * 4.0
					* (sine[PN0_0][wave] * siniPhi[2] + cosine[PN0_0][wave] * cosiPhi[2]);

		}
	case PN0_5: {
		REAL8 deltamPerM = sqrt(1.0 - 4.0 * params->eta);
		for (UINT2 wave = HP; wave < WAVEFORM; wave++) {
			REAL8 coefSine3Phi = -9.0 * (sine[PN0_0][wave] * e2[Z] - cosine[PN0_0][wave] * e1[Z]);
			REAL8 coefCosine3Phi = -9.0 * (cosine[PN0_0][wave] * e2[Z] + sine[PN0_0][wave] * e1[Z]);
			h[wave] -= vPowi_3[3] * deltamPerM
					* (cosine[PN0_5][wave] * cosiPhi[1] - sine[PN0_5][wave] * siniPhi[1]
							+ coefSine3Phi * siniPhi[3] + coefCosine3Phi * cosiPhi[3]) / 2.0;
		}
	}
	case PN1_5:
		break;
	default:
		XLALPrintError("XLAL Error - %s: Amp. corrections not known "
				"to PN order %d, highest is %d\n", __func__, orderOfAmplitude, PN1_5);
		XLAL_ERROR(XLAL_EINVAL);
		break;
	}
	REAL8 amp = -params->totalMass * params->eta * LAL_MRSUN_SI / distance;
	*hp *= amp;
	*hc *= amp;
	return XLAL_SUCCESS;
}

inline static REAL8 innerProduct(const REAL8 left[], const REAL8 right[]) {
	REAL8 product = 0.0;
	for (UINT2 dim = X; dim < DIMENSIONS; dim++) {
		product += left[dim] * right[dim];
	}
	return (product);
}

static int XLALDerivator(double t, const double values[], double dvalues[], void *mparams) {
	Parameters*params = mparams;
	UNUSED(t);
	memset(dvalues, 0, EVOLVING_VARIABLES * sizeof(REAL8));
	params->coeff.omegaPowi_3[0] = 1.;
	params->coeff.omegaPowi_3[1] = cbrt(values[OMEGA]);
	for (UINT2 i = 2; i < PNORDER; i++) {
		params->coeff.omegaPowi_3[i] = params->coeff.omegaPowi_3[i - 1]
				* params->coeff.omegaPowi_3[1];
	}
	for (UINT2 i = PN0_0; i <= params->orderOfPhase; i++) {
		dvalues[OMEGA] += params->coeff.omega[i] * params->coeff.omegaPowi_3[i];
	}
	if (params->orderOfPhase >= PN3_0) {
		dvalues[OMEGA] -= 856.0 / 105.0 * log(params->coeff.omegaPowi_3[2]);
	}
	if ((params->chiAmp[0] || params->chiAmp[1]) && params->interactionFlags) {
		REAL8 lnhatCrossChih[2][DIMENSIONS];
		const REAL8 *chi_p[2] = { values + CHI1X, values + CHI2X };
		params->coeff.chih1chih2 = innerProduct(chi_p[0], chi_p[1]);
		for (UINT2 current = 0; current < 2; current++) {
			params->coeff.lnhchih[current] = innerProduct(values + LNHATX, chi_p[current]);
			crossProduct(lnhatCrossChih[current], values + LNHATX, chi_p[current]);
		}
		switch (params->orderOfPhase) {
		case PNALL:
		case PN4_0:
		case PN3_5:
		case PN3_0:
		case PN2_5:
		case PN2_0: {
			REAL8 temp = 0.0;
			if (isSet(params->interactionFlags, LAL_SIM_INSPIRAL_INTERACTION_QUAD_MONO_2PN)) {
				for (UINT2 current = 0; current < 2; current++) {
					temp += params->coeff.omegaQM[current]
							* (3.0 * sq(params->coeff.lnhchih[current]) - 1.0);
					for (UINT2 dim = X; dim < DIMENSIONS; dim++) {
						dvalues[CHI1X + DIMENSIONS * current + dim] += params->coeff.chihQM[current]
								* params->coeff.lnhchih[current] * lnhatCrossChih[current][dim]
								* sq(values[OMEGA]);
					}
				}
			}
			if (isSet(params->interactionFlags, LAL_SIM_INSPIRAL_INTERACTION_SPIN_SPIN_2PN)) {
				temp += (params->coeff.omegaSS[0] * params->coeff.lnhchih[0]
						* params->coeff.lnhchih[1]
						+ params->coeff.omegaSS[1] * params->coeff.chih1chih2);
				REAL8 chihjCrossChihi[DIMENSIONS];
				for (UINT2 current = 0; current < 2; current++) {
					crossProduct(chihjCrossChihi, chi_p[!current], chi_p[current]);
					for (UINT2 dim = X; dim < DIMENSIONS; dim++) {
						dvalues[CHI1X + DIMENSIONS * current + dim] += params->coeff.chihSS[current]
								* (chihjCrossChihi[dim]
										- 3.0 * params->coeff.lnhchih[!current]
												* lnhatCrossChih[current][dim]) * sq(values[OMEGA]);
					}
				}

			}
			if (isSet(params->interactionFlags, LAL_SIM_INSPIRAL_INTERACTION_SPIN_SPIN_SELF_2PN)) {
				for (UINT2 current = 0; current < 2; current++) {
					temp += params->coeff.omegaSELF[current]
							* (7 - sq(params->coeff.lnhchih[current]));
				}
			}
			dvalues[OMEGA] += temp * params->coeff.omegaPowi_3[4];
		}
		case PN1_5:
			if (isSet(params->interactionFlags, LAL_SIM_INSPIRAL_INTERACTION_SPIN_ORBIT_15PN)) {
				for (UINT2 current = 0; current < 2; current++) {
					dvalues[OMEGA] += params->coeff.omegaSO[current]
							* params->coeff.lnhchih[current] * values[OMEGA];
					for (UINT2 dim = X; dim < DIMENSIONS; dim++) {
						dvalues[CHI1X + DIMENSIONS * current + dim] += params->coeff.chihSO[current]
								* lnhatCrossChih[current][dim] * params->coeff.omegaPowi_3[5];
					}
				}
			}
			for (UINT2 dim = X; dim < DIMENSIONS; dim++) {
				dvalues[LNHATX + dim] += (params->coeff.lnhat[0] * dvalues[CHI1X + dim]
						+ params->coeff.lnhat[1] * dvalues[CHI2X + dim])
						* params->coeff.omegaPowi_3[1];
			}
			REAL8 omegaLNhat[DIMENSIONS], omegaE[DIMENSIONS], omegaLNhatLNhat;
			crossProduct(omegaLNhat, values + LNHATX, dvalues + LNHATX);
			omegaLNhatLNhat = innerProduct(omegaLNhat, values + LNHATX);
			for (UINT2 dim = X; dim < DIMENSIONS; dim++) {
				omegaE[dim] = omegaLNhat[dim] - omegaLNhatLNhat * values[LNHATX];
			}
			crossProduct(dvalues + E1X, omegaE, values + E1X);
		case PN1_0:
		case PN0_5:
		case PN0_0:
			break;
		default:
			break;
		}
	}
	dvalues[OMEGA] *= params->coeff.omegaGlobal * params->coeff.omegaPowi_3[7]
			* params->coeff.omegaPowi_3[4];
	dvalues[PHASE] = values[OMEGA];
	return GSL_SUCCESS;
}

static int XLALStop(double t, const double values[], double dvalues[], void *mparams) {
	UNUSED(t);
	Parameters *params = mparams;
	REAL8 meco = params->coeff.meco[0] / params->coeff.omegaPowi_3[1];
	UINT2 i;
	for (i = PN0_0 + 2; i <= params->orderOfPhase; i += 2) {
		meco += params->coeff.meco[i] * params->coeff.omegaPowi_3[i - 1];
	}
	if ((params->chiAmp[0] || params->chiAmp[1]) && params->interactionFlags) {
		switch (params->orderOfPhase) {
		case PNALL:
		case PN4_0:
		case PN3_5:
		case PN3_0:
		case PN2_5:
		case PN2_0: {
			REAL8 temp = 0.0;
			if (isSet(params->interactionFlags, LAL_SIM_INSPIRAL_INTERACTION_QUAD_MONO_2PN)) {
				for (UINT2 current = 0; current < 2; current++) {
					temp += params->coeff.omegaQM[current]
							* (3.0 * sq(params->coeff.lnhchih[current]) - 1.0);
				}
				meco += params->coeff.mecoQM * temp * values[OMEGA];
			}
			if (isSet(params->interactionFlags, LAL_SIM_INSPIRAL_INTERACTION_SPIN_SPIN_2PN)) {
				meco += params->coeff.mecoSS * values[OMEGA]
						* (params->coeff.chih1chih2
								- 3.0 * params->coeff.lnhchih[0] * params->coeff.lnhchih[1]);
			}
		}
		case PN1_5:
			if (isSet(params->interactionFlags, LAL_SIM_INSPIRAL_INTERACTION_SPIN_ORBIT_15PN)) {
				for (UINT2 current = 0; current < 2; current++) {
					meco += params->coeff.mecoSO[current] * params->coeff.lnhchih[current]
							* params->coeff.omegaPowi_3[2];
				}
			}
		case PN1_0:
		case PN0_5:
		case PN0_0:
			break;
		default:
			break;
		}
	}
	if (meco > 0.0) {
		return TEST_ENERGY;
	}
	if (dvalues[OMEGA] < 0.0) {
		return TEST_OMEGADOT;
	}
	if (values[OMEGA] / (params->totalMass * LAL_PI) > 0.5 / params->samplingTime) {
		return TEST_NYQUIST;
	}
	if (isnan(values[OMEGA])) {
		return TEST_OMEGANAN;
	}
	return GSL_SUCCESS;
}
