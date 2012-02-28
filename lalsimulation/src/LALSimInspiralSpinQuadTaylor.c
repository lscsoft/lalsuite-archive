/**	@file   LALSimInspiralSpinQuadTaylor.c
 *	@author László Veréb
 *	@date   21 Feb 2012
 *	@brief  
 */

#include <math.h>
#include <string.h>

#include <lal/AVFactories.h>
#include <lal/LALAdaptiveRungeKutta4.h>
#include <lal/LALConstants.h>
#include <lal/XLALError.h>

#include "LALSimInspiralSpinQuadTaylor.h"

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
	REAL8 MECO[PNORDER]; ///< for MECO
	REAL8 MECO_SO[2]; ///< MECO spin-orbit contribution
	REAL8 MECO_SS; ///< MECO spin-spin contribution
	REAL8 MECO_QM; ///< MECO quad-mono contribution
	REAL8 omegaPowi_3[PNORDER]; ///< current powers of the PN-parameter \f$v^{(i/3)}\f$
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
		Parameters *params, INT4 orderOfAmplitude);

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
		REAL8 e1y, REAL8 e1z, REAL8 endPhase, REAL8 initialFrequency, REAL8 samplingTime,
		INT4 orderOfPhase, INT4 orderOfAmplitude, LALSimInspiralInteraction interactionFlags) {
	Parameters parameter;
	memset(&parameter, 0, sizeof(Parameters));
	XLALCalculateConstantParameters(&parameter, mass1, mass2, qm1, qm2, chi1x, chi1y, chi1z, chi2x,
			chi2y, chi2z, lnhatx, lnhaty, lnhatz, e1x, e1y, e1z, samplingTime, orderOfPhase,
			interactionFlags);
	REAL8 length = (5.0 / 256.0) * pow(LAL_PI, -8.0 / 3.0)
			* pow(parameter.chirpMass * initialFrequency, -5.0 / 3.0) / initialFrequency;
	REAL8Array *out = NULL;
	XLALEvolveParameters(&out, &length, &parameter, initialFrequency);
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
			chi2y, chi2z, lnhatx, lnhaty, lnhatz, e1x, e1y, e1z, samplingTime, orderOfPhase,
			interactionFlags);
	REAL8 length = (5.0 / 256.0) * pow(LAL_PI, -8.0 / 3.0)
			* pow(parameter.chirpMass * initialFrequency, -5.0 / 3.0) / initialFrequency;
	REAL8Array *out = NULL;
	XLALEvolveParameters(&out, &length, &parameter, initialFrequency);
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
			chi2y, chi2z, lnhatx, lnhaty, lnhatz, e1x, e1y, e1z, samplingTime, orderOfPhase,
			interactionFlags);
	REAL8 length = (5.0 / 256.0) * pow(LAL_PI, -8.0 / 3.0)
			* pow(parameter.chirpMass * initialFrequency, -5.0 / 3.0) / initialFrequency;
	REAL8Array *out = NULL;
	XLALEvolveParameters(&out, &length, &parameter, initialFrequency);
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
	for (UINT2 order = PN0_0; order <= param->orderOfPhase; order += 2) {
		param->coeff.MECO[order] = -param->eta * (REAL8) (order + 2) / 6.0;
	}
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
		param->coeff.MECO[PN3_0] *= (34445.0 / 576.0 - piPow2 * 205.0 / 96.0) * param->eta
				- etaPow3 * 35.0 / 5184.0 - etaPow2 * 155.0 / 96.0 - 675.0 / 64.0;
	case PN2_5:
		param->coeff.omega[PN2_5] = -LAL_PI * (4159.0 + 15876.0 * param->eta) / 672.0;
	case PN2_0:
		param->coeff.omega[PN2_0] = (34103.0 + 122949.0 * param->eta + 59472.0 * etaPow2) / 18144.0;
		param->coeff.MECO[PN2_0] *= (-81.0 + 57.0 * param->eta - etaPow2) / 24.0;
		if (isSet(param->interactionFlags, LAL_SIM_INSPIRAL_INTERACTION_SPIN_SPIN_2PN)) {
			param->coeff.chihSS[0] = spin_MPow2[1] / 2.0;
			param->coeff.chihSS[1] = spin_MPow2[0] / 2.0;
			param->coeff.omegaSS[0] = 721.0 * param->eta * param->chiAmp[0] * param->chiAmp[1]
					/ 48.0;
			param->coeff.omegaSS[1] = -247.0 * param->coeff.omegaSS[0] / 721.0;
			param->coeff.MECO_SS = -spin_MPow2[0] * spin_MPow2[1];
		}
		if (isSet(param->interactionFlags, LAL_SIM_INSPIRAL_INTERACTION_QUAD_MONO_2PN)) {
			for (UINT2 i = 0; i < 2; i++) {
				param->coeff.omegaQM[i] = spin_MPow2[i] * param->chiAmp[i] * param->qm[i] * 2.5;
				param->coeff.chihQM[i] = -param->qm[i] * param->eta * param->chiAmp[i] * 1.5;
			}
			param->coeff.MECO_QM = param->eta / 5.0;
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
				param->coeff.MECO_SO[i] = -spin_MPow2[i] * 5.0 * param->eta * (4.0 + 3.0 * mj_mi[i])
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
		param->coeff.MECO[PN1_0] *= -(9.0 + param->eta) / 12.0;
	case PN0_5:
	case PN0_0:
		param->coeff.omega[PN0_0] = 1.0;
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

inline static REAL8 innerProduct(const REAL8 left[], const REAL8 right[]) {
	REAL8 product = 0.0;
	for (UINT2 dim = X; dim < DIMENSIONS; dim++) {
		product += left[dim] * right[dim];
	}
	return (product);
}

inline static void crossProduct(REAL8 out[], const REAL8 left[], const REAL8 right[]) {
	out[0] = left[1] * right[2] - left[2] * right[1];
	out[1] = left[2] * right[0] - left[0] * right[2];
	out[2] = left[0] * right[1] - left[1] * right[0];
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
	dvalues[OMEGA] *= params->coeff.omegaGlobal * params->coeff.omegaPowi_3[7]
			* params->coeff.omegaPowi_3[4];
	dvalues[PHASE] = values[OMEGA];
	return GSL_SUCCESS;
}

static int XLALStop(double t, const double values[], double dvalues[], void *mparams) {
	UNUSED(t);
	Parameters *params = mparams;
	REAL8 meco = params->coeff.MECO[0] / params->coeff.omegaPowi_3[1];
	UINT2 i;
	for (i = PN0_0 + 2; i <= params->orderOfPhase; i += 2) {
		meco += params->coeff.MECO[i] * params->coeff.omegaPowi_3[i - 1];
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
