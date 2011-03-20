/**
 * @file LALSQTPNWaveform.c
 *		Contains the function definition to create GWforms.
 * @author László Veréb
 * @date 2010.05.21.
 */

#include <lal/LALSQTPNWaveform.h>
#include <lal/LALSQTPNIntegrator.h>
#include <lal/LALSQTPNWaveformInterface.h>

NRCSID (LALSQTPNWAVEFORMC, "$Id LALSQTPN_Waveform.c$");

/**		The macro function calculates the scalar product of two vectors.
 * @param[in]  a1	: the left vector
 * @param[in]  a2	: the right vector
 * @return	the product
 */
#define SCALAR_PRODUCT3(a1, a2) \
	((a1)[0] * (a2)[0] + (a1)[1] * (a2)[1] + (a1)[2] * (a2)[2]);

/**		The macro function calculates the vector product of two vectors.
 * @param[in]  left		: the left vector
 * @param[in]  right	: the right vector
 * @param[out] product	: the vector product
 */
#define VECTOR_PRODUCT3(left, right, product)\
	(product)[0] = ((left)[1] * (right)[2] - (left)[2] * (right)[1]);\
	(product)[1] = ((left)[2] * (right)[0] - (left)[0] * (right)[2]);\
	(product)[2] = ((left)[0] * (right)[1] - (left)[1] * (right)[0]);

void XLALSQTPNFillCoefficients(LALSQTPNWaveformParams * const params) {

	// variable declaration and initialization
	REAL8 thetahat = 1039. / 4620.;
	REAL8 spin_MPow2[2];
	REAL8 m_m[2] = { params->mass[1] / params->mass[0], params->mass[0] / params->mass[1] };
	REAL8 piPow2 = SQT_SQR(LAL_PI);
	REAL8 etaPow2 = SQT_SQR(params->eta);
	REAL8 etaPow3 = etaPow2 * params->eta;
	INT2 i;
	for (i = 0; i < 2; i++) {
		spin_MPow2[i] = params->chiAmp[i] * SQT_SQR(params->mass[i]) / SQT_SQR(params->totalMass);
	}

	// calculating the coefficients
	params->coeff.domegaGlobal = params->eta * 96. / 5.;
	for (i = LAL_PNORDER_NEWTONIAN; i < LAL_PNORDER_PSEUDO_FOUR; i += 2) {
		params->coeff.meco[i] = -0.5 * params->eta * (REAL8) (i + 2) / 3.;
	}
	switch (params->order) {
	case LAL_PNORDER_THREE_POINT_FIVE:
		params->coeff.domega[LAL_PNORDER_THREE_POINT_FIVE] = (-4415. / 4032. + params->eta
				* 358675. / 6048. + etaPow2 * 91495. / 1512.) * LAL_PI;
	case LAL_PNORDER_THREE:
		params->coeff.domega[LAL_PNORDER_THREE] = (16447322263. / 139708800. - LAL_GAMMA * 1712.
				/ 105. + piPow2 * 16. / 3.) + (-273811877. / 1088640. + piPow2 * 451. / 48.
				- thetahat * 88. / 3.) * params->eta + etaPow2 * 541. / 896. - etaPow3 * 5605.
				/ 2592.;
		params->coeff.domegaLN = -856. / 105.;
		params->coeff.meco[LAL_PNORDER_THREE] *= -675. / 64. + (209323. / 4032. - 205. * piPow2
				/ 96. + (110. / 9.) * (1987. / 3080.)) * params->eta - 155. * etaPow2 / 96. - 35.
				* etaPow3 / 5184.;
	case LAL_PNORDER_TWO_POINT_FIVE:
		params->coeff.domega[LAL_PNORDER_TWO_POINT_FIVE] = -(1.0 / 672.0) * (4159. + 15876.
				* params->eta) * LAL_PI;
	case LAL_PNORDER_TWO:
		params->coeff.domega[LAL_PNORDER_TWO] = 34103. / 18144. + params->eta * 13661. / 2016.
				+ etaPow2 * 59. / 18.;
		params->coeff.domegaSSselfConst = 0.;
		params->coeff.domegaQMConst = 0.;
		if ((params->spinInteraction & LAL_SSInter) == LAL_SSInter) {
			params->coeff.dchihSS[0] = spin_MPow2[1] / 2.;
			params->coeff.dchihSS[1] = spin_MPow2[0] / 2.;
			params->coeff.domegaSS[0] = 721. * params->eta * params->chiAmp[0] * params->chiAmp[1]
					/ 48.;
			params->coeff.domegaSS[1] = -247. * params->coeff.domegaSS[0] / 721.;
			params->coeff.mecoSS = -spin_MPow2[0] * spin_MPow2[1];
		}
		if ((params->spinInteraction & LAL_SSselfInter) == LAL_SSselfInter) {
			for (i = 0; i < 2; i++) {
				params->coeff.domegaSSself[i] = -spin_MPow2[i] * params->chiAmp[i] / 96.;
				params->coeff.domegaSSselfConst -= 7. * params->coeff.domegaSSself[i];
			}
		}
		if ((params->spinInteraction & LAL_QMInter) == LAL_QMInter) {
			for (i = 0; i < 2; i++) {
				params->coeff.domegaQM[i] = spin_MPow2[i] * params->chiAmp[i]
						* params->qmParameter[i] * 7.5;
				params->coeff.domegaQMConst -= params->coeff.domegaQM[i] / 3.;
				params->coeff.dchihQM[i] = -params->qmParameter[i] * params->eta
						* params->chiAmp[i] * 3. / 2.;
			}
			params->coeff.mecoQM = 2. * params->eta;
		}
		params->coeff.meco[LAL_PNORDER_TWO] *= (-81. + 57. * params->eta - etaPow2) / 24.;
	case LAL_PNORDER_ONE_POINT_FIVE:
		params->coeff.domega[LAL_PNORDER_ONE_POINT_FIVE] = 4. * LAL_PI;
		if ((params->spinInteraction & LAL_SOInter) == LAL_SOInter) {
			for (short j = 0; j < 2; j++) {
				params->coeff.domegaSO[j] = -spin_MPow2[j] * (113.0 + 75.0 * m_m[j]) / 12.0;
				params->coeff.mecoSO[j] = -spin_MPow2[j] * 5.0 * params->eta * (4.0 + 3.0 * m_m[j])
						/ 9.;
				params->coeff.dchihSO[j] = (4.0 + 3.0 * m_m[j]) * params->eta / 2.0;
			}
		}
		if (params->spinInteraction != 0) {
			for (i = 0; i < 2; i++) {
				params->coeff.dLNh[i] = -spin_MPow2[i] / params->eta;
			}
		}
	case LAL_PNORDER_ONE:
		params->coeff.domega[LAL_PNORDER_ONE] = -(1.0 / 336.0) * (743. + 924. * params->eta);
		params->coeff.meco[LAL_PNORDER_ONE] *= -(9. + params->eta) / 12.;
	case LAL_PNORDER_HALF:
		params->coeff.domega[LAL_PNORDER_HALF] = 0.;
	case LAL_PNORDER_NEWTONIAN:
		params->coeff.domega[LAL_PNORDER_NEWTONIAN] = 1.;
	default:
		break;
	}
}

int LALSQTPNDerivator(REAL8 t, const REAL8 values[], REAL8 dvalues[], void * param) {

	// variable declaration and initialization
	LALSQTPNWaveformParams *params = param;
	UNUSED(t);
	const REAL8 *chi_p[2] = { values + LALSQTPN_CHIH1_1, values + LALSQTPN_CHIH2_1 };
	UINT2 i; // indexes
	memset(dvalues, 0, LALSQTPN_NUM_OF_VAR * sizeof(REAL8));
	REAL8 omegaPowi_3[8];
	omegaPowi_3[0] = 1.;
	omegaPowi_3[1] = cbrt(values[LALSQTPN_OMEGA]);
	omegaPowi_3[1] = pow(values[LALSQTPN_OMEGA], 1. / 3.);
	for (i = 2; i < 8; i++) {
		omegaPowi_3[i] = omegaPowi_3[i - 1] * omegaPowi_3[1];
	}
	REAL8 SS_Omega, SSself_Omega, QM_Omega;
	SS_Omega = SSself_Omega = QM_Omega = 0.;
	REAL8 chih1chih2, LNhchih[2], LNhxchih[2][3];
	chih1chih2 = SCALAR_PRODUCT3(chi_p[0], chi_p[1]);
	for (i = 0; i < 2; i++) {
		LNhchih[i] = SCALAR_PRODUCT3(values + LALSQTPN_LNH_1, chi_p[i]);
		VECTOR_PRODUCT3(values + LALSQTPN_LNH_1, chi_p[i], LNhxchih[i]);
	}

	// calculating domega and MECO without the spin components
	for (i = LAL_PNORDER_NEWTONIAN; i <= params->order; i++) {
		dvalues[LALSQTPN_OMEGA] += params->coeff.domega[i] * omegaPowi_3[i];
	}
	dvalues[LALSQTPN_MECO] += params->coeff.meco[0] / omegaPowi_3[1];
	for (i = LAL_PNORDER_NEWTONIAN + 2; i <= params->order; i += 2) {
		dvalues[LALSQTPN_MECO] += params->coeff.meco[i] * omegaPowi_3[i - 1];
	}

	// calculating the other derivatives and the domega and MECO with spin
	// components

	switch (params->order) {
	case LAL_PNORDER_THREE_POINT_FIVE:
	case LAL_PNORDER_THREE:
		dvalues[LALSQTPN_OMEGA] += params->coeff.domegaLN * log(16. * omegaPowi_3[2])
				* omegaPowi_3[LAL_PNORDER_THREE];
	case LAL_PNORDER_TWO_POINT_FIVE:
	case LAL_PNORDER_TWO:
		if ((params->spinInteraction & LAL_SSInter) == LAL_SSInter) {
			XLALSQTPNAddSSContributions(params, values, dvalues);
		}
		if ((params->spinInteraction & LAL_SSselfInter) == LAL_SSselfInter) {
			XLALSQTPNAddSelfContributions(params, values, dvalues);
		}
		if ((params->spinInteraction & LAL_QMInter) == LAL_QMInter) {
			XLALSQTPNAddQMContributions(params, values, dvalues);
		}
	case LAL_PNORDER_ONE_POINT_FIVE:
		if ((params->spinInteraction & LAL_SOInter) == LAL_SOInter) {
			XLALSQTPNAddSOContributions(params, values, dvalues);
		}
		if (params->spinInteraction) {
			for (i = 0; i < 3; i++) {
				dvalues[LALSQTPN_LNH_1 + i] += (params->coeff.dLNh[0] * dvalues[LALSQTPN_CHIH1_1
						+ i] + params->coeff.dLNh[1] * dvalues[LALSQTPN_CHIH2_1 + i])
						* omegaPowi_3[1];
			}
		}
	case LAL_PNORDER_ONE:
	case LAL_PNORDER_HALF:
	case LAL_PNORDER_NEWTONIAN:
	default:
		break;
	}
	dvalues[LALSQTPN_OMEGA] *= params->coeff.domegaGlobal * omegaPowi_3[7] * omegaPowi_3[4];
	dvalues[LALSQTPN_PHASE] = values[LALSQTPN_OMEGA] + values[LALSQTPN_LNH_3]
			* (values[LALSQTPN_LNH_2] * dvalues[LALSQTPN_LNH_1] - values[LALSQTPN_LNH_1]
					* dvalues[LALSQTPN_LNH_2]) / (SQT_SQR(values[LALSQTPN_LNH_1]) + SQT_SQR(
			values[LALSQTPN_LNH_2]));
	return GSL_SUCCESS;
}

void LALSQTPNGenerator(LALStatus *status, LALSQTPNWave *waveform, LALSQTPNWaveformParams *params) {
	INITSTATUS(status, "LALSQTPNGenerator", LALSQTPNWAVEFORMC);
	ATTATCHSTATUSPTR(status);
	ASSERT(params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
	ASSERT(waveform, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
	// variable declaration and initialization
	UINT4 i = 0; // index
	REAL8 time = 0.;
	REAL8 LNhztol = 1.0e-8;
	REAL8 alpha, amp;
	const REAL8 geometrized_m_total = params->totalMass * LAL_MTSUN_SI;
	const REAL8 freq_Step = geometrized_m_total * LAL_PI;
	const REAL8 step = params->samplingTime / geometrized_m_total;
	REAL8 values[LALSQTPN_NUM_OF_VAR], dvalues[LALSQTPN_NUM_OF_VAR];
	LALSQTPNIntegratorSystem integrator;
	xlalErrno = 0;
	if (XLALSQTPNIntegratorInit(&integrator, LALSQTPN_NUM_OF_VAR, params, LALSQTPNDerivator)) {
		if (XLAL_ENOMEM == XLALClearErrno()) {
			ABORT(status, LALINSPIRALH_EMEM, LALINSPIRALH_MSGEMEM);
		} else {
			ABORTXLAL(status);
		}
	}

	// initializing the dynamic variables
	values[LALSQTPN_PHASE] = params->phi;
	values[LALSQTPN_OMEGA] = params->lowerFreq * freq_Step;
	values[LALSQTPN_LNH_1] = sin(params->inclination); ///< \f$\hat{L_N}=\sin\iota\f$
	values[LALSQTPN_LNH_2] = 0.; ///< \f$\hat{L_N}=0\f$
	values[LALSQTPN_LNH_3] = cos(params->inclination); ///< \f$\hat{L_N}=\cos\iota\f$
	values[LALSQTPN_MECO] = 0.;
	for (i = 0; i < 3; i++) {
		values[LALSQTPN_CHIH1_1 + i] = params->chih[0][i];
		values[LALSQTPN_CHIH2_1 + i] = params->chih[1][i];
	}

	// filling the LALSQTPNCoefficients
	xlalErrno = 0;
	XLALSQTPNFillCoefficients(params);
	if (xlalErrno) {
		ABORTXLAL(status);
	}
	LALSQTPNDerivator(time, values, dvalues, params);
	dvalues[LALSQTPN_MECO] = -1.; // to be able to start the loop
	i = 0;
	do {
		alpha = atan2(values[LALSQTPN_LNH_2], values[LALSQTPN_LNH_1]);
		amp = params->signalAmp * pow(values[LALSQTPN_OMEGA], 2. / 3.);

		// calculating the waveform components
		if (waveform->waveform->h) {
			XLALSQTPNCalculateHPHC2(params, values, &(waveform->waveform->h->data->data[2 * i]));
		}
		if (waveform->waveform->a) {
			waveform->waveform->a->data->data[2 * i] = -amp * 0.5 * (1. + values[LALSQTPN_LNH_3]
					* values[LALSQTPN_LNH_3]);
			waveform->waveform->a->data->data[2 * i + 1] = -amp * values[LALSQTPN_LNH_3];
			waveform->waveform->phi->data->data[i] = 2. * (values[LALSQTPN_PHASE] - params->phi);
			waveform->waveform->shift->data->data[i] = 2. * alpha;
			waveform->waveform->f->data->data[i] = values[LALSQTPN_OMEGA] / freq_Step;
		}

		// evolving
		time = i++ * params->samplingTime;
		xlalErrno = 0;
		if (XLALSQTPNIntegratorFunc(values, &integrator, step)) {
			ABORTXLAL(status);
		}
		// if one of the variables is nan, the PN approximation braked down
		if (isnan(values[LALSQTPN_PHASE]) || isnan(values[LALSQTPN_OMEGA]) || isnan(
				values[LALSQTPN_LNH_1]) || isnan(values[LALSQTPN_LNH_2]) || isnan(
				values[LALSQTPN_LNH_3]) || isnan(values[LALSQTPN_CHIH1_1]) || isnan(
				values[LALSQTPN_CHIH1_2]) || isnan(values[LALSQTPN_CHIH1_3]) || isnan(
				values[LALSQTPN_CHIH2_1]) || isnan(values[LALSQTPN_CHIH2_2]) || isnan(
				values[LALSQTPN_CHIH2_3])) {
			break;
		}
		LALSQTPNDerivator(time, values, dvalues, params);
		if ((waveform->waveform && i == waveform->waveform->f->data->length)) {
			XLALSQTPNIntegratorFree(&integrator);
			ABORT(status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
		}
	} while (dvalues[LALSQTPN_MECO] < 0. && dvalues[LALSQTPN_OMEGA] > 0.0 && SQT_SQR(
			values[LALSQTPN_LNH_3]) < 1. - LNhztol && values[LALSQTPN_OMEGA] / freq_Step
			< params->samplingFreq / 2./* && values[LALSQTPN_OMEGA] / freq_Step < params->finalFreq*/);
	if (waveform->waveform->h) {
		params->finalFreq = values[LALSQTPN_OMEGA] / (LAL_PI * geometrized_m_total);
		params->finalFreq = waveform->waveform->f->data->data[i - 1];
		params->coalescenceTime = time;
	}
	if (waveform->waveform->a) {
		params->finalFreq = waveform->waveform->f->data->data[i - 1];
	}

	waveform->length = i;
	XLALSQTPNIntegratorFree(&integrator);
	DETATCHSTATUSPTR(status);
	RETURN(status);
}

void XLALSQTPNCalculateHPHC(LALSQTPNWaveformParams *params, REAL8 values[], REAL4 *h) {
	REAL8 deltaM = (params->mass[0] - params->mass[1]);
	REAL8 amp = -2. * params->totalMass * params->eta * LAL_MRSUN_SI / params->distance * pow(
			values[LALSQTPN_OMEGA], 2. / 3.);
	REAL8 alpha = atan2(values[LALSQTPN_LNH_2], values[LALSQTPN_LNH_1]);
	REAL8 twoAlpha = 2. * alpha;
	REAL8 cos_Iota = values[LALSQTPN_LNH_3];
	REAL8 sin_Iota = sin(acos(cos_Iota));
	REAL8 sinPow2_Iota = SQT_SQR(sin_Iota);
	REAL8 cos_iPhi[5], sin_iPhi[5];

	REAL8 contribution[3][2];
	REAL8 cosine_Part[3][2];
	REAL8 sine_Part[3][2];
	REAL8 K[2];
	UINT2 i;
	for (i = 0; i < 5; i++) {
		cos_iPhi[i] = cos((double) i * values[LALSQTPN_PHASE]);
		sin_iPhi[i] = sin((double) i * values[LALSQTPN_PHASE]);
	}
	if ((params->amplitudeContribution & LALSQTPN_1_0) == LALSQTPN_1_0) {
		REAL8 DELTAX = params->totalMass * (params->chi[1][0] * params->mass[1] - params->chi[0][0]
				* params->mass[0]);
		REAL8 DELTAY = params->totalMass * (params->chi[1][1] * params->mass[1] - params->chi[0][1]
				* params->mass[0]);
		cosine_Part[LAL_PNORDER_ONE][LALSQTPN_PLUS] = -(DELTAY * sin(alpha) - DELTAX * cos(alpha))
				/ SQT_SQR(params->totalMass);
		cosine_Part[LAL_PNORDER_ONE][LALSQTPN_CROSS] = (DELTAY * cos(alpha) + DELTAX * sin(alpha))
				/ SQT_SQR(params->totalMass);
		REAL8 cosine = cos_Iota * cos(alpha);
		REAL8 sine = cos_Iota * sin(alpha);
		sine_Part[LAL_PNORDER_ONE][LALSQTPN_PLUS] = -(DELTAY * cosine + DELTAX * sine)
				/ SQT_SQR(params->totalMass);
		sine_Part[LAL_PNORDER_ONE][LALSQTPN_CROSS] = (DELTAY * sine + DELTAX * cosine)
				/ SQT_SQR(params->totalMass);
	}
	if ((params->amplitudeContribution & LALSQTPN_0_5) == LALSQTPN_0_5
			|| (params->amplitudeContribution & LALSQTPN_1_0) == LALSQTPN_1_0) {
		K[LALSQTPN_PLUS] = -1. / 2. * cos(twoAlpha) * sinPow2_Iota;
		K[LALSQTPN_CROSS] = -1. / 2. * sin(twoAlpha) * sinPow2_Iota;
	}
	//if ((params->amplitudeContribution & LAL_PNORDER_NEWTONIAN) == LAL_PNORDER_NEWTONIAN) {
	cosine_Part[LAL_PNORDER_NEWTONIAN][LALSQTPN_PLUS] = -1. / 2. * (1. + SQT_SQR(cos_Iota)) * cos(
			twoAlpha);
	cosine_Part[LAL_PNORDER_NEWTONIAN][LALSQTPN_CROSS] = -1. / 2. * (1. + SQT_SQR(cos_Iota)) * sin(
			twoAlpha);
	sine_Part[LAL_PNORDER_NEWTONIAN][LALSQTPN_PLUS] = cos_Iota * sin(twoAlpha);
	sine_Part[LAL_PNORDER_NEWTONIAN][LALSQTPN_CROSS] = -cos_Iota * cos(twoAlpha);
	//}
	if ((params->amplitudeContribution & LALSQTPN_0_0) == LALSQTPN_0_0
			|| (params->amplitudeContribution & LALSQTPN_1_0) == LALSQTPN_1_0) {
		for (i = LALSQTPN_PLUS; i <= LALSQTPN_CROSS; i++) {
			contribution[LAL_PNORDER_NEWTONIAN][i] = -2. * (cosine_Part[LAL_PNORDER_NEWTONIAN][i]
					* cos_iPhi[2] + sine_Part[LAL_PNORDER_NEWTONIAN][i] * sin_iPhi[2]);
		}
	}
	if ((params->amplitudeContribution & LALSQTPN_0_5) == LALSQTPN_0_5) {
		for (i = LALSQTPN_PLUS; i <= LALSQTPN_CROSS; i++) {
			contribution[LAL_PNORDER_HALF][i] = 1. / 4. * deltaM / params->totalMass * (3.
					* cosine_Part[LAL_PNORDER_NEWTONIAN][i] * sin_Iota * (3. * cos_iPhi[3]
					- cos_iPhi[1]) + 3. * sine_Part[LAL_PNORDER_NEWTONIAN][i] * sin_Iota * (3.
					* sin_iPhi[3] - sin_iPhi[1]) - 2. * K[i] * sin_Iota * cos_iPhi[1]);
		}
	}
	if ((params->amplitudeContribution & LALSQTPN_1_0) == LALSQTPN_1_0) {
		for (i = LALSQTPN_PLUS; i <= LALSQTPN_CROSS; i++) {
			contribution[LAL_PNORDER_ONE][i] = -8. / 3. * (1. - 3. * params->eta) * sinPow2_Iota
					* (cosine_Part[LAL_PNORDER_NEWTONIAN][i] * cos_iPhi[4]
							+ sine_Part[LAL_PNORDER_NEWTONIAN][i] * sin_iPhi[4])
					+ cosine_Part[LAL_PNORDER_ONE][i] * cos_iPhi[1] + sine_Part[LAL_PNORDER_ONE][i]
					* sin_iPhi[1] + 1. / 6. * (4. * (1. - 3. * params->eta) * sinPow2_Iota
					* cos_iPhi[2] * K[i] - (4. * (1. - 3. * params->eta) * sinPow2_Iota + (19. - 3.
					* params->eta)) * contribution[LAL_PNORDER_NEWTONIAN][i]);
		}
	}
	for (i = LALSQTPN_PLUS; i <= LALSQTPN_CROSS; i++) {
		h[i] = 0.0;
		if ((params->amplitudeContribution & LALSQTPN_1_0) == LALSQTPN_1_0) {
			h[i] += contribution[LAL_PNORDER_ONE][i] * pow(values[LALSQTPN_OMEGA], 2. / 3.);
		}
		if ((params->amplitudeContribution & LALSQTPN_0_5) == LALSQTPN_0_5) {
			h[i] += contribution[LAL_PNORDER_HALF][i] * pow(values[LALSQTPN_OMEGA], 1. / 3.);
		}
		if ((params->amplitudeContribution & LALSQTPN_0_0) == LALSQTPN_0_0) {
			h[i] += contribution[LAL_PNORDER_NEWTONIAN][i];
		}
		h[i] *= amp;
	}
}

void XLALSQTPNCalculateCoefficients0_0order(REAL8 values[], REAL8 cosine[], REAL8 sine[]) {
	REAL8 alpha = atan2(values[LALSQTPN_LNH_2], values[LALSQTPN_LNH_1]);
	REAL8 twoAlpha = 2. * alpha;
	REAL8 cos_Iota = values[LALSQTPN_LNH_3];
	cosine[LALSQTPN_PLUS] = -1. / 2. * (1. + SQT_SQR(cos_Iota)) * cos(twoAlpha);
	cosine[LALSQTPN_CROSS] = -1. / 2. * (1. + SQT_SQR(cos_Iota)) * sin(twoAlpha);
	sine[LALSQTPN_PLUS] = cos_Iota * sin(twoAlpha);
	sine[LALSQTPN_CROSS] = -cos_Iota * cos(twoAlpha);
}

void XLALSQTPNCalculateCoefficients0_5order(REAL8 values[], REAL8 cosine[]) {
	REAL8 alpha = atan2(values[LALSQTPN_LNH_2], values[LALSQTPN_LNH_1]);
	REAL8 twoAlpha = 2. * alpha;
	REAL8 cos_Iota = values[LALSQTPN_LNH_3];
	REAL8 sin_Iota = sin(acos(cos_Iota));
	REAL8 sinPow2_Iota = SQT_SQR(sin_Iota);
	cosine[LALSQTPN_PLUS] = -1. / 2. * cos(twoAlpha) * sinPow2_Iota;
	cosine[LALSQTPN_CROSS] = -1. / 2. * sin(twoAlpha) * sinPow2_Iota;
}

void XLALSQTPNCalculateCoefficients1_0order(LALSQTPNWaveformParams *params, REAL8 values[],
		REAL8 cosine[], REAL8 sine[]) {
	REAL8 alpha = atan2(values[LALSQTPN_LNH_2], values[LALSQTPN_LNH_1]);
	REAL8 cos_Iota = values[LALSQTPN_LNH_3];
	REAL8 DELTAX = params->totalMass * (params->chi[1][0] * params->mass[1] - params->chi[0][0]
			* params->mass[0]);
	REAL8 DELTAY = params->totalMass * (params->chi[1][1] * params->mass[1] - params->chi[0][1]
			* params->mass[0]);
	cosine[LALSQTPN_PLUS] = -(DELTAY * sin(alpha) - DELTAX * cos(alpha))
			/ SQT_SQR(params->totalMass);
	cosine[LALSQTPN_CROSS] = (DELTAY * cos(alpha) + DELTAX * sin(alpha))
			/ SQT_SQR(params->totalMass);
	REAL8 cosineX = cos_Iota * cos(alpha);
	REAL8 sineX = cos_Iota * sin(alpha);
	sine[LALSQTPN_PLUS] = -(DELTAY * cosineX + DELTAX * sineX) / SQT_SQR(params->totalMass);
	sine[LALSQTPN_CROSS] = (DELTAY * sineX + DELTAX * cosineX) / SQT_SQR(params->totalMass);
}

void XLALSQTPNCalculateAmplitudeContribution0_0(REAL8 values[], REAL8 contribution[]) {
	REAL8 cosine_Part[2];
	REAL8 sine_Part[2];
	REAL8 sin_2Phi = sin(2 * values[LALSQTPN_PHASE]);
	REAL8 cos_2Phi = cos(2 * values[LALSQTPN_PHASE]);
	XLALSQTPNCalculateCoefficients0_0order(values, cosine_Part, sine_Part);
	for (short i = LALSQTPN_PLUS; i <= LALSQTPN_CROSS; i++) {
		contribution[i] = -2. * (cosine_Part[i] * cos_2Phi + sine_Part[i] * sin_2Phi);
	}
}

void XLALSQTPNCalculateAmplitudeContribution0_5(LALSQTPNWaveformParams *params, REAL8 values[],
		REAL8 contribution[]) {
	REAL8 K[2];
	REAL8 cosine_Part[2];
	REAL8 sine_Part[2];
	REAL8 cos_Iota = values[LALSQTPN_LNH_3];
	REAL8 sin_Iota = sin(acos(cos_Iota));
	REAL8 deltaM = (params->mass[0] - params->mass[1]);
	REAL8 sin_Phi = sin(values[LALSQTPN_PHASE]);
	REAL8 cos_Phi = cos(values[LALSQTPN_PHASE]);
	REAL8 sin_3Phi = sin(3 * values[LALSQTPN_PHASE]);
	REAL8 cos_3Phi = cos(3 * values[LALSQTPN_PHASE]);
	XLALSQTPNCalculateCoefficients0_0order(values, cosine_Part, sine_Part);
	XLALSQTPNCalculateCoefficients0_5order(values, K);
	for (short i = LALSQTPN_PLUS; i <= LALSQTPN_CROSS; i++) {
		contribution[i] = 1. / 4. * deltaM / params->totalMass * (3. * cosine_Part[i] * sin_Iota
				* (3. * cos_3Phi - cos_Phi) + 3. * sine_Part[i] * sin_Iota * (3. * sin_3Phi
				- sin_Phi) - 2. * K[i] * sin_Iota * cos_Phi);
	}
}

void XLALSQTPNCalculateAmplitudeContribution1_0(LALSQTPNWaveformParams *params, REAL8 values[],
		REAL8 contribution[]) {

	REAL8 cos_Iota = values[LALSQTPN_LNH_3];
	REAL8 sin_Iota = sin(acos(cos_Iota));
	REAL8 sinPow2_Iota = SQT_SQR(sin_Iota);
	REAL8 contribution0[2];
	REAL8 cosine_Part[3][2];
	REAL8 sine_Part[3][2];
	REAL8 sin_Phi = sin(values[LALSQTPN_PHASE]);
	REAL8 cos_Phi = cos(values[LALSQTPN_PHASE]);
	REAL8 cos_2Phi = cos(2 * values[LALSQTPN_PHASE]);
	REAL8 sin_4Phi = sin(4 * values[LALSQTPN_PHASE]);
	REAL8 cos_4Phi = cos(4 * values[LALSQTPN_PHASE]);
	XLALSQTPNCalculateCoefficients0_0order(values, cosine_Part[LAL_PNORDER_NEWTONIAN],
			sine_Part[LAL_PNORDER_NEWTONIAN]);
	XLALSQTPNCalculateCoefficients0_5order(values, cosine_Part[LAL_PNORDER_HALF]);
	XLALSQTPNCalculateCoefficients1_0order(params, values, cosine_Part[LAL_PNORDER_ONE],
			sine_Part[LAL_PNORDER_ONE]);
	XLALSQTPNCalculateAmplitudeContribution0_0(values, contribution0);
	for (short i = LALSQTPN_PLUS; i <= LALSQTPN_CROSS; i++) {
		contribution[i] = -8. / 3. * (1. - 3. * params->eta) * sinPow2_Iota
				* (cosine_Part[LAL_PNORDER_NEWTONIAN][i] * cos_4Phi
						+ sine_Part[LAL_PNORDER_NEWTONIAN][i] * sin_4Phi)
				+ cosine_Part[LAL_PNORDER_ONE][i] * cos_Phi + sine_Part[LAL_PNORDER_ONE][i]
				* sin_Phi + 1. / 6. * (4. * (1. - 3. * params->eta) * sinPow2_Iota * cos_2Phi
				* cosine_Part[LAL_PNORDER_HALF][i] - (4. * (1. - 3. * params->eta) * sinPow2_Iota
				+ (19. - 3. * params->eta)) * contribution0[i]);
	}
}

void XLALSQTPNCalculateHPHC2(LALSQTPNWaveformParams *params, REAL8 values[], REAL4 *h) {
	REAL8 contribution[3][2];
	h[LALSQTPN_PLUS] = h[LALSQTPN_CROSS] = 0.0;
	REAL8 amp = -2. * params->totalMass * params->eta * LAL_MRSUN_SI / params->distance * pow(
			values[LALSQTPN_OMEGA], 2. / 3.);
	if ((params->amplitudeContribution & LALSQTPN_1_0) == LALSQTPN_1_0) {
		XLALSQTPNCalculateAmplitudeContribution1_0(params, values, contribution[LAL_PNORDER_ONE]);
		for (short i = LALSQTPN_PLUS; i <= LALSQTPN_CROSS; i++) {
			h[i] += contribution[LAL_PNORDER_ONE][i] * pow(values[LALSQTPN_OMEGA], 2. / 3.);
		}
	}
	if ((params->amplitudeContribution & LALSQTPN_0_5) == LALSQTPN_0_5) {
		XLALSQTPNCalculateAmplitudeContribution0_5(params, values, contribution[LAL_PNORDER_HALF]);
		for (short i = LALSQTPN_PLUS; i <= LALSQTPN_CROSS; i++) {
			h[i] += contribution[LAL_PNORDER_HALF][i] * pow(values[LALSQTPN_OMEGA], 1. / 3.);
		}
	}
	if ((params->amplitudeContribution & LALSQTPN_0_0) == LALSQTPN_0_0) {
		XLALSQTPNCalculateAmplitudeContribution0_0(values, contribution[LAL_PNORDER_NEWTONIAN]);
		for (short i = LALSQTPN_PLUS; i <= LALSQTPN_CROSS; i++) {
			h[i] += contribution[LAL_PNORDER_NEWTONIAN][i];
		}
	}
	for (short i = LALSQTPN_PLUS; i <= LALSQTPN_CROSS; i++) {
		h[i] *= amp;
	}
}

void XLALSQTPNAddQMContributions(LALSQTPNWaveformParams *params, const REAL8 values[],
		REAL8 dvalues[]) {
	REAL8 temp = params->coeff.domegaQMConst;
	REAL8 LNhchih[2], LNhxchih[2][3];
	const REAL8 *chi_p[2] = { values + LALSQTPN_CHIH1_1, values + LALSQTPN_CHIH2_1 };
	for (short i = 0; i < 2; i++) {
		LNhchih[i] = SCALAR_PRODUCT3(values + LALSQTPN_LNH_1, chi_p[i]);
		VECTOR_PRODUCT3(values + LALSQTPN_LNH_1, chi_p[i], LNhxchih[i]);
	}
	for (short i = 0; i < 2; i++) {
		temp += params->coeff.domegaQM[i] * SQT_SQR(LNhchih[i]);
		for (short j = 0; j < 3; j++) {
			dvalues[LALSQTPN_CHIH1_1 + 3 * i + j] += params->coeff.dchihQM[i] * LNhchih[i]
					* LNhxchih[i][j] * SQT_SQR(values[LALSQTPN_OMEGA]);
		}
	}
	dvalues[LALSQTPN_OMEGA] += temp * pow(values[LALSQTPN_OMEGA], 4. / 3.);
	dvalues[LALSQTPN_MECO] += params->coeff.mecoQM * temp * values[LALSQTPN_OMEGA];
}

void XLALSQTPNAddSSContributions(LALSQTPNWaveformParams *params, const REAL8 values[],
		REAL8 dvalues[]) {
	REAL8 chih1chih2, LNhchih[2], LNhxchih[2][3], chih1xchih2[2][3];
	const REAL8 *chi_p[2] = { values + LALSQTPN_CHIH1_1, values + LALSQTPN_CHIH2_1 };
	chih1chih2 = SCALAR_PRODUCT3(chi_p[0], chi_p[1]);
	for (short i = 0; i < 2; i++) {
		LNhchih[i] = SCALAR_PRODUCT3(values + LALSQTPN_LNH_1, chi_p[i]);
		VECTOR_PRODUCT3(values + LALSQTPN_LNH_1, chi_p[i], LNhxchih[i]);
	}
	for (short i = 0; i < 2; i++) {
		short k = (i + 1) % 2;
		VECTOR_PRODUCT3(chi_p[k], chi_p[i], chih1xchih2[i]);
		for (short j = 0; j < 3; j++) {
			dvalues[LALSQTPN_CHIH1_1 + 3 * i + j] += params->coeff.dchihSS[i] * (chih1xchih2[i][j]
					- 3. * LNhchih[k] * LNhxchih[i][j]) * SQT_SQR(values[LALSQTPN_OMEGA]);
		}
	}
	dvalues[LALSQTPN_MECO] += params->coeff.mecoSS * (chih1chih2 - 3 * LNhchih[0] * LNhchih[1])
			* values[LALSQTPN_OMEGA];
	dvalues[LALSQTPN_OMEGA] += (params->coeff.domegaSS[0] * LNhchih[0] * LNhchih[1]
			+ params->coeff.domegaSS[1] * chih1chih2) * pow(values[LALSQTPN_OMEGA], 4. / 3.);
}

inline void XLALSQTPNAddSelfContributions(LALSQTPNWaveformParams *params, const REAL8 values[],
		REAL8 dvalues[]) {
	REAL8 temp = params->coeff.domegaSSselfConst;
	REAL8 LNhchih;
	const REAL8 *chi_p[2] = { values + LALSQTPN_CHIH1_1, values + LALSQTPN_CHIH2_1 };
	for (short i = 0; i < 2; i++) {
		LNhchih = SCALAR_PRODUCT3(values + LALSQTPN_LNH_1, chi_p[i]);
		temp += params->coeff.domegaSSself[i] * SQT_SQR(LNhchih);
	}
	dvalues[LALSQTPN_OMEGA] += temp * pow(values[LALSQTPN_OMEGA], 4. / 3.);
}


inline void XLALSQTPNAddSOContributions(LALSQTPNWaveformParams *params, const REAL8 values[],
		REAL8 dvalues[]) {
	REAL8 LNhchih, LNhxchih[2][3];
	REAL8 omegaPow5_3 = pow(values[LALSQTPN_OMEGA], 5.0 / 3.0);
	const REAL8 *chi_p[2] = { values + LALSQTPN_CHIH1_1, values + LALSQTPN_CHIH2_1 };
	for (short i = 0; i < 2; i++) {
		LNhchih = SCALAR_PRODUCT3(values + LALSQTPN_LNH_1, chi_p[i]);
		VECTOR_PRODUCT3(values + LALSQTPN_LNH_1, chi_p[i], LNhxchih[i]);
		dvalues[LALSQTPN_OMEGA] += params->coeff.domegaSO[i] * LNhchih * values[LALSQTPN_OMEGA];
		dvalues[LALSQTPN_MECO] += params->coeff.mecoSO[i] * LNhchih * pow(values[LALSQTPN_OMEGA],
				2.0 / 3.0);
	}
	for (short i = 0; i < 3; i++) {
		dvalues[LALSQTPN_CHIH1_1 + i] += params->coeff.dchihSO[0] * LNhxchih[0][i] * omegaPow5_3;
		dvalues[LALSQTPN_CHIH2_1 + i] += params->coeff.dchihSO[1] * LNhxchih[1][i] * omegaPow5_3;
	}
}
