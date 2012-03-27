/**
 * @file LALSQTPNWaveformInterface.c
 *		Contains function definitions to integrate the SpinQuadTaylor code into the other parts of the LALSuit.
 *	If you want to run the program use the \ref LALSQTPNWaveformTest.c file int the
 *	test directory.
 * @author László Veréb
 * @date 2011.07.11.
 */

#include <lal/LALSQTPNWaveformInterface.h>
#include <lal/LALSQTPNWaveform.h>

extern int switchMode;

/**		The function calculates the parameters from the InspiralTemplate
 * structure. <em>The used parameters are:</em>
 * <ul>
 *	<li>masses of the BHs (or NSs) \f$m_i\f$ in \f$M_\odot\f$</li>
 *	<li>the spin components \f$\chi_{ij}\f$, the values of \f$\sqrt{\sum_j\chi_{ij}}\f$, are between 0 and 1</li>
 *	<li>the quadrupole parameters \f$w_i\in(4,8)\f$ for NSs [1] and \f$w_i=1\f$ for BHs[2] are 1 (default 1)</li>
 *	<li>the inclination (angle between the line of sight and Newtonian orbital angular momentum) \f$\iota\f$ in \f$rad\f$
 *	<li>the initial frequency \f$f_L\f$ in \f$Hz\f$</li>
 *	<li>the distance \f$d\f$ in \f$Mpc\f$</li>
 *	<li>the sampling time \f$t_s\f$ in \f$s\f$</li>
 *	<li>the PN order, see #LALPNOrder</li>
 *	<li>level of accuracy in including spin and quadrupole contributions, see
 *	#LALSQTPNSpinInteraction</li>
 * <ul><br />
 * <em>The calculated parameters:</em>
 *	\f{center}{
 *	\begin{gather}
 *		\displaystyle M_{in}=m_1+m_2,\quad
 *		\mu=\frac{m_1m_2}{M_{in}},\quad
 *		\eta=\frac{\mu}{M_{in}},\\
 *		\chi_i=\sqrt{\sum_{j}\chi_{ij}^2},\quad
 *		\hat{\chi}_{ij}=\dfrac{\chi_{ij}}{\chi_i},\\
 *		f_s=t_s^{-1}\\
 *		A=\frac{4\cdot\eta M_{in}M_\odot\displaystyle\frac{G}{c^2}}{d\cdot3.0856775807\cdot10^{16}\cdot10^6}
 *	\end{gather}
 *	\f}
 * and the initial phase \f$\phi=0\f$
 * Assuming that:
 * <ul>
 * <li>masses are positive</li>
 * <li>eta is positive</li>
 * <li>sampling frequency is positive</li>
 * <li>distance is positive</li>
 * </ul>
 * @param[out]	wave		: the used parameters
 * @param[in]		params	: the inspiral parameters
 */
static void XLALSQTPNFillParams(LALSQTPNWaveformParams *wave, InspiralTemplate *params) {
	wave->mass[0] = params->mass1;
	wave->mass[1] = params->mass2;
	wave->totalMass = wave->mass[0] + wave->mass[1];
	wave->mu = wave->mass[0] * wave->mass[1] / wave->totalMass;
	wave->eta = wave->mu / wave->totalMass;
	wave->chirpMass = wave->totalMass * pow(wave->eta, 3.0 / 5.0);
	wave->deltam_M = sqrt(1.0 - 4.0 * wave->eta);
	wave->chiAmp[0] = wave->chiAmp[1] = 0.;
	UINT2 i;
	for (i = 0; i < 3; i++) {
		wave->chi[0][i] = params->spin1[i];
		wave->chi[1][i] = params->spin2[i];
		wave->chiAmp[0] += SQT_SQR(wave->chi[0][i]);
		wave->chiAmp[1] += SQT_SQR(wave->chi[1][i]);
	}
	wave->chiAmp[0] = sqrt(wave->chiAmp[0]);
	wave->chiAmp[1] = sqrt(wave->chiAmp[1]);
	for (i = 0; i < 3; i++) {
		if (wave->chiAmp[0] != 0.) {
			wave->chih[0][i] = wave->chi[0][i] / wave->chiAmp[0];
		} else {
			wave->chih[0][i] = 0.;
		}
		if (wave->chiAmp[1] != 0.) {
			wave->chih[1][i] = wave->chi[1][i] / wave->chiAmp[1];
		} else {
			wave->chih[1][i] = 0.;
		}
	}
	wave->qmParameter[0] = 1.0; //params->qmParameter[0];
	wave->qmParameter[1] = 1.0; //params->qmParameter[1];
	wave->distance = params->distance;
	wave->inclination = params->inclination;
	wave->lowerFreq = params->fLower;
	wave->finalFreq = (
		params->fFinal < params->fLower ? params->fCutoff :
			(params->fCutoff < params->fFinal ? params->fCutoff : params->fFinal));
	wave->samplingFreq = params->tSampling;
	wave->samplingTime = 1.0 / wave->samplingFreq;
	wave->phi = 0.;
	wave->signalAmp = 4.0 * wave->totalMass * wave->eta * LAL_MRSUN_SI / wave->distance;
	wave->order = params->order;
	wave->spinInteraction = params->interaction;
	wave->amplitudeContribution = params->ampOrder;
}

void LALSQTPNWaveformTemplates(LALStatus *status, REAL4Vector *signalvec1, REAL4Vector *signalvec2,
	InspiralTemplate *params) {

	XLALPrintDeprecationWarning("LALSQTPNWaveformTemplates", "XLALSQTPNWaveformTemplates");
	INITSTATUS(status);
	ATTATCHSTATUSPTR(status);

	if (XLALSQTPNWaveformTemplates(signalvec1, signalvec2, params))
		ABORTXLAL(status);

	DETATCHSTATUSPTR(status);
	RETURN(status);
}

int XLALSQTPNWaveformTemplates(REAL4Vector *signalvec1, REAL4Vector *signalvec2,
	InspiralTemplate *params) {

	// Check the relevant pointers
	if (!signalvec1 || !signalvec1->data || !signalvec2 || !signalvec2->data || !params)
		XLAL_ERROR(XLAL_EFAULT);

	// Check the parameters are sane
	if (params->nStartPad < 0 || params->nEndPad < 0 || params->fLower <= 0
		|| params->tSampling <= 0 || params->totalMass <= 0.)
		XLAL_ERROR(XLAL_EINVAL);

	InspiralInit paramsInit;
	LALSQTPNWaveformParams wave_Params;
	LALSQTPNWave wave;

	XLALInspiralInit(params, &paramsInit);
	if (xlalErrno)
		XLAL_ERROR(XLAL_EFUNC);

	memset(signalvec1->data, 0, signalvec1->length * sizeof(REAL4));
	memset(signalvec2->data, 0, signalvec2->length * sizeof(REAL4));
	XLALSQTPNFillParams(&wave_Params, params);

	wave.waveform = NULL;
	wave.hp = signalvec1;
	wave.hc = signalvec2;

	/* Call the engine function */
	switch (switchMode) {
	case LALSQTPN_FIXED:
		if (XLALSQTPNGeneratorFixed(&wave, &wave_Params)) {
			XLALSQTPNDestroyCoherentGW(wave.waveform);
			XLAL_ERROR(XLAL_EFUNC);
		}
		break;
	case LALSQTPN_PRECESSING:
		if (XLALSQTPNGenerator(&wave, &wave_Params)) {
			XLALSQTPNDestroyCoherentGW(wave.waveform);
			XLAL_ERROR(XLAL_EFUNC);
		}
		break;
	default:
		break;
	}
	return XLAL_SUCCESS;
}

void LALSQTPNWaveform(LALStatus *status, REAL4Vector *signalvec, InspiralTemplate *params) {

	XLALPrintDeprecationWarning("LALSQTPNWaveform", "XLALSQTPNWaveform");
	INITSTATUS(status);
	ATTATCHSTATUSPTR(status);

	if (XLALSQTPNWaveform(signalvec, params))
		ABORTXLAL(status);

	DETATCHSTATUSPTR(status);
	RETURN(status);
}

int XLALSQTPNWaveform(REAL4Vector *signalvec, InspiralTemplate *params) {

	// Check the relevant pointers
	if (!signalvec || !signalvec->data || !params)
		XLAL_ERROR(XLAL_EFAULT);

	// Check the parameters are sane
	if (params->nStartPad < 0 || params->nEndPad < 0 || params->fLower <= 0
		|| params->tSampling <= 0 || params->totalMass <= 0.)
		XLAL_ERROR(XLAL_EINVAL);

	InspiralInit paramsInit;
	LALSQTPNWaveformParams wave_Params;
	LALSQTPNWave wave;
	memset(&wave, 0, sizeof(LALSQTPNWave));
	wave.waveform = NULL;
	wave.hp = signalvec;
	wave.hc = NULL;

	XLALInspiralSetup(&(paramsInit.ak), params);
	if (xlalErrno)
		XLAL_ERROR(XLAL_EFUNC);
	if (XLALInspiralChooseModel(&(paramsInit.func), &(paramsInit.ak), params))
		XLAL_ERROR(XLAL_EFUNC);

	XLALSQTPNFillParams(&wave_Params, params);
	//wave_Params.distance *= LAL_PC_SI * 1.e6;
	//wave_Params.signalAmp /= LAL_PC_SI * 1.e6;

	/* Call the engine function */
	switch (switchMode) {
	case LALSQTPN_FIXED:
		if (XLALSQTPNGeneratorFixed(&wave, &wave_Params)) {
			XLALSQTPNDestroyCoherentGW(wave.waveform);
			XLAL_ERROR(XLAL_EFUNC);
		}
		break;
	case LALSQTPN_PRECESSING:
		if (XLALSQTPNGenerator(&wave, &wave_Params)) {
			XLALSQTPNDestroyCoherentGW(wave.waveform);
			XLAL_ERROR(XLAL_EFUNC);
		}
		break;
	default:
		break;
	}
	//params->tC = wave_Params.coalescenceTime;

	return XLAL_SUCCESS;
}

void LALSQTPNWaveformForInjection(LALStatus *status, CoherentGW *waveform, InspiralTemplate *params,
	PPNParamStruc *ppnParams) {

	XLALPrintDeprecationWarning("LALSQTPNWaveformForInjection", "XLALSQTPNWaveformForInjection");
	INITSTATUS(status);
	ATTATCHSTATUSPTR(status);

	if (XLALSQTPNWaveformForInjection(waveform, params, ppnParams))
		ABORTXLAL(status);

	DETATCHSTATUSPTR(status);
	RETURN(status);
}

int XLALSQTPNWaveformForInjection(CoherentGW *waveform, InspiralTemplate *params,
	PPNParamStruc *ppnParams) {

	// Check the relevant pointers
	if (!waveform || !params || waveform->a || waveform->f || waveform->phi || waveform->shift)
		XLAL_ERROR(XLAL_EFAULT);

	// variable declaration and initialization
	InspiralInit paramsInit;

	// Compute some parameters
	XLALInspiralInit(params, &paramsInit);
	if (xlalErrno)
		XLAL_ERROR(XLAL_EFUNC);
	if (paramsInit.nbins == 0) {
		XLALPrintWarning(
			"Warning! Waveform of zero length requested in %s. Returning empty waveform.\n ",
			__func__);
		return XLAL_SUCCESS;
	}

	// Allocate the waveform structures.
	xlalErrno = 0;
	XLALSQTPNAllocateHplusAndHcrossInCoherentGW(waveform, paramsInit.nbins);

	LALSQTPNWave wave;
	wave.waveform = waveform;
	wave.hp = wave.hc = NULL;
	LALSQTPNWaveformParams wave_Params;

	// filling the parameters
	XLALSQTPNFillParams(&wave_Params, params);
	if (xlalErrno)
		XLAL_ERROR(XLAL_EFUNC);

	// calling the engine function
	switch (switchMode) {
	case LALSQTPN_FIXED:
		if (XLALSQTPNGeneratorFixed(&wave, &wave_Params)) {
			XLALSQTPNDestroyCoherentGW(wave.waveform);
			XLAL_ERROR(XLAL_EFUNC);
		}
		break;
	case LALSQTPN_PRECESSING:
		if (XLALSQTPNGenerator(&wave, &wave_Params)) {
			XLALSQTPNDestroyCoherentGW(wave.waveform);
			XLAL_ERROR(XLAL_EFUNC);
		}
		break;
	default:
		break;
	}
	params->fFinal = wave_Params.finalFreq;
	{
		// --- fill some output ---
		ppnParams->tc = (REAL8) (wave.length - 1) / params->tSampling;
		ppnParams->length = wave.length;
		ppnParams->fStop = params->fFinal;
		ppnParams->termCode = GENERATEPPNINSPIRALH_EFSTOP;
		ppnParams->termDescription = GENERATEPPNINSPIRALH_MSGEFSTOP;
		ppnParams->fStart = ppnParams->fStartIn;
	} // end phase condition

	return XLAL_SUCCESS;
}

int XLALSQTPNAllocateHplusAndHcrossInCoherentGW(CoherentGW *wave, UINT4 length) {
	if (!wave) {
		XLAL_ERROR(XLAL_EFAULT);
	}
	if (length <= 0) {
		XLAL_ERROR(XLAL_EBADLEN);
	}
	if (wave->h) {
		XLAL_ERROR(XLAL_EFAULT);
	}
	wave->h = (REAL4TimeVectorSeries *) LALMalloc(sizeof(REAL4TimeVectorSeries));
	if (!wave->h) {
		XLALSQTPNDestroyCoherentGW(wave);
		XLAL_ERROR(XLAL_ENOMEM);
	}
	xlalErrno = 0;
	wave->h->data = XLALCreateREAL4VectorSequence(length, 2);
	memset(wave->h->data->data, 0, 2 * length * sizeof(REAL4));
	if (!wave->h->data) {
		XLALSQTPNDestroyCoherentGW(wave);
		XLAL_ERROR(XLAL_ENOMEM);
	}
	return XLAL_SUCCESS;
}

int XLALSQTPNAllocateAmplitudeInCoherentGW(CoherentGW *wave, UINT4 length) {
	if (!wave) {
		XLAL_ERROR(XLAL_EFAULT);
	}
	if (length <= 0) {
		XLAL_ERROR(XLAL_EBADLEN);
	}
	if (wave->a || wave->f || wave->phi || wave->shift) {
		XLAL_ERROR(XLAL_EFAULT);
	}
	wave->a = (REAL4TimeVectorSeries *) LALMalloc(sizeof(REAL4TimeVectorSeries));
	wave->f = (REAL4TimeSeries *) LALMalloc(sizeof(REAL4TimeSeries));
	wave->phi = (REAL8TimeSeries *) LALMalloc(sizeof(REAL8TimeSeries));
	wave->shift = (REAL4TimeSeries *) LALMalloc(sizeof(REAL4TimeSeries));
	if (!(wave->a && wave->f && wave->phi && wave->shift)) {
		XLALSQTPNDestroyCoherentGW(wave);
		XLAL_ERROR(XLAL_ENOMEM);
	}
	xlalErrno = 0;
	wave->a->data = XLALCreateREAL4VectorSequence(length, 2);
	wave->f->data = XLALCreateREAL4Vector(length);
	wave->phi->data = XLALCreateREAL8Vector(length);
	wave->shift->data = XLALCreateREAL4Vector(length);
	if (!(wave->a->data && wave->f->data && wave->phi->data && wave->shift->data)) {
		XLALSQTPNDestroyCoherentGW(wave);
		XLAL_ERROR(XLAL_ENOMEM);
	}
	return XLAL_SUCCESS;
}

void XLALSQTPNDestroyCoherentGW(CoherentGW *wave) {
	if (wave->h) {
		if (wave->h->data) {
			XLALDestroyREAL4VectorSequence(wave->h->data);
		}
		XLALFree(wave->h);
		wave->h = NULL;
	}
	if (wave->a) {
		if (wave->a->data) {
			XLALDestroyREAL4VectorSequence(wave->a->data);
		}
		XLALFree(wave->a);
		wave->a = NULL;
	}
	if (wave->f) {
		if (wave->f->data) {
			XLALDestroyREAL4Vector(wave->f->data);
		}
		XLALFree(wave->f);
		wave->f = NULL;
	}
	if (wave->phi) {
		if (wave->phi->data) {
			XLALDestroyREAL8Vector(wave->phi->data);
		}
		XLALFree(wave->phi);
		wave->phi = NULL;
	}
	if (wave->shift) {
		if (wave->shift->data) {
			XLALDestroyREAL4Vector(wave->shift->data);
		}
		XLALFree(wave->shift);
		wave->shift = NULL;
	}
}
