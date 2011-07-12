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

NRCSID (LALSQTPNWAVEFORMINTERFACEC, "$Id LALSQTPNWaveformInterface.c$");

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
	wave->chirpMass = wave->totalMass * pow(wave->eta, 3. / 5.);
	wave->deltam_M = sqrt(1.0 - 4.0 * wave->eta);
	wave->chiAmp[0] = wave->chiAmp[1] = 0.;
	INT2 i;
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
	wave->qmParameter[0] = 1.;//params->qmParameter[0];
	wave->qmParameter[1] = 1.;//params->qmParameter[1];
	wave->distance = params->distance;
	wave->inclination = params->inclination;
	wave->lowerFreq = params->fLower;
	wave->finalFreq = (params->fFinal < params->fLower ? params->fCutoff : (params->fCutoff
			< params->fFinal ? params->fCutoff : params->fFinal));
	wave->samplingFreq = params->tSampling;
	wave->samplingTime = 1. / wave->samplingFreq;
	wave->phi = 0.;
	wave->signalAmp = 4. * wave->totalMass * wave->eta * LAL_MRSUN_SI / wave->distance;
	wave->order = params->order;
	wave->spinInteraction = params->spinInteraction;
	if (wave->spinInteraction) {
		wave->spinInteraction |= LAL_SOInter;
	}
	wave->amplitudeContribution = params->ampOrder;
}

void LALSQTPNWaveformTemplates(LALStatus *status, REAL4Vector *signalvec1, REAL4Vector *signalvec2,
		InspiralTemplate *params) {

	InspiralInit paramsInit;

	INITSTATUS(status, "LALSTPNWaveform", LALSQTPNWAVEFORMINTERFACEC);
	ATTATCHSTATUSPTR(status);

	ASSERT(signalvec1, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
	ASSERT(signalvec1->data, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
	ASSERT(signalvec2, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
	ASSERT(signalvec2->data, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
	ASSERT(params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
	ASSERT(params->nStartPad >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
	ASSERT(params->nEndPad >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
	ASSERT(params->fLower > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
	ASSERT(params->tSampling > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
	ASSERT(params->totalMass > 0., status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
	LALSQTPNWaveformParams wave_Params;
	LALSQTPNWave wave;

	TRY(LALInspiralInit(status->statusPtr, params, &paramsInit), status);

	memset(signalvec1->data, 0, signalvec1->length * sizeof(REAL4));
	memset(signalvec2->data, 0, signalvec2->length * sizeof(REAL4));
	XLALSQTPNFillParams(&wave_Params, params);

	wave.waveform = NULL;
	wave.hp = signalvec1;
	wave.hc = signalvec2;

	/* Call the engine function */
	switch (switchMode) {
		case LALSQTPN_ADAPTIVE:
			LALSQTPNGeneratorAdaptive(status->statusPtr, &wave, &wave_Params);
			break;
		case LALSQTPN_PRECESSING:
			LALSQTPNGenerator(status->statusPtr, &wave, &wave_Params);
			break;
		default:
			break;
	}
	CHECKSTATUSPTR(status);

	DETATCHSTATUSPTR(status);
	RETURN(status);
}

void LALSQTPNWaveform(LALStatus *status, REAL4Vector *signalvec, InspiralTemplate *params) {
	INITSTATUS(status, "LALSQTPNWaveform", LALSQTPNWAVEFORMINTERFACEC);
	ATTATCHSTATUSPTR(status);
	InspiralInit paramsInit;
	LALSQTPNWaveformParams wave_Params;
	LALSQTPNWave wave;
	memset(&wave, 0, sizeof(LALSQTPNWave));
	wave.waveform = NULL;
	wave.hp = signalvec;
	wave.hc = NULL;

	ASSERT(signalvec, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
	ASSERT(signalvec->data, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
	ASSERT(params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
	ASSERT(params->nStartPad >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
	ASSERT(params->nEndPad >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
	ASSERT(params->fLower > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
	ASSERT(params->tSampling > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
	ASSERT(params->totalMass > 0., status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
	LALInspiralSetup(status->statusPtr, &(paramsInit.ak), params);
	CHECKSTATUSPTR(status);
	LALInspiralChooseModel(status->statusPtr, &(paramsInit.func), &(paramsInit.ak), params);
	CHECKSTATUSPTR(status);
	XLALSQTPNFillParams(&wave_Params, params);
	//wave_Params.distance *= LAL_PC_SI * 1.e6;
	//wave_Params.signalAmp /= LAL_PC_SI * 1.e6;

	/* Call the engine function */
	switch (switchMode) {
		case LALSQTPN_ADAPTIVE:
			LALSQTPNGeneratorAdaptive(status->statusPtr, &wave, &wave_Params);
			break;
		case LALSQTPN_PRECESSING:
			LALSQTPNGenerator(status->statusPtr, &wave, &wave_Params);
			break;
		default:
			break;
	}
	//params->tC = wave_Params.coalescenceTime;
	CHECKSTATUSPTR(status);
	DETATCHSTATUSPTR(status);
	RETURN(status);
}

void LALSQTPNWaveformForInjection(LALStatus *status, CoherentGW *waveform,
		InspiralTemplate *params, PPNParamStruc *ppnParams) {
	INITSTATUS(status, "LALSQTPNWaveformInterface", LALSQTPNWAVEFORMINTERFACEC);
	ATTATCHSTATUSPTR(status);
	// variable declaration and initialization
	UINT4 i;
	InspiralInit paramsInit;

	ASSERT(params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
	ASSERT(waveform, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
	ASSERT(!(waveform->a), status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
	ASSERT(!(waveform->f), status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
	ASSERT(!(waveform->phi), status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
	ASSERT(!(waveform->shift), status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
	// Compute some parameters
	LALInspiralInit(status->statusPtr, params, &paramsInit);
	CHECKSTATUSPTR(status);
	if (paramsInit.nbins == 0) {
		DETATCHSTATUSPTR(status);
		RETURN(status);
	}

	// Allocate the waveform structures.
	xlalErrno = 0;
	XLALSQTPNAllocateCoherentGW(waveform, paramsInit.nbins);

	LALSQTPNWave wave;
	wave.waveform = waveform;
	wave.hp = wave.hc = NULL;
	LALSQTPNWaveformParams wave_Params;

	// filling the parameters
	XLALSQTPNFillParams(&wave_Params, params);

	// calling the engine function
	switch (switchMode) {
		case LALSQTPN_ADAPTIVE:
			puts("Adaptive.");
			LALSQTPNGeneratorAdaptive(status->statusPtr, &wave, &wave_Params);
			break;
		case LALSQTPN_PRECESSING:
			puts("Precessing.");
			LALSQTPNGenerator(status->statusPtr, &wave, &wave_Params);
			break;
		default:
			break;
	}
	BEGINFAIL(status) {
		XLALSQTPNDestroyCoherentGW(waveform);
	}ENDFAIL(status);
	params->fFinal = wave_Params.finalFreq;
	for (i = 0; i < wave.length; i++) {
		if (waveform->phi->data->data[i] != 0.) {
			break;
		}
		if (i == wave.length - 1) {
			XLALSQTPNDestroyCoherentGW(waveform);
			DETATCHSTATUSPTR(status);
			RETURN(status);
		}
	}

	{
		if (waveform->a != NULL) {
			waveform->f->data->length = waveform->phi->data->length = waveform->shift->data->length
					= wave.length;
			waveform->a->data->length = 2 * wave.length;
			waveform->a->deltaT = waveform->f->deltaT = waveform->phi->deltaT
					= waveform->shift->deltaT = 1. / params->tSampling;

			waveform->a->sampleUnits = lalStrainUnit;
			waveform->f->sampleUnits = lalHertzUnit;
			waveform->phi->sampleUnits = lalDimensionlessUnit;
			waveform->shift->sampleUnits = lalDimensionlessUnit;

			waveform->position = ppnParams->position;
			waveform->psi = ppnParams->psi;

			snprintf(waveform->a->name, LALNameLength, "STPN inspiral amplitudes");
			snprintf(waveform->f->name, LALNameLength, "STPN inspiral frequency");
			snprintf(waveform->phi->name, LALNameLength, "STPN inspiral phase");
			snprintf(waveform->shift->name, LALNameLength, "STPN inspiral polshift");
		}
		// --- fill some output ---
		ppnParams->tc = (REAL8)(wave.length - 1) / params->tSampling;
		ppnParams->length = wave.length;
		ppnParams->dfdt = ((REAL4)(waveform->f->data->data[wave.length - 1]
				- waveform->f->data->data[wave.length - 2])) * ppnParams->deltaT;
		ppnParams->fStop = params->fFinal;
		ppnParams->termCode = GENERATEPPNINSPIRALH_EFSTOP;
		ppnParams->termDescription = GENERATEPPNINSPIRALH_MSGEFSTOP;

		ppnParams->fStart = ppnParams->fStartIn;
	} // end phase condition
	DETATCHSTATUSPTR(status);
	RETURN(status);
}

int XLALSQTPNAllocateCoherentGW(CoherentGW *wave, UINT4 length) {
	if (!wave) {
		XLAL_ERROR(__func__, XLAL_EFAULT);
	}
	if (length <= 0) {
		XLAL_ERROR(__func__, XLAL_EBADLEN);
	}
	if (wave->a || wave->f || wave->phi || wave->shift) {
		XLAL_ERROR(__func__, XLAL_EFAULT);
	}
	wave->a = (REAL4TimeVectorSeries *)LALMalloc(sizeof(REAL4TimeVectorSeries));
	wave->f = (REAL4TimeSeries *)LALMalloc(sizeof(REAL4TimeSeries));
	wave->phi = (REAL8TimeSeries *)LALMalloc(sizeof(REAL8TimeSeries));
	wave->shift = (REAL4TimeSeries *)LALMalloc(sizeof(REAL4TimeSeries));
	if (!(wave->a && wave->f && wave->phi && wave->shift)) {
		XLALSQTPNDestroyCoherentGW(wave);
		XLAL_ERROR(__func__, XLAL_ENOMEM);
	}
	xlalErrno = 0;
	wave->a->data = XLALCreateREAL4VectorSequence(length, 2);
	wave->f->data = XLALCreateREAL4Vector(length);
	wave->phi->data = XLALCreateREAL8Vector(length);
	wave->shift->data = XLALCreateREAL4Vector(length);
	if (!(wave->a->data && wave->f->data && wave->phi->data && wave->shift->data)) {
		XLALSQTPNDestroyCoherentGW(wave);
		XLAL_ERROR(__func__, XLAL_ENOMEM);
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
