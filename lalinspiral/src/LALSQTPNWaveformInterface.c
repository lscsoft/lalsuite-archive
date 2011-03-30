/**
 * @file LALSQTPNWaveformInterface.c
 *		Contains function definitions to integrate the SpinQuadTaylor code into the other parts of the LALSuit.
 *	If you want to run the program use the LALSQTPNWaveformTest.c file int the
 *	test directory.
 * @author László Veréb
 * @date 2010.06.27.
 */

#include <lal/LALSQTPNWaveformInterface.h>
#include <lal/LALSQTPNWaveform.h>

NRCSID (LALSQTPNWAVEFORMINTERFACEC, "$Id LALSQTPNWaveformInterface.c$");

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
	LALSQTPNGenerator_Old(status->statusPtr, &wave, &wave_Params);
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
	LALSQTPNGenerator_Old(status->statusPtr, &wave, &wave_Params);
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
	LALSQTPNGenerator_Old(status->statusPtr, &wave, &wave_Params);
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
			for (i = 0; i < wave.length; i++) {
				// (PPNParamStruct)ppnParams->phi === (InspiralTemplate)params->startPhase === (SimInspiralTable)injparams/this_event->coa_phase it is set to 0 in LALSQTPNWaveformTest.c at line 83.
				waveform->phi->data->data[i] = waveform->phi->data->data[i]
				      - waveform->phi->data->data[wave.length - 1] + ppnParams->phi;
			}
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
	static const char *func = "LALSQTPNAllocateCoherentGW";
	if (!wave) {
		XLAL_ERROR(func, XLAL_EFAULT);
	}
	if (length <= 0) {
		XLAL_ERROR(func, XLAL_EBADLEN);
	}
	if (wave->a || wave->f || wave->phi || wave->shift) {
		XLAL_ERROR(func, XLAL_EFAULT);
	}
	wave->h = (REAL4TimeVectorSeries *)LALMalloc(sizeof(REAL4TimeVectorSeries));
	wave->a = (REAL4TimeVectorSeries *)LALMalloc(sizeof(REAL4TimeVectorSeries));
	wave->f = (REAL4TimeSeries *)LALMalloc(sizeof(REAL4TimeSeries));
	wave->phi = (REAL8TimeSeries *)LALMalloc(sizeof(REAL8TimeSeries));
	wave->shift = (REAL4TimeSeries *)LALMalloc(sizeof(REAL4TimeSeries));
	if (!(wave->h && wave->a && wave->f && wave->phi && wave->shift)) {
		XLALSQTPNDestroyCoherentGW(wave);
		XLAL_ERROR(func, XLAL_ENOMEM);
	}
	xlalErrno = 0;
	wave->h->data = XLALCreateREAL4VectorSequence(length, 2);
	wave->a->data = XLALCreateREAL4VectorSequence(length, 2);
	wave->f->data = XLALCreateREAL4Vector(length);
	wave->phi->data = XLALCreateREAL8Vector(length);
	wave->shift->data = XLALCreateREAL4Vector(length);
	if (!(wave->h->data && wave->a->data && wave->f->data && wave->phi->data && wave->shift->data)) {
		XLALSQTPNDestroyCoherentGW(wave);
		XLAL_ERROR(func, XLAL_ENOMEM);
	}
	return XLAL_SUCCESS;
}

void XLALSQTPNDestroyCoherentGW(CoherentGW *wave) {
	//static const char *func = "LALSQTPNDestroyCoherentGW";
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

void XLALSQTPNFillParams(LALSQTPNWaveformParams *wave, InspiralTemplate *params) {
	wave->mass[0] = params->mass1;
	wave->mass[1] = params->mass2;
	wave->totalMass = wave->mass[0] + wave->mass[1];
	wave->mu = wave->mass[0] * wave->mass[1] / wave->totalMass;
	wave->eta = wave->mu / wave->totalMass;
	wave->chirpMass = wave->totalMass * pow(wave->eta, 3. / 5.);
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

void XLALSQTPNCalcHpHcFromCoherentGWat(INT4 ind, INT2 amplitude,CoherentGW *wave, REAL8 mass[], REAL8 distance, REAL4 h[]) {
	REAL8 contribution[3][2];
	h[LALSQTPN_PLUS] = h[LALSQTPN_CROSS] = 0.0;
	LALSQTPNWaveformParams params;
	params.mass[0] = mass[0];
	params.mass[1] = mass[1];
	params.totalMass = params.mass[0] + params.mass[1];
	params.eta = params.mass[0] * params.mass[1] / (params.totalMass*params.totalMass);
	params.distance = distance * LAL_PC_SI * 1e6;
	REAL8 values[LALSQTPN_NUM_OF_VAR];
	const REAL8 freq_Step = params.totalMass * LAL_MTSUN_SI * LAL_PI;
	values[LALSQTPN_OMEGA] = wave->f->data->data[ind] * freq_Step;
	REAL8 amp1 = 4. * params.totalMass * params.eta * LAL_MRSUN_SI / params.distance * pow(values[LALSQTPN_OMEGA], 2.0 / 3.0);
	values[LALSQTPN_PHASE] = (wave->phi->data->data[ind]-wave->phi->data->data[0]) / 2.0;
	values[LALSQTPN_LNH_3] = -wave->a->data->data[2*ind+1] / amp1;
	REAL8 twoAlpha = wave->shift->data->data[ind];
	INT2 i;
	REAL8 amp2 = -2.0 * pow(values[LALSQTPN_OMEGA], 2.0 / 3.0) * params.totalMass * params.eta
			* LAL_MRSUN_SI / params.distance;
	if ((amplitude & LALSQTPN_1_0) == LALSQTPN_1_0) {
		XLALSQTPNCalculateAmplitudeContribution1_0(&params, values, twoAlpha, contribution[LAL_PNORDER_ONE]);
		for (i = LALSQTPN_PLUS; i <= LALSQTPN_CROSS; i++) {
			h[i] += contribution[LAL_PNORDER_ONE][i] * pow(values[LALSQTPN_OMEGA], 2.0 / 3.0);
		}
	}
	if ((amplitude & LALSQTPN_0_5) == LALSQTPN_0_5) {
		XLALSQTPNCalculateAmplitudeContribution0_5(&params, values, twoAlpha, contribution[LAL_PNORDER_HALF]);
		for (i = LALSQTPN_PLUS; i <= LALSQTPN_CROSS; i++) {
			h[i] += contribution[LAL_PNORDER_HALF][i] * pow(values[LALSQTPN_OMEGA], 1.0 / 3.0);
		}
	}
	if ((amplitude & LALSQTPN_0_0) == LALSQTPN_0_0) {
		XLALSQTPNCalculateAmplitudeContribution0_0(values, twoAlpha, contribution[LAL_PNORDER_NEWTONIAN]);
		for (i = LALSQTPN_PLUS; i <= LALSQTPN_CROSS; i++) {
			h[i] += contribution[LAL_PNORDER_NEWTONIAN][i];
		}
	}
	for (i = LALSQTPN_PLUS; i <= LALSQTPN_CROSS; i++) {
		h[i] *= amp2;
	}
}
