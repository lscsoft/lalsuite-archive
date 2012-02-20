/**
 * @file LALSQTPNWaveformInterface.h
 *		Contains function declarations to integrate the SpinQuadTaylor code into the other parts of the LALSuit.
 * @author László Veréb
 * @date 2011.07.11.
 */

#ifndef LALSQTPNWAVEFORMINTERFACE_H
#define LALSQTPNWAVEFORMINTERFACE_H

#include <lal/Units.h>
#include <lal/LALInspiral.h>
#include <lal/SeqFactories.h>
#include <lal/LALStatusMacros.h>
#include <lal/GenerateInspiral.h>
#include <lal/LALSQTPNWaveform.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

#define LALSQTPN_MSGPPNPARAMS "the PPNParamsStruct structure is null"
#define LALSQTPN_MSGINSPIRALTEMPLATE "the InspiralTemplate structure is null"
#define LALSQTPN_ZEROLENGTH 3
#define LALSQTPN_MSGZEROLENGTH "the given length is not positive"

typedef enum tagLALSQTPNSwitchMode {
	LALSQTPN_FIXED, LALSQTPN_PRECESSING,
} LALSQTPNSwitchMode;

int XLALSQTPNWaveformTemplates(REAL4Vector *signalvec1, REAL4Vector *signalvec2, InspiralTemplate *params);

// LAL wrapper to XLAL function
void LALSQTPNWaveformTemplates(LALStatus *status, REAL4Vector *signalvec1, REAL4Vector *signalvec2, InspiralTemplate *params);

/**		The function returns the generated waveform.
 * @param[in,out]	status	: LAL universal status structure
 * @param[out]	signalvec	: array containing the waveform \f$(h_+, h_\times)\f$
 * @param[in]	params		: structure containing the inspiral parameters
 */
int XLALSQTPNWaveform (REAL4Vector *signalvec, InspiralTemplate *params);

// LAL wrapper to XLAL function
void LALSQTPNWaveform (LALStatus *status, REAL4Vector *signalvec, InspiralTemplate *params);

/**		The function returns the generated waveform for injection.
 * @param[in,out]	status	: LAL universal status structure
 * @param[out]	wave_out	: structure containing the waveform \f$(a_1,
 * a_2, \Phi, \alpha)\f$
 * @param[in]	params		: structure containing the inspiral parameters
 * @param[in]	ppnParams	: parameters for restricted post-Newtonian waveform
 */
int XLALSQTPNWaveformForInjection(CoherentGW *wave_out, InspiralTemplate *params, PPNParamStruc *ppnParams);

// LAL wrapper to XLAL function
void LALSQTPNWaveformForInjection(LALStatus *status, CoherentGW *wave_out,
		InspiralTemplate *params, PPNParamStruc *ppnParams);

/**		The function allocates memory for the waveform's \f$h_+\f$, \f$h_\times\f$.
 * @param[out] wave	  : pointer to the allocated waveform
 * @param[in]  length : the length of the waveform
 * @return
 */
int XLALSQTPNAllocateHplusAndHcrossInCoherentGW(CoherentGW *wave, UINT4 length);

/**		The function allocates memory for the waveform's \f$a_1\f$, \f$a_2\f$,
 * \f$\Phi\f$ and \f$\alpha\f$.
 * @param[out]		waveform	: pointer to the allocated waveform
 * @param[in]		length		: the length of the waveform
 */
int XLALSQTPNAllocateAmplitudeInCoherentGW(CoherentGW *waveform, UINT4 length);

/**		The function deallocates memory of the waveform.
 * @param[out]		waveform	: pointer to the allocated waveform
 */
void XLALSQTPNDestroyCoherentGW(CoherentGW *waveform);

#ifdef __cplusplus
}
#endif

#endif /* LALSQTPNWAVEFOMRINTERFACE_H */
