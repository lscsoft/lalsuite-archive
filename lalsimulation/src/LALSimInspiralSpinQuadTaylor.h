/**	@file   LALSimInspiralSpinQuadTaylor.h
 *	@author László Veréb
 *	@date   21 Feb 2012
 *	@brief  Interfaceses to generate the SpinQuadTaylor waveform.
 */

#ifndef LALSIMINSPIRALSPINQUADTAYLOR_H_
#define LALSIMINSPIRALSPINQUADTAYLOR_H_

#include <lal/LALSimInspiral.h>
#include <lal/TimeSeries.h>

/**	Computes the \f$h_+\f$ and \f$h_\cross\f$ polarized waveforms. The spin parameter is
 * \f$\chi_i=S_i/m_i^2\f$.
 *
 * @param[out] hp				: \f$h_+\f$ polarized waveform
 * @param[out] hc				: \f$h_\cross\f$ polarized waveform
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
 * @param[in]  orderOfAmplitude	: twice amplitude post-Newtonian order
 * @param[in]  interactionFlags : flag to control spin effects
 * @return
 */
int XLALSimInspiralSpinQuadTaylorEvolveWaveform(REAL8TimeSeries **hp, REAL8TimeSeries **hc,
		REAL8 mass1, REAL8 mass2, REAL8 qm1, REAL8 qm2, REAL8 chi1x, REAL8 chi1y, REAL8 chi1z,
		REAL8 chi2x, REAL8 chi2y, REAL8 chi2z, REAL8 lnhatx, REAL8 lnhaty, REAL8 lnhatz, REAL8 e1x,
		REAL8 e1y, REAL8 e1z, REAL8 endPhase, REAL8 initialFrequency, REAL8 samplingTime,
		INT4 orderOfPhase, INT4 orderOfAmplitude, LALSimInspiralInteraction interactionFlags);

/**	Computes the orbital equations. The spin parameter is \f$\chi_i=S_i/m_i^2\f$.
 *
 * @param[out] V				: post-Newtonian parameter
 * @param[out] Phi				: orbital phase
 * @param[out] S1x				: first spin vector x component
 * @param[out] S1y				: first spin vector y component
 * @param[out] S1z				: first spin vector z component
 * @param[out] S2x				: second spin vector x component
 * @param[out] S2y				: second spin vector y component
 * @param[out] S2z				: second spin vector z component
 * @param[out] LNhatx			: unit orbital angular momentum x component
 * @param[out] LNhaty			: unit orbital angular momentum y component
 * @param[out] LNhatz			: unit orbital angular momentum z component
 * @param[out] E1x				: orbital plane basis vector x component
 * @param[out] E1y				: orbital plane basis vector y component
 * @param[out] E1z				: orbital plane basis vector z component
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
int XLALSimInspiralSpinQuadTaylorEvolveOrbit(REAL8TimeSeries **V, REAL8TimeSeries **Phi,
		REAL8TimeSeries **S1x, REAL8TimeSeries **S1y, REAL8TimeSeries **S1z, REAL8TimeSeries **S2x,
		REAL8TimeSeries **S2y, REAL8TimeSeries **S2z, REAL8TimeSeries **LNhatx,
		REAL8TimeSeries **LNhaty, REAL8TimeSeries **LNhatz, REAL8TimeSeries **E1x,
		REAL8TimeSeries **E1y, REAL8TimeSeries **E1z, REAL8 mass1, REAL8 mass2, REAL8 qm1,
		REAL8 qm2, REAL8 chi1x, REAL8 chi1y, REAL8 chi1z, REAL8 chi2x, REAL8 chi2y, REAL8 chi2z,
		REAL8 lnhatx, REAL8 lnhaty, REAL8 lnhatz, REAL8 e1x, REAL8 e1y, REAL8 e1z, REAL8 endPhase,
		REAL8 initialFrequency, REAL8 samplingTime, INT4 orderOfPhase,
		LALSimInspiralInteraction interactionFlags);

/** Computes the \f$h_+\f$ and \f$h_\cross\f$ polarized waveforms with the orbital equations.
 * The spin parameter is \f$\chi_i=S_i/m_i^2\f$.
 *
 * @param[out] hp				: \f$h_+\f$ polarized waveform
 * @param[out] hc				: \f$h_\cross\f$ polarized waveform
 * @param[out] V				: post-Newtonian parameter
 * @param[out] Phi				: orbital phase
 * @param[out] S1x				: first spin vector x component
 * @param[out] S1y				: first spin vector y component
 * @param[out] S1z				: first spin vector z component
 * @param[out] S2x				: second spin vector x component
 * @param[out] S2y				: second spin vector y component
 * @param[out] S2z				: second spin vector z component
 * @param[out] LNhatx			: unit orbital angular momentum x component
 * @param[out] LNhaty			: unit orbital angular momentum y component
 * @param[out] LNhatz			: unit orbital angular momentum z component
 * @param[out] E1x				: orbital plane basis vector x component
 * @param[out] E1y				: orbital plane basis vector y component
 * @param[out] E1z				: orbital plane basis vector z component
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
 * @param[in]  orderOfAmplitude	: twice amplitude post-Newtonian order
 * @param[in]  interactionFlags : flag to control spin effects
 * @return
 */
int XLALSimInspiralSpinQuadTaylorEvolveAll(REAL8TimeSeries **hp, REAL8TimeSeries **hc,
		REAL8TimeSeries **V, REAL8TimeSeries **Phi, REAL8TimeSeries **S1x, REAL8TimeSeries **S1y,
		REAL8TimeSeries **S1z, REAL8TimeSeries **S2x, REAL8TimeSeries **S2y, REAL8TimeSeries **S2z,
		REAL8TimeSeries **LNhatx, REAL8TimeSeries **LNhaty, REAL8TimeSeries **LNhatz,
		REAL8TimeSeries **E1x, REAL8TimeSeries **E1y, REAL8TimeSeries **E1z, REAL8 mass1,
		REAL8 mass2, REAL8 qm1, REAL8 qm2, REAL8 chi1x, REAL8 chi1y, REAL8 chi1z, REAL8 chi2x,
		REAL8 chi2y, REAL8 chi2z, REAL8 lnhatx, REAL8 lnhaty, REAL8 lnhatz, REAL8 e1x, REAL8 e1y,
		REAL8 e1z, REAL8 endPhase, REAL8 initialFrequency, REAL8 samplingTime, INT4 orderOfPhase,
		INT4 orderOfAmplitude, LALSimInspiralInteraction interactionFlags);

#endif /* LALSIMINSPIRALSPINQUADTAYLOR_H_ */
