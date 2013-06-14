#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_integration.h>

#include <lal/LALDatatypes.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/TimeSeries.h>
#include <lal/Units.h>
#include <lal/SphericalHarmonics.h>

#include <lal/LALSimBlackHoleRingdownTiger.h>
#define EPS LAL_REAL4_EPS
#define TINY LAL_REAL8_MIN
#define MAXITER 16384

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

// TODO: Add a top-level function to get combination of all modes


/**
 * Computes the waveform for the ringdown of a black hole
 * quasinormal mode (l,m).
 * The amplitudes and QNM frequen are taken
 * according to Gossan et. al [arXiv: 1111.5819]
 *
 * \note The dimensionless spin assumes values between -1 and 1.
 *
 * \todo Extend so that overtones can be computed too.
 */
int XLALSimBlackHoleRingdown(
	REAL8TimeSeries **hplus,	/**< plus-polarization waveform [returned] */
	REAL8TimeSeries **hcross,	/**< cross-polarization waveform [returned] */
	const LIGOTimeGPS *t0,		/**< start time of ringdown */
	double phi0,			/**< initial phase of ringdown (rad) */
	double deltaT,			/**< sampling interval (s) */
	double mass,			/**< black hole mass (kg) */
	double a,	/**< black hole dimensionless spin parameter */
//	double fractional_mass_loss,	/**< fraction of mass radiated in this mode */
	double distance,		/**< distance to source (m) */
	double inclination,		/**< inclination of source's spin axis (rad) */
	int l,				/**< polar mode number */
	int m				/**< azimuthal mode number */
)
{	
	const int s = -2; /* spin weight for gravitational radiation */
	double cosi = cos(inclination);
	COMPLEX16 A, omega, Q, tau;
	COMPLEX16 Yplus, Ycross;
	size_t length;
	size_t j;
	int errnum;

	

/*
	if (XLALSimBlackHoleRingdownModeEigenvaluesLeaver(&A, &omega, a, l, m, s) < 0)
		XLAL_ERROR(XLAL_EFUNC);
	XLAL_TRY(sphwf1 = XLALSimBlackHoleRingdownSpheroidalWaveFunctionLeaver(mu, a, l, m, s, A, omega), errnum);
	XLAL_TRY(sphwf2 = XLALSimBlackHoleRingdownSpheroidalWaveFunctionLeaver(-mu, a, l, m, s, A, omega), errnum);
	if (errnum)
		XLAL_ERROR(XLAL_EFUNC);
	omega *= 0.5; /* convert from Leaver's convention 2M = 1 */
	
	/* compute length of waveform to compute */
/*	length = ceil(log(LAL_REAL8_EPS) * LAL_G_SI * mass / (pow(LAL_C_SI, 3.0) * cimag(omega) * deltaT));
	if (length < 1)
		XLAL_ERROR(XLAL_EBADLEN);
*/
	/* compute the amplitude factors for the +m and -m modes */
/*	A1 = A2 = -4.0 * (LAL_G_SI*mass/(pow(LAL_C_SI, 2.0)*distance))
		* sqrt(-0.5*cimag(omega)*fractional_mass_loss)/cabs(omega)
 		* cexp(I*m*phi0);
	A1 *= sphwf1;
	A2 = conj(A2*sphwf2);

	omega_dt = pow(LAL_C_SI, 3.0)*omega*deltaT/(LAL_G_SI*mass);

	*hplus = XLALCreateREAL8TimeSeries("H_PLUS", t0, 0.0, deltaT, &lalStrainUnit, length);
	*hcross = XLALCreateREAL8TimeSeries("H_CROSS", t0, 0.0, deltaT, &lalStrainUnit, length);
	if (hplus == NULL || hcross == NULL)
		XLAL_ERROR(XLAL_EFUNC);
*/

	/* compute the waveforms */
	for (j = 0; j < length; ++j) {
		COMPLEX16 h;
		h = A1*cexp(-I*omega_dt*j) + A2*cexp(I*conj(omega_dt)*j);
		(*hplus)->data->data[j] = creal(h);
		(*hcross)->data->data[j] = -cimag(h);
	}

	return 0;
}

REAL8 XLALSimSphericalHarmonicPlus(int l, int m, float iota){
	COMPLEX16 Yplus = 0.0;
	Yplus = XLALSpinWeightedSphericalHarmonic(iota, 0.0, -2, l, m) + (2.0*(l % 2) - 1.0)*XLALSpinWeightedSphericalHarmonics(iota, 0.0, -2, l, -m);
	return Yplus;
}

REAL8 XLALSimSphericalHarmonicCross(int l, int m, float iota){
	COMPLEX16 Ycross = 0.0;
	Ycross = XLALSpinWeightedSphericalHarmonic(iota, 0.0, -2, l, m) - (2.0*(l % 2) - 1.0)*XLALSpinWeightedSphericalHarmonics(iota, 0.0, -2, l, -m);
	return Ycross;
}

REAL8 XLALSimRingdownModeAmplitudes(INT4 l, INT4 m, REAL8 eta){
	REAL8 A = 0.864*eta;
	if (l==2 && m==2){ A *= 1.0 }
	else if (l==2 && m==1){ A *= 0.52*pow(1.0 - 4*eta, 0.71) }
	else if (l==3 && m==3){ A *= 0.44*pow(1.0 - 4*eta, 0.45) }
	else if (l==4 && m==4){ A *= 5.4*((eta - 0.22)*(eta - 0.22) + 0.04)}
	else A = 0.0;
	return A;
	}

