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
 * The amplitudes and QNM frequen are calculated
 * according to Gossan et. al 2012 [arXiv: 1111.5819]
 *
 * \note The dimensionless spin assumes values between -1 and 1.
 *
 * \todo Extend so that overtones can be computed too.
 */
int XLALSimBlackHoleRingdownTiger(
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


/**
 * Calculates the amplitudes for the QNM l, m, n=0 for
 * a given symmetric mass-ratio eta of the initial binary.
 * Based on an interpolation for NON-SPINNING binary black
 * holes, derived in Gossan et al (2012) [arXiv: 1111.5819]
 * 
 * TO DO: rewrite this properly, add initial (effective) spin dependence
 **/
REAL8 XLALSimRingdownQNMAmplitudes(INT4 l, INT4 m, REAL8 eta, REAL8 spin1[3], REAL8 spin2[3]){
	REAL8 A = 0.8639*eta;
	if (l==2 && m==2){ A *= 1.0 }
	else if (l==2 && abs(m)==1){ A *= 0.52*pow(1.0 - 4.*eta, 0.71); }	// - 0.23 * chiEff
	else if (l==3 && abs(m)==3){ A *= 0.44*pow(1.0 - 4.*eta, 0.45); }
	else if (l==3 && abs(m)==2){ A *= 3.69*(eta - 0.2)*(eta - 0.2) + 0.053; }
	else if (l==4 && abs(m)==4){ A *= 5.41*((eta - 0.22)*(eta - 0.22) + 0.04); }
	else A = 0.0;
	return A;
	}
	
/**
 * Computes the complex dimensionless eigen-frequency for the
 * quasinormal mode (l,m), given a dimensionless spin a \in [-1,1].
 * This is based on the interpolation formula 
 *
 * \note The dimensionless spin assumes values between -1 and 1.
 *
 * \todo Extend so that overtones can be computed too.
 */
COMPLEX16 XLALSimRingdownInterpOmega(UINT4 l, INT4 m, UINT4 n, REAL8 a){
	COMPLEX16 omega=0.0;
	if (abs(a) > 1.){
		fprintf(stderr, "ERROR: Dimensionless spin magnitude larger than 1! Aborting...\n");
		exit(-1);
	}
	switch (l){
	case 2:
		switch (abs(m)){
		case 2:
			omega = (1.5251 - 1.1568*pow(1.e0-a,0.1292));
			omega += I*(2*(0.7000 + 1.4187*pow(1.e0-a,-0.4990)))/omega;
			break;
		case 1:
			omega = (0.6000 - 0.2339*pow(1.e0-a,0.4175));
			omega += I*(2*(-0.3000 + 2.3561*pow(1.e0-a,-0.2277)))/omega;
			break;
		default:
			break;
		}
		break;
	case 3:
		switch (abs(m)){
		case 3:
			omega = (1.8956 - 1.3043*pow(1.e0-a,0.1818));
			omega += I*(2*(0.9000 + 2.3430*pow(1.e0-a,-0.4810)))/omega;
			break;
		case 2:
			omega = (1.1481 - 0.5552*pow(1.e0-a,0.3002));
			omega += I*(2*(0.8313 + 2.3773*pow(1.e0-a,-0.3655)))/omega;
			break;
		default:
			break;
		}
		break;
	case 4:
		switch (abs(m)){
		case 4:
			omega = (2.3000 - 1.5056*pow(1.e0-a,0.2244));
			omega += I*(2*(1.1929 + 3.1191*pow(1.e0-a,-0.4825)))/omega;
			break;
		default:
			break;
		}
		break;
	default:
		break;
	}
	if (omega = 0.0){
		fprintf(stderr, "ERROR: Invalid mode l = %u, m = %d, n = %u\n", l, m, n);
		exit(-1);
	}
	return omega;
}

/** 
 * Frequency in rad/sec given dimensionless complex frequency and M
 **/
REAL8 XLALQNMFreqOfOmega(COMPLEX16 omega, REAL8 mtot){
	return (creal(omega)/mtot);
}

/** 
 * Damping time in sec given dimensionless complex frequency and M
 **/
REAL8 XLALQNMTauOfOmega(COMPLEX16 omega, REAL8 mtot){
	REAL8 tau = 0;
	tau = mtot/cimag(omega);
	return tau;
}

