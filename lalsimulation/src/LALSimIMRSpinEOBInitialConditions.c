#include <lal/LALSimInspiral.h>
#include <lal/LALSimIMR.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_deriv.h>

#include "LALSimIMRSpinEOB.h"
#include "LALSimIMRSpinEOBHcapNumericalDerivative.c"
#include "LALSimIMRSpinEOBHcapNumericalDerivativePrec.c"
#include "LALSimIMRSpinEOBHamiltonian.c"
#include "LALSimIMRSpinEOBHamiltonianPrec.c"
#include "LALSimIMREOBFactorizedWaveform.c"
#include "LALSimIMRSpinEOBInitialConditionsPrec.c"

#ifndef _LALSIMIMRSPINEOBINITIALCONDITIONS_C
#define _LALSIMIMRSPINEOBINITIALCONDITIONS_C

static int
XLALFindSphericalOrbitV2(const gsl_vector * x,	/**<< Parameters requested by gsl root finder */
			 void *params,	/**<< Spin EOB parameters */
			 gsl_vector * f	/**<< Function values for the given parameters */
)
{
	SEOBRootParams *rootParams = (SEOBRootParams *) params;

	REAL8		py      , pz, r, ptheta, pphi;

	/* Numerical derivative of Hamiltonian wrt given value */
	REAL8		dHdx    , dHdpy, dHdpz;
	REAL8		dHdr    , dHdptheta, dHdpphi;

	/* Populate the appropriate values */
	/* In the special theta=pi/2 phi=0 case, r is x */
	rootParams->values[0] = r = gsl_vector_get(x, 0);
	rootParams->values[4] = py = gsl_vector_get(x, 1);
	rootParams->values[5] = pz = gsl_vector_get(x, 2);

	//printf("Values r = %.16e, py = %.16e, pz = %.16e\n", r, py, pz);

	ptheta = -r * pz;
	pphi = r * py;

	/* dHdR */
	dHdx = XLALSpinPrecHcapNumDerivWRTParam(0, rootParams->values, rootParams->params);
	if (XLAL_IS_REAL8_FAIL_NAN(dHdx)) {
		XLAL_ERROR(XLAL_EFUNC);
	}
	//printf("dHdx = %.16e\n", dHdx);

	/* dHdPphi (I think we can use dHdPy in this coord system) */
	/* TODO: Check this is okay */
	dHdpy = XLALSpinPrecHcapNumDerivWRTParam(4, rootParams->values, rootParams->params);
	if (XLAL_IS_REAL8_FAIL_NAN(dHdpy)) {
		XLAL_ERROR(XLAL_EFUNC);
	}
	/* dHdPtheta (I think we can use dHdPz in this coord system) */
	/* TODO: Check this is okay */
	dHdpz = XLALSpinPrecHcapNumDerivWRTParam(5, rootParams->values, rootParams->params);
	if (XLAL_IS_REAL8_FAIL_NAN(dHdpz)) {
		XLAL_ERROR(XLAL_EFUNC);
	}
	/* Now convert to spherical polars */
	dHdr = dHdx - dHdpy * pphi / (r * r) + dHdpz * ptheta / (r * r);
	dHdptheta = -dHdpz / r;
	dHdpphi = dHdpy / r;

	/* populate the function vector */
	gsl_vector_set(f, 0, dHdr);
	gsl_vector_set(f, 1, dHdptheta);
	gsl_vector_set(f, 2, dHdpphi - rootParams->omega);

	//printf("Current funcvals = %.16e %.16e %.16e\n", gsl_vector_get(f, 0), gsl_vector_get(f, 1),
		 //gsl_vector_get(f, 2) /* dHdpphi */ );

	return XLAL_SUCCESS;
}


/**
 * Main function for calculating the spinning EOB initial conditions, following the
 * quasi-spherical initial conditions described in Sec. IV A of
 * Buonanno, Chen & Damour PRD 74, 104005 (2006).
 * All equation numbers in the comments of this function refer to this paper.
 * 
 * Inputs are binary parameters (masses, spin vectors and inclination),
 * EOB model parameters and initial frequency.
 * 
 * Outputs are initial dynamical variables in a vector
 * (x, y, z, px, py, pz, S1x, S1y, S1z, S2x, S2y, S2z).
 * 
 * In the paper, the initial conditions are solved for a given initial radius,
 * while in this function, they are solved for a given inital orbital frequency.
 * 
 * This difference is reflected in solving Eq. (4.8).
 * The calculation is done in 5 steps:
 * STEP 1) Rotate to LNhat0 along z-axis and N0 along x-axis frame, where
 * LNhat0 and N0 are initial normal to orbital plane and initial orbital separation;
 * STEP 2) After rotation in STEP 1, in spherical coordinates,
 * phi0 and theta0 are given directly in Eq. (4.7),
 * r0, pr0, ptheta0 and pphi0 are obtained by solving Eqs. (4.8) and (4.9)
 * (using gsl_multiroot_fsolver).
 * At this step, we find initial conditions for a spherical orbit without
 * radiation reaction.
 * STEP 3) Rotate to L0 along z-axis and N0 along x-axis frame, where
 * L0 is the initial orbital angular momentum and
 * L0 is calculated using initial position and linear momentum obtained in STEP 2.
 * STEP 4) In the L0-N0 frame, we calculate (dE/dr)|sph using Eq. (4.14),
 * then initial dr/dt using Eq. (4.10), and finally pr0 using Eq. (4.15).
 * STEP 5) Rotate back to the original inertial frame by inverting the rotation of STEP 3 and
 * then inverting the rotation of STEP 1.
 */

static int UNUSED 
XLALSimIMRSpinEOBInitialConditionsV2(
				     REAL8Vector * initConds,	/**<< OUTPUT, Initial dynamical variables */
				     const REAL8 mass1,	/**<< mass 1 */
				     const REAL8 mass2,	/**<< mass 2 */
				     const REAL8 fMin,	/**<< Initial frequency (given) */
				     const REAL8 inc,	/**<< Inclination */
				     const REAL8 spin1[],	/**<< Initial spin vector 1 */
				     const REAL8 spin2[],	/**<< Initial spin vector 2 */
				     SpinEOBParams * params	/**<< Spin EOB parameters */
)
{

#ifndef LAL_NDEBUG
	if (!initConds) {
		XLAL_ERROR(XLAL_EINVAL);
	}
#endif


	static const int UNUSED lMax = 8;

	int		i;

	/* Variable to keep track of whether the user requested the tortoise */
	int		tmpTortoise;

	UINT4		SpinAlignedEOBversion;

	REAL8		mTotal;
	REAL8		eta;
	REAL8		omega   , v0;	/* Initial velocity and angular
					 * frequency */

	REAL8		ham;	/* Hamiltonian */

	REAL8		LnHat    [3];	/* Initial orientation of angular
					 * momentum */
	REAL8		rHat     [3];	/* Initial orientation of radial
					 * vector */
	REAL8		vHat     [3];	/* Initial orientation of velocity
					 * vector */
	REAL8		Lhat     [3];	/* Direction of relativistic ang mom */
	REAL8		qHat     [3];
	REAL8		pHat     [3];

	/* q and p vectors in Cartesian and spherical coords */
	REAL8		qCart    [3], pCart[3];
	REAL8		qSph     [3], pSph[3];

	/* We will need to manipulate the spin vectors */
	/* We will use temporary vectors to do this */
	REAL8		tmpS1    [3];
	REAL8		tmpS2    [3];
	REAL8		tmpS1Norm[3];
	REAL8		tmpS2Norm[3];

	REAL8Vector	qCartVec, pCartVec;
	REAL8Vector	s1Vec, s2Vec, s1VecNorm, s2VecNorm;
	REAL8Vector	sKerr, sStar;
	REAL8		sKerrData[3], sStarData[3];
	REAL8		a = 0.;
	//, chiS, chiA;
	//REAL8 chi1, chi2;

	/*
	 * We will need a full values vector for calculating derivs of
	 * Hamiltonian
	 */
	REAL8		sphValues[12];
	REAL8		cartValues[12];

	/* Matrices for rotating to the new basis set. */
	/* It is more convenient to calculate the ICs in a simpler basis */
	gsl_matrix     *rotMatrix = NULL;
	gsl_matrix     *invMatrix = NULL;
	gsl_matrix     *rotMatrix2 = NULL;
	gsl_matrix     *invMatrix2 = NULL;

	/* Root finding stuff for finding the spherical orbit */
	SEOBRootParams	rootParams;
	const gsl_multiroot_fsolver_type *T = gsl_multiroot_fsolver_hybrids;
	gsl_multiroot_fsolver *rootSolver = NULL;

	gsl_multiroot_function rootFunction;
	gsl_vector     *initValues = NULL;
	gsl_vector     *finalValues = NULL;
	int		gslStatus;
	const int	maxIter = 100;

	memset(&rootParams, 0, sizeof(rootParams));

	mTotal = mass1 + mass2;
	eta = mass1 * mass2 / (mTotal * mTotal);
	memcpy(tmpS1, spin1, sizeof(tmpS1));
	memcpy(tmpS2, spin2, sizeof(tmpS2));
	memcpy(tmpS1Norm, spin1, sizeof(tmpS1Norm));
	memcpy(tmpS2Norm, spin2, sizeof(tmpS2Norm));
	for (i = 0; i < 3; i++) {
		tmpS1Norm[i] /= mTotal * mTotal;
		tmpS2Norm[i] /= mTotal * mTotal;
	}
	SpinAlignedEOBversion = params->seobCoeffs->SpinAlignedEOBversion;
	/* We compute the ICs for the non-tortoise p, and convert at the end */
	tmpTortoise = params->tortoise;
	params->tortoise = 0;

	EOBNonQCCoeffs *nqcCoeffs = NULL;
	nqcCoeffs = params->nqcCoeffs;

	/*
	 * STEP 1) Rotate to LNhat0 along z-axis and N0 along x-axis frame,
	 * where LNhat0 and N0 are initial normal to orbital plane and
	 * initial orbital separation;
	 */

	/* Set the initial orbital ang mom direction. Taken from STPN code */
	LnHat[0] = sin(inc);
	LnHat[1] = 0.;
	LnHat[2] = cos(inc);

	/*
	 * Set the radial direction - need to take care to avoid singularity
	 * if L is along z axis
	 */
	if (LnHat[2] > 0.9999) {
		rHat[0] = 1.;
		rHat[1] = rHat[2] = 0.;
	} else {
		REAL8		theta0 = atan(-LnHat[2] / LnHat[0]);	/* theta0 is between 0
									 * and Pi */
		rHat[0] = sin(theta0);
		rHat[1] = 0;
		rHat[2] = cos(theta0);
	}

	/* Now we can complete the triad */
	vHat[0] = CalculateCrossProductPrec(0, LnHat, rHat);
	vHat[1] = CalculateCrossProductPrec(1, LnHat, rHat);
	vHat[2] = CalculateCrossProductPrec(2, LnHat, rHat);

	NormalizeVectorPrec(vHat);

	/* XXX Test code XXX */
	/*
	 * for ( i = 0; i < 3; i++ ) { printf ( " LnHat[%d] = %.16e, rHat[%d]
	 * = %.16e, vHat[%d] = %.16e\n", i, LnHat[i], i, rHat[i], i, vHat[i]
	 * ); }
	 * 
	 * printf("\n\n" ); for ( i = 0; i < 3; i++ ) { printf ( " s1[%d] =
	 * %.16e, s2[%d] = %.16e\n", i, tmpS1[i], i, tmpS2[i] ); }
	 */

	/* Allocate and compute the rotation matrices */
	XLAL_CALLGSL(rotMatrix = gsl_matrix_alloc(3, 3));
	XLAL_CALLGSL(invMatrix = gsl_matrix_alloc(3, 3));
	if (!rotMatrix || !invMatrix) {
		if (rotMatrix)
			gsl_matrix_free(rotMatrix);
		if (invMatrix)
			gsl_matrix_free(invMatrix);
		XLAL_ERROR(XLAL_ENOMEM);
	}
	if (CalculateRotationMatrixPrec(rotMatrix, invMatrix, rHat, vHat, LnHat) == XLAL_FAILURE) {
		gsl_matrix_free(rotMatrix);
		gsl_matrix_free(invMatrix);
		XLAL_ERROR(XLAL_ENOMEM);
	}
	/* Rotate the orbital vectors and spins */
	ApplyRotationMatrixPrec(rotMatrix, rHat);
	ApplyRotationMatrixPrec(rotMatrix, vHat);
	ApplyRotationMatrixPrec(rotMatrix, LnHat);
	ApplyRotationMatrixPrec(rotMatrix, tmpS1);
	ApplyRotationMatrixPrec(rotMatrix, tmpS2);
	ApplyRotationMatrixPrec(rotMatrix, tmpS1Norm);
	ApplyRotationMatrixPrec(rotMatrix, tmpS2Norm);

	/* XXX Test code XXX */
	/*
	 * printf( "\nAfter applying rotation matrix:\n\n" ); for ( i = 0; i
	 * < 3; i++ ) { printf ( " LnHat[%d] = %.16e, rHat[%d] = %.16e,
	 * vHat[%d] = %.16e\n", i, LnHat[i], i, rHat[i], i, vHat[i] ); }
	 * 
	 * printf("\n\n" ); for ( i = 0; i < 3; i++ ) { printf ( " s1[%d] =
	 * %.16e, s2[%d] = %.16e\n", i, tmpS1[i], i, tmpS2[i] ); }
	 */

	/*
	 * STEP 2) After rotation in STEP 1, in spherical coordinates, phi0
	 * and theta0 are given directly in Eq. (4.7), r0, pr0, ptheta0 and
	 * pphi0 are obtained by solving Eqs. (4.8) and (4.9) (using
	 * gsl_multiroot_fsolver). At this step, we find initial conditions
	 * for a spherical orbit without radiation reaction.
	 */

	/* Calculate the initial velocity from the given initial frequency */
	omega = LAL_PI * mTotal * LAL_MTSUN_SI * fMin;
	v0 = cbrt(omega);

	/* Given this, we can start to calculate the initial conditions */
	/* for spherical coords in the new basis */
	rootParams.omega = omega;
	rootParams.params = params;

	/* To start with, we will just assign Newtonian-ish ICs to the system */
	rootParams.values[0] = 1. / (v0 * v0);	/* Initial r */
	rootParams.values[4] = v0;	/* Initial p */
	memcpy(rootParams.values + 6, tmpS1, sizeof(tmpS1));
	memcpy(rootParams.values + 9, tmpS2, sizeof(tmpS2));

	//printf("ICs guess: r = %.16e, p = %.16e\n", rootParams.values[0], rootParams.values[4]);

	/* Initialise the gsl stuff */
	XLAL_CALLGSL(rootSolver = gsl_multiroot_fsolver_alloc(T, 3));
	if (!rootSolver) {
		gsl_matrix_free(rotMatrix);
		gsl_matrix_free(invMatrix);
		XLAL_ERROR(XLAL_ENOMEM);
	}
	XLAL_CALLGSL(initValues = gsl_vector_calloc(3));
	if (!initValues) {
		gsl_multiroot_fsolver_free(rootSolver);
		gsl_matrix_free(rotMatrix);
		gsl_matrix_free(invMatrix);
		XLAL_ERROR(XLAL_ENOMEM);
	}
	gsl_vector_set(initValues, 0, rootParams.values[0]);
	gsl_vector_set(initValues, 1, rootParams.values[4]);

	rootFunction.f = XLALFindSphericalOrbitV2;
	rootFunction.n = 3;
	rootFunction.params = &rootParams;

	gsl_multiroot_fsolver_set(rootSolver, &rootFunction, initValues);

	/* We are now ready to iterate to find the solution */
	i = 0;

	do {
		XLAL_CALLGSL(gslStatus = gsl_multiroot_fsolver_iterate(rootSolver));
		if (gslStatus != GSL_SUCCESS) {
			XLALPrintError("Error in GSL iteration function!\n");
			gsl_multiroot_fsolver_free(rootSolver);
			gsl_vector_free(initValues);
			gsl_matrix_free(rotMatrix);
			gsl_matrix_free(invMatrix);
			XLAL_ERROR(XLAL_EFUNC);
		}
		XLAL_CALLGSL(gslStatus = gsl_multiroot_test_residual(rootSolver->f, 1.0e-10));
		i++;
	}
	while (gslStatus == GSL_CONTINUE && i <= maxIter);

	if (i > maxIter && gslStatus != GSL_SUCCESS) {
		gsl_multiroot_fsolver_free(rootSolver);
		gsl_vector_free(initValues);
		gsl_matrix_free(rotMatrix);
		gsl_matrix_free(invMatrix);
		XLAL_ERROR(XLAL_EMAXITER);
	}
	finalValues = gsl_multiroot_fsolver_root(rootSolver);

	/*
	 * printf( "Spherical orbit conditions here given by the
	 * following:\n" ); printf( " x = %.16e, py = %.16e, pz = %.16e\n",
	 * gsl_vector_get( finalValues, 0 ), gsl_vector_get( finalValues, 1
	 * ), gsl_vector_get( finalValues, 2 ) );
	 */

	memset(qCart, 0, sizeof(qCart));
	memset(pCart, 0, sizeof(pCart));

	qCart[0] = gsl_vector_get(finalValues, 0);
	pCart[1] = gsl_vector_get(finalValues, 1);
	pCart[2] = gsl_vector_get(finalValues, 2);

	/* Free the GSL root finder, since we're done with it */
	gsl_multiroot_fsolver_free(rootSolver);
	gsl_vector_free(initValues);

	/*
	 * STEP 3) Rotate to L0 along z-axis and N0 along x-axis frame, where
	 * L0 is the initial orbital angular momentum and L0 is calculated
	 * using initial position and linear momentum obtained in STEP 2.
	 */

	/* Now we can calculate the relativistic L and rotate to a new basis */
	memcpy(qHat, qCart, sizeof(qCart));
	memcpy(pHat, pCart, sizeof(pCart));

	NormalizeVectorPrec(qHat);
	NormalizeVectorPrec(pHat);

	Lhat[0] = CalculateCrossProductPrec(0, qHat, pHat);
	Lhat[1] = CalculateCrossProductPrec(1, qHat, pHat);
	Lhat[2] = CalculateCrossProductPrec(2, qHat, pHat);

	NormalizeVectorPrec(Lhat);

	XLAL_CALLGSL(rotMatrix2 = gsl_matrix_alloc(3, 3));
	XLAL_CALLGSL(invMatrix2 = gsl_matrix_alloc(3, 3));

	if (CalculateRotationMatrixPrec(rotMatrix2, invMatrix2, qHat, pHat, Lhat) == XLAL_FAILURE) {
		gsl_matrix_free(rotMatrix);
		gsl_matrix_free(rotMatrix2);
		gsl_matrix_free(invMatrix);
		gsl_matrix_free(invMatrix2);
		XLAL_ERROR(XLAL_ENOMEM);
	}
	ApplyRotationMatrixPrec(rotMatrix2, rHat);
	ApplyRotationMatrixPrec(rotMatrix2, vHat);
	ApplyRotationMatrixPrec(rotMatrix2, LnHat);
	ApplyRotationMatrixPrec(rotMatrix2, tmpS1);
	ApplyRotationMatrixPrec(rotMatrix2, tmpS2);
	ApplyRotationMatrixPrec(rotMatrix2, tmpS1Norm);
	ApplyRotationMatrixPrec(rotMatrix2, tmpS2Norm);
	ApplyRotationMatrixPrec(rotMatrix2, qCart);
	ApplyRotationMatrixPrec(rotMatrix2, pCart);
        
        gsl_matrix_free(rotMatrix);
        gsl_matrix_free(rotMatrix2);

	/*
	 * STEP 4) In the L0-N0 frame, we calculate (dE/dr)|sph using Eq.
	 * (4.14), then initial dr/dt using Eq. (4.10), and finally pr0 using
	 * Eq. (4.15).
	 */

	/* Now we can calculate the flux. Change to spherical co-ords */
	CartesianToSphericalPrec(qSph, pSph, qCart, pCart);
	memcpy(sphValues, qSph, sizeof(qSph));
	memcpy(sphValues + 3, pSph, sizeof(pSph));
	memcpy(sphValues + 6, tmpS1, sizeof(tmpS1));
	memcpy(sphValues + 9, tmpS2, sizeof(tmpS2));

	memcpy(cartValues, qCart, sizeof(qCart));
	memcpy(cartValues + 3, pCart, sizeof(pCart));
	memcpy(cartValues + 6, tmpS1, sizeof(tmpS1));
	memcpy(cartValues + 9, tmpS2, sizeof(tmpS2));

	REAL8		dHdpphi , d2Hdr2, d2Hdrdpphi;
	REAL8		rDot    , dHdpr, flux, dEdr;

	d2Hdr2 = XLALCalculateSphHamiltonianDeriv2Prec(0, 0, sphValues, params);
	d2Hdrdpphi = XLALCalculateSphHamiltonianDeriv2Prec(0, 5, sphValues, params);

	//printf("d2Hdr2 = %.16e, d2Hdrdpphi = %.16e\n", d2Hdr2, d2Hdrdpphi);

	dHdpphi = XLALSpinPrecHcapNumDerivWRTParam(4, cartValues, params) / sphValues[0];

	dEdr = -dHdpphi * d2Hdr2 / d2Hdrdpphi;

	//printf("d2Hdr2 = %.16e d2Hdrdpphi = %.16e dHdpphi = %.16e\n", d2Hdr2, d2Hdrdpphi, dHdpphi);

	if (d2Hdr2 != 0.0) {
		/* We will need to calculate the Hamiltonian to get the flux */
		s1Vec.length = s2Vec.length = s1VecNorm.length = s2VecNorm.length = sKerr.length = sStar.length = 3;
		s1Vec.data = tmpS1;
		s2Vec.data = tmpS2;
		s1VecNorm.data = tmpS1Norm;
		s2VecNorm.data = tmpS2Norm;
		sKerr.data = sKerrData;
		sStar.data = sStarData;

		qCartVec.length = pCartVec.length = 3;
		qCartVec.data = qCart;
		pCartVec.data = pCart;

		//chi1 = tmpS1[0] * LnHat[0] + tmpS1[1] * LnHat[1] + tmpS1[2] * LnHat[2];
		//chi2 = tmpS2[0] * LnHat[0] + tmpS2[1] * LnHat[1] + tmpS2[2] * LnHat[2];

		//printf("magS1 = %.16e, magS2 = %.16e\n", chi1, chi2);

		//chiS = 0.5 * (chi1 / (mass1 * mass1) + chi2 / (mass2 * mass2));
		//chiA = 0.5 * (chi1 / (mass1 * mass1) - chi2 / (mass2 * mass2));

		XLALSimIMRSpinEOBCalculateSigmaKerr(&sKerr, mass1, mass2, &s1Vec, &s2Vec);
		XLALSimIMRSpinEOBCalculateSigmaStar(&sStar, mass1, mass2, &s1Vec, &s2Vec);

		/*
		 * The a in the flux has been set to zero, but not in the
		 * Hamiltonian
		 */
		a = sqrt(sKerr.data[0] * sKerr.data[0] + sKerr.data[1] * sKerr.data[1] + sKerr.data[2] * sKerr.data[2]);
		//XLALSimIMREOBCalcSpinFacWaveformCoefficients(params->eobParams->hCoeffs, mass1, mass2, eta, /* a */ 0.0, chiS, chiA);
		//XLALSimIMRCalculateSpinPrecEOBHCoeffs(params->seobCoeffs, eta, a);
		ham = XLALSimIMRSpinPrecEOBHamiltonian(eta, &qCartVec, &pCartVec, &s1VecNorm, &s2VecNorm, &sKerr, &sStar, params->tortoise, params->seobCoeffs);

		//printf("hamiltonian at this point is %.16e\n", ham);

		/* And now, finally, the flux */
		REAL8Vector	polarDynamics;
		REAL8		polarData[4];

		polarDynamics.length = 4;
		polarDynamics.data = polarData;

		polarData[0] = qSph[0];
		polarData[1] = 0.;
		polarData[2] = pSph[0];
		polarData[3] = pSph[2];

		flux = XLALInspiralSpinFactorizedFlux(&polarDynamics, nqcCoeffs, omega, params, ham, lMax, SpinAlignedEOBversion);
		flux = flux / eta;

		rDot = -flux / dEdr;

		/*
		 * We now need dHdpr - we take it that it is safely linear up
		 * to a pr of 1.0e-3
		 */
		cartValues[3] = 1.0e-3;
        REAL8		csi = sqrt(XLALSimIMRSpinPrecEOBHamiltonianDeltaT(params->seobCoeffs, qSph[0], eta, a)*XLALSimIMRSpinPrecEOBHamiltonianDeltaR(params->seobCoeffs, qSph[0], eta, a)) / (qSph[0] * qSph[0] + a * a);

        dHdpr = csi*csi*XLALSpinPrecHcapNumDerivWRTParam(3, cartValues, params);

		
//		  printf( "Ingredients going into prDot:\n" ); printf( "flux= %.16e, dEdr = %.16e, dHdpr = %.16e\n", flux, dEdr, dHdpr);
		 

		/*
		 * We can now calculate what pr should be taking into account
		 * the flux
		 */
		pSph[0] = rDot / (dHdpr / cartValues[3]);
	} else {
		/*
		 * Since d2Hdr2 has evaluated to zero, we cannot do the
		 * above. Just set pr to zero
		 */
		//printf("d2Hdr2 is zero!\n");
		pSph[0] = 0;
	}

	/* Now we are done - convert back to cartesian coordinates ) */
	SphericalToCartesianPrec(qCart, pCart, qSph, pSph);

	/*
	 * STEP 5) Rotate back to the original inertial frame by inverting
	 * the rotation of STEP 3 and then inverting the rotation of STEP 1.
	 */

	/* Undo rotations to get back to the original basis */
	/* Second rotation */
	ApplyRotationMatrixPrec(invMatrix2, rHat);
	ApplyRotationMatrixPrec(invMatrix2, vHat);
	ApplyRotationMatrixPrec(invMatrix2, LnHat);
	ApplyRotationMatrixPrec(invMatrix2, tmpS1);
	ApplyRotationMatrixPrec(invMatrix2, tmpS2);
	ApplyRotationMatrixPrec(invMatrix2, tmpS1Norm);
	ApplyRotationMatrixPrec(invMatrix2, tmpS2Norm);
	ApplyRotationMatrixPrec(invMatrix2, qCart);
	ApplyRotationMatrixPrec(invMatrix2, pCart);

	/* First rotation */
	ApplyRotationMatrixPrec(invMatrix, rHat);
	ApplyRotationMatrixPrec(invMatrix, vHat);
	ApplyRotationMatrixPrec(invMatrix, LnHat);
	ApplyRotationMatrixPrec(invMatrix, tmpS1);
	ApplyRotationMatrixPrec(invMatrix, tmpS2);
	ApplyRotationMatrixPrec(invMatrix, tmpS1Norm);
	ApplyRotationMatrixPrec(invMatrix, tmpS2Norm);
	ApplyRotationMatrixPrec(invMatrix, qCart);
	ApplyRotationMatrixPrec(invMatrix, pCart);
        
        gsl_matrix_free(invMatrix);
        gsl_matrix_free(invMatrix2);

	/* If required, apply the tortoise transform */
	if (tmpTortoise) {
		REAL8		r = sqrt(qCart[0] * qCart[0] + qCart[1] * qCart[1] + qCart[2] * qCart[2]);
		REAL8		deltaR = XLALSimIMRSpinPrecEOBHamiltonianDeltaR(params->seobCoeffs, r, eta, a);
		REAL8		deltaT = XLALSimIMRSpinPrecEOBHamiltonianDeltaT(params->seobCoeffs, r, eta, a);
		REAL8		csi = sqrt(deltaT * deltaR) / (r * r + a * a);

		REAL8		pr = (qCart[0] * pCart[0] + qCart[1] * pCart[1] + qCart[2] * pCart[2]) / r;

		params->tortoise = tmpTortoise;

		//printf("Applying the tortoise to p (csi = %.26e)\n", csi);

		for (i = 0; i < 3; i++) {
			pCart[i] = pCart[i] + qCart[i] * pr * (csi - 1.) / r;
		}
	}
	/* Now copy the initial conditions back to the return vector */
	memcpy(initConds->data, qCart, sizeof(qCart));
	memcpy(initConds->data + 3, pCart, sizeof(pCart));
	memcpy(initConds->data + 6, tmpS1Norm, sizeof(tmpS1Norm));
	memcpy(initConds->data + 9, tmpS2Norm, sizeof(tmpS2Norm));
    for (i=0; i<12; i++) {
        if (fabs(initConds->data[i]) <=1.0e-15) {
            initConds->data[i] = 0.;
        }
    }
        //gsl_multiroot_fsolver_free(rootSolver);
	//printf("THE FINAL INITIAL CONDITIONS:\n");
	
//printf( " %.16e %.16e %.16e\n%.16e %.16e %.16e\n%.16e %.16e %.16e\n%.16e %.16e %.16e\n", initConds->data[0],initConds->data[1], initConds->data[2], initConds->data[3],initConds->data[4], initConds->data[5], initConds->data[6],initConds->data[7], initConds->data[8], initConds->data[9],initConds->data[10], initConds->data[11] );
	 

	return XLAL_SUCCESS;
}

#endif				/* _LALSIMIMRSPINEOBINITIALCONDITIONS_C */
