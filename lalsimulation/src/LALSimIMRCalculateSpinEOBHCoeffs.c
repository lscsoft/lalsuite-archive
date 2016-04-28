#ifndef _LALSIMIMRCALCULATESPINEOBHCOEFFS_C
#define _LALSIMIMRCALCULATESPINEOBHCOEFFS_C

/*------------------------------------------------------------------------------------------
 *
 *          Prototypes of functions defined in this code.
 *
 *------------------------------------------------------------------------------------------
 */

static int XLALSimIMRCalculateSpinEOBHCoeffs(
        SpinEOBHCoeffs *coeffs,
        const REAL8    eta,
        const REAL8    a,
        const UINT4    SpinAlignedEOBversion
        );


/*------------------------------------------------------------------------------------------
 *
 *          Defintions of functions.
 *
 *------------------------------------------------------------------------------------------
 */

/**
 *
 * This function is used to calculate some coefficients which will be used in the
 * spinning EOB Hamiltonian. It takes the following inputs:
 *
 * coeffs - a (non-null) pointer to a SpinEOBParams structure. This will be populated
 * with the output.
 * eta - the symmetric mass ratio.
 * sigmaKerr - the spin of the effective Kerr background (a combination of the individual spins).
 *
 * If all goes well, the function will return XLAL_SUCCESS. Otherwise, XLAL_FAILURE is returned.
 */
static int XLALSimIMRCalculateSpinEOBHCoeffs(
        SpinEOBHCoeffs *coeffs, /**<< OUTPUT, EOB parameters including pre-computed coefficients */
        const REAL8    eta,     /**<< symmetric mass ratio */
        const REAL8    a,       /**<< Normalized deformed Kerr spin */
        const UINT4    SpinAlignedEOBversion  /**<< 1 for SEOBNRv1; 2 for SEOBNRv2 */
        )
{

  REAL8 KK, k0, k1, k2, k3, k4, k5, k5l, k1p2, k1p3;
  REAL8 m1PlusEtaKK;

  coeffs->SpinAlignedEOBversion = SpinAlignedEOBversion;

  /* Constants are fits taken from Eq. 37 */
  static const REAL8 c0  = 1.4467; /* needed to get the correct self-force results */
  static const REAL8 c1  = -1.7152360250654402;
  static const REAL8 c2  = -3.246255899738242;

  static const REAL8 c20  = 1.712;
  static const REAL8 c21  = -1.803949138004582;
  static const REAL8 c22  = -39.77229225266885;
  static const REAL8 c23  = 103.16588921239249;

  if ( !coeffs )
  {
    XLAL_ERROR( XLAL_EINVAL );
  }


  coeffs->b3  = 0.;
  coeffs->bb3 = 0.;
  coeffs->KK = KK  = c0 + c1*eta + c2*eta*eta;
  if ( SpinAlignedEOBversion == 2)
  {
     coeffs->KK = KK = c20 + c21*eta + c22*eta*eta + c23*eta*eta*eta;
  }

  m1PlusEtaKK = -1. + eta*KK;
  /* Eqs. 5.77 - 5.81 of BB1 */
  coeffs->k0 = k0 = KK*(m1PlusEtaKK - 1.);
  coeffs->k1 = k1 = - 2.*(k0 + KK)*m1PlusEtaKK;
  k1p2= k1*k1;
  k1p3= k1*k1p2;
  coeffs->k2 = k2 = (k1 * (k1 - 4.*m1PlusEtaKK)) / 2. - a*a*k0*m1PlusEtaKK*m1PlusEtaKK;
  coeffs->k3 = k3 = -k1*k1*k1/3. + k1*k2 + k1*k1*m1PlusEtaKK - 2.*(k2 - m1PlusEtaKK)*m1PlusEtaKK - a*a*k1*m1PlusEtaKK*m1PlusEtaKK;
  coeffs->k4 = k4 = (24.*k1*k1*k1*k1 - 96.*k1*k1*k2 + 48.*k2*k2 - 64.*k1*k1*k1*m1PlusEtaKK
      + 48.*a*a*(k1*k1 - 2.*k2)*m1PlusEtaKK*m1PlusEtaKK +
      96.*k1*(k3 + 2.*k2*m1PlusEtaKK) - m1PlusEtaKK*(192.*k3 + m1PlusEtaKK*(-3008. + 123.*LAL_PI*LAL_PI)))/96.;
  coeffs->k5 = k5 = 0.0;
  coeffs->k5l= k5l= 0.0;
  if ( SpinAlignedEOBversion == 2 )
  {
    coeffs->k5 = k5 = m1PlusEtaKK*m1PlusEtaKK
	       * (-4237./60.+128./5.*LAL_GAMMA+2275.*LAL_PI*LAL_PI/512.
	       - 1./3.*a*a*(k1p3-3.*k1*k2+3.*k3)
	       - (k1p3*k1p2-5.*k1p3*k2+5.*k1*k2*k2+5.*k1p2*k3-5.*k2*k3-5.*k1*k4)/5./m1PlusEtaKK/m1PlusEtaKK
	       + (k1p2*k1p2-4.*k1p2*k2+2.*k2*k2+4.*k1*k3-4.*k4)/2./m1PlusEtaKK+256./5.*log(2.));
    coeffs->k5l = k5l = m1PlusEtaKK*m1PlusEtaKK * 64./5.;
  }

  /*printf( "a = %.16e, k0 = %.16e, k1 = %.16e, k2 = %.16e, k3 = %.16e, k4 = %.16e, b3 = %.16e, bb3 = %.16e, KK = %.16e\n",
            a, coeffs->k0, coeffs->k1, coeffs->k2, coeffs->k3, coeffs->k4, coeffs->b3, coeffs->bb3, coeffs->KK );
  */

  /* Now calibrated parameters for spin models */
  coeffs->d1 = coeffs->d1v2 = 0.0;
  coeffs->dheffSS = coeffs->dheffSSv2 = 0.0;
  switch ( SpinAlignedEOBversion )
  {
     case 1:
       coeffs->d1 = -69.5;
       coeffs->dheffSS = 2.75;
       break;
     case 2:
       coeffs->d1v2 = -74.71 - 156.*eta + 627.5*eta*eta;
       coeffs->dheffSSv2 = 8.127 - 154.2*eta + 830.8*eta*eta;
       break;
     default:
       XLALPrintError( "XLAL Error - %s: wrong SpinAlignedEOBversion value, must be 1 or 2!\n", __func__ );
       XLAL_ERROR( XLAL_EINVAL );
       break;
  }

  return XLAL_SUCCESS;
}

#endif // _LALSIMIMRCALCULATESPINEOBHCOEFFS_C
