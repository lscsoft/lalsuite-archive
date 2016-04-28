#ifndef _LALSIMIMRSPINPRECEOBEULERANGLES_V3OPT_C
#define _LALSIMIMRSPINPRECEOBEULERANGLES_V3OPT_C

/**
 * Computes RHS of ODE for gamma. Eq. 10 of PRD 89, 084006 (2014)
 */
static double f_alphadotcosi( double x, void * inparams )
{
	PrecEulerAnglesIntegration* params = (PrecEulerAnglesIntegration*) inparams;

	REAL8 alphadot = gsl_spline_eval_deriv( params->alpha_spline, x, params->alpha_acc );
	REAL8 beta = gsl_spline_eval( params->beta_spline, x, params->beta_acc );

	return -1. * alphadot * cos(beta);

}

/**
 * Given the trajectory in an inertial frame, this computes Euler angles
 * of the roation from the inertial frame to the minimal-rotation frame
 * that co-precesses with LN(t) = rvec(t) x rdotvec(t)
 */
static int EulerAnglesI2P(REAL8Vector *Alpha, /**<< output: alpha Euler angle */
                 REAL8Vector *Beta, /**<< output: beta Euler angle */
                 REAL8Vector *Gamma, /**<< output: gamma Euler angle */
                 INT4 *phaseCounterA, /**<< output: counter for unwrapping of alpha */
                 INT4 *phaseCounterB, /**<< output: counter for unwrapping of beta */
                 const REAL8Vector tVec, /**<< time series */
                 const REAL8Vector posVecx, /**<< x time series */
                 const REAL8Vector posVecy, /**<< y time series */
                 const REAL8Vector posVecz, /**<< z time series */
                 const UINT4 retLen, /**<< Array length of the trajectory */
                 const REAL8 InitialAlpha, /**<< Initial alpha (used only if flag_highSR=1) */
                 const REAL8 InitialGamma, /**<< Initial gamma */
                 UINT4 flag_highSR, /**<< Flag to indicate whether one is analyzing the high SR trajectory */
                 UNUSED REAL8 * dynamics, /**<< Dynamical variables from ODE solution */ /*OPTV3: For the ODE solution derivative*/
                 UNUSED SpinEOBParams * seobParams /**<< Spin EOB Parameters */ /*OPTV3: For the ODE solution derivative*/
                 ){
    UINT4 i = 0;
    REAL8Vector *LN_x = NULL, *LN_y = NULL, *LN_z = NULL;
    REAL8 tmpR[3], tmpRdot[3], magLN;
    REAL8 precEulerresult = 0, precEulererror = 0;
    gsl_integration_workspace * precEulerw = gsl_integration_workspace_alloc (1000);
    gsl_function precEulerF;
    PrecEulerAnglesIntegration precEulerparams;
    REAL8 inGamma = InitialGamma;

    LN_x = XLALCreateREAL8Vector( retLen );
    LN_y = XLALCreateREAL8Vector( retLen );
    LN_z = XLALCreateREAL8Vector( retLen );

    gsl_spline *x_spline = gsl_spline_alloc( gsl_interp_cspline, retLen );
    gsl_spline *y_spline = gsl_spline_alloc( gsl_interp_cspline, retLen );
    gsl_spline *z_spline = gsl_spline_alloc( gsl_interp_cspline, retLen );

    gsl_interp_accel *x_acc    = gsl_interp_accel_alloc();
    gsl_interp_accel *y_acc    = gsl_interp_accel_alloc();
    gsl_interp_accel *z_acc    = gsl_interp_accel_alloc();
//    if(seobParams->almostAlignedSpins){ //OPTv3
        gsl_spline_init( x_spline, tVec.data, posVecx.data, retLen );
        gsl_spline_init( y_spline, tVec.data, posVecy.data, retLen );
        gsl_spline_init( z_spline, tVec.data, posVecz.data, retLen );

    for( i=0; i < retLen; i++ )
    {
        tmpR[0] = posVecx.data[i]; tmpR[1] = posVecy.data[i]; tmpR[2] = posVecz.data[i];

/* OPTV3: We would prefer to use the ODE solution directly (stored in dynamics[]),
    as the overhead from computing gsl_spline-based derivatives can be significant
    (*and* it is less accurate!), but when the spins are nearly aligned, the ODE
    solution is transformed into (we think) a spin-aligned system
    [see Step 2 of XLALSimIMRSpinEOBWaveformAll_opt()]. So it seems we must revert
    to spline-based derivatives in this case.

//            for(int ijk=0; ijk<3; ijk++)
//                tmpRdot[ijk]=dynamics[i+(ijk+15)*retLen];
*/

        tmpRdot[0] = gsl_spline_eval_deriv( x_spline, tVec.data[i], x_acc );
        tmpRdot[1] = gsl_spline_eval_deriv( y_spline, tVec.data[i], y_acc );
        tmpRdot[2] = gsl_spline_eval_deriv( z_spline, tVec.data[i], z_acc );

        LN_x->data[i] = tmpR[1] * tmpRdot[2] - tmpR[2] * tmpRdot[1];
        LN_y->data[i] = tmpR[2] * tmpRdot[0] - tmpR[0] * tmpRdot[2];
        LN_z->data[i] = tmpR[0] * tmpRdot[1] - tmpR[1] * tmpRdot[0];

        magLN = sqrt(LN_x->data[i] * LN_x->data[i] + LN_y->data[i] * LN_y->data[i]
                     + LN_z->data[i] * LN_z->data[i]);
        LN_x->data[i] /= magLN; LN_y->data[i] /= magLN; LN_z->data[i] /= magLN;

        /*  Eq. 19 of PRD 89, 084006 (2014) */
        /*  Also unwrap the two angles */
        if (fabs(LN_x->data[i]) <= 1.e-10 && fabs(LN_y->data[i]) <=1.e-10){
            Alpha->data[i] = 0.0;
            inGamma = 0.0;
        } else {
            Alpha->data[i] = atan2( LN_y->data[i], LN_x->data[i] )
                             +  *phaseCounterA * LAL_TWOPI;
            if (i==0 && flag_highSR != 1){
                inGamma = -Alpha->data[i];
            }
        }

        if( i>0 && Alpha->data[i] - Alpha->data[i-1] > 5. ) {
            *phaseCounterA = *phaseCounterA - 1;
            Alpha->data[i] -= LAL_TWOPI;
        }
        else if ( i && Alpha->data[i] - Alpha->data[i-1] < -5. ) {
            *phaseCounterA = *phaseCounterA + 1;
            Alpha->data[i] += LAL_TWOPI;
        }
        if (LN_z->data[i] >1.) {
            LN_z->data[i] = 1.;
        }
        if (LN_z->data[i] <-1.) {
            LN_z->data[i] = -1.;
        }
        if ( flag_highSR == 1) {
            Alpha->data[i] -= (Alpha->data[0] - InitialAlpha);
        }

        if (fabs(1.0 - LN_z->data[i]) < 1.e-12){
            REAL8 LN_xy;
            LN_xy = sqrt(LN_x->data[i]*LN_x->data[i] +
			 LN_y->data[i]*LN_y->data[i]);
            //LN_z->data[i] = sqrt(1.0 - LN_xy*LN_xy);
            Beta->data[i] = atan2(LN_xy, LN_z->data[i]);
            //printf("here   ");
        }else{
            Beta->data[i] = acos( LN_z->data[i] );
        }
        if( i>0 && Beta->data[i] > Beta->data[i-1] ) {
            *phaseCounterB = *phaseCounterB - 1;
        }
    }
    /* Integrate \dot{\alpha} \cos{\beta} to get the final Euler angle
     Eq. 20 of PRD 89, 084006 (2014) */
    gsl_spline_init( x_spline, tVec.data, Alpha->data, retLen );
    gsl_spline_init( y_spline, tVec.data, Beta->data, retLen );

    precEulerparams.alpha_spline = x_spline;
    precEulerparams.alpha_acc    = x_acc;
    precEulerparams.beta_spline  = y_spline;
    precEulerparams.beta_acc     = y_acc;

    precEulerF.function = &f_alphadotcosi;
    precEulerF.params   = &precEulerparams;

    for( i = 0; i < retLen; i++ )
    {
        //if( i==0 ) { Gamma->data[i] = InitialGamma; }
        if( i==0 ) { Gamma->data[i] = inGamma; }
        else
        {
            gsl_integration_qags (&precEulerF, tVec.data[i-1], tVec.data[i], 1e-9, 1e-9, 1000, precEulerw, &precEulerresult, &precEulererror);
            Gamma->data[i] = Gamma->data[i-1] + precEulerresult;
        }
    }
    gsl_integration_workspace_free( precEulerw );
    gsl_spline_free( x_spline );
    gsl_spline_free( y_spline );
    gsl_spline_free( z_spline );
    gsl_interp_accel_free( x_acc );
    gsl_interp_accel_free( y_acc );
    gsl_interp_accel_free( z_acc );
    XLALDestroyREAL8Vector( LN_x );
    XLALDestroyREAL8Vector( LN_y );
    XLALDestroyREAL8Vector( LN_z );
    return XLAL_SUCCESS;
}

/**
 * Given Euler angles to go from initial inertial frame to precessing frama
 * and the LNhat vector, this functions computes the Euler angles to
 * go from the precessing frame to the frame of the total angular
 * momentum
 */
static void EulerAnglesP2J(
                REAL8 *aP2J, /**<< alpha Euler angle from precessing to final-J frame */
                REAL8 *bP2J, /**<< beta Euler angle from precessing to final-J frame */
                REAL8 *gP2J, /**<< gamma Euler angle from precessing to final-J frame */
                const REAL8 aI2P, /**<< alpha Euler angle from inertial to precessing frame */
                const REAL8 bI2P, /**<< beta Euler angle from inertial to precessing frame */
                const REAL8 gI2P, /**<< gamma Euler angle from inertial to precessing frame */
                const REAL8 LNhx, /**<< x component of LNhat */
                const REAL8 LNhy, /**<< y component of LNhat */
                const REAL8 LNhz, /**<< z component of LNhat */
                const REAL8 JframeEx[], /**<< x-axis of the total-angular-momentum frame */
                const REAL8 JframeEy[], /**<< y-axis of the total-angular-momentum frame */
                const REAL8 JframeEz[] /**<< z-axis of the total-angular-momentum frame */
) {
    REAL8 LframeEx[3] = {0,0,0}, LframeEy[3] = {0,0,0}, LframeEz[3] = {0,0,0};
    LframeEx[0] =  cos(aI2P)*cos(bI2P)*cos(gI2P) - sin(aI2P)*sin(gI2P);
    LframeEx[1] =  sin(aI2P)*cos(bI2P)*cos(gI2P) + cos(aI2P)*sin(gI2P);
    LframeEx[2] = -sin(bI2P)*cos(gI2P);
    LframeEy[0] = -cos(aI2P)*cos(bI2P)*sin(gI2P) - sin(aI2P)*cos(gI2P);
    LframeEy[1] = -sin(aI2P)*cos(bI2P)*sin(gI2P) + cos(aI2P)*cos(gI2P);
    LframeEy[2] =  sin(bI2P)*sin(gI2P);
    LframeEz[0] =  LNhx;
    LframeEz[1] =  LNhy;
    LframeEz[2] =  LNhz;
    REAL8 normJ, normLz;
    normJ = JframeEz[0]*JframeEz[0]+JframeEz[1]*JframeEz[1]+JframeEz[2]*JframeEz[2];
    normLz = LframeEz[0]*LframeEz[0]+LframeEz[1]*LframeEz[1]+LframeEz[2]*LframeEz[2];
    *aP2J = atan2(JframeEz[0]*LframeEy[0]+JframeEz[1]*LframeEy[1]+JframeEz[2]*LframeEy[2],
                 JframeEz[0]*LframeEx[0]+JframeEz[1]*LframeEx[1]+JframeEz[2]*LframeEx[2]);
    *bP2J = acos( JframeEz[0]*LframeEz[0]+JframeEz[1]*LframeEz[1]+JframeEz[2]*LframeEz[2]);
    if (*bP2J < 1.e-4){
        *bP2J = acos((JframeEz[0]*LframeEz[0]+JframeEz[1]*LframeEz[1]+JframeEz[2]*LframeEz[2])/sqrt(normJ*normLz));

    }

    *gP2J = atan2(  JframeEy[0]*LframeEz[0]+JframeEy[1]*LframeEz[1]+JframeEy[2]*LframeEz[2],
                 -(JframeEx[0]*LframeEz[0]+JframeEx[1]*LframeEz[1]+JframeEx[2]*LframeEz[2]));

    /* I2P Euler angles are stored only for debugging purposes */
    if ( fabs(*bP2J-LAL_PI) < 1.e-10){
        *gP2J = 0.0;
        *aP2J = atan2( JframeEx[1], JframeEx[0]);
    }

    if ( fabs(*bP2J) < 1.e-10){
        *gP2J = 0.0;
        *aP2J = atan2( JframeEx[1], JframeEx[0]);
    }
}
#endif // _LALSIMIMRSPINPRECEOBEULERANGLES_C
