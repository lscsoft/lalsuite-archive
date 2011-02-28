/**
 * \file LALAdaptiveRungeKutta4.h
 * \ingroup LALAdaptiveRungeKutta4
 * \author M.Vallisneri
 *
 * \breif Functions associated with the XLALAdaptiveRungeKutta4.
 *
 * Prototypes:
 *   - integrator
 *     Integration structure (quasi-class). Created using XLALAdaptiveRungeKutta4Init().
 *   - ...
 *   .
 *
 * Description
 *
 * The code LALAdaptiveRungeKutta4.c evolves a system of \f$n\f$ coupled first--order differential equations.
 * Internally, it uses GSL routines to perform adaptive-step evolution, and then interpolates the resulting
 * trajectories to a fixed step size.
 *
 * Prior to evolving a system using XLALAdaptiveRungeKutta4(), it is necessary to create an integrator structure using
 * XLALAdaptiveRungeKutta4Init(). Once you are done with the integrator, free it with XLALAdaptiveRungeKutta4Free().
 *
 * Algorithm
 *
 * TBF.
 *
 * Uses
 *
 * For updated SpinTaylor waveforms.
 *
 * Notes
 *
 * None so far...
 */

#ifndef _LALADAPTIVERUNGEKUTTA4_H
#define _LALADAPTIVERUNGEKUTTA4_H

#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_spline.h>

#include <lal/LALGSL.h>
#include <lal/SeqFactories.h>

#ifdef  __cplusplus
extern "C" {
#pragma }
#endif

#define XLAL_BEGINGSL \
	{ \
		gsl_error_handler_t *saveGSLErrorHandler_; \
		XLALGSL_PTHREAD_MUTEX_LOCK; \
		saveGSLErrorHandler_ = gsl_set_error_handler_off();
#define XLAL_ENDGSL \
		gsl_set_error_handler( saveGSLErrorHandler_ ); \
		XLALGSL_PTHREAD_MUTEX_UNLOCK; \
	}

/**
 * The structure for the XLALAdaptiveRungeKutta4 integrator
 */
typedef struct
tagark4GSLIntegrator
{
	gsl_odeiv_step    *step;
	gsl_odeiv_control *control;
	gsl_odeiv_evolve  *evolve;

	gsl_odeiv_system  *sys;   

	int (* dydt) (double t, const double y[], double dydt[], void * params);
	int (* stop) (double t, const double y[], double dydt[], void * params);

	int retries;         /* retries with smaller step when derivatives encounter singularity */
	int stopontestonly;  /* stop only on test, use tend to size buffers only */

	int returncode;
} ark4GSLIntegrator;

/**
 * The initialization function for the XLALAdaptiveRungeKutta4 integrator
 */
ark4GSLIntegrator *XLALAdaptiveRungeKutta4Init(
		int dim,  /**< dimensionality of the integration functions */
		int (* dydt) (double t, const double y[], double dydt[], void * params),  /**< derivatives functions */ /* These are XLAL functions! */
		int (* stop) (double t, const double y[], double dydt[], void * params),  /**< stopping comditions */
		double eps_abs,  /**< abosolute error tolerance */ /* FIXME: is this right? */
		double eps_rel  /**< relative error tolerance */ /* FIXME: is this right? */
		);

/**
 * The function that frees a XLALAdaptiveRungeKutta4 integrator
 */
void XLALAdaptiveRungeKutta4Free(
		ark4GSLIntegrator *integrator /**< pointer to the integrator */
		);

/**
 * The XLALAdaptiveRungeKutta4 integrator itself
 */
unsigned int XLALAdaptiveRungeKutta4(
		/* FIXME: are these descriptions correct? */
		ark4GSLIntegrator *integrator,  /**< pointer to integrator */
		void *params,  /**< integrator parameters */
		REAL8 *yinit,  /**< initial parameters */
		REAL8 tinit,  /**< initial time */
		REAL8 tend,  /**< end time */
		REAL8 deltat,  /**< fixed deltat between sample to interpolate onto */
		REAL8Array **yout  /**< output array */
		);

#ifdef  __cplusplus
#pragma {
}
#endif

#endif /* _LALADAPTIVERUNGEKUTTA4_H */
