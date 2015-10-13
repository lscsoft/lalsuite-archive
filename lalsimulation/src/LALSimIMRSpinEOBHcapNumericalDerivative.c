/*
 * Copyright (C) 2010 Craig Robinson, Yi Pan
 * 
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation; either version 2 of the License, or (at your option)
 * any later version.
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 * 
 * You should have received a copy of the GNU General Public License along with
 * with program; see the file COPYING. If not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */


/**
 * \author Craig Robinson, Yi Pan
 *
 * \brief In newer versions of the EOBNR approximant, we
 * do not have an analytic expression for the derivative of the waveform.
 * As such, it is necessary to calculate the derivatives numerically. This
 * function provides the means to do just that.
 *
 */

#ifndef _LALSIMIMRSPINEOBHCAPNUMERICALDERIVATIVE_C
#define _LALSIMIMRSPINEOBHCAPNUMERICALDERIVATIVE_C

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

/* #include "LALSimIMRSpinEOBHcapNumericalDerivative.h" */

#include <unistd.h>
#include <math.h>
#include <gsl/gsl_deriv.h>
#include <lal/LALSimInspiral.h>
#include <lal/LALSimIMR.h>

#include "LALSimIMRSpinEOB.h"

#include "LALSimIMRSpinEOBAuxFuncs.c"
#include "LALSimIMRSpinEOBHamiltonian.c"
#include "LALSimIMRSpinEOBHamiltonianPrec.c"
#include "LALSimIMRSpinEOBFactorizedFlux.c"
#include "LALSimIMRSpinEOBFactorizedWaveformCoefficients.c"
#include "LALSimIMRSpinEOBFactorizedWaveformCoefficientsPrec.c"
#include "LALSimIMREOBFactorizedWaveform.c"


#endif				/* _LALSIMIMRSPINEOBHCAPNUMERICALDERIVATIVE_C */