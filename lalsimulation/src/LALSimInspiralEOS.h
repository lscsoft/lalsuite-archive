/*
 *  Copyright (C) 2012 Walter Del Pozzo, Tjonnie Li, Michalis Agathos
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with with program; see the file COPYING. If not, write to the
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 *  MA  02111-1307  USA
 */

#ifndef _LALSIMINSPIRALEOS_H  /* Double-include protection. */
#define _LALSIMINSPIRALEOS_H

#ifdef  __cplusplus   /* C++ protection. */
extern "C" {
#endif

#include <lal/LALDatatypes.h>

typedef enum
{   
    LAL_SIM_INSPIRAL_EOS_NONE = 0,
    LAL_SIM_INSPIRAL_EOS_MS1 = 1,
    LAL_SIM_INSPIRAL_EOS_H4 = 2,
    LAL_SIM_INSPIRAL_EOS_SQM3 = 3,
    LAL_SIM_INSPIRAL_EOS_MPA1 = 4,
    LAL_SIM_INSPIRAL_EOS_GNH3 = 5
} LALEquationOfState;

LALEquationOfState XLALSimEOSfromString(char eos_name[]);

REAL8 XLALSimInspiralEOSLambda(LALEquationOfState eos_type, REAL8 m_intr_msun);


#ifdef  __cplusplus
}
#endif  /* C++ protection. */
#endif  /* Double-include protection. */