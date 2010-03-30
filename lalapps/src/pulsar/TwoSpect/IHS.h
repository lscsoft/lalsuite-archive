/*
*  Copyright (C) 2010 Evan Goetz
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

#ifndef __IHS_H__
#define __IHS_H__

#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include "TwoSpectTypes.h"
#include "TwoSpect.h"
#include "candidates.h"

ihsMaximaStruct *new_ihsMaxima(ffdataStruct *ffdata, INT4 columns);
void free_ihsMaxima(ihsMaximaStruct *data);
void runIHS(ihsMaximaStruct *out, ffdataStruct *in, INT4 columns);

ihsVals * new_ihsVals(void);
void free_ihsVals(ihsVals *ihsvals);
void incHarmSum(ihsVals *out, REAL8Vector *in);

ihsfarStruct * new_ihsfarStruct(INT4 columns);
void free_ihsfarStruct(ihsfarStruct *ihsfarstruct);
void genIhsFar(ihsfarStruct *out, ffdataStruct *ffdata, INT4 columns, REAL8 threshold);

void findIHScandidates(candidate *candlist[], INT4 *numofcandidates, ihsfarStruct *ihsfarstruct, REAL8Vector *aveFFnoise, inputParamsStruct *inputParams, ffdataStruct *ffdata, ihsMaximaStruct *ihsmaxima, REAL4 ra, REAL4 dec);

REAL8Vector * ihsSums(REAL8Vector *ihss, INT4 cols);

REAL8 ihsFOM(REAL8Vector *ihss, INT4Vector *locs, REAL8Vector *expect, REAL8Vector *sigma);
REAL8 ihsLoc(REAL8Vector *ihss, INT4Vector *locs, REAL8Vector *expect, REAL8Vector *sigma);


#endif



