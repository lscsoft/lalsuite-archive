/* Copyright (C) 2012 Walter Del pozzo
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

#ifndef _LALSIMINSPIRALWAVEFORMGRFLAGS_H
#define _LALSIMINSPIRALWAVEFORMGRFLAGS_H

#include <stdio.h> 
#include <stdlib.h> 
#include <string.h>
#include <lal/LALMalloc.h>
#include <lal/LALError.h>

/** Default values for all enumerated flags */ 
#define LAL_SIM_INSPIRAL_GR_TEST_DEFAULT LAL_SIM_INSPIRAL_GR_TEST_NONE

/***************************************
* Enumerator that contains all flags to turn on testing each PN cefficient
***************************************/

typedef enum {
    LAL_SIM_INSPIRAL_GR_TEST_NONE = 1, /**< No GR test */
    LAL_SIM_INSPIRAL_GR_TEST_0PN = 0, /**< Leading order */
    LAL_SIM_INSPIRAL_GR_TEST_05PN = 0, /**< 0.5PN */
    LAL_SIM_INSPIRAL_GR_TEST_1PN = 0, /**< 1PN */
    LAL_SIM_INSPIRAL_GR_TEST_15PN = 0, /**< 1.5PN */
    LAL_SIM_INSPIRAL_GR_TEST_2PN = 0, /**< 2PN */
    LAL_SIM_INSPIRAL_GR_TEST_25PN = 0, /**< 2.5PN */
    LAL_SIM_INSPIRAL_GR_TEST_25PNL = 0, /**< 2.5PN logarithmic term */
    LAL_SIM_INSPIRAL_GR_TEST_3PN = 0, /**< 3PN */    
    LAL_SIM_INSPIRAL_GR_TEST_3PNL = 0, /**< 3PN logarithmic term*/   
    LAL_SIM_INSPIRAL_GR_TEST_35PN = 0, /**< 3.5PN */
} LALSimInspiralGRTestFlags;

/***************************************
* Linked list node for the testing GR parameters
***************************************/

typedef struct tagLALSimGRTestParamData
{
    char        name[32];
    double    value;
} LALSimInspiralGRTestParamData;

typedef struct tagLALSimGRTestParam
{
    struct tagLALSimGRTestParamData *data;
    struct tagLALSimGRTestParam *next;
}  LALSimGRTestParam;

LALSimGRTestParam *XLALSimCreateGRParam(const char *name, 
                                                        double value);
    
void XLALSimAddGRParam(LALSimGRTestParam *parameter,
                                                    const char *name, 
                                                    double value);

void XLALSimSetGRParamValue(LALSimGRTestParam *parameter, 
                                                        const char *name, 
                                                        double value);
  
double XLALSimGetGRParamValue(LALSimGRTestParam *parameter, 
                                                                const char *name);

int XLALSimGRParamExists(LALSimGRTestParam *parameter, 
                                                    const char *name);

//void XLALSimCopyGRPara(LALSimGRTestParam **parameterOutPtr, 
//                                                   LALSimGRTestParam *parameterIn);

void XLALSimPrintGRParamStruct(FILE *fp, LALSimGRTestParam *parameter);
                                                        
void XLALSimFreeGRParam(LALSimGRTestParam *parameter);

#endif /* _LALSIMINSPIRALWAVEFORMGRFLAGS_H */