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
#include <stdbool.h>
#include <stdlib.h> 
#include <string.h>
#include <lal/LALMalloc.h>
#include <lal/LALError.h>

/**
 * Linked list node for the testing GR parameters
 */
typedef struct tagLALSimGRTestParamData
{
    char name[32]; 	/**< Name of the variable */
    double value;  	/**< Value of the variable */
} LALSimInspiralGRTestParamData;

/**
 * Linked list of any number of parameters for testing GR
 */
typedef struct tagLALSimGRTestParam
{
    struct tagLALSimGRTestParamData *data; /**< Current variable */
    struct tagLALSimGRTestParam *next; /**< The next variable in linked list */
}  LALSimGRTestParam;

/**
 * Function that creates the head node of the GR parameters linked list.
 * It is initialized with a single parameter with given name and value
 */
LALSimGRTestParam *XLALSimCreateGRParam(
        const char *name, /**< Name of first parameter in new linked list */
        double value 	 /**< Value of first parameter in new linked list */
        );

/**
 * Function that adds a prameter to the GR parameters linked list. If the
 * parameter already exists, it prints a warning and does nothing.
*/
void XLALSimAddGRParam(
        LALSimGRTestParam *parameter, 	/**< Linked list of parameters */
        const char *name, 		/**< Parameter name */
        double value 			/**< Parameter value */
        );

/**
 * Function that sets the value of the desired parameter in the GR parameters
 * linked list to 'value'.  Throws an error if the parameter is not found
 */
int XLALSimSetGRParamValue(
        LALSimGRTestParam *parameter, 	/**< Linked list to be modified */
        const char *name, 		/**< Name of parameter to be modified */
        double value 			/**< New value for parameter */
        );

/**
 * Function that returns the value of the desired parameters in the
 * GR parameters linked list.  Aborts if the parameter is not found
 */
double XLALSimGetGRParamValue(
        LALSimGRTestParam *parameter, 	/**< Linked list to retrieve from */
        const char *name 	   /**< Name of parameter to be retrieved */
        );

/**
 * Function that checks whether the requested parameter exists within the
 * GR parameters linked list.  Returns true (1) or false (0) accordingly
 */
bool XLALSimGRParamExists(
        LALSimGRTestParam *parameter, 	/**< Linked list to check */
        const char *name 		/**< Parameter name to check for */
        );

/** Function that prints off the whole list */
void XLALSimPrintGRParamStruct(
        FILE *fp, 			/** FILE pointer to write to */
        LALSimGRTestParam *parameter 	/**< Linked list to print */
        );

/** Function that destroys the list */
void XLALSimDestroyGRParam(
        LALSimGRTestParam *parameter 	/**< Linked list to destroy */
        );

#endif /* _LALSIMINSPIRALWAVEFORMGRFLAGS_H */
