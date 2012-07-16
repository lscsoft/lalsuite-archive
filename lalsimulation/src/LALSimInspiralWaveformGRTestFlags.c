/* Copyright (C) 2012 Walter Del Pozzo
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
 
#include  <lal/LALSimInspiralWaveformGRTestFlags.h>

/* Function that creates the head node of the GR parameters linked list. */

LALSimGRTestParam *XLALSimCreateGRParam(const char *name, double value)
{
        LALSimGRTestParam *parameter = (LALSimGRTestParam *)malloc(sizeof(LALSimGRTestParam));
        parameter->data =  (LALSimInspiralGRTestParamData *)malloc(sizeof(LALSimInspiralGRTestParamData));
        memcpy(parameter->data->name, name, 32);
        parameter->data->value = value;
        parameter->next=NULL;
        return parameter;
}

/* Function that adds a prameter to the GR parameters linked list. */

void XLALSimAddGRParam(LALSimGRTestParam *parameter,
                                                    const char *name, 
                                                    double value)
{
    if (!XLALSimGRParamExists(parameter, name))
    {
        LALSimGRTestParam *newParam = (LALSimGRTestParam *)malloc(sizeof(LALSimGRTestParam));
        newParam->data =  (LALSimInspiralGRTestParamData *)malloc(sizeof(LALSimInspiralGRTestParamData));
        memcpy(newParam->data->name, name, 32);
        newParam->data->value = value;
        newParam->next = parameter->next;
        parameter->next = newParam;
    }
    else 
    {
        fprintf(stderr,"WARNING! Parameter '%s' exists already! Not added to the structure\n",name);
    }

}

/* Function that checks whether the requested parameter exists within the GR parameters linked list.
* Returns 1 for yes, 0 for no */

int XLALSimGRParamExists(LALSimGRTestParam *parameter, const char *name)
{
  while(parameter) {if(!strcmp(parameter->data->name, name)) return 1; else parameter=parameter->next;}
  return 0;
}

/* Function that returns the value of the desired parameters in the GR parameters linked list.
* Aborts if the parameter is not found */

double XLALSimGetGRParamValue(LALSimGRTestParam *parameter, const char *name)
{
    if (XLALSimGRParamExists(parameter, name)) 
        {
            while(parameter) 
            {
                if(!strcmp(parameter->data->name, name)) return parameter->data->value;
                parameter=parameter->next;
            }
        }
    else 
    {
        fprintf( stderr,"ERROR: parameter '%s' unknown!\n",name);
        exit(-1);
    }
    return 0.0;
}

/* Function that sets the value of the desired parameters in the GR parameters linked list to value.
* Aborts if the parameter is not found */

void XLALSimSetGRParamValue(LALSimGRTestParam *parameter, const char *name, double value)
{
    if (XLALSimGRParamExists(parameter, name)) 
    {
        while(parameter) 
        {
            if(!strcmp(parameter->data->name, name)) parameter->data->value = value;
            parameter=parameter->next;
        }
    }
    else 
    {
        fprintf( stderr,"ERROR: parameter '%s' unknown!\n",name );
        exit(-1);
    }
}

/* Function that prints off the whole list */

void XLALSimPrintGRParamStruct(FILE *fp, LALSimGRTestParam *parameter)
{
        while(parameter) 
        {
            fprintf(fp,"%s %10.5f\n",parameter->data->name,parameter->data->value);
            parameter=parameter->next;
        }
}

/* Function that frees the list */

void XLALSimFreeGRParam(LALSimGRTestParam *parameter)
{
    if( parameter!= NULL ) 
    { 
        XLALSimFreeGRParam(parameter->next);
        parameter->next = NULL;
    }
    free(parameter);
}