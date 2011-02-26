/*
*  Copyright (C) 2007 Jolien Creighton
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

/* <lalVerbatim> */
#include <stdlib.h>
#include <lal/LALETNullStream.h>

NRCSID( LALETNULLSTREAMTESTC, "$Id$" );

#define USAGE "Usage: %s [--GPSStart INT] [--Nsegs INT] [--SegDur REAL8] [--Srate REAL8]\n"

extern int lalDebugLevel;

/*******************************************************************************
 *
 *  GLOBAL VARIABLES
 * 
 ******************************************************************************/

// GLOBAL VARIABLES USED IN ARGUMENT PARSING - SETS THE DEFAULT VALUE
REAL8 duration_global=1000.0;
LIGOTimeGPS GPSStart_global; // NB: NO INFRASTRUCTURE FOR NANOSECOND INPUT!
REAL8 SampleRate_global=2048.0; // JUSTIFY CHOICE WDP: the data have been generated with this Srate

/*******************************************************************************
 *
 *  FUNCTION PROTOTYPES
 * 
 ******************************************************************************/

void initialise(INT4 argc, char *argv[]);

/*******************************************************************************
 *
 *  MAIN FUNCTION
 * 
 ******************************************************************************/

int main( int argc, char *argv[] )
{
    /***************************************************************************
     *
     *  TEST FUNCTION TO COMPUTE THE ET NULLSTREAM
     * 
     **************************************************************************/
    
    /***************************************************************************
     *
     *  ERROR HANDLING
     * 
     **************************************************************************/
    
    static LALStatus status;
    lalDebugLevel = 0;
    
    /***************************************************************************
     *
     *  DECLARE VARIABLES
     * 
     **************************************************************************/
    
    UINT4 j;
    REAL8TimeSeries *NullStream;
    REAL8FrequencySeries *inverse_spectrum;
    REAL8FFTPlan *fwdplan = NULL;
    REAL8FFTPlan *revplan = NULL;
    REAL8Window  *windowplan = NULL;
    REAL8 deltaF=0.0;
    UINT4 nSegs=100;
    REAL8 segDur,strideDur;
    UINT4 seglen,stride;
    GPSStart_global.gpsSeconds=802338816;
    NullStream = (REAL8TimeSeries *) LALMalloc(sizeof(REAL8TimeSeries));
    segDur = duration_global/(REAL8)nSegs;
    seglen=(UINT4)(segDur*SampleRate_global);
    segDur = seglen/SampleRate_global;
    nSegs =(INT4)floor(duration_global/segDur);

    fprintf(stderr,"Choosing %i segments length %i, (%f s)\n",nSegs,seglen,segDur);

    stride = seglen; /* Overlap the padding */
    strideDur = stride / SampleRate_global;
    
    /***************************************************************************
     *
     *  COMPUTE NULLSTREAM
     * 
     **************************************************************************/
    
    /* Get the arguments and act on them */
    initialise(argc,argv); 
    NullStream = LALETNullStream(&GPSStart_global,duration_global);
    FILE *nsout;
    nsout=fopen("temp_ns.dat","w");
    uint k;
    for (k=0;k<NullStream->data->length;k++){
	fprintf(nsout,"%g %g\n",k*NullStream->deltaT,NullStream->data->data[k]);
	}
    fclose(nsout); 
    // RESAMPLE TIMESERIES (WHY IS THIS DONE?)
    //XLALResampleREAL8TimeSeries(NullStream,1.0/SampleRate_global);
    deltaF=(REAL8)SampleRate_global/seglen;//(REAL8)1.0/((REAL8) NullStream->data->length*NullStream->deltaT);
    printf("deltaF : %g\n",deltaF);    
    /***************************************************************************
     *
     *  COMPUTE PSD (PUT IN FUNCTION INSTEAD OF IN THE TEST FUNCTION)
     *      + OUTPUT TO FILE
     * 
     **************************************************************************/
    
    windowplan = XLALCreateTukeyREAL8Window(seglen,0.1*(REAL8)8.0*SampleRate_global/(REAL8)seglen);
    int check=0;
    fwdplan = XLALCreateForwardREAL8FFTPlan( seglen, 0 );
    revplan = XLALCreateReverseREAL8FFTPlan( seglen, 0 );
    fprintf(stderr,"Shrinking... (lost %d samples from end)\n",NullStream->data->length-(seglen*nSegs));
    NullStream=(REAL8TimeSeries *)XLALShrinkREAL8TimeSeries(NullStream,(size_t) 0, (size_t) seglen*nSegs);
    fprintf(stderr,"Computing power spectrum, seglen %i stride %i\n",seglen,stride);
    inverse_spectrum = (REAL8FrequencySeries *)XLALCreateREAL8FrequencySeries("inverse spectrum",&(NullStream->epoch),0.0,deltaF,&lalDimensionlessUnit,seglen/2+1);
    check=XLALREAL8AverageSpectrumMedian(inverse_spectrum,NullStream,(UINT4)seglen,(UINT4)stride,windowplan,fwdplan);
    if (check) {fprintf(stderr,"Failed! \n");exit(-1);}
    check|=XLALREAL8SpectrumInvertTruncate(inverse_spectrum, 10.0, seglen, (seglen-stride)/4, fwdplan, revplan ); 
    if (check) {fprintf(stderr,"Failed computing spectrum! Exiting...\n");exit(-1);}    
    // PSD OUTPUT FILE
    FILE *psdout;
	char psdname[100];
	sprintf(psdname,"nullstream_%g.dat",SampleRate_global);
	psdout=fopen(psdname,"w");
    
    // PRINT FILE
	for(j=0;j<inverse_spectrum->data->length;j++) {
       fprintf(psdout,"%10.10lf %10.10e\n",j*deltaF,inverse_spectrum->data->data[j]); 
    }
	fclose(psdout);
    
    /***************************************************************************
     *
     *  CLEAN AND EXIT
     * 
     **************************************************************************/
      
    XLALDestroyREAL8TimeSeries(NullStream);
    XLALDestroyREAL8FrequencySeries(inverse_spectrum);
    XLALDestroyREAL8FFTPlan(fwdplan);
    XLALDestroyREAL8Window(windowplan);
    
    REPORTSTATUS( &status );
    return status.statusCode;
}

/*******************************************************************************
 *
 *  FUNCTION
 * 
 ******************************************************************************/

void initialise(INT4 argc, char *argv[])
{
    
	/***************************************************************************
     *
     *  Argument Parsing
     * 
     **************************************************************************/
	
    // ARGUMENT LOOP INTEGER
    INT4 arg = 0;
    
    
    // GET INPUT FROM --GPSSTART FLAG (OPTIONAL, DEFAULT DEFINED AS GLOBAL)
	for (arg=1; arg<argc; arg++)
	{
		if (!strcmp(argv[arg],"--GPSstart")) 
		{
			if(argc<=arg+1)
			{
				LALPrintError( USAGE, *argv );
				exit(-1);
			}
			else if( (argv[arg+1][0]=='-') )
			{
				LALPrintError( USAGE, *argv );
				exit(-1);
			}
			else if( argc>arg+2 )
			{
				if(argv[arg+2][0]!='-')
				{
					LALPrintError( USAGE, *argv );
					exit(-1);
				}
				else
				{
					GPSStart_global.gpsSeconds = atoi(argv[++arg]);
				}
			}
			else
			{
				GPSStart_global.gpsSeconds = atoi(argv[++arg]);
			}
        }
        
        // GET INPUT FROM --DURATION FLAG (OPTIONAL, DEFAULT DEFINED AS GLOBAL)
		else if (!strcmp(argv[arg],"--duration")) 
		{
			if(argc<=arg+1)
			{
				LALPrintError( USAGE, *argv );
				exit(-1);
			}
			else if( (argv[arg+1][0]=='-') )
			{
				LALPrintError( USAGE, *argv );
				exit(-1);
			}
			else if( argc>arg+2 )
			{
				if(argv[arg+2][0]!='-')
				{
					LALPrintError( USAGE, *argv );
					exit(-1);
				}
				else
				{
					duration_global = atof(argv[++arg]);
				}
			}
			else
			{
				duration_global = atof(argv[++arg]);
			}
		}
        
        // GET INPUT FROM --SRATE FLAG (OPTIONAL, DEFAULT DEFINED AS GLOBAL)
		else if (!strcmp(argv[arg],"--Srate")) 
		{
			if(argc<=arg+1)
			{
				LALPrintError( USAGE, *argv );
				exit(-1);
			}
			else if( (argv[arg+1][0]=='-') )
			{
				LALPrintError( USAGE, *argv );
				exit(-1);
			}
			else if( argc>arg+2 )
			{
				if(argv[arg+2][0]!='-')
				{
					LALPrintError( USAGE, *argv );
					exit(-1);
				}
				else
				{
					SampleRate_global = atof(argv[++arg]);
				}
			}
			else
			{
				SampleRate_global = atof(argv[++arg]);
			}
		}
    }
    
    
}

/* </lalVerbatim> */
