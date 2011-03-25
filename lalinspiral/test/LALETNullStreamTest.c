/*
*  Copyright (C) 2011 Walter Del Pozzo, Tjonnie Li
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
#include <getopt.h>
#include <lal/LALETNullStream.h>

NRCSID( LALETNULLSTREAMTESTC, "$Id$" );

#define USAGE "Usage: LALETNullStreamTest [--duration REAL8] [--GPSstart INT] [--Nsegs INT] [--o output] \n \
	 [OPTIONAL] [--Srate REAL8 [DEFAULT:8192]]\n"

extern int lalDebugLevel;

/*******************************************************************************
 *
 *  GLOBAL VARIABLES
 * 
 ******************************************************************************/

// GLOBAL VARIABLES USED IN ARGUMENT PARSING - SETS THE DEFAULT VALUE


CHAR ifo[LALNameLength] = "NS";
CHAR channel[LALNameLength]="NS:STRAIN";

CHAR outfile[FILENAME_MAX];
REAL8 duration_global=100.0;
LIGOTimeGPS GPSStart_global; // NB: NO INFRASTRUCTURE FOR NANOSECOND INPUT!
UINT4 SampleRate_global=8192; // JUSTIFY CHOICE WDP: the data have been generated with this Srate
UINT4 Nsegs_global=10; 
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
    UINT4 nSegs=0;
    REAL8 segDur,strideDur;
    UINT4 seglen,stride;
    REAL8 end_freq=10.0; /* cutoff frequency */ 
    FrOutPar frSerie = {ifo, channel, ADCDataChannel, 1, 0, 0 };
    initialise(argc,argv);	
    nSegs=Nsegs_global;
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
    NullStream = LALETNullStream(&GPSStart_global,duration_global);
    
    /* write the null stream to a frame file */
    printf("writing the NullStream time serie to file...\n");
    LALFrWriteREAL8TimeSeries(&status,NullStream,&frSerie);
    printf("done\n");
    /*
    FILE *nsout;
    nsout=fopen("temp_ns.dat","w");
    uint k;
    for (k=0;k<NullStream->data->length;k++){
	fprintf(nsout,"%g %g\n",k*NullStream->deltaT,NullStream->data->data[k]);
	}
    fclose(nsout);
    */ 
    // RESAMPLE TIMESERIES (WHY IS THIS DONE?)
    if (SampleRate_global!=8192) {
	fprintf(stderr,"Sample rate %d, resampling...\r",SampleRate_global);
	XLALResampleREAL8TimeSeries(NullStream,1.0/SampleRate_global);
	fprintf(stderr,"done\n");
	}
    deltaF=(REAL8)SampleRate_global/seglen;
    fprintf(stderr,"deltaF : %g\n",deltaF);    
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
    check|=XLALREAL8SpectrumInvertTruncate(inverse_spectrum, end_freq, seglen, (seglen-stride)/4, fwdplan, revplan ); 
    if (check) {fprintf(stderr,"Failed computing spectrum! Exiting...\n");exit(-1);}
    /* Normalize */
    for (j=0;j<inverse_spectrum->data->length;j++) {
    	inverse_spectrum->data->data[j]*=1.0L/3.0L;
    }    
    // PSD OUTPUT FILE
    FILE *psdout;
	//char psdname[100];
	//sprintf(psdname,"nullstream_%i.dat",SampleRate_global);
	psdout=fopen(outfile,"w");
    
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
    
    return status.statusCode;
}

/*******************************************************************************
 *
 *  FUNCTION
 * 
 ******************************************************************************/

void initialise(int argc, char *argv[])
{
	int i;
	double GPS;
	/*	sprintf(outfile,"default.dat"); */
	/* Sets up global variables from the command line */
	static struct option long_options[]=
	{	{"duration",required_argument,0,'a'},
		{"GPSstart",required_argument,0,'b'},
		{"Nsegs",required_argument,0,'c'},
		{"Srate",required_argument,0,'d'},
                {"out",required_argument,0,'e'},
		{0,0,0,0}};
	if(argc<=1) {fprintf(stderr,USAGE); exit(-1);}
	while((i=getopt_long(argc,argv,"a:b:c:d:e",long_options,&i))!=-1)
		{ 
		switch(i) {
			case 'a':
				duration_global=atof(optarg);
				break;	
			case 'b':
				GPS=atof(optarg);
				XLALGPSSetREAL8(&GPSStart_global,GPS);
				break;
			case 'c':
				Nsegs_global=atoi(optarg);
				break;
			case 'd':
				SampleRate_global=atof(optarg);
				fprintf(stderr,"Setting the sample rate to %i. This might not be what you want!\n",SampleRate_global);
				break;
                        case 'e':
				strcpy(outfile,optarg);
                                break;  
			default:
				fprintf(stdout,USAGE); exit(0);
				break;
			}
		}
}

