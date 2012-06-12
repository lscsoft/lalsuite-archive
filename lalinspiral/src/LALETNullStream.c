/*
*  Copyright (C) 2010 Walter Del Pozzo, Tjonnie Li, Chris Van Den Broeck
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

#include <lal/LALETNullStream.h>

//NRCSID (LALETNULLSTREAMC,"$Id$");

REAL8TimeSeries * LALETNullStream (LIGOTimeGPS *GPSStart, REAL8 duration )
{
    /***************************************************************************
     *
     *  FUNCTION TO COMPUTE THE ETNULLSTREAM FROM A 3 CHANNEL CACHE FILES
     *
     **************************************************************************/
    
    /***************************************************************************
     *
     *  DECLARE VARIABLES
     *
     **************************************************************************/    
    
    UINT4 i=0;
    UINT4 j=0;
    REAL8TimeSeries *RawData[3]; // INDIVIDUAL DETECTOR STRAINS
    REAL8TimeSeries *NullStream; // NULLSTREAM STRAIN
    
    /***************************************************************************
     *
     *  INPUT DATA SPECIFICS - SET BY HAND CURRENTLY, MIGHT WANT TO ADD TO ARGUMENT
     *      PARSING
     *
     **************************************************************************/  
    
    const CHAR *ChannelNames[3] = {"E1:STRAIN", "E2:STRAIN", "E3:STRAIN"};
    const CHAR *CacheFileNames[3] = {"/home/tania/cache/E1GW.cache","/home/tania/cache/E2GW.cache","/home/tania/cache/E3GW.cache"};
    //const CHAR *CacheFileNames[3] = {"/home/tania/cache/E1new.cache","/home/tania/cache/E2new.cache","/home/tania/cache/E3new.cache"};
    //const CHAR *CacheFileNames[3] = {"/home/tania/cache/E1_GWOnlyNewton.cache","/home/tania/cache/E2_GWOnlyNewton.cache","/home/tania/cache/E3_GWOnlyNewton.cache"};
    
    /************************************************************************** 
     *
     *  READ DATA FROM CACHE FILE
     *
     **************************************************************************/  
    
    /* Read the data from disk into a vector (RawData) */
    for (i=0; i<3; i++) {
        RawData[i] = ReadTimeSerieFromCache(CacheFileNames[i], ChannelNames[i],GPSStart,duration);
        if(&RawData[i]==NULL){fprintf(stderr,"Error opening %s in %s\n",ChannelNames[i],CacheFileNames[i]); exit(-1);}
    }
    
    /***************************************************************************
     *
     *  COMPUTE NULL STREAM
     *		- DEFINED AS THE SUM OF THE STRAIN TIME SERIES OF THE THREE DETECTORS
     *		- SEE EQN 32 OF THE MDC PAPER
     *
     **************************************************************************/  
    
    /* Initialise the NullStream */
    NullStream = XLALCreateREAL8TimeSeries("ETNullStream", GPSStart, 0.0, RawData[1]->deltaT, &lalStrainUnit, RawData[1]->data->length);
    
    // SUM ALL DETECTOR RESPONSES
    for (i=0; i<3; i++) { 
			for (j=0; j<NullStream->data->length; j++) {
					NullStream->data->data[j] += RawData[i]->data->data[j];
			}
		}
   
   // PRINT NULL STREAM TO FILE (TESTING PURPOSES)
	FILE *ns_file;
	char NS_name[124];
sprintf(NS_name,"%s%d%s","ns_",GPSStart->gpsSeconds,".dat");
                        ns_file = fopen(NS_name, "w");

for (j=0; j<NullStream->data->length/512; j++) {
fprintf(ns_file,"%10.10lf %10.10e %10.10e %10.10e %10.10e\n",j*512*NullStream->deltaT,NullStream->data->data[j*512], RawData[0]->data->data[j*512], RawData[1]->data->data[j*512], RawData[2]->data->data[j*512]);
                        }

fclose(ns_file);
 
    /***************************************************************************
     *
     *  CLEANING UP
     *
     **************************************************************************/ 
    
    for (i=0; i<3; i++) {
        XLALDestroyREAL8TimeSeries(RawData[i]);
    }
    
    /***************************************************************************
     *
     *  EXIT FUNCTION
     *
     **************************************************************************/ 
    
	return NullStream;
}

REAL8TimeSeries *ReadTimeSerieFromCache(const CHAR *cachefile, const CHAR *channel, LIGOTimeGPS *start, REAL8 duration)
{
    LALStatus status;
    FrCache *cache = NULL;
    FrStream *stream = NULL;
    REAL8TimeSeries *output = NULL;
    
    cache = XLALFrImportCache(cachefile);
    if (cache==NULL) {
        fprintf(stderr,"ERROR: unable to import cache file %s\n",cachefile);
        exit(-1);
    }
    stream = XLALFrCacheOpen(cache);
    if (stream==NULL) {
        fprintf(stderr,"ERROR: unable to open stream from frame cache\n");
        exit(-1);
    }
    output = XLALFrInputREAL8TimeSeries(stream,channel,start,duration,0);
    if (output==NULL) {
        fprintf(stderr,"ERROR: unable to read channel %s from %s at time %i \n", cachefile, channel, start->gpsSeconds);
        fprintf(stderr,"Check that the specified data exist or that the duration is not too long\n");
        exit(-1);
    }
    LALDestroyFrCache(&status,&cache);
    LALFrClose(&status,&stream);
    return output;
}

void PopulateNullStream(REAL8TimeSeries *NullStream, REAL8TimeSeries *RawData)
{
    
    /***************************************************************************
     *
     *  FUNCTION TO POPULATE THE NULLSTREAM USING THE 3 DETECTOR RAW INPUT
     *
     **************************************************************************/
    
    UINT4 i=0;
    REAL8 onethird=1.0L/3.0;
    for (i=0; i<NullStream->data->length; i++) {
        NullStream->data->data[i]+=onethird*RawData->data->data[i];
    }
    
    // (TJONNIE) SUPERFLUOUS, MIGHT WANT TO DELETE
    NullStream->deltaT=RawData->deltaT;
    NullStream->sampleUnits=RawData->sampleUnits;
    NullStream->f0=RawData->f0;
}

REAL8FrequencySeries * ComputeSingleDetectorInvPSDfromNullStream(LIGOTimeGPS *GPSStart, REAL8 duration, UINT4 SampleRate, UINT4 nSegs)
{
    /***************************************************************************
     *
     *  FUNCTION TO COMPUTE THE ET INVPSD FROM THE NULLSTREAM
     *		- PSD OF THE NULL STREAM ACCORDING TO EQN 34 OF THE MDC PAPER
     *
     **************************************************************************/
    
    /***************************************************************************
     *
     *  DECLARE VARIABLES
     * 
     **************************************************************************/
    
    int check = 0;
    UINT4 j=0;
    
    REAL8 segDur = duration/(REAL8)nSegs;
    UINT4 seglen = (UINT4)(segDur*SampleRate); // GET THE LENGTH OF A SEGMENT
    segDur = seglen/SampleRate; // UPDATE WITH THE TRUNCATED SEGLEN
  	nSegs = (INT4)floor(duration/segDur); // UPDATE NSEGS TO ACCOMODATE SAMPLERATE
    UINT4 stride = seglen; /* Overlap the padding */
    //REAL8 strideDur = stride / SampleRate;
    REAL8 deltaF=(REAL8)SampleRate/seglen;
    REAL8 end_freq=10.0; /* cutoff frequency */ 
    
    REAL8Window  *windowplan = XLALCreateTukeyREAL8Window(seglen,0.1*(REAL8)8.0*SampleRate/(REAL8)seglen);;
    REAL8FFTPlan *fwdplan = XLALCreateForwardREAL8FFTPlan( seglen, 1 );
    REAL8FFTPlan *revplan = XLALCreateReverseREAL8FFTPlan( seglen, 1 );
    REAL8FrequencySeries *inverse_spectrum = XLALCreateREAL8FrequencySeries("inverse spectrum",GPSStart,0.0,deltaF,&lalDimensionlessUnit,seglen/2+1);
    
    /***************************************************************************
     *
     *  COMPUTING THE NULLSTREAM
     * 
     **************************************************************************/
     
     REAL8TimeSeries * NullStream = LALETNullStream(GPSStart,duration);
     
    // RESAMPLE TIMESERIES - OPTIONAL
    if (SampleRate!=((UINT4)(1.0/NullStream->deltaT))) {
			fprintf(stdout,"... Sample rate %d, resampling...\r",SampleRate);
			XLALResampleREAL8TimeSeries(NullStream,1.0/SampleRate);
		}
     
    /***************************************************************************
     *
     *  COMPUTING THE INVERSE PSD FROM THE NULLSTREAM
     * 
     **************************************************************************/
    
    // SHRINKING NULLSTREAM INTO INTEGER SEGMENTS
    fprintf(stdout,"... Shrinking - (lost %d samples from end)\n",NullStream->data->length-(seglen*nSegs));
    NullStream=(REAL8TimeSeries *)XLALShrinkREAL8TimeSeries(NullStream,(size_t) 0, (size_t) seglen*nSegs);
    fprintf(stdout,"... Computing power spectrum, seglen %i stride %i\n",seglen,stride);
    
    // CALCULATE THE INVERSE SPECTRUM 
		check = XLALREAL8AverageSpectrumMedian(inverse_spectrum,NullStream,(UINT4)seglen,(UINT4)stride,windowplan,fwdplan);
    if (check) {fprintf(stderr,"Failed computing the XLALREAL8AverageSpectrumMedian \n");exit(-1);}
    
    // TRUNCATE INVERSE SPECTRUM
    check = XLALREAL8SpectrumInvertTruncate(inverse_spectrum, end_freq, seglen, (seglen-stride)/4, fwdplan, revplan ); 
    if (check) {fprintf(stderr,"Failed computing XLALREAL8SpectrumInvertTruncate \n");exit(-1);}
    
    /* Normalize, ACCORDING TO EQN 35 FROM ET MDC PAPER - NB: EQN35 IS FOR THE SPECTRAL DENSITY, BUT THE INVERSE SPECTRAL DENSITY IS CALCULATED HERE */
    for (j=0;j<inverse_spectrum->data->length;j++) {
    	inverse_spectrum->data->data[j]*=3.0L;
    }    
    
    /***************************************************************************
     *
     *  PRINT INVPSD TO FILE - OPTIONAL
     * 
     **************************************************************************/
    /*
    FILE *psdout;
		//char psdname[100];
		//sprintf(psdname,"nullstream_%i.dat",SampleRate_global);
		psdout=fopen(outfile,"w");
    
    // PRINT FILE
		for(j=0;j<inverse_spectrum->data->length;j++) {
       fprintf(psdout,"%10.10lf %10.10e\n",j*deltaF,inverse_spectrum->data->data[j]); 
    }
    fclose(psdout);
    */
    
    /* write the inverse spectrum from a frequency serie to a frame file */     
    //printf("writing the inverse spectrum to a frame file\n");
    //LALFrWriteREAL8FrequencySeries(&status,inverse_spectrum,&frSerie_freq,1);
    
    return inverse_spectrum;
}
