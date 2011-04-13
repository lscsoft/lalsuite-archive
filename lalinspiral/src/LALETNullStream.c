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

NRCSID (LALETNULLSTREAMC,"$Id$");

REAL8TimeSeries * LALETNullStream (LIGOTimeGPS *GPSStart, REAL8 duration )
{
    /***************************************************************************
     *
     *  FUNCTION TO COMPUTE THE ETNULLSTREAM FROM A 3 CHANNEL CACHE FILES
     *
     **************************************************************************/
    
    /***************************************************************************
     *
     *  INITIALISE STATUSPOINTER
     *
     **************************************************************************/
    
    //INITSTATUS( status, "LALETNullStream", LALETNULLSTREAMC);
	//ATTATCHSTATUSPTR( status );
    
    /***************************************************************************
     *
     *  DECLARE VARIABLES
     *
     **************************************************************************/    
    
    UINT4 i=0;
    REAL8TimeSeries *RawData[3];
    REAL8TimeSeries *NullStream;
    
    /***************************************************************************
     *
     *  INPUT DATA - SET BY HAND CURRENTLY, MIGHT WANT TO ADD TO ARGUMENT
     *      PARSING
     *
     **************************************************************************/  
    
    const CHAR *ChannelNames[3] = {"E1:STRAIN", "E2:STRAIN", "E3:STRAIN"};
    const CHAR *CacheFileNames[3] = {"/atlas/user/scr01/tania/data/cache/E1.cache","/atlas/user/scr01/tania/data/cache/E2.cache","/atlas/user/scr01/tania/data/cache/E3.cache"};
    
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
     *  COMPUTE NULLSTREAM
     *
     **************************************************************************/  
    
    /* Initialise the NullStream */
    NullStream = XLALCreateREAL8TimeSeries("ETNullStream", GPSStart, 0.0, RawData[1]->deltaT, &lalStrainUnit, RawData[1]->data->length);
    
    for (i=0; i<3; i++) {
        /* Sum the strain to the null stream */
        PopulateNullStream(NullStream,RawData[i]);
    }
    
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

void PopulateHplus(REAL8TimeSeries *hplus, REAL8TimeSeries *RawData)
{
    /***************************************************************************
     *
     *  FUNCTION TO POPULATE THE HPLUS USING THE 3 DETECTOR RAW INPUT
     *      NB: NOT CODED YET
     *
     **************************************************************************/
    hplus=NULL;
    RawData = NULL;

    return;
}

void PopulateHcross(REAL8TimeSeries *hcross, REAL8TimeSeries *RawData)
{
    /***************************************************************************
     *
     *  FUNCTION TO POPULATE THE HCROSS USING THE 3 DETECTOR RAW INPUT
     *      NB: NOT CODED YET
     *
     **************************************************************************/
    
    hcross = NULL;
    RawData = NULL;
    
    return;
}

void ComputePSDfromNullStream(REAL8FrequencySeries *psd, REAL8TimeSeries *RawData)
{
    /***************************************************************************
     *
     *  FUNCTION TO COMPUTE THE ET PSD FROM THE NULLSTREAM
     *
     **************************************************************************/
    
    psd = NULL;
    RawData = NULL;
    
    return;
}
