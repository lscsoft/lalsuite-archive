/*
*  Copyright (C) 2007 Duncan Brown, Jolien Creighton, Lisa M. Goggin
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

#include <math.h>
#include <string.h>
#include <limits.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/LALConstants.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/AVFactories.h>
#include <lal/VectorOps.h>
#include <lal/Units.h>
#include <lal/TimeFreqFFT.h>
#include <lal/RingUtils.h>
#include <lal/Date.h>
#include <lal/LALFrameL.h>

#include "lalapps.h"
#include "lalappsfrutils.h"
#include "errutil.h"
#include "gpstime.h"
#include "ring.h"

RCSID( "$Id$" );

static REAL8 compute_template_variance(
    COMPLEX8FrequencySeries  *stilde,
    COMPLEX8FrequencySeries  *stildeSine,
    REAL4FrequencySeries     *invspec,
    REAL8                     dynRangeFac,
    int                       writeCData
    );

static int filter_segment_template(
    REAL4TimeSeries          *result,
    COMPLEX8FrequencySeries  *rtilde,
    COMPLEX8FrequencySeries  *stilde,
    COMPLEX8FrequencySeries  *segment,
    REAL4FFTPlan             *plan
    );

static int filter_segment_templateSine(
    REAL4TimeSeries          *resultSine,
    COMPLEX8FrequencySeries  *rtildeSine,
    COMPLEX8FrequencySeries  *stildeSine,
    COMPLEX8FrequencySeries  *segment,
    REAL4FFTPlan             *plan
    );

static SnglRingdownTable * find_events(
    SnglRingdownTable        *events,
    UINT4                    *numEvents,
    REAL4TimeSeries          *result,
    REAL4TimeSeries          *resultSine,
    REAL4                     tmpltSigma,
    REAL4                     tmpltFrequency,
    REAL4                     tmpltQuality,
    struct ring_params       *params
    );

typedef struct
tagCDataNode
{
  CHAR cdataStrNode[LALNameLength];
  struct tagCDataNode *next;
}
CDataNode;

int XLALAddCData( CDataNode **cdataStrCat, CHAR *id );

int XLALAddCData( CDataNode **cdataStrCat, CHAR *id )
{
      int notPresent = 1;
      int addedCData = 0;
      //CDataNode *prevCData = NULL;
      CDataNode *nextCData = NULL;
      CDataNode *thisCData = NULL;

      thisCData = *cdataStrCat;
      nextCData = (*cdataStrCat)->next;
      *cdataStrCat = NULL;

      while ( thisCData ) {
	if ( strcmp(thisCData->cdataStrNode, id ) )  {
	  notPresent *= 1;
	}
	else {
	  notPresent *= 0;
	}

	if ( ! *cdataStrCat ) {
	  *cdataStrCat = thisCData;
	}

	//prevCData = thisCData;
	thisCData = thisCData->next;
	//if ( thisCData )
	//  nextCData = thisCData->next;
      }

      verbose( "Entered XLALAddCData loop; notPresent is %d\n", notPresent);
      verbose( "Here addedCData is %d\n", addedCData); 
      verbose( "Here cdataStrCat is %s\n", (*cdataStrCat)->cdataStrNode );
      verbose( "Here cdataStr is %s\n", &id); 

      if ( notPresent ) {
	(*cdataStrCat)->next = (CDataNode *) LALCalloc( 1, sizeof(CDataNode) );
	strcpy( (*cdataStrCat)->next->cdataStrNode , id );
	addedCData = 1;
      }

      verbose( "Here addedCData AFTER check is %d\n", addedCData);
      verbose( "Here cdataStrCat is %s\n", (*cdataStrCat)->next->cdataStrNode );

      return( addedCData );
}

void XLALFindChirpCreateCoherentInput(
     COMPLEX8TimeSeries         **coherentInputData,
     COMPLEX8TimeSeries         *input,
     SnglRingdownTable          *templt,
     REAL4                      coherentSegmentLength,/*in seconds*/
     INT4                       corruptedDataLength /*in timepoints */
);

SnglRingdownTable * ring_filter(
    RingDataSegments         *segments,
    RingTemplateBank         *bank,
    REAL4FrequencySeries     *invSpectrum,
    REAL4FFTPlan             *fwdPlan,
    REAL4FFTPlan             *revPlan,
    struct ring_params       *params
    )
{
  SnglRingdownTable       *events = NULL; /* head of linked list of events */
  REAL4TimeSeries          signalvec;
  REAL4TimeSeries          signalSine;
  REAL4TimeSeries          result;
  REAL4TimeSeries          resultSine;
  COMPLEX8FrequencySeries  stilde;
  COMPLEX8FrequencySeries  stildeSine;
  COMPLEX8FrequencySeries  rtilde;
  COMPLEX8FrequencySeries  rtildeSine;

  UINT4 segmentLength;
  UINT4 sgmnt;
  UINT4 tmplt;

  if ( ! params->doFilter )
    return NULL;

  segmentLength = floor( params->segmentDuration * params->sampleRate + 0.5 );

  memset( &signalvec, 0, sizeof( signalvec ) );
  memset( &signalSine, 0, sizeof( signalSine ) );
  memset( &result, 0, sizeof( result ) );
  memset( &resultSine, 0, sizeof( resultSine ) );
  memset( &stilde, 0, sizeof( stilde ) );
  memset( &stildeSine, 0, sizeof( stildeSine ) );
  memset( &rtilde, 0, sizeof( rtilde ) );
  memset( &rtildeSine, 0, sizeof( rtildeSine ) );

  signalvec.deltaT = 1.0/params->sampleRate;
  signalSine.deltaT = 1.0/params->sampleRate;
  signalvec.sampleUnits = lalStrainUnit;
  signalSine.sampleUnits = lalStrainUnit;
  rtilde.deltaF = 1.0 / params->segmentDuration;
  rtildeSine.deltaF = 1.0 / params->segmentDuration;
  signalvec.data = XLALCreateREAL4Vector( segmentLength );
  signalSine.data = XLALCreateREAL4Vector( segmentLength );
  stilde.data = XLALCreateCOMPLEX8Vector( segmentLength/2 + 1 );
  stildeSine.data = XLALCreateCOMPLEX8Vector( segmentLength/2 + 1 );
  rtilde.data = XLALCreateCOMPLEX8Vector( segmentLength/2 + 1 );
  rtildeSine.data = XLALCreateCOMPLEX8Vector( segmentLength/2 + 1 );
  resultSine.data = XLALCreateREAL4Vector( segmentLength );
  result.data = XLALCreateREAL4Vector( segmentLength );

  /* loop over all elements in the template bank */
  for ( tmplt = 0; tmplt < bank->numTmplt; ++tmplt )
  {
    SnglRingdownTable *thisTmplt = bank->tmplt + tmplt;
    UINT4 numEvents = 0;
    REAL8 sigma;

    verbose( "creating template %d\n", tmplt );

    /* make template and fft it */
    XLALComputeRingTemplate( &signalvec, thisTmplt );
    /* write template time series if requested */
    if ( params->writeTemplateTimeSeries )
    {
     snprintf( signalvec.name, sizeof(signalvec.name), "TMPLT_%u", tmplt );
     write_REAL4TimeSeries( &signalvec );
    }
    XLALREAL4TimeFreqFFT( &stilde, &signalvec, fwdPlan );
    /* write template fft if requested */
    if ( params->writeTemplateFFT )
    {
     snprintf( stilde.name, sizeof(stilde.name), "TMPLT_%u_FFT", tmplt );
     write_COMPLEX8FrequencySeries( &stilde );
    }

    /* make Sinetemplate and fft it */
    XLALComputeRingTemplateSine( &signalSine, thisTmplt );
    /* write Sinetemplate time series if requested */
    if ( params->writeTemplateTimeSeries )
    {
     snprintf( signalSine.name, sizeof(signalSine.name), "TMPLT_Sine_%u", tmplt );
     write_REAL4TimeSeries( &signalSine );
    }
    XLALREAL4TimeFreqFFT( &stildeSine, &signalSine, fwdPlan );
    /* write Sinetemplate fft if requested */
    if ( params->writeTemplateFFT )
    {
     snprintf( stildeSine.name, sizeof(stildeSine.name), "TMPLT_Sine_%u_FFT", tmplt );
     write_COMPLEX8FrequencySeries( &stildeSine );
    }


    /* compute sigma for this template */
    sigma = sqrt( compute_template_variance( &stilde, &stildeSine, invSpectrum,
          params->dynRangeFac, params->writeCData ) );
    /* loop over segments */
    for ( sgmnt = 0; sgmnt < segments->numSgmnt; ++sgmnt )
    {
      verbose( "  filtering segment %d against template %d\n", sgmnt, tmplt );

      /* filter the segment with the template */
      filter_segment_template( &result, &rtilde, &stilde,
          &segments->sgmnt[sgmnt], revPlan );

      if ( params->writeCData ) {
        /* filter the segment with the Sinetemplate */
        filter_segment_templateSine( &resultSine, &rtildeSine, &stildeSine,
          &segments->sgmnt[sgmnt], revPlan );
      }

      /* search through results for threshold crossings and record events */
      events = find_events( events, &numEvents, &result, &resultSine, sigma,
          thisTmplt->frequency, thisTmplt->quality, params );

      /* write filter output if requested */
      if ( params->writeFilterOutput )
      { /* guess we better normalize it so it is SNR-like... */
        REAL4 snrFactor = 2 * params->dynRangeFac / sigma;
        UINT4 k;
        for ( k = 0; k < result.data->length; ++k )
          result.data->data[k] *= snrFactor;
        write_REAL4TimeSeries( &result );

        if ( params->writeCData ) {
          for ( k = 0; k < resultSine.data->length; ++k )
            resultSine.data->data[k] *= snrFactor;
          write_REAL4TimeSeries( &resultSine );
        }
      }
    }
    params->numEvents += numEvents;
    verbose( "found %u event%s in template %u\n", numEvents,
        numEvents == 1 ? "" : "s", tmplt );
  }

  verbose( "found %u event%s\n", params->numEvents,
      params->numEvents == 1 ? "" : "s" );

  /* cleanup */
  XLALDestroyREAL4Vector( result.data );
  XLALDestroyREAL4Vector( resultSine.data );
  XLALDestroyCOMPLEX8Vector( rtilde.data );
  XLALDestroyCOMPLEX8Vector( rtildeSine.data );
  XLALDestroyCOMPLEX8Vector( stilde.data );
  XLALDestroyCOMPLEX8Vector( stildeSine.data );
  XLALDestroyREAL4Vector( signalvec.data );
  XLALDestroyREAL4Vector( signalSine.data );

  return events;
}

SnglRingdownTable * ring_filter_cdata(
    RingDataSegments         *segments,
    SnglRingdownTable        *bankHead,
    REAL4FrequencySeries     *invSpectrum,
    REAL4FFTPlan             *fwdPlan,
    REAL4FFTPlan             *revPlan,
    FrameHNode              **coherentFrames,
    struct ring_params       *params
    )
{
  SnglRingdownTable       *events = NULL; /* head of linked list of events */
  SnglRingdownTable       *bank = NULL;
  REAL4TimeSeries          signalvec;
  REAL4TimeSeries          signalSine;
  REAL4TimeSeries          result;
  REAL4TimeSeries          resultSine;
  COMPLEX8FrequencySeries  stilde;
  COMPLEX8FrequencySeries  stildeSine;
  COMPLEX8FrequencySeries  rtilde;
  COMPLEX8FrequencySeries  rtildeSine;
  COMPLEX8TimeSeries       cdataSeries;

  UINT4 segmentLength;
  UINT4 sgmnt;
  UINT4 tmplt;
  REAL8 trigTime = 0.0;
  REAL8 lowerBound = 0.0;
  REAL8 upperBound = 0.0;
  CHAR  cdataStr[LALNameLength];
  INT4  CDataAdded = 0;
  CDataNode  *cdataStrCat = NULL;
  CDataNode  *thisCDataStr = NULL;
  COMPLEX8TimeSeries           *coherentInputData = NULL;
  FrameHNode *thisCoherentFrame = NULL;
  struct FrameH *outFrame   = NULL;

  if ( ! params->doFilter )
    return NULL;

  segmentLength = floor( params->segmentDuration * params->sampleRate + 0.5 );

  memset( &signalvec, 0, sizeof( signalvec ) );
  memset( &signalSine, 0, sizeof( signalSine ) );
  memset( &result, 0, sizeof( result ) );
  memset( &resultSine, 0, sizeof( resultSine ) );
  memset( &stilde, 0, sizeof( stilde ) );
  memset( &stildeSine, 0, sizeof( stildeSine ) );
  memset( &rtilde, 0, sizeof( rtilde ) );
  memset( &rtildeSine, 0, sizeof( rtildeSine ) );
  memset( &cdataSeries, 0, sizeof( cdataSeries ) );

  signalvec.deltaT = 1.0/params->sampleRate;
  signalSine.deltaT = 1.0/params->sampleRate;
  signalvec.sampleUnits = lalStrainUnit;
  signalSine.sampleUnits = lalStrainUnit;
  rtilde.deltaF = 1.0 / params->segmentDuration;
  rtildeSine.deltaF = 1.0 / params->segmentDuration;
  signalvec.data = XLALCreateREAL4Vector( segmentLength );
  signalSine.data = XLALCreateREAL4Vector( segmentLength );
  stilde.data = XLALCreateCOMPLEX8Vector( segmentLength/2 + 1 );
  stildeSine.data = XLALCreateCOMPLEX8Vector( segmentLength/2 + 1 );
  rtilde.data = XLALCreateCOMPLEX8Vector( segmentLength/2 + 1 );
  rtildeSine.data = XLALCreateCOMPLEX8Vector( segmentLength/2 + 1 );
  resultSine.data = XLALCreateREAL4Vector( segmentLength );
  result.data = XLALCreateREAL4Vector( segmentLength );
  cdataSeries.data = XLALCreateCOMPLEX8Vector( segmentLength );
  cdataSeries.deltaT = 1.0/params->sampleRate;

  /* loop over all elements in the template bank */
  for ( bank = bankHead,  tmplt = 0; bank ; bank = bank->next, ++tmplt )
  {
    UINT4 numEvents = 0;
    REAL8 sigma;

    /* make template and fft it */
    XLALComputeRingTemplate( &signalvec, bank );
    XLALREAL4TimeFreqFFT( &stilde, &signalvec, fwdPlan );

    /* make Sinetemplate and fft it */
    XLALComputeRingTemplateSine( &signalSine, bank );
    XLALREAL4TimeFreqFFT( &stildeSine, &signalSine, fwdPlan );

    /* compute sigma for this template */
    sigma = sqrt( compute_template_variance( &stilde, &stildeSine, invSpectrum,
          params->dynRangeFac, params->writeCData ) );
    /* loop over segments */
    for ( sgmnt = 0; sgmnt < segments->numSgmnt; ++sgmnt )
    {

      if ( params->writeCData ) {
        INT8 fcSegStartTimeNS;
        INT8 fcSegEndTimeNS;
        int temp_frequency = 0;
        int temp_quality = 0;
        REAL8 buffer = 0.0; /*CHECK: 64.0;*/

        fcSegStartTimeNS = XLALGPSToINT8NS( &(segments->sgmnt[sgmnt].epoch) );
        fcSegEndTimeNS = fcSegStartTimeNS + (INT8) ( 1e9 * 256 );

        memcpy( &(cdataSeries.epoch), &(segments->sgmnt[sgmnt].epoch), sizeof(LIGOTimeGPS) );

        trigTime = bank->start_time.gpsSeconds + 1e-9 *
                     bank->start_time.gpsNanoSeconds;

        lowerBound = rint(fcSegStartTimeNS/1000000000L) + buffer;
        upperBound = rint(fcSegEndTimeNS/1000000000L) - buffer;

        verbose( "  filtering segment %d against template %d\n", sgmnt, tmplt );

        verbose( "GPS end time of bankCurrent in s and ns are %d and %d; trigtime in (s) is %12.3f; lower and upper bounds in (s) is %12.3f, %12.3f\n",
                      bank->start_time.gpsSeconds,bank->start_time.gpsNanoSeconds, trigTime, lowerBound, upperBound);

        if ( trigTime >= lowerBound && trigTime <= upperBound ) {
          verbose("ENTERED cdata LOOP\n");

          /* filter the segment with the template */
          filter_segment_template( &result, &rtilde, &stilde,
            &segments->sgmnt[sgmnt], revPlan );
          /* filter the segment with the Sinetemplate */
          filter_segment_templateSine( &resultSine, &rtildeSine, &stildeSine,
            &segments->sgmnt[sgmnt], revPlan );

          /* search through results for threshold crossings and record events */
          events = find_events( events, &numEvents, &result, &resultSine, sigma,
            bank->frequency, bank->quality, params );

          temp_frequency = floor( bank->frequency * 1000.0 );
          temp_quality = floor( bank->quality * 1000.0 );
          snprintf( cdataStr, LALNameLength*sizeof(CHAR),
            "%d_%d_%d_%d", bank->start_time.gpsSeconds,
            (bank->start_time.gpsNanoSeconds - (bank->start_time.gpsNanoSeconds % 1000000))/1000000,
            temp_frequency, temp_quality );

          verbose( "cdataStr IS %s\n", cdataStr );

          /* guess we better normalize it so it is SNR-like... */
          REAL4 snrFactor = 2 * params->dynRangeFac / sigma;
          UINT4 k;
          for ( k = 0; k < result.data->length; ++k ) {
            result.data->data[k] *= snrFactor;
            resultSine.data->data[k] *= snrFactor;
            cdataSeries.data->data[k].re = result.data->data[k];
            cdataSeries.data->data[k].im = resultSine.data->data[k];
          }
	  if ( !cdataStrCat ) {
	    XLALFindChirpCreateCoherentInput(
                        &coherentInputData, &cdataSeries,
			bank, params->CDataLength/2,
			(INT4) rint(0.5*params->CDataLength*params->sampleRate) );

	    thisCDataStr = cdataStrCat = (CDataNode *) LALCalloc( 1, sizeof(CDataNode) );
	    strcpy( thisCDataStr->cdataStrNode , cdataStr );
	    verbose("tempcdataStrNode is %s\n", thisCDataStr->cdataStrNode);
	    CDataAdded = 1;
	  }
	  else {
	    CDataAdded = 0;
            verbose( "The event id in ELSE is %" LAL_UINT8_FORMAT "\n",
                                        bank->event_id->id);
	    CDataAdded = XLALAddCData( &cdataStrCat, cdataStr );
	    if ( CDataAdded ) {
	      XLALFindChirpCreateCoherentInput(
                          &coherentInputData, &cdataSeries,
                          bank, params->CDataLength/2,
			  (INT4) rint(0.5*params->CDataLength*params->sampleRate) );
	    }
	  }

	  if ( CDataAdded && coherentInputData )
	    {
              verbose( "coherentInputData Name string is %s\n",
				    coherentInputData->name);

	      snprintf( coherentInputData->name,
			LALNameLength*sizeof(CHAR),
			"%s:CBC-CData", params->ifoName );

	      if ( ! *coherentFrames )
		{
		  thisCoherentFrame = *coherentFrames = (FrameHNode *)
		    LALCalloc( 1, sizeof(FrameHNode) );
		}
	      else
		{
		  thisCoherentFrame = thisCoherentFrame->next =
		    (FrameHNode *) LALCalloc( 1, sizeof(FrameHNode) );
		}

	      thisCoherentFrame->frHeader = fr_add_proc_COMPLEX8TimeSeries(
		  outFrame, coherentInputData, "none", cdataStr );

              verbose( "GPS end time in s and ns are: %d and %d\n",
                            coherentInputData->epoch.gpsSeconds, coherentInputData->epoch.gpsNanoSeconds);
              verbose( "Event ID used for C Data is %" LAL_UINT8_FORMAT "\n",
                            bank->event_id->id);
              verbose( "C Data string is %s\n",
			    cdataStr);
              verbose( "coherentInputData Name string is %s\n",
			    coherentInputData->name);

	      XLALDestroyCOMPLEX8Vector( coherentInputData->data );
	      LALFree( coherentInputData );
	      coherentInputData = NULL;
	    }

        } /*Closes writeCData */
      }
      else {
        /* filter the segment with the template */
        filter_segment_template( &result, &rtilde, &stilde,
          &segments->sgmnt[sgmnt], revPlan );
        /* filter the segment with the Sinetemplate */
        filter_segment_templateSine( &resultSine, &rtildeSine, &stildeSine,
          &segments->sgmnt[sgmnt], revPlan );

        /* search through results for threshold crossings and record events */
        events = find_events( events, &numEvents, &result, &resultSine, sigma,
          bank->frequency, bank->quality, params );

        /* write filter output if requested */
        if ( params->writeFilterOutput )
        { /* guess we better normalize it so it is SNR-like... */
          REAL4 snrFactor = 2 * params->dynRangeFac / sigma;
          UINT4 k;
          for ( k = 0; k < result.data->length; ++k ) {
            result.data->data[k] *= snrFactor;
            resultSine.data->data[k] *= snrFactor;
            cdataSeries.data->data[k].re = result.data->data[k];
            cdataSeries.data->data[k].im = resultSine.data->data[k];
          }
        }
      }
    }

    params->numEvents += numEvents;
    verbose( "found %u event%s in template %u\n", numEvents,
        numEvents == 1 ? "" : "s", tmplt );
  }

  verbose( "found %u event%s\n", params->numEvents,
      params->numEvents == 1 ? "" : "s" );

  /* cleanup */
  XLALDestroyREAL4Vector( result.data );
  XLALDestroyREAL4Vector( resultSine.data );
  XLALDestroyCOMPLEX8Vector( rtilde.data );
  XLALDestroyCOMPLEX8Vector( rtildeSine.data );
  XLALDestroyCOMPLEX8Vector( stilde.data );
  XLALDestroyCOMPLEX8Vector( stildeSine.data );
  XLALDestroyREAL4Vector( signalvec.data );
  XLALDestroyREAL4Vector( signalSine.data );
  XLALDestroyCOMPLEX8Vector( cdataSeries.data );


  return events;
}

static REAL8 compute_template_variance(
    COMPLEX8FrequencySeries  *stilde,
    COMPLEX8FrequencySeries  *stildeSine,
    REAL4FrequencySeries     *invspec,
    REAL8                     dynRangeFac,
    int                       writeCData
    )
{
  UINT4 k;
  REAL8 var;

  var = 0;
  for ( k = 0; k < stilde->data->length; ++k )
  {
    REAL8 re = stilde->data->data[k].re;
    REAL8 im = stilde->data->data[k].im;
    if ( writeCData ) {
      REAL8 reSine = stildeSine->data->data[k].re;
      REAL8 imSine = stildeSine->data->data[k].im;
      var += invspec->data->data[k] * (re*re + im*im + reSine*reSine + imSine*imSine);
    }
    else {
      var += invspec->data->data[k] * (re*re + im*im);
    }
  }

  var *= 4.0 * dynRangeFac * dynRangeFac * stilde->deltaF;
  return var;
}

static int filter_segment_template(
    REAL4TimeSeries          *result,
    COMPLEX8FrequencySeries  *rtilde,
    COMPLEX8FrequencySeries  *stilde,
    COMPLEX8FrequencySeries  *segment,
    REAL4FFTPlan             *plan
    )
{
  char *s;

  /* name of rtilde */
  snprintf( rtilde->name, sizeof( rtilde->name ), "%s_%s",
      segment->name, stilde->name );
  /* name of result is the same but without the _FFT */
  strncpy( result->name, rtilde->name, sizeof( result->name ) - 1 );
  /* make sure that the string ends with _FFT */
  s = result->name + strlen( result->name ) - strlen( "_FFT" );
  if ( strcmp( s, "_FFT" ) == 0 )
    *s = 0; /* it does: terminate here */

  /* multiply segment by filter and store in fft of result */
  XLALCCVectorMultiplyConjugate( rtilde->data, segment->data, stilde->data );
  XLALUnitMultiply( &rtilde->sampleUnits, &segment->sampleUnits, &stilde->sampleUnits );
  rtilde->epoch = segment->epoch;
/* write_COMPLEX8FrequencySeries( rtilde ); */

  /* inverse fft to obtain result */
  XLALREAL4FreqTimeFFT( result, rtilde, plan );
/* write_REAL4TimeSeries( result ); */
  result->epoch = rtilde->epoch;
  return 0;
}

static int filter_segment_templateSine(
    REAL4TimeSeries          *resultSine,
    COMPLEX8FrequencySeries  *rtildeSine,
    COMPLEX8FrequencySeries  *stildeSine,
    COMPLEX8FrequencySeries  *segment,
    REAL4FFTPlan             *plan
    )
{
  char *s;

  /* name of rtildeSine */
  snprintf( rtildeSine->name, sizeof( rtildeSine->name ), "%s_%s",
      segment->name, stildeSine->name );
  /* name of resultSine is the same but without the _FFT */
  strncpy( resultSine->name, rtildeSine->name, sizeof( resultSine->name ) - 1 );
  /* make sure that the string ends with _FFT */
  s = resultSine->name + strlen( resultSine->name ) - strlen( "_FFT" );
  if ( strcmp( s, "_FFT" ) == 0 )
    *s = 0; /* it does: terminate here */

  /* multiply segment by Sinefilter and store in fft of resultSine */
  XLALCCVectorMultiplyConjugate( rtildeSine->data, segment->data, stildeSine->data );
  XLALUnitMultiply( &rtildeSine->sampleUnits, &segment->sampleUnits, &stildeSine->sampleUnits );
  rtildeSine->epoch = segment->epoch;
/* write_COMPLEX8FrequencySeries( rtildeSine ); */

  /* inverse fft to obtain result */
  XLALREAL4FreqTimeFFT( resultSine, rtildeSine, plan );
/* write_REAL4TimeSeries( result ); */
  resultSine->epoch = rtildeSine->epoch;
  return 0;
}

/* NOTE: numEvents must be number of events _FROM_CURRENT_TEMPLATE_ so far. */
/* It must be set to zero when filtering against a new template is started. */
static SnglRingdownTable * find_events(
    SnglRingdownTable  *events,
    UINT4              *numEvents,
    REAL4TimeSeries    *result,
    REAL4TimeSeries    *resultSine,
    REAL4               tmpltSigma,
    REAL4               tmpltFrequency,
    REAL4               tmpltQuality,
    struct ring_params *params
    )
{
  const REAL4 efolds = 10.0; /* number of efolds of ringdown in template */
  SnglRingdownTable *thisEvent = events; /* the current event */
  REAL4 snrFactor; /* factor to convert from filter result to snr */
  REAL4 threshold; /* modified threshold on filter result (rather than snr) */
  REAL4 thresholdStat = 0.0; /* statistic on which the threshold is applied */
  REAL4 filterDuration;
  INT8  lastTimeNS;
  INT8  gapTimeNS;
  UINT4 segmentStride;
  UINT4 eventCount = 0;
  UINT4 jmin;
  UINT4 jmax;
  UINT4 j;
  INT8  t0 = 0;
  INT8  tpeak = 0;
  INT8  tmp;
  /* LALMSTUnitsAndAcc     gmstUnits = { MST_HRS, LALLEAPSEC_STRICT };*/

  /* compute filter duration: sum of rindown duration and spec trunc duration */
  filterDuration  = efolds * tmpltQuality / (LAL_PI * tmpltFrequency);
  filterDuration += params->truncateDuration;

  /* gap time: filter duration unless explicitly specified */
  if ( params->maximizeEventDuration < 0 ) /* maximize based on duration*/
    gapTimeNS = sec_to_ns( filterDuration );
  else /* maximize with specified duration */
    gapTimeNS = sec_to_ns( params->maximizeEventDuration );

  /* determine if we are maximizing over current event or not */
  /* set lastTimeNS to be time of last event or before any possible */
  /* event in this segment depending on whether there was a previous */
  /* event from this template */
  if ( *numEvents ) /* continue maximizing over the most recent event */
    lastTimeNS = epoch_to_ns( &thisEvent->start_time );
  else /* set last time to a time before any possible event in this segment */
    lastTimeNS = epoch_to_ns( &result->epoch ) - gapTimeNS - (INT8)(1000000000);

  /* compute modified threshold on filter output rather than snr */
  snrFactor = 2 * params->dynRangeFac / tmpltSigma;
  threshold = params->threshold / snrFactor;

  /* compute start and stop index for scanning time series */
  segmentStride = floor( params->strideDuration / result->deltaT + 0.5 );
  jmin = segmentStride/2;
  jmax = jmin + segmentStride;

  for ( j = jmin; j < jmax; ++j ) {
    REAL4 phase;

    if ( params->writeCData ) {
      thresholdStat = sqrt( pow(result->data->data[j],2)
                              + pow(resultSine->data->data[j],2) );
      phase = atan2( resultSine->data->data[j], result->data->data[j]);
    }
    else {
      thresholdStat = fabs( result->data->data[j] );
      phase = 0.0;
    }

    if ( thresholdStat > threshold ) /* threshold crossing */
    {
      REAL4 snr;
      INT8  timeNS;
      REAL4 sigma;
      REAL4 amp;

      snr     = thresholdStat * snrFactor;
      timeNS  = epoch_to_ns( &result->epoch );
      timeNS += sec_to_ns( j * result->deltaT );

      if ( timeNS > lastTimeNS + gapTimeNS ) /* new distinct event */
      {
        /* create a new node on the linked list */
        thisEvent = LALCalloc( 1, sizeof( *thisEvent ) );
        thisEvent->next = events;
        events = thisEvent;
        t0 = timeNS;
        tpeak = timeNS;

        /* copy general information about the filter */
        strncpy( thisEvent->ifo, params->ifoName, sizeof( thisEvent->ifo ) );
        strncpy( thisEvent->channel, strchr( params->channel, ':') + 1,
            sizeof( thisEvent->channel ) );
/*       LAL_CALL( LALGPStoGMST1( &status, &(thisEvent->start_time_gmst),
                        &(thisEvent->start_time), &gmstUnits ), &status); */
        /*this isnt working properly, will just leave awhile, as it is not
         * needed */

        thisEvent->frequency = tmpltFrequency;
        thisEvent->quality = tmpltQuality;
        thisEvent->mass = XLALBlackHoleRingMass( thisEvent->frequency, thisEvent->quality);
        thisEvent->spin = XLALBlackHoleRingSpin( thisEvent->quality);

        /* specific information about this threshold crossing */
        ns_to_epoch( &thisEvent->start_time, timeNS );
        thisEvent->snr = snr;

        amp = XLALBlackHoleRingAmplitude( thisEvent->frequency, thisEvent->quality, 0.5, 0.01 );
        sigma=tmpltSigma * amp;
        thisEvent->sigma_sq = pow(sigma, 2.0);
        thisEvent->eff_dist = sigma / thisEvent->snr;
        thisEvent->amplitude=amp;
        thisEvent->epsilon = 0.0;  /* borrowing these this column */
        thisEvent->phase = phase;    /* and this one */
        ++eventCount;
      }
      else if ( snr > thisEvent->snr ) /* maximize within a set of crossings */
      {
        /* update to specific information about this threshold crossing */
        ns_to_epoch( &thisEvent->start_time, timeNS );
        tpeak=timeNS;
        thisEvent->snr = snr;

        amp = XLALBlackHoleRingAmplitude( thisEvent->frequency, thisEvent->quality, 0.5, 0.01 );
        sigma=tmpltSigma * amp;
        thisEvent->eff_dist = sigma / thisEvent->snr;
        thisEvent->amplitude = amp;
        /* time clustered before the peak snr*/
        tmp = timeNS - t0;
        thisEvent->epsilon =(REAL4)tmp /1000000000;
      }

      /* update last threshold crossing time */
      lastTimeNS = timeNS;
      /* time clustered after the peak snr*/
      tmp = timeNS - tpeak;
      thisEvent->phase = phase;
    }

  } /* Closes loop over j */

  *numEvents += eventCount;
  return events;
}

void XLALFindChirpCreateCoherentInput(
     COMPLEX8TimeSeries         **coherentInputData,
     COMPLEX8TimeSeries         *input,
     SnglRingdownTable          *templt,
     REAL4                      coherentSegmentLength,/*in seconds*/
     INT4                       corruptedDataLength /*in timepoints */
     )
{
  COMPLEX8TimeSeries      *cohInputData = NULL;
  LIGOTimeGPS              start_time;
  INT4                     numPoints = 0;
  REAL4                    cohSegLength = 0.0;
  INT4                     inputEpochSeconds = 0;
  INT4                     inputEpochNanoSeconds = 0;
  REAL8                    deltaT = 0.0;
  REAL8                    intpart = 0.0;
  REAL8                    fracpart = 0.0;
  REAL8                    tempTime = 0.0;
  INT4                     crupDataLength = 0;
  INT4                     cohSegStart = 0;
  INT4                     nonCorruptEnd = 0;
  INT4                     nonCorruptStart = 0;
  INT4                     fullCohSegLength = 0;
  INT4                     eventTimePoint = 0;

  /* Ensure that arguments are reasonable */

  /* check that the output handle exists, but points to a null pointer */
  start_time = templt->start_time;
  numPoints = input->data->length;
  cohSegLength = coherentSegmentLength;
  inputEpochSeconds = input->epoch.gpsSeconds;
  inputEpochNanoSeconds = input->epoch.gpsNanoSeconds;
  deltaT = input->deltaT;
  crupDataLength = corruptedDataLength;

  /* Now determine if the event lies in corrupted data */

  nonCorruptStart = crupDataLength;
  nonCorruptEnd = numPoints - crupDataLength;
  cohSegStart = (INT4) rint( ((start_time.gpsSeconds + start_time.gpsNanoSeconds*1.0e-9) - (inputEpochSeconds + inputEpochNanoSeconds*1.0e-9)-cohSegLength)/deltaT - 1.0);

  eventTimePoint = (INT4) rint( ((start_time.gpsSeconds + start_time.gpsNanoSeconds*1.0e-9) - (inputEpochSeconds + inputEpochNanoSeconds*1.0e-9))/deltaT - 1.0);

  verbose( "templt->start_time, start_time.gpsNanoSeconds, numPoints, cohSegLength, deltaT are %d, %d, %d, %f, %f\n", start_time.gpsSeconds, start_time.gpsNanoSeconds, numPoints, cohSegLength, deltaT );

  verbose( "cohSegStart, nonCorruptStart, eventTimePoint, nonCorruptEnd are %d, %d, %d, %d\n", cohSegStart, nonCorruptStart, eventTimePoint, nonCorruptEnd );

  if( eventTimePoint < nonCorruptEnd && eventTimePoint >= nonCorruptStart )
    {
      /* The coherent chunk is outside of the corrupted data */

      cohInputData = *coherentInputData = (COMPLEX8TimeSeries *)
	LALCalloc(1, sizeof(COMPLEX8TimeSeries) );

      fullCohSegLength = 2*cohSegLength/deltaT;

      cohInputData->data = XLALCreateCOMPLEX8Vector( fullCohSegLength );

      /* store C-data snippet and its associated info */
      memcpy(cohInputData->name, input->name, LALNameLength * sizeof(CHAR) );
      memcpy(cohInputData->data->data, &(input->data->data[cohSegStart]), fullCohSegLength * sizeof(COMPLEX8) );
      cohInputData->deltaT = deltaT;
      /*CHECK: Below is a temporary fix for communicating sigmasq
	from inspiral to coherent_inspiral via the f0 member of
	the COMPLEX8TimeSeries structure*/
      cohInputData->f0 = (REAL8) templt->sigma_sq;
      tempTime = inputEpochSeconds + inputEpochNanoSeconds*1.0e-9 + cohSegStart * deltaT;
      fracpart = modf(tempTime, &intpart);
      cohInputData->epoch.gpsSeconds = (INT4) intpart;
      cohInputData->epoch.gpsNanoSeconds = (INT4) (fracpart*1.0e9);
    }
  else if ( eventTimePoint <= numPoints && eventTimePoint >= 0 )
    {
      /* return cohInputData with length 2*cohSegLength, 
	 with strings of zeros at the beginning and end */
      int j = 0;
      INT4 eventTimeFromSegEnd = 0;
      INT4 cohSegTimePointsBy2 = 0;
      INT4 cohInputDataStart = 0;
      INT4 cohInputDataLength = 0;

      cohSegTimePointsBy2 = floor(cohSegLength/deltaT);
      eventTimeFromSegEnd = numPoints - eventTimePoint;

      cohInputData = *coherentInputData = (COMPLEX8TimeSeries *)
	LALCalloc(1, sizeof(COMPLEX8TimeSeries) );

      fullCohSegLength = 2*cohSegLength/deltaT;

      cohInputData->data = XLALCreateCOMPLEX8Vector( fullCohSegLength );

      /* store C-data snippet and its associated info */
      memcpy(cohInputData->name, input->name, LALNameLength * sizeof(CHAR) );

      verbose( "Initializing coherent input COMPLEX8 vector");

      for(j=0;j<fullCohSegLength;j++)
	{
	  cohInputData->data->data[j].re = 0.0;
	  cohInputData->data->data[j].im = 0.0;
	}

      if ( ( eventTimeFromSegEnd <= cohSegTimePointsBy2 ) )
	{
          verbose( "Deriving parameters for c-data snippet");
	  cohSegStart = eventTimePoint - eventTimeFromSegEnd;
	  cohInputDataStart = cohSegTimePointsBy2 - eventTimeFromSegEnd;
	  cohInputDataLength = 2*eventTimeFromSegEnd;

	  verbose( "When eventTimeFromSegEnd <= cohSegTimePointsBy2, then cohInputDataStart, cohSegStart, cohInputDataLength are %d, %d, %d\n", cohInputDataStart, cohSegStart, cohInputDataLength );

	  memcpy(cohInputData->data->data + cohInputDataStart, &(input->data->data[cohSegStart]), cohInputDataLength * sizeof(COMPLEX8) );

	}
      else if ( ( eventTimePoint <= cohSegTimePointsBy2 ) )
	{
	  cohSegStart = 0;
	  cohInputDataStart = cohSegTimePointsBy2 - eventTimePoint;
	  cohInputDataLength = 2*eventTimePoint;

	  verbose( "When eventTimePoint <= cohSegTimePointsBy2, then cohInputDataStart, eventTimePoint, cohInputDataLength are %d, %d, %d\n", cohInputDataStart, eventTimePoint, cohInputDataLength );

	  if ( ( eventTimePoint <= cohSegTimePointsBy2 ) )
	    memcpy(cohInputData->data->data + cohInputDataStart, &(input->data->data[cohSegStart]), cohInputDataLength * sizeof(COMPLEX8) );

	}

      cohInputData->deltaT = deltaT;
      /*CHECK: Below is a temporary fix for communicating sigmasq
	from inspiral to coherent_inspiral via the f0 member of
	the COMPLEX8TimeSeries structure*/
      cohInputData->f0 = (REAL8) templt->sigma_sq;
      tempTime = inputEpochSeconds + inputEpochNanoSeconds*1.0e-9 + cohSegStart * deltaT;
      fracpart = modf(tempTime, &intpart);
      cohInputData->epoch.gpsSeconds = (INT4) intpart;
      cohInputData->epoch.gpsNanoSeconds = (INT4) (fracpart*1.0e9);

    }
  else
    {
      /* return a null pointer cohInputData */
    }
}
