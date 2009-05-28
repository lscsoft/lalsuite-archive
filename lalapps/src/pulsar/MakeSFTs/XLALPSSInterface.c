#include <string.h>
#include "XLALPSSInterface.h"
#include "lal/XLALError.h"
#include "lal/Units.h"
#include "lal/Date.h"

PSSEventParams *XLALCreatePSSEventParams(UINT4 length) { 
  PSSEventParams*ep = NULL;
  if ( length == 0 )
    XLAL_ERROR_NULL( "XLALCreatePSSEventParams", XLAL_EINVAL ); 
  ep = crea_evenparam(length);
  if ( ep == NULL )
    XLAL_ERROR_NULL( "XLALCreatePSSEventParams", XLAL_EFAULT ); 
  return ep;
} 

void XLALDestroyPSSEventParams(PSSEventParams *ts) {
  if ( !ts )
    XLAL_ERROR_VOID( "XLALCreatePSSTimeseries", XLAL_EFAULT );
  if(ts->xamed)
    free(ts->xamed);
  if(ts->xastd)
    free(ts->xastd);
  if(ts->begin)
    free(ts->begin);
  if(ts->duration)
    free(ts->duration);
  if(ts->imax)
    free(ts->imax);
  if(ts->crmax)
    free(ts->crmax);
  if(ts->ener)
    free(ts->ener);
  free(ts);
}

REAL4TimeSeries *XLALConvertREAL8TimeSeriesToREAL4TimeSeries(REAL4TimeSeries *r4ts, REAL8TimeSeries *r8ts) {
  UINT4 i;
  if (r4ts->data->length != r8ts->data->length)
    XLAL_ERROR_NULL( "XLALCreatePSSEventParams", XLAL_EINVAL ); 
  for(i = 0; i < r8ts->data->length; i++)
    r4ts->data->data[i] = r8ts->data->data[i];
  strncpy(r4ts->name, r8ts->name, LALNameLength);
  r4ts->epoch       = r8ts->epoch;
  r4ts->deltaT      = r8ts->deltaT;
  r4ts->f0          = r8ts->f0;
  r4ts->sampleUnits = r8ts->sampleUnits;
  return r4ts;
}

PSSTimeseries *XLALCreatePSSTimeseries(UINT4 length) {
  PSSTimeseries *pssGD;
  if ( length == 0 )
    XLAL_ERROR_NULL( "XLALCreatePSSTimeseries", XLAL_EINVAL );
  pssGD = crea_gd(length,0,0,"byXLALCreatePSSTimeseries");
  if( !pssGD )
    XLAL_ERROR_NULL( "XLALCreatePSSTimeseries", XLAL_EFAULT );
  return(pssGD);
}

void XLALDestroyPSSTimeseries(PSSTimeseries *ts) {
  if ( !ts )
    XLAL_ERROR_VOID( "XLALDestroyPSSTimeseries", XLAL_EFAULT );
  if(ts->x)
    free(ts->x);
  if(ts->y)
    free(ts->y);
  if(ts->capt)
    free(ts->capt);
  free(ts);
}

REAL8TimeSeries *XLALConvertPSSTimeseriesToREAL8Timeseries(REAL8TimeSeries *ts, PSSTimeseries *tsPSS) {
  UINT4 i;

  if ( !tsPSS || !ts )
    XLAL_ERROR_NULL( "XLALConvertPSSTimeseriesToREAL8Timeseries", XLAL_EFAULT );
  if (tsPSS->n != (int)ts->data->length)
    XLAL_ERROR_NULL( "XLALConvertPSSTimeseriesToREAL8Timeseries", XLAL_EINVAL ); 
  for(i = 0; i < ts->data->length; i++)
    ts->data->data[i] = tsPSS->y[i];
  strncpy(ts->name, tsPSS->name, LALNameLength);
  /* XLALFloatToGPS(&(ts->epoch), tsPSS->ini); /* FIXME: needs newer LAL */
  ts->deltaT = tsPSS->dx;
  ts->f0     = 0.0; /* no heterodyning */
  XLALParseUnitString( &(ts->sampleUnits), tsPSS->capt );
  return ts;
}

PSSTimeseries *XLALConvertREAL8TimeseriesToPSSTimeseries(PSSTimeseries *tsPSS, REAL8TimeSeries *ts) {
  UINT4 i;
  char unit[LALNameLength];

  if ( !tsPSS || !ts )
    XLAL_ERROR_NULL( "XLALConvertREAL8TimeseriesToPSSTimeseries", XLAL_EFAULT );
  if (tsPSS->nall != (int)ts->data->length)
    XLAL_ERROR_NULL( "XLALConvertREAL8TimeseriesToPSSTimeseries", XLAL_EINVAL ); 
  if( XLALUnitAsString( unit, LALNameLength, &(ts->sampleUnits) ) ) {
    tsPSS->capt = (char*)malloc(strlen(unit)+1);
    if( !tsPSS->capt )
      XLAL_ERROR_NULL( "XLALConvertREAL8TimeseriesToPSSTimeseries", XLAL_ENOMEM );
    strncpy(tsPSS->capt,unit,strlen(unit)+1);
  } else {
    tsPSS->capt = NULL;
  }
  for(i = 0; i < ts->data->length; i++)
    tsPSS->y[i] = ts->data->data[i];
  tsPSS->n = ts->data->length;
  strncpy(tsPSS->name, ts->name, sizeof(tsPSS->name)); /* sizeof(tsPSS->name) = 20 */
  /* XLALFloatToGPS(&(ts->epoch), tsPSS->ini); /* FIXME: needs newer LAL */
  tsPSS->dx = ts->deltaT;

  return tsPSS;
}

PSSTimeseries *XLALPSSHighpassData(PSSTimeseries *tsout, PSSTimeseries *tsin, PSSHeaderParams* hp, REAL4 f) {
  int ret;
  if ( !tsout || !tsin || !hp )
    XLAL_ERROR_NULL( "XLALCreatePSSTimeseries", XLAL_EFAULT );
  ret = highpass_data_bil(tsout,tsin,hp,f);
  if ( ret != 0 )
    XLAL_ERROR_NULL( "XLALCreatePSSTimeseries", XLAL_EFAULT );
  return tsout;
}

PSSEventParams *XLALIdentifyPSSCleaningEvents(PSSEventParams *events, PSSTimeseries *ts, PSSHeaderParams* hp) {
  int ret;
  if ( !events || !ts )
    XLAL_ERROR_NULL( "XLALIdentifyPSSCleaningEvents", XLAL_EFAULT );
  ret = sn_medsig(ts,events,hp);
  if ( ret != 0 )
    XLAL_ERROR_NULL( "XLALIdentifyPSSCleaningEvents", XLAL_EFUNC );
  ret = even_anst(ts,events);
  if ( ret != 0 )
    XLAL_ERROR_NULL( "XLALIdentifyPSSCleaningEvents", XLAL_EFUNC );
  return events;
}

PSSTimeseries *XLALSubstractPSSCleaningEvents(PSSTimeseries *tsout, PSSTimeseries *tsin,
					      PSSTimeseries *tshp, PSSEventParams *events,
					      PSSHeaderParams* hp) {
  int ret;
  if ( !tsout || !tsin || !tshp || !events )
    XLAL_ERROR_NULL( "XLALSubstractPSSCleaningEvents", XLAL_EFAULT );
  ret = purge_data_subtract(tsout,tsin,tshp,events,hp);
  if ( ret != 0 )
    XLAL_ERROR_NULL( "XLALSubstractPSSCleaningEvents", XLAL_EFUNC );
  return tsout;
}
