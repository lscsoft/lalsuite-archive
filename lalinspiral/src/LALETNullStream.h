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

#ifndef _LALETNULLSTREAM_H    /* Protect against double-inclusion */
#define _LALETNULLSTREAM_H

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConfig.h>
#include <lal/TimeFreqFFT.h>
#include <lal/LALConstants.h>
#include <lal/DetectorSite.h>
#include <lal/LALDetectors.h>
#include <lal/LALStatusMacros.h>
#include <lal/LALNoiseModels.h>
#include <lal/LALComplex.h>
#include <lal/Units.h>
#include <lal/Date.h>
#include <lal/AVFactories.h>
#include <lal/VectorOps.h>
#include <lal/TimeSeries.h>
#include <lal/FrequencySeries.h>
#include <lal/LALCalibration.h>
#include <lal/FrameCache.h>
#include <lal/FrameStream.h>
#include <lal/ResampleTimeSeries.h>

#ifdef  __cplusplus
extern "C" {
#endif
	
//NRCSID (LALETNULLSTREAMH,"$Id$");
    
    REAL8TimeSeries * LALETNullStream (LIGOTimeGPS *GPSStart, REAL8 duration );

    REAL8TimeSeries *ReadTimeSerieFromCache(const CHAR *cachefile, const CHAR *channel, LIGOTimeGPS *start, REAL8 duration);
    
    void PopulateNullStream(REAL8TimeSeries *NullStream, REAL8TimeSeries *RawData);
    
   REAL8FrequencySeries * ComputeSingleDetectorInvPSDfromNullStream(LIGOTimeGPS *GPSStart, REAL8 duration, UINT4 SampleRate, UINT4 nSegs); 

#ifdef  __cplusplus
}
#endif      /* Close C++ protection */

#endif      /* Close double-include protection */
