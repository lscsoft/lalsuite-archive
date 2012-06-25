# Copyright (C) 2012 Duncan M. Macleod
# 
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#   
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
# 
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

"""
This module provides a bunch of user-friendly wrappers to the SWIG-bound REAL8TimeSeries and REAL8FrequencySeries objects and their associated functions.
"""

# =============================================================================
# Preamble
# =============================================================================

from __future__ import division
import os,re,tempfile,warnings

# import swig LALSuite packages
import swiglal as lal
import swiglalframe as lalframe

# import gluecache
from glue import lal as _gluecache

# type code dict
_typestr = {\
            lal.LAL_I2_TYPE_CODE: 'INT2',\
            lal.LAL_I4_TYPE_CODE: 'INT4',\
            lal.LAL_I8_TYPE_CODE: 'INT8',\
            lal.LAL_U2_TYPE_CODE: 'UINT2',\
            lal.LAL_U4_TYPE_CODE: 'UINT4',\
            lal.LAL_U8_TYPE_CODE: 'UINT8',\
            lal.LAL_S_TYPE_CODE:  'REAL4',\
            lal.LAL_D_TYPE_CODE:  'REAL8',\
            lal.LAL_C_TYPE_CODE:  'COMPLEX8',\
            lal.LAL_Z_TYPE_CODE:  'COMPLEX16',\
            }

# =============================================================================
# Data conversion
# =============================================================================

def fromarray(array, name="", epoch=lal.LIGOTimeGPS(), deltaT=1,\
              f0=0, sampleUnits=lal.lalDimensionlessUnit,\
              frequencyseries=False):
    """
    Convert numpy.array to REAL8TimeSeries. Use frequencyseries to return
    REAL8FrequencySeries.

    Arguments:

        array : numpy.array
            input data array

    Keyword arguments:

        name : str
            name of data source
        epoch : lal.LIGOTimeGPS
            GPS start time for array
        deltaT : float
            sampling time for data (or frequency spacing for FrequencySeries)
        f0 : float
            lower frequency limit for data
        sampleUnits : lal.LALUnit
            amplitude unit for array
    """
  
    if frequencyseries:
      series = lal.XLALCreateREAL8FrequencySeries(name, epoch, f0, deltaT,\
                                                  sampleUnits, array.size)
    else:
      series = lal.XLALCreateREAL8TimeSeries(name, epoch, f0, deltaT,\
                                             sampleUnits, array.size)
    series.data.data = array.astype(float)
    return series

# =============================================================================
# Frame reading helpers
# =============================================================================

def fromNDS(chname, start, duration, server="nds.ligo.caltech.edu",\
            port=31200):
    """
    Read data from NDS into a REAL?TimeSeries
    """

    import nds

    connection = nds.daq(server, port)
    data       = connection.fetch(start, start+duration, chname)[0]
    epoch      = lal.LIGOTimeGPS(start)
    deltaT     = duration/data.size
    return fromarray(data, name=chname, epoch=epoch, deltaT=deltaT)

def fromframefile(framefile, chname, start=-1, duration=1,\
                  datatype=-1, verbose=False):
    """
    Read data from the GWF format framefile into a REAL?TimeSeries.
    Restrict data with the gpsstart and duration parameters.
 
    Arguments:

        framefile : str
            path to input frame file
        chname : str
            name of channel to read

    Keyword arguments:

        start : lal.LIGOTimeGPS
            GPS start time of requested data
        duration : float
            length of requested data series in seconds
        datatype : int
            LAL enum for requested datatype, -1 == read from frame metadata
        verbose : [ True | False ]
            verbose output
    """

    # open frame
    framefile = os.path.abspath(framefile)
    stream = lalframe.XLALFrOpen('', framefile)

    # return
    return fromFrStream(stream, chname, start=start, duration=duration,\
                        datatype=datatype, verbose=verbose)

def fromlalcache(cache, chname, start=-1, duration=1, datatype=-1,\
                 verbose=False):
    """
    Read data from a cache object into a REAL?TimeSeries.
    Restrict data with the gpsstart and duration parameters.
 
    Arguments:

        cache : [ filepath | file | glue.lal.Cache | lalframe.FrCache ]
            object to be read to lalframe.FrCache
        chname : str
            name of channel to read

    Keyword arguments:

        start : lal.LIGOTimeGPS
            GPS start time of requested data
        duration : float
            length of requested data series in seconds
        datatype : int
            LAL enum for requested datatype, -1 == read from frame metadata
        verbose : [ True | False ]
            verbose output
    """

    # read cache from file object into gluecache for the time being
    if hasattr(cache, 'readline'):
        cache = gluecache.Cache.fromfile(cache)

    # read cache from glue.lal.Cache, or lalframe.FrCache, or file
    if cache and isinstance(cache, _gluecache.Cache):
        cache = gluecache_to_FrCache(cache)
    elif isinstance(cache, str):
        cache = lalframe.XLALFrImportCache(cache)

    # open cache to stream
    stream = lalframe.XLALFrCacheOpen(cache)

    # return
    return fromFrStream(stream, chname, start=start, duration=duration,\
                        datatype=datatype, verbose=verbose)

def fromFrStream(stream, chname, start=-1, duration=1, datatype=-1,\
                 verbose=False):
    """
    Read data from the lalframe.FrStream object stream into a REAL?TimeSeries.
    Restrict data with the gpsstart and duration parameters.
 
    Arguments:

        stream : lalframe.FrStream
            frame stream to read
        chname : str
            name of channel to read

    Keyword arguments:

        start : lal.LIGOTimeGPS
            GPS start time of requested data
        duration : float
            length of requested data series in seconds
        datatype : int
            LAL enum for requested datatype, -1 == read from frame metadata
        verbose : [ True | False ]
            verbose output
    """

    # set mode
    if verbose:  mode = lalframe.LAL_FR_VERBOSE_MODE
    else:        mode = lalframe.LAL_FR_DEFAULT_MODE
    lalframe.XLALFrSetMode(stream, mode)

    # set time
    if int(start) == -1:  start = stream.epoch
    start = lal.LIGOTimeGPS(start)
    duration = float(duration)
    
    # get series type
    if datatype == -1:
        datatype = lalframe.XLALFrGetTimeSeriesType(chname, stream)
      
    # read to series
    func = getattr(lalframe, 'XLALFrRead%sTimeSeries' % _typestr[datatype])
    series = func(stream, chname, start, duration, 0)

    # return
    return series

def resample(series, rate, inplace=True):
    """
    Resample a REAL?TimeSeries to the new rate (Hertz). By default resampling
    is done in-place and a copy of the series is returned, but if requested
    the original series is duplicated before resampling and so is unaffected.

    Arguments:

        series : swiglal.REAL?TimeSeries
            input timeseries to resample
        rate : float
            sampling rate (Hertz) to resample to

    Keyword arguments:

        inplace : [ True | False ]
            perform resampling inplace on input series, default: True.
            If False, input series is duplicated and so is left in the original
            state.
    """

    datatype = _series_type_code(type(series))
    TYPESTR  = _typestr[datatype]

    func = getattr(lal, 'XLALResample%sTimeSeries' % TYPESTR)
    if inplace:
      func(series, 1/rate)
      return series
    else:
      series2 = duplicate(series)
      func(series2, 1/rate)
      return series2

def highpass(series, frequency, amplitude=0.9, filtorder=8, inplace=True):
    """
    Apply Butterworth high pass filter to REAL?TimeSeries.

    Arguments:

        series : swiglal.REAL?TimeSeries
            input series to highpass
        frequency : float
            frequency below which to attenuate

    Keyword arguments:

        amplitude : float
            desired amplitude response of filter
        filtorder : int
            desired order of filter
        inplace : [ True | False ]
            perform resampling inplace on input series, default: True.
            If False, input series is duplicated and so is left in the original
            state.
    """

    datatype = _series_type_code(type(series))
    TYPESTR  = _typestr[datatype]

    func = getattr(lal, 'XLALHighPass%sTimeSeries' % TYPESTR)
    if inplace:
      func(series, frequency, amplitude, filtorder)
      return series
    else:
      series2 = duplicate(series)
      func(series2, frequency, amplitude, filtorder)
      return series2

def lowpass(series, frequency, amplitude=0.9, filtorder=8, inplace=True):
    """
    Apply Butterworth low pass filter to REAL?TimeSeries.

    Arguments:

        series : swiglal.REAL?TimeSeries
            input series to lowpass
        frequency : float
            frequency above which to attenuate

    Keyword arguments:

        amplitude : float
            desired amplitude response of filter
        filtorder : int
            desired order of filter
        inplace : [ True | False ]
            perform resampling inplace on input series, default: True.
            If False, input series is duplicated and so is left in the original
            state.
    """

    datatype = _series_type_code(type(series))
    TYPESTR  = _typestr[datatype]

    func = getattr(lal, 'XLALLowPass%sTimeSeries' % TYPESTR)
    if inplace:
      func(series, frequency, amplitude, filtorder)
      return series
    else:
      series2 = duplicate(series) 
      func(series2, frequency, amplitude, filtorder)
      return series2

def bandpass(series, fmin, fmax, amplitude=0.9, filtorder=8, inplace=True):
    """
    Apply Butterworth band pass filter to REAL?TimeSeries.
    
    Arguments:

        series : swiglal.REAL?TimeSeries
            input series to lowpass
        fmin : float
            frequency below which to attenuate
        fmax : float
            frequency above which to attenuate

    Keyword arguments:

        amplitude : float
          desired amplitude response of filter
        filtorder : int
           desired order of filter
    """

    series = highpass(series, fmin, amplitude, filtorder, inplace) 
    series = lowpass(series, fmax, amplitude, filtorder, inplace) 
    return series

def duplicate(series):
    """
    Duplicate a REAL?TimeSeries.

    Arguments:

        series : swiglal.REAL?TimeSeries
            input series to duplicate
    """
    datatype = _series_type_code(type(series))
    TYPESTR  = _typestr[datatype]
    func     = getattr(lal, 'XLALCreate%sTimeSeries' % TYPESTR)
    out      = func(series.name, series.epoch, series.f0, series.deltaT,\
                    series.sampleUnits, series.data.length)
    out.data.data = series.data.data
    return out


# =============================================================================
# Average spectrum
# =============================================================================

def compute_average_spectrum(series, seglen, stride, window=None, plan=None,\
                             average='medianmean', unit=lal.lalStrainUnit):
    """
    Calculate the average (power) spectrum of the given REAL?TimeSeries

    Arguments:

        seglen : int
            length of FFT (in samples)
        stride : int
            gap between FFTs (in samples)

    Keyword arguments:

        window : lal.REAL8Window
            swiglal window
        plan : lal.REAL8FFTPlan
            plan for FFT
        spectrum : [ 'median' | 'medianmean' | 'welch' ]
            averaging method for spectrum, default: 'medianmean'
        unit : lal.lalUnit
            LAL unit for data
    """

    datatype = _series_type_code(type(series))
    TYPESTR  = _typestr[datatype]

    seglen = int(seglen)
    stride = int(stride)

    # check data
    reclen  = series.data.length
    numseg  = 1 + int((reclen - seglen)/stride)
    worklen = (numseg - 1)*stride + seglen
    if worklen != reclen:
        warnings.warn("Data is not long enough to be covered completely, the "+\
                      "trailing %d samples will not be used for spectrum "\
                      "calculation" % (reclen-worklen))
        series = duplicate(series)
        func = getattr(lal, 'XLALResize%sTimeSeries' % TYPESTR)
        func(series, 0, worklen)
        reclen  = series.data.length
        numseg  = 1 + int((reclen - seglen)/stride)
    if numseg % 2 == 1:
        warnings.warn("Data is not long enough to be covered completely, the "+\
                      "trailing %d samples will not be used for spectrum "\
                      "calculation" % (seglen-stride))
        worklen = reclen - (seglen-stride)
        series = duplicate(series)
        func = getattr(lal, 'XLALResize%sTimeSeries' % TYPESTR)
        func(series, 0, worklen)

    # generate window
    if not window:
        func   = getattr(lal, "XLALCreateKaiser%sWindow" % TYPESTR)
        window = func(seglen, 24)
    # generate FFT plan
    if not plan:
        func = getattr(lal, "XLALCreateForward%sFFTPlan" % TYPESTR)
        plan = func(seglen, 1)

    # generate output spectrum
    f0       = (1/series.deltaT) * (1/seglen)
    deltaF   = (1/seglen)
    func     = getattr(lal, "XLALCreate%sFrequencySeries" % TYPESTR)
    spectrum = func(series.name, series.epoch, f0, deltaF, lal.lalStrainUnit,\
                    seglen//2+1)

    # calculate medianmean spectrum
    if re.match('median?mean\Z', average, re.I):
        func = getattr(lal, "XLAL%sAverageSpectrumMedianMean" % TYPESTR)
    elif re.match('median\Z', average, re.I):
        func = getattr(lal, "XLAL%sAverageSpectrumMedian" % TYPESTR)
    elif re.match('welch\Z', average, re.I):
        func = getattr(lal, "XLAL%sAverageSpectrumWelch" % TYPESTR)
    else:
        raise NotImplementedError("Sorry, only 'median' and 'medianmean' "+\
                                  " and 'welch' average methods are available.")
    func(spectrum, series, seglen, stride, window, plan)

    # dereference
    del window
    del plan

    # return
    return spectrum

# =============================================================================
# Get data type from series
# =============================================================================

def _series_type_code(t):
    """
    Return LAL type code for this type string.
    """
    if re.search('REAL4', str(t)):   return lal.LAL_S_TYPE_CODE
    elif re.search('REAL8', str(t)): return lal.LAL_D_TYPE_CODE
    else: raise TypeError("Cannot find type code for %s" % str(t))

# =============================================================================
# Helper functions
# =============================================================================

def gluecache_to_FrCache(cache):
    """
    Convert a glue.lal.Cache object to a lalframe.FrCache object.
    Writes cache to temporary file and reads to FrCache.

    Arguments:

        cache : glue.lal.Cache
            LAL cache object from GLUE to convert
    """
    t = tempfile.NamedTemporaryFile(delete=False)
    cache.tofile(t)
    frcache = lalframe.XLALFrImportCache(t.name)
    t.close()
    os.remove(t.name)
    return frcache

