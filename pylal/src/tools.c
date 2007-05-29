#include <Python.h>
#include <lal/XLALError.h>
#include <lal/LALStdio.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/CoincInspiralEllipsoid.h>
#include <lal/TimeDelay.h>
#include <lal/DetectorSite.h>
#include <lal/DetResponse.h>
#include <lal/Date.h>

SnglInspiralTable *PySnglInspiral2CSnglInspiral(PyObject *row) {
    /* Convert a Python SnglInspiral (row) to a C SnglInspiralTable.
    Used in function PyCalculateEThincaParameter and
      PyThincaParameterForInjection.
    Calls LAL functions LALSnprintf and LALCalloc. */
    
    SnglInspiralTable *event; /* Return value */
    PyObject *temp; /* Holds each datum to copy */
    
    /* allocate new memory for row */
    event = (SnglInspiralTable *) LALCalloc(1, sizeof(SnglInspiralTable));
    
    /* copy to C SnglInspiral row */
    /* Procedure for each variable:
       - Extract; this increases the Python object's reference count
       - Copy
       - Decrement the Python object's reference count; omission => mem leak
     */
    temp = PyObject_GetAttrString(row, "ifo");
    LALSnprintf( event->ifo, LIGOMETA_IFO_MAX  * sizeof(CHAR),
                 "%s", PyString_AsString(temp));
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row,  "search");
    LALSnprintf( event->search, LIGOMETA_SEARCH_MAX  * sizeof(CHAR),
                 "%s", PyString_AsString(temp));
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "channel");
    LALSnprintf( event->channel, LIGOMETA_CHANNEL_MAX  * sizeof(CHAR),
                 "%s", PyString_AsString(temp));
    Py_XDECREF(temp);
    
    temp = PyObject_GetAttrString(row, "end_time");
    event->end_time.gpsSeconds = PyInt_AsLong(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "end_time_ns");
    event->end_time.gpsNanoSeconds = PyInt_AsLong(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "end_time_gmst");
    event->end_time_gmst = PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "impulse_time");
    event->impulse_time.gpsSeconds = PyInt_AsLong(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "impulse_time_ns");
    event->impulse_time.gpsNanoSeconds = PyInt_AsLong(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "template_duration");
    event->template_duration = PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "event_duration");
    event->event_duration = PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "amplitude");
    event->amplitude = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "eff_distance");
    event->eff_distance = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "coa_phase");
    event->coa_phase = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "mass1");
    event->mass1 = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "mass2");
    event->mass2 = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "mchirp");
    event->mchirp = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "mtotal");
    event->mtotal = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "eta");
    event->eta = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "tau0");
    event->tau0 = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "tau2");
    event->tau2 = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "tau3");
    event->tau3 = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "tau4");
    event->tau4 = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "tau5");
    event->tau5 = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "ttotal");
    event->ttotal = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "psi0");
    event->psi0 = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "psi3");
    event->psi3 = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "alpha");
    event->alpha = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "alpha1");
    event->alpha1 = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "alpha2");
    event->alpha2 = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "alpha3");
    event->alpha3 = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "alpha4");
    event->alpha4 = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "alpha5");
    event->alpha5 = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "alpha6");
    event->alpha6 = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "beta");
    event->beta = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "f_final");
    event->f_final = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "snr");
    event->snr = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "chisq");
    event->chisq = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "chisq_dof");
    event->chisq_dof = PyInt_AsLong(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "sigmasq");
    event->sigmasq = PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "rsqveto_duration");
    event->rsqveto_duration = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    
    temp = PyObject_GetAttrString(row, "Gamma0");
    event->Gamma[0] = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "Gamma1");
    event->Gamma[1] = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "Gamma2");
    event->Gamma[2] = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "Gamma3");
    event->Gamma[3] = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "Gamma4");
    event->Gamma[4] = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "Gamma5");
    event->Gamma[5] = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "Gamma6");
    event->Gamma[6] = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "Gamma7");
    event->Gamma[7] = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "Gamma8");
    event->Gamma[8] = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "Gamma9");
    event->Gamma[9] = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    
    event->event_id = (EventIDColumn *) LALCalloc(1, sizeof(EventIDColumn));
    temp = PyObject_GetAttrString(row, "event_id");
    event->event_id->id = PyLong_AsLongLong(temp);
    Py_XDECREF(temp);
    
    event->next = NULL;
    
    return event;
}

SimInspiralTable *PySimInspiral2CSimInspiral(PyObject *row) {
    /* Convert a Python SimInspiral (row) to a C SimInspiralTable.
    Used in function PyThincaParameterForInjection.
    Calls LAL functions LALSnprintf and LALCalloc. */
    
    SimInspiralTable *event; /* Return value */
    PyObject *temp; /* Holds each datum for copy and refcount decrement */
    
    /* allocate new memory for row */
    event = (SimInspiralTable *) LALCalloc(1, sizeof(SimInspiralTable));
    
    /* copy to C SimInspiral row */
    /* Procedure for each variable:
       - Extract; this increases the Python object's reference count
       - Copy
       - Decrement the Python object's reference count; omission => mem leak
     */
    temp = PyObject_GetAttrString(row, "waveform");
    LALSnprintf( event->waveform, LIGOMETA_IFO_MAX  * sizeof(CHAR),
                 "%s", PyString_AsString(temp));
    Py_XDECREF(temp);
    
    temp = PyObject_GetAttrString(row, "geocent_end_time");
    event->geocent_end_time.gpsSeconds = PyInt_AsLong(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "geocent_end_time_ns");
    event->geocent_end_time.gpsNanoSeconds = PyInt_AsLong(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "h_end_time");
    event->h_end_time.gpsSeconds = PyInt_AsLong(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "h_end_time_ns");
    event->h_end_time.gpsNanoSeconds = PyInt_AsLong(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "l_end_time");
    event->l_end_time.gpsSeconds = PyInt_AsLong(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "l_end_time_ns");
    event->l_end_time.gpsNanoSeconds = PyInt_AsLong(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "g_end_time");
    event->g_end_time.gpsSeconds = PyInt_AsLong(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "g_end_time_ns");
    event->g_end_time.gpsNanoSeconds = PyInt_AsLong(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "t_end_time");
    event->t_end_time.gpsSeconds = PyInt_AsLong(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "t_end_time_ns");
    event->t_end_time.gpsNanoSeconds = PyInt_AsLong(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "v_end_time");
    event->v_end_time.gpsSeconds = PyInt_AsLong(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "v_end_time_ns");
    event->v_end_time.gpsNanoSeconds = PyInt_AsLong(temp);
    Py_XDECREF(temp);
    
    temp = PyObject_GetAttrString(row, "end_time_gmst");
    event->end_time_gmst = PyFloat_AsDouble(temp);
    
    temp = PyObject_GetAttrString(row, "source");
    LALSnprintf( event->source, LIGOMETA_IFO_MAX  * sizeof(CHAR),
                 "%s", PyString_AsString(temp));
    Py_XDECREF(temp);
    
    temp = PyObject_GetAttrString(row, "mass1");
    event->mass1 = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "mass2");
    event->mass2 =  (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "eta");
    event->eta = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "distance");
    event->distance = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "longitude");
    event->longitude = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "latitude");
    event->latitude = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "inclination");
    event->inclination = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "coa_phase");
    event->coa_phase = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "polarization");
    event->polarization = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "psi0");
    event->psi0 = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "psi3");
    event->psi3 = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "alpha");
    event->alpha = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "alpha1");
    event->alpha1 = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "alpha2");
    event->alpha2 = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "alpha3");
    event->alpha3 = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "alpha4");
    event->alpha4 = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "alpha5");
    event->alpha5 = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "alpha6");
    event->alpha6 = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "beta");
    event->beta = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "spin1x");
    event->spin1x = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "spin1y");
    event->spin1y = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "spin1z");
    event->spin1z = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "spin2x");
    event->spin2x = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "spin2y");
    event->spin2y = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "spin2z");
    event->spin2z = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "theta0");
    event->theta0 = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "phi0");
    event->phi0 = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "f_lower");
    event->f_lower = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "f_final");
    event->f_final = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "mchirp");
    event->mchirp = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "eff_dist_h");
    event->eff_dist_h = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "eff_dist_l");
    event->eff_dist_l = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "eff_dist_g");
    event->eff_dist_g = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "eff_dist_t");
    event->eff_dist_t = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "eff_dist_v");
    event->eff_dist_v = (float)PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    
    event->event_id = (EventIDColumn *) LALCalloc(1, sizeof(EventIDColumn));
    temp = PyObject_GetAttrString(row, "simulation_id");
    LALSnprintf( event->event_id->textId, LIGOMETA_IFO_MAX  * sizeof(CHAR),
                 "%s", PyString_AsString(temp));
    Py_XDECREF(temp);
    
    event->next = NULL;
    
    return event;
}

static PyObject *PyCalculateEThincaParameter(PyObject *self, PyObject *args) {
    /* Take two Python SnglInspiral values (rows of SnglInspiralTable) and
    call XLALCalculateEThincaParameter on their contents. */
    
    double result;
    PyObject *py_row1, *py_row2;
    SnglInspiralTable *c_row1, *c_row2;
    InspiralAccuracyList* accuracyParams=NULL;
    
    if (! PyArg_ParseTuple(args, "OO", &py_row1, &py_row2))
        return NULL;
    
    /* Get rows into a format suitable for the LAL call */
    c_row1 = PySnglInspiral2CSnglInspiral(py_row1);
    c_row2 = PySnglInspiral2CSnglInspiral(py_row2);
    
    accuracyParams=(InspiralAccuracyList*)LALMalloc( sizeof( InspiralAccuracyList) );

    XLALPopulateAccuracyParams( accuracyParams );

    /* This is the main call */
    result = (double) XLALCalculateEThincaParameter(c_row1, c_row2, accuracyParams);
    if (XLAL_IS_REAL8_FAIL_NAN((REAL8) result)) {
        /* convert XLAL exception to Python exception */
        XLALClearErrno();
        PyErr_SetString(PyExc_ValueError, "SnglInspiral triggers are not coincident.");
        return NULL;
    }
    
    /* Free temporary memory */
    LALFree(c_row1);
    LALFree(c_row2);
    LALFree(accuracyParams);
    
    return Py_BuildValue("d", result);
}

static PyObject *PyCalculateEThincaParameterExt(PyObject *self, PyObject *args) {
    /* Take two Python SnglInspiral values (rows of SnglInspiralTable) and
    call XLALCalculateEThincaParameter on their contents. */
    
    double result;
    double ra_deg, dec_deg;
    long gps;
    LIGOTimeGPS *gpstime;
    PyObject    *py_row1, *py_row2;
    SnglInspiralTable *c_row1, *c_row2;
    InspiralAccuracyList* accuracyParams;
    
    if (! PyArg_ParseTuple(args, "OOldd", &py_row1, &py_row2, &gps, &ra_deg, &dec_deg ))
        return NULL;
    
    /* check the values */
    if (ra_deg<0 || ra_deg > 360) 
    {
      XLALPrintError("Right ascension value outside [0; 360]. Value given: %f\n", ra_deg);
      return NULL;
    }
    if (dec_deg<-90 || dec_deg>90) 
    {
      XLALPrintError("Declination value outside [-90; 90]. Value given: %f\n", dec_deg);
      return NULL;
    }

    /* Get rows into a format suitable for the LAL call */
    c_row1 = PySnglInspiral2CSnglInspiral(py_row1);
    c_row2 = PySnglInspiral2CSnglInspiral(py_row2);
    
    gpstime=(LIGOTimeGPS*)LALMalloc( sizeof(LIGOTimeGPS) );
    gpstime->gpsSeconds=gps;
    gpstime->gpsNanoSeconds=0;

    accuracyParams=(InspiralAccuracyList*)LALMalloc( sizeof( InspiralAccuracyList) );
    XLALPopulateAccuracyParamsExt( accuracyParams, gpstime, ra_deg, dec_deg );

    /* This is the main call */    
    result = (double) XLALCalculateEThincaParameter(c_row1, c_row2, accuracyParams);
    if (XLAL_IS_REAL8_FAIL_NAN((REAL8) result)) {
        /* convert XLAL exception to Python exception */
        XLALClearErrno();
        PyErr_SetString(PyExc_ValueError, "SnglInspiral triggers are not coincident.");
        return NULL;
    }
    
    /* Free temporary memory */
    LALFree(c_row1);
    LALFree(c_row2);
    LALFree(accuracyParams);
    LALFree(gpstime);
    
    return Py_BuildValue("d", result);
}

static PyObject *PyEThincaParameterForInjection(PyObject *self, PyObject *args) {
    /* Take a Python SimInspiral value and a Python SnglInspiral value
    (rows of SimInspiralTable and SnglInspiralTable, respectively) and
    call XLALEThincaParameterForInjection on their contents. */
    
    double result;
    PyObject *py_row1, *py_row2;
    SimInspiralTable *c_row1;
    SnglInspiralTable *c_row2;
    
    if (! PyArg_ParseTuple(args, "OO", &py_row1, &py_row2))
        return NULL;
    
    /* Get rows into a format suitable for the LAL call */
    c_row1 = PySimInspiral2CSimInspiral(py_row1);
    c_row2 = PySnglInspiral2CSnglInspiral(py_row2);
    
    /* This is the main call */
    result = (double) XLALEThincaParameterForInjection(c_row1, c_row2);
    
    /* Free temporary memory */
    LALFree(c_row1);
    LALFree(c_row2);
    
    return Py_BuildValue("d", result);
}


static PyObject *PyExtTrigTool(PyObject *self, PyObject *args) {
    /* Takes a bunch of informations (like time and sky location)
       and calculates the antenna responses for a given detector,
       angles given in radians! */

  static LALStatus status;
  double gps, alpha, delta, incl, pol;
  double fplus, fcross, fave, qvalue;
  double ci, cc;
  char* detName=NULL;
  LALDetector detector;
  LIGOTimeGPS  timeGPS;
  LALGPSandAcc timeAcc;    
  SkyPosition skysource;
  LALSource* source=NULL;
  LALDetAndSource* detSource=NULL;
  LALDetAMResponse* pResponse=NULL;
  
  /* parsing arguments */
  if (! PyArg_ParseTuple(args, "ddddds", &gps, &alpha, &delta, &incl,  &pol, &detName))
    return NULL;

  source    = (LALSource*)LALMallocShort( sizeof(LALSource));
  detSource = (LALDetAndSource*)LALMallocShort( sizeof(LALDetAndSource));
  pResponse = (LALDetAMResponse*)LALMallocShort( sizeof( LALDetAMResponse ));
  
  /* get detector */
  if ( !strcmp(detName, "H1") ) {
    detector=lalCachedDetectors[LALDetectorIndexLHODIFF];
  } else if ( !strcmp(detName, "H2") ) {
    detector=lalCachedDetectors[LALDetectorIndexLHODIFF]; 
  } else if ( !strcmp(detName, "L1") ) {
    detector=lalCachedDetectors[LALDetectorIndexLLODIFF];
  } else if ( !strcmp(detName, "G1") ) {
    detector=lalCachedDetectors[LALDetectorIndexGEO600DIFF];
  } else if ( !strcmp(detName, "V1") ) {
    detector=lalCachedDetectors[LALDetectorIndexVIRGODIFF];
  }

  /* set the time */
  XLALFloatToGPS( &timeGPS, gps);
  timeAcc.gps      = timeGPS;
  timeAcc.accuracy = LALLEAPSEC_STRICT;

  /* set the sky position */
  skysource.longitude = alpha*LAL_PI_180;
  skysource.latitude  = delta*LAL_PI_180;
  skysource.system    = COORDINATESYSTEM_EQUATORIAL;
 
  /* fill LALSource structure */
  source->equatorialCoords = skysource;
  source->orientation      = pol*LAL_PI_180;

  /* fill LALDetAndSource */
  detSource->pDetector = &detector;
  detSource->pSource   = source;

  /* call the C-function itself */
  LALComputeDetAMResponse( &status, pResponse, detSource, &timeAcc);

  /* calculate the average response and the q-value */
  fplus=pResponse->plus;
  fcross=pResponse->cross;
  fave=sqrt( (fplus*fplus + fcross*fcross)/2.0 );
  ci=cos( incl*LAL_PI_180 );
  cc=ci*ci;
  qvalue=sqrt( fplus*fplus*(1+cc)*(1+cc)/4.0 + fcross*fcross*cc ); /* ref: Duncans PhD, eq. (4.3) on page 57 */

  /* set the return value */
  return Py_BuildValue("dddd", fplus, fcross, fave, qvalue );
}

static PyObject *PyTimeDelay(PyObject *self, PyObject *args) {
    /* Takes a bunch of informations (like time and sky location)
       and calculates the time delay between two detectors */
    
  static LALStatus status;
  double gps, alpha, delta;
  char* detName1=NULL;
  char* detName2=NULL;
  LALDetector detector1, detector2;
  LIGOTimeGPS  timeGPS;
  LALGPSandAcc timeAcc;    
  LALPlaceAndGPS* site1;
  LALPlaceAndGPS* site2;
  SkyPosition skysource;
  TwoDetsTimeAndASource* sourceDets;
  double delay;

  /* parsing arguments */
  if (! PyArg_ParseTuple(args, "dddss", &gps, &alpha, &delta, &detName1, &detName2))
    return NULL;
  
  site1     = (LALPlaceAndGPS*)LALMallocShort( sizeof(LALPlaceAndGPS ));
  site2     = (LALPlaceAndGPS*)LALMallocShort( sizeof(LALPlaceAndGPS ));
  sourceDets= (TwoDetsTimeAndASource*)LALMallocShort( sizeof(TwoDetsTimeAndASource));
  
  /* get detector */
  if ( !strcmp(detName1, "H1") ) {
    detector1=lalCachedDetectors[LALDetectorIndexLHODIFF];
  } else if ( !strcmp(detName1, "H2") ) {
    detector1=lalCachedDetectors[LALDetectorIndexLHODIFF]; 
  } else if ( !strcmp(detName1, "L1") ) {
    detector1=lalCachedDetectors[LALDetectorIndexLLODIFF];
  } else if ( !strcmp(detName1, "G1") ) {
    detector1=lalCachedDetectors[LALDetectorIndexGEO600DIFF];
  } else if ( !strcmp(detName1, "V1") ) {
    detector1=lalCachedDetectors[LALDetectorIndexVIRGODIFF];
  }
  if ( !strcmp(detName2, "H1") ) {
    detector2=lalCachedDetectors[LALDetectorIndexLHODIFF];
  } else if ( !strcmp(detName2, "H2") ) {
    detector2=lalCachedDetectors[LALDetectorIndexLHODIFF]; 
  } else if ( !strcmp(detName2, "L1") ) {
    detector2=lalCachedDetectors[LALDetectorIndexLLODIFF];
  } else if ( !strcmp(detName2, "G1") ) {
    detector2=lalCachedDetectors[LALDetectorIndexGEO600DIFF];
  } else if ( !strcmp(detName2, "V1") ) {
    detector2=lalCachedDetectors[LALDetectorIndexVIRGODIFF];
  }

  /* set the time */
  XLALFloatToGPS( &timeGPS, gps);
  timeAcc.gps      = timeGPS;
  timeAcc.accuracy = LALLEAPSEC_STRICT;  

  /* set the detectors and the time */
  site1->p_detector = &detector1;
  site1->p_gps      = &timeGPS;
  site2->p_detector = &detector2;
  site2->p_gps      = &timeGPS;

  /* set the sky position */
  skysource.longitude = alpha*LAL_PI_180;
  skysource.latitude  = delta*LAL_PI_180;
  skysource.system    = COORDINATESYSTEM_EQUATORIAL;
  
  /* put all together */
  sourceDets->p_det_and_time1 = site1;
  sourceDets->p_det_and_time2 = site2;
  sourceDets->p_source        = &skysource;  
 
  /* now finally after all this chunk: do the calculation... */
  LALTimeDelay( &status, &delay, sourceDets );

  return Py_BuildValue("d", delay);
}

static struct PyMethodDef tools_methods[] = {
    {"XLALCalculateEThincaParameter", PyCalculateEThincaParameter,
     METH_VARARGS,
     "XLALCalculateEThincaParameter(SnglInspiral1, SnglInspiral2)\n"
     "\n"
     "Takes two SnglInspiral objects (rows of a SnglInspiralTable) and\n"
     "calculates the overlap factor between them."},
    {"XLALCalculateEThincaParameterExt", PyCalculateEThincaParameterExt,
     METH_VARARGS,
     "XLALCalculateEThincaParameterExt(SnglInspiral1, SnglInspiral2, gps, ra_deg, dec_deg)\n"
     "\n"
     "Takes two SnglInspiral objects (rows of a SnglInspiralTable) and\n"
     "calculates the overlap factor between them, using the known time delay \n"
     "between the IFO's at a given time for a given sky location (given in degrees)."},
     {"XLALEThincaParameterForInjection", PyEThincaParameterForInjection,
      METH_VARARGS, \
      "XLALEThincaParameterForInjection(SimInspiral, SnglInspiral)\n"
      "\n"
      "Takes a SimInspiral and a SnglInspiral object (rows of\n"
      "SimInspiralTable and SnglInspiralTable, respectively) and\n"
      "calculates the ethinca parameter required to put the SimInspiral\n"
      "parameters inside the SnglInspiral template's ellipse."},
    {"XLALExtTrigTool", PyExtTrigTool,
      METH_VARARGS, \
     "XLALExtTrigTool(gps, alpha, delta, iota, psi, det)\n"
     "\n"
     "Calculates the antenna factors at a given gps time (as double)\n"
     "and for a given sky location (alpha, delta) [degree]\n"
     "for a source with given inclination iota and polarization psi\n"
     "(all given in degrees) for a given detector (name, like 'H1').\n" 
     "The returned values are (f-plus, f-cross, f-average, q-value). \n"
     "Example: XLALExtTrigTool(854378604.780, 11.089, 42.308, 0, 0, 'H1' )"}, 
    {"XLALTimeDelay", PyTimeDelay,
     METH_VARARGS,					\
     "XLALTimeDelay(gps, alpha, delta, det1, det2)\n"
     "\n"
     "Calculates the time delay for two given detectors (names, like 'H1')\n"
     "and for a given sky location (alpha, delta) [degree] \n"
     "and a given gps time (as double). The returned value is the\n " 
     "time delay in seconds.\n"
     "Example: XLALTimeDelay(854378604.780, 11.089, 42.308, 'H1','L1' )"    },   
    {NULL, NULL, 0}
};

void inittools (void) {
    (void) Py_InitModule("pylal.tools", tools_methods);
}
