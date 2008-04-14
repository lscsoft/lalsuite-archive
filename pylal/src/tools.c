#include <Python.h>
#include <lal/XLALError.h>
#include <lal/LALStdio.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/CoincInspiralEllipsoid.h>
#include <lal/CoincInspiralEllipsoid.h>
#include <lal/TimeDelay.h>
#include <lal/DetectorSite.h>
#include <lal/DetResponse.h>
#include <lal/Date.h>
#include <lal/EllipsoidOverlapTools.h>

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
    if ( temp == Py_None ) { Py_DECREF(temp); temp = PyString_FromString(""); }
    LALSnprintf( event->ifo, LIGOMETA_IFO_MAX  * sizeof(CHAR),
                 "%s", PyString_AsString(temp));
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row,  "search");
    if ( temp == Py_None ) { Py_DECREF(temp); temp = PyString_FromString(""); }
    LALSnprintf( event->search, LIGOMETA_SEARCH_MAX  * sizeof(CHAR),
                 "%s", PyString_AsString(temp));
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "channel");
    if ( temp == Py_None ) { Py_DECREF(temp); temp = PyString_FromString(""); }
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
    PyObject *temp2;
    
    /* allocate new memory for row */
    event = (SimInspiralTable *) LALCalloc(1, sizeof(SimInspiralTable));
    
    /* copy to C SimInspiral row */
    /* Procedure for each variable:
       - Extract; this increases the Python object's reference count
       - Copy
       - Decrement the Python object's reference count; omission => mem leak
     */
    temp = PyObject_GetAttrString(row, "waveform");
    if ( temp == Py_None ) { Py_DECREF(temp); temp = PyString_FromString(""); }
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
    if ( temp == Py_None ) { Py_DECREF(temp); temp = PyString_FromString(""); }
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
    temp2 = PyObject_Str(temp);
    LALSnprintf( event->event_id->textId, LIGOMETA_IFO_MAX  * sizeof(CHAR),
                 "%s", PyString_AsString(temp2));
    Py_XDECREF(temp);
    Py_XDECREF(temp2);
    
    event->next = NULL;
    
    return event;
}

static PyObject *PyCalculateEThincaParameter(PyObject *self, PyObject *args) {
    /* Take two Python SnglInspiral values (rows of SnglInspiralTable) and
    call XLALCalculateEThincaParameter on their contents. */
    
    double result;
    PyObject *py_row1, *py_row2;
    SnglInspiralTable *c_row1, *c_row2;
    InspiralAccuracyList accuracyParams;
    
    if (! PyArg_ParseTuple(args, "OO", &py_row1, &py_row2))
        return NULL;
    
    /* Get rows into a format suitable for the LAL call */
    c_row1 = PySnglInspiral2CSnglInspiral(py_row1);
    c_row2 = PySnglInspiral2CSnglInspiral(py_row2);

    memset(&accuracyParams, 0, sizeof(accuracyParams));
    XLALPopulateAccuracyParams( &accuracyParams );

    /* This is the main call */
    result = (double) XLALCalculateEThincaParameter(c_row1, c_row2, &accuracyParams);
    
    /* Free temporary memory */
    LALFree(c_row1->event_id);
    LALFree(c_row1);
    LALFree(c_row2->event_id);
    LALFree(c_row2);

    if (XLAL_IS_REAL8_FAIL_NAN((REAL8) result)) {
        /* convert XLAL exception to Python exception */
        XLALClearErrno();
        PyErr_SetString(PyExc_ValueError, "SnglInspiral triggers are not coincident.");
        return NULL;
    }
    
    return PyFloat_FromDouble(result);
}

static PyObject *PyCalculateEThincaParameterExt(PyObject *self, PyObject *args) {
    /* Take two Python SnglInspiral values (rows of SnglInspiralTable) and
    call XLALCalculateEThincaParameter on their contents. */
    
    double result;
    double ra_deg, dec_deg;
    long gps;
    LIGOTimeGPS gpstime;
    PyObject    *py_row1, *py_row2;
    SnglInspiralTable *c_row1, *c_row2;
    InspiralAccuracyList accuracyParams;
    
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

    XLALGPSSet(&gpstime, gps, 0);
    memset(&accuracyParams, 0, sizeof(accuracyParams));
    XLALPopulateAccuracyParamsExt( &accuracyParams, &gpstime, ra_deg, dec_deg );

    /* This is the main call */    
    result = (double) XLALCalculateEThincaParameter(c_row1, c_row2, &accuracyParams);
    
    /* Free temporary memory */
    LALFree(c_row1->event_id);
    LALFree(c_row1);
    LALFree(c_row2->event_id);
    LALFree(c_row2);

    if (XLAL_IS_REAL8_FAIL_NAN((REAL8) result)) {
        /* convert XLAL exception to Python exception */
        XLALClearErrno();
        PyErr_SetString(PyExc_ValueError, "SnglInspiral triggers are not coincident.");
        return NULL;
    }
    
    return PyFloat_FromDouble(result);
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
    LALFree(c_row1->event_id);
    LALFree(c_row1);
    LALFree(c_row2->event_id);
    LALFree(c_row2);
    
    return PyFloat_FromDouble(result);
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
      "calculates the ethinca parameter required to put the SimInspiral\n"},   
    {NULL, NULL, 0}
};

void inittools (void) {
    (void) Py_InitModule("pylal.tools", tools_methods);
}
