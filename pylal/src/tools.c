#include <Python.h>
#include <lal/LALstdio.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/CoincInspiralEllipsoid.h>

SnglInspiralTable *PySnglInspiral2CSnglInspiral(PyObject *row) {
    /* Convert a Python SnglInspiral (row) to a C SnglInspiralTable.
    Used in function PyCalculateEThincaParameter.
    Calls LAL functions LALSnprintf and LALCalloc. */
    
    SnglInspiralTable *event; /* Return value */
    PyObject *temp; /* Holds each datum for copy and refcount decrement */

    /* allocate new memory for row */
    event = (SnglInspiralTable *) LALCalloc( 1, sizeof(SnglInspiralTable) );
    
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
    event->amplitude = PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "eff_distance");
    event->eff_distance = PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "coa_phase");
    event->coa_phase = PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "mass1");
    event->mass1  = PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "mass2");
    event->mass2  =  PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "mchirp");
    event->mchirp =  PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "mtotal");
    event->mtotal =  PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "eta");
    event->eta    = PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "tau0");
    event->tau0   = PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "tau2");
    event->tau2   = PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "tau3");
    event->tau3   = PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "tau4");
    event->tau4   = PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "tau5");
    event->tau5   = PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "ttotal");
    event->ttotal = PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "psi0");
    event->psi0   = PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "psi3");
    event->psi3   = PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "alpha");
    event->alpha  = PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "alpha1");
    event->alpha1 = PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "alpha2");
    event->alpha2 = PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "alpha3");
    event->alpha3 = PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "alpha4");
    event->alpha4 = PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "alpha5");
    event->alpha5 = PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "alpha6");
    event->alpha6 = PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "beta");
    event->beta   = PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "f_final");
    event->f_final= PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "snr");
    event->snr    = PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "chisq");
    event->chisq  = PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "chisq_dof");
    event->chisq_dof = PyInt_AsLong(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "sigmasq");
    event->sigmasq = PyFloat_AsDouble(temp);
    Py_XDECREF(temp);
    temp = PyObject_GetAttrString(row, "rsqveto_duration");
    event->rsqveto_duration = PyFloat_AsDouble(temp);
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
    
    event->event_id = (EventIDColumn *) LALCalloc( 1, sizeof(EventIDColumn) );
    temp = PyObject_GetAttrString(row, "event_id");
    event->event_id->id = PyLong_AsLongLong(temp);
    Py_XDECREF(temp);
    
    event->next = NULL; /* just in case */
    
    return event;
}

static PyObject *PyCalculateEThincaParameter(PyObject *self, PyObject *args) {
    /* Take two Python SnglInspiral values (rows of SnglInspiralTable) and
    call XLALCalculateEThincaParameter on their contents. */
    
    double result;
    PyObject *py_row1, *py_row2;
    SnglInspiralTable *c_row1, *c_row2, *temp;
    
    if (! PyArg_ParseTuple(args, "OO", &py_row1, &py_row2))
        return NULL;
    
    /* Get rows into a format suitable for the LAL call */
    c_row1 = PySnglInspiral2CSnglInspiral(py_row1);
    c_row2 = PySnglInspiral2CSnglInspiral(py_row2);
	
    /* This is the main call */
    result = XLALCalculateEThincaParameter(c_row1, c_row2);
	
    /* Free rows (being unduly paranoid, since these are single rows) */
    while(c_row1) {
      temp = c_row1;
      c_row1 = c_row1->next;
      LALFree(temp);
    }
    while(c_row2) {
      temp = c_row2;
      c_row2 = c_row2->next;
      LALFree(temp);
    }
    
    return Py_BuildValue("d", result);
}

static struct PyMethodDef tools_methods[] = {
    {"XLALCalculateEThincaParameter", PyCalculateEThincaParameter, 1},
    {NULL, NULL, 0}
};

void inittools (void) {
    (void) Py_InitModule("pylal.tools", tools_methods);
}