#include <Python.h>
#include <structmember.h>

#include <numpy/arrayobject.h>

#include <lal/LALInference.h>
#include <lal/Sequence.h>
#include <lal/LALBarycenter.h>
#include <lal/LALDetectors.h>

#include <complex16frequencyseries.h>
#include <complex16timeseries.h>
#include <real8fftplan.h>
#include <real8frequencyseries.h>
#include <real8timeseries.h>
#include <ligotimegps.h>
#include <real8window.h>
#include <tools.h>

#include "LALVariables.h"
#include "LALIFOData.h"
#include "BarycenterInput.h"
#include "EphemerisData.h"

#define MODULE_NAME "pylal._lalinference"

const char LIDocString[] =
"This module provides data types and function wrappers for"
"LALInference.";

int setLALIFODataFromData(LALIFOData* internal,PyObject* old,PyObject* new){
    /* require li_LALIFOData*/
    if (new == NULL) {
        PyErr_SetString(PyExc_TypeError, "Cannot delete the REAL8Window attribute");
        return -1;
    }
    if(!PyObject_TypeCheck(new, &li_LALIFOData_Type)){
        PyErr_SetObject(PyExc_TypeError, new);
        return -1;
    }

    if(old) Py_DECREF(old);

    Py_INCREF(new);
    old=new;

    li_LALIFOData* newdata=(li_LALIFOData*)new;
    internal=newdata->data;
    return 0;
}

PyObject* getLALIFODataFromData(LALIFOData* internal,PyObject* owner){
    if(!owner){
        Py_INCREF(Py_None);
        return Py_None;
    }else{
        
        return li_LALIFOData_new(internal,owner);
    }
}

int setBarycenterInputFromData(BarycenterInput* internal,PyObject* old,PyObject* new){
    /* require li_BarycenterInput*/
    if (new == NULL) {
        PyErr_SetString(PyExc_TypeError, "Cannot delete the BarycenterInput attribute");
        return -1;
    }
    if(!PyObject_TypeCheck(new, &li_BarycenterInput_Type)){
        PyErr_SetObject(PyExc_TypeError, new);
        return -1;
    }

    if(old) Py_DECREF(old);

    Py_INCREF(new);
    old=new;

    li_BarycenterInput* newbary=(li_BarycenterInput*)new;
    internal=newbary->data;
    return 0;
}

PyObject* getBarycenterInputFromData(BarycenterInput* internal,PyObject* owner){
    if(!owner){
        Py_INCREF(Py_None);
        return Py_None;
    }else{
        
        return li_BarycenterInput_new(internal,owner);
    }
}

int setREAL8WindowFromData(REAL8Window* internal,PyObject* old,PyObject* new){
    /* require pylal_REAL8Window*/
    if (new == NULL) {
        PyErr_SetString(PyExc_TypeError, "Cannot delete the REAL8Window attribute");
        return -1;
    }
    if(!PyObject_TypeCheck(new, &pylal_REAL8Window_Type)){
        PyErr_SetObject(PyExc_TypeError, new);
        return -1;
    }

    if(old) Py_DECREF(old);

    Py_INCREF(new);
    old=new;

    pylal_REAL8Window* newwin=(pylal_REAL8Window*)new;
    internal=newwin->window;
    return 0;
}

PyObject* getREAL8WindowFromData(REAL8Window* internal,PyObject* owner){
    if(!owner){
        Py_INCREF(Py_None);
        return Py_None;
    }else{
        
        return pylal_REAL8Window_new(internal,owner);
    }
}

int setREAL8FFTPlanFromData(REAL8FFTPlan* internal,PyObject* old,PyObject* new){
    /* require pylal_REAL8FFTPlan*/
    if (new == NULL) {
        PyErr_SetString(PyExc_TypeError, "Cannot delete the REAL8FFTPlan attribute");
        return -1;
    }
    if(!PyObject_TypeCheck(new, &pylal_REAL8FFTPlan_Type)){
        PyErr_SetObject(PyExc_TypeError, new);
        return -1;
    }

    if(old) Py_DECREF(old);

    Py_INCREF(new);
    old=new;

    pylal_REAL8FFTPlan* newvec=(pylal_REAL8FFTPlan*)new;
    internal=newvec->plan;
    return 0;
}

PyObject* getREAL8FFTPlanFromData(REAL8FFTPlan* internal,PyObject* owner){
    if(!owner){
        Py_INCREF(Py_None);
        return Py_None;
    }else{
        
        return pylal_REAL8FFTPlan_new(internal,owner);
    }
}

int setCOMPLEX16TimeSeriesFromData(COMPLEX16TimeSeries* internal,PyObject* old,PyObject* new){
    /* require pylal_COMPLEX16TimeSeries*/
    if (new == NULL) {
        PyErr_SetString(PyExc_TypeError, "Cannot delete the COMPLEX16TimeSeries attribute");
        return -1;
    }
    if(!PyObject_TypeCheck(new, &pylal_COMPLEX16TimeSeries_Type)){
        PyErr_SetObject(PyExc_TypeError, new);
        return -1;
    }

    if(old) Py_DECREF(old);

    Py_INCREF(new);
    old=new;

    pylal_COMPLEX16TimeSeries* newvec=(pylal_COMPLEX16TimeSeries*)new;
    internal=newvec->series;
    return 0;
}
PyObject* getCOMPLEX16TimeSeriesFromData(COMPLEX16TimeSeries* internal,PyObject* owner){
    
    if(!owner){
        Py_INCREF(Py_None);
        return Py_None;
    }else{
        
        return pylal_COMPLEX16TimeSeries_new(internal,owner);
    }
}

int setCOMPLEX16FrequencySeriesFromData(COMPLEX16FrequencySeries* internal,PyObject* old,PyObject* new){
    /* require pylal_COMPLEX16FrequencySeries*/
    if (new == NULL) {
        PyErr_SetString(PyExc_TypeError, "Cannot delete the COMPLEX16FrequencySeries attribute");
        return -1;
    }
    if(!PyObject_TypeCheck(new, &pylal_COMPLEX16FrequencySeries_Type)){
        PyErr_SetObject(PyExc_TypeError, new);
        return -1;
    }

    if(old) Py_DECREF(old);

    Py_INCREF(new);
    old=new;

    pylal_COMPLEX16FrequencySeries* newvec=(pylal_COMPLEX16FrequencySeries*)new;
    internal=newvec->series;
    return 0;
}

PyObject* getCOMPLEX16FrequencySeriesFromData(COMPLEX16FrequencySeries* internal,PyObject* owner){
    if(!owner){
        Py_INCREF(Py_None);
        return Py_None;
    }else{
        
        return pylal_COMPLEX16FrequencySeries_new(internal,owner);
    }
}

int setREAL8FrequencySeriesFromData(REAL8FrequencySeries* internal,PyObject* old,PyObject* new){
    /* require pylal_REAL8FrequencySeries*/
    if (new == NULL) {
        PyErr_SetString(PyExc_TypeError, "Cannot delete the REAL8FrequencySeries attribute");
        return -1;
    }
    if(!PyObject_TypeCheck(new, &pylal_REAL8FrequencySeries_Type)){
        PyErr_SetObject(PyExc_TypeError, new);
        return -1;
    }

    if(old) Py_DECREF(old);

    Py_INCREF(new);
    old=new;

    pylal_REAL8FrequencySeries* newvec=(pylal_REAL8FrequencySeries*)new;
    internal=newvec->series;
    return 0;
}

PyObject* getREAL8FrequencySeriesFromData(REAL8FrequencySeries* internal,PyObject* owner){
    if(!owner){
        Py_INCREF(Py_None);
        return Py_None;
    }else{
        
        return pylal_REAL8FrequencySeries_new(internal,owner);
    }
}

int setREAL8TimeSeriesFromData(REAL8TimeSeries* internal,PyObject* old,PyObject* new){
    /* require pylal_REAL8TimeSeries */
    if (new == NULL) {
        PyErr_SetString(PyExc_TypeError, "Cannot delete the REAL8TimeSeries attribute");
        return -1;
    }
    if(!PyObject_TypeCheck(new, &pylal_REAL8TimeSeries_Type)){
        PyErr_SetObject(PyExc_TypeError, new);
        return -1;
    }

    if(old) Py_DECREF(old);

    Py_INCREF(new);
    old=new;

    pylal_REAL8TimeSeries* newvec=(pylal_REAL8TimeSeries*)new;
    internal=newvec->series;
    return 0;
}

PyObject* getREAL8TimeSeriesFromData(REAL8TimeSeries* internal,PyObject* owner){
    if(!owner){
        Py_INCREF(Py_None);
        return Py_None;
    }else{
        
        return pylal_REAL8TimeSeries_new(internal,owner);
    }
}

int setLALDetectorFromData(LALDetector* target,PyObject* old,PyObject* origin){
    if (origin == NULL) {
        PyErr_SetString(PyExc_TypeError, "Cannot delete the LALDetector attribute");
        return -1;
    }
    
    if(!PyObject_TypeCheck(origin, &pylal_LALDetector_Type)){
        PyErr_SetObject(PyExc_TypeError, origin);
        return -1;
    }
    if(old) Py_DECREF(old);
    
    pylal_LALDetector* origindet=(pylal_LALDetector*)origin;
    Py_INCREF(origin);
    old=origin;
    memcpy((void*)target,(void*)&origindet->detector,sizeof(LALDetector));
    
    return 0;
}

PyObject* getLALDetectorFromData(LALDetector internal){
    PyObject *empty_tuple = PyTuple_New(0);
    pylal_LALDetector *obj = (pylal_LALDetector *) PyType_GenericNew(&pylal_LALDetector_Type,empty_tuple,NULL);
    Py_DECREF(empty_tuple);

    if(!obj) {
        return NULL;
    }
    
    memcpy((void*)&obj->detector,(void*)&internal,sizeof(LALDetector));
    return (PyObject *) obj;
}

int setLIGOTimeGPSFromData(LIGOTimeGPS target,PyObject* origin){
    if(!PyObject_TypeCheck(origin, &pylal_LIGOTimeGPS_Type)){
        PyErr_SetObject(PyExc_TypeError, origin);
        return -1;
    }
    
    pylal_LIGOTimeGPS* origingpstime=(pylal_LIGOTimeGPS*)origin;
    XLALGPSSet(&target, origingpstime->gps.gpsSeconds, origingpstime->gps.gpsNanoSeconds);
    
    return 0;
}


/*
 * ============================================================================
 *
 *                            Module methods
 *
 * ============================================================================
 */

/*
 * ============================================================================
 *
 *                            Module Registration
 *
 * ============================================================================
 */
 
static PyMethodDef module_methods[] = {
    {NULL}  /* Sentinel */
};

PyMODINIT_FUNC
init_lalinference(void)
{
    PyObject *m;

    import_array();
    
    pylal_complex16frequencyseries_import();
    pylal_complex16timeseries_import();
    pylal_real8fftplan_import();
    pylal_real8frequencyseries_import();
    pylal_real8timeseries_import();
    pylal_ligotimegps_import();
    pylal_real8window_import();

    LALVariables_import();
    BarycenterInput_import();
    EmpherisData_import();
    m = Py_InitModule3(MODULE_NAME,module_methods,LIDocString);  
}
