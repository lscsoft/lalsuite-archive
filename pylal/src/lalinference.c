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

#include "lalinference/LALIFOData.h"
#include "lalinference/BarycenterInput.h"
#include "lalinference/EphemerisData.h"
#include "lalinference/LALVariables.h"
#include "lalinference/PosVelAcc.h"
#include "lalinference/LALInferenceRunState.h"

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

int setEphemerisDataFromData(EphemerisData* internal,PyObject* old,PyObject* new){
    /* require li_BarycenterInput*/
    if (new == NULL) {
        PyErr_SetString(PyExc_TypeError, "Cannot delete the BarycenterInput attribute");
        return -1;
    }
    if(!PyObject_TypeCheck(new, &li_EphemerisData_Type)){
        PyErr_SetObject(PyExc_TypeError, new);
        return -1;
    }

    if(old) Py_DECREF(old);

    Py_INCREF(new);
    old=new;

    li_EphemerisData* newephem=(li_EphemerisData*)new;
    internal=newephem->data;
    return 0;
}

PyObject* getEphemerisDataFromData(EphemerisData* internal,PyObject* owner){
    if(!owner){
        Py_INCREF(Py_None);
        return Py_None;
    }else{
        
        return li_EphemerisData_new(internal,owner);
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

int setLALDetectorFromData(LALDetector target,PyObject* old,PyObject* origin){
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
    target=origindet->detector;
    
    return 0;
}

PyObject* getLALDetectorFromData(LALDetector internal,PyObject* owner){
    if(!owner){
        Py_INCREF(Py_None);
        return Py_None;
    }else{
        
        return pylal_LALDetector_new(internal,owner);
    }
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

PyObject* getLIGOTimeGPSFromData(LIGOTimeGPS target){
    
    return pylal_LIGOTimeGPS_new(target);
}

/*
 * ============================================================================
 *
 *                            EphemerisData
 *
 * ============================================================================
 */

/*Methods*/

 /* Destructor for LALIFOData */
static void EphemerisData_dealloc(li_EphemerisData *self)
{
    self->ob_type->tp_free((PyObject *)self);
}

static int EphemerisData__init__(li_EphemerisData *self, PyObject *args, PyObject *kwds)
{
    self->data=(EphemerisData*)malloc(sizeof(EphemerisData));
    memset((void*)self->data,0,sizeof(EphemerisData));
    return 0;
}

static PyObject* EphemerisData__new__(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    
    li_EphemerisData *obj = (li_EphemerisData*) PyType_GenericNew(type, args, kwds);
    if(!obj)
        return NULL;
    obj->data=(EphemerisData*)malloc(sizeof(EphemerisData));
    memset((void*)obj->data,0,sizeof(EphemerisData));
    obj->ephemE=NULL;
    obj->ephemS=NULL;
    obj->ephiles=NULL;
    return (PyObject*)obj;
}

int setPosVelAccFromData(PosVelAcc* internal,PyObject* old,PyObject* new){
    if (new == NULL) {
        PyErr_SetString(PyExc_TypeError, "Cannot delete the PosVelAcc attribute");
        return -1;
    }
    
    if(!PyObject_TypeCheck(new, &li_PosVelAcc_Type)){
        PyErr_SetObject(PyExc_TypeError, new);
        return -1;
    }

    if(old) Py_DECREF(old);
    Py_INCREF(new);
    old=new;

    li_PosVelAcc* newpva=(li_PosVelAcc*)new;
    internal=newpva->data;
    
    return 0;
}

PyObject* getPosVelAccFromData(PosVelAcc* internal,PyObject* pva){
    if(pva==NULL){
        Py_INCREF(Py_None);
        return Py_None;
    }else{
        return li_PosVelAcc_new(internal,pva);
    }
}

int setTwoTupleOfCharsFromData(EphemerisFilenames internal,PyObject* old,PyObject* new){
    
    return 0;
}

PyObject* getTwoTupleOfCharsFromData(PyObject* pva){
    if(pva==NULL){
        Py_INCREF(Py_None);
        return Py_None;
    }else{
        Py_INCREF(pva);
        return pva;
    }
}


/*ephemE*/
static PyObject* EphemerisData_getephemE(li_EphemerisData *self, void *closure){return getPosVelAccFromData(self->data->ephemE,self->ephemE);}
static int EphemerisData_setephemE(li_EphemerisData *self, PyObject *value, void *closure){return setPosVelAccFromData(self->data->ephemE,self->ephemE,value);}

/*ephemS*/
static PyObject* EphemerisData_getephemS(li_EphemerisData *self, void *closure){return getPosVelAccFromData(self->data->ephemS,self->ephemS);}
static int EphemerisData_setephemS(li_EphemerisData *self, PyObject *value, void *closure){return setPosVelAccFromData(self->data->ephemS,self->ephemS,value);}

/*ephemS*/
static PyObject* EphemerisData_getephiles(li_EphemerisData *self, void *closure){return getTwoTupleOfCharsFromData(self->ephiles);}
static int EphemerisData_setephiles(li_EphemerisData *self, PyObject *value, void *closure){return setTwoTupleOfCharsFromData(self->data->ephiles,self->ephiles,value);}


static PyMethodDef EphemerisData_methods[]= {
    {NULL} /* Sentinel */
};

static PyGetSetDef EphemerisData_getseters[] = {
    {"ephemE",(getter)EphemerisData_getephemE,(setter)EphemerisData_setephemE,"ephemE",NULL},
    {"ephemS",(getter)EphemerisData_getephemS,(setter)EphemerisData_setephemS,"ephemS",NULL},
    {"ephiles",(getter)EphemerisData_getephiles,(setter)EphemerisData_setephiles,"ephiles",NULL},
    
    {NULL}  /* Sentinel */
};

static struct PyMemberDef EphemerisData_members[] = {
    
    {"nentriesE", T_INT, offsetof(li_EphemerisData, data)+offsetof(EphemerisData,nentriesE), 0, "nentriesE"},
    {"nentriesS", T_INT, offsetof(li_EphemerisData, data)+offsetof(EphemerisData,nentriesS), 0, "nentriesS"},
    {"dtEtable", T_DOUBLE, offsetof(li_EphemerisData, data)+offsetof(EphemerisData,dtEtable), 0, "dtEtable"},
    {"dtStable", T_DOUBLE, offsetof(li_EphemerisData, data)+offsetof(EphemerisData,dtStable), 0, "dtStable"},
    
    {NULL,}
};

static PyTypeObject li_ephemerisdata_type = {
    PyObject_HEAD_INIT(NULL)
    0,              /* obj_size - unused (must be 0) */
    "lalinference.EphemerisData",    /* tp_name, name of type */
    sizeof(li_EphemerisData),  /* tp_basicsize */
    0,              /* tp_itemsize, need to check */
    (destructor)EphemerisData_dealloc,  /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,        /*tp_flags*/
    "", /* tp_doc */
    0,                     /* tp_traverse */
    0,                     /* tp_clear */
    0,                     /* tp_richcompare */
    0,                     /* tp_weaklistoffset */
    0,                     /* tp_iter */
    0,                     /* tp_iternext */
    EphemerisData_methods,             /* tp_methods */
    EphemerisData_members,             /* tp_members */
    EphemerisData_getseters,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)EphemerisData__init__,      /* tp_init */
    0,                         /* tp_alloc */
    EphemerisData__new__,                 /* tp_new */
};


/*
 * ============================================================================
 *
 *                            BarycenterInput
 *
 * ============================================================================
 */

/*Methods*/

 /* Destructor for LALIFOData */
static void BarycenterInput_dealloc(li_BarycenterInput *self)
{
    self->ob_type->tp_free((PyObject *)self);
}

static int BarycenterInput__init__(li_BarycenterInput *self, PyObject *args, PyObject *kwds)
{
    self->data=(BarycenterInput*)malloc(sizeof(BarycenterInput));
    memset((void*)self->data,0,sizeof(BarycenterInput));
    return 0;
}

static PyObject* BarycenterInput__new__(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    
    li_BarycenterInput *obj = (li_BarycenterInput*) PyType_GenericNew(type, args, kwds);
    if(!obj)
        return NULL;
    obj->data=(BarycenterInput*)malloc(sizeof(BarycenterInput));
    memset((void*)obj->data,0,sizeof(BarycenterInput));
    obj->tgps=NULL;
    obj->site=NULL;
    return (PyObject*)obj;
}

/*tgps*/
static PyObject* BarycenterInput_gettgps(li_BarycenterInput *self, void *closure){return pylal_LIGOTimeGPS_new(self->data->tgps);}
static int BarycenterInput_settgps(li_BarycenterInput *self, PyObject *value, void *closure){return setLIGOTimeGPSFromData(self->data->tgps,value);}

/*site*/
static PyObject* BarycenterInput_getsite(li_BarycenterInput *self, void *closure){return getLALDetectorFromData(self->data->site,self->site);}
static int BarycenterInput_setsite(li_BarycenterInput *self, PyObject *value, void *closure){return setLALDetectorFromData(self->data->site,self->site,value);}


static PyMethodDef BarycenterInput_methods[]= {
    {NULL} /* Sentinel */
};

static PyGetSetDef BarycenterInput_getseters[] = {
    {"tgps",(getter)BarycenterInput_gettgps,(setter)BarycenterInput_settgps,"name",NULL},
    {"site",(getter)BarycenterInput_getsite,(setter)BarycenterInput_setsite,"name",NULL},
    
    {NULL}  /* Sentinel */
};

static struct PyMemberDef BarycenterInput_members[] = {
    //REAL8's
    {"alpha", T_DOUBLE, offsetof(li_BarycenterInput, data)+offsetof(BarycenterInput,alpha), 0, "alpha"},
    {"delta", T_DOUBLE, offsetof(li_BarycenterInput, data)+offsetof(BarycenterInput,delta), 0, "delta"},
    {"dInv", T_DOUBLE, offsetof(li_BarycenterInput, data)+offsetof(BarycenterInput,dInv), 0, "dInv"},
    
    {NULL,}
};

static PyTypeObject li_barycenterinput_type = {
    PyObject_HEAD_INIT(NULL)
    0,              /* obj_size - unused (must be 0) */
    "lalinference.BarycenterInput",    /* tp_name, name of type */
    sizeof(li_BarycenterInput),  /* tp_basicsize */
    0,              /* tp_itemsize, need to check */
    (destructor)BarycenterInput_dealloc,  /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,        /*tp_flags*/
    "", /* tp_doc */
    0,                     /* tp_traverse */
    0,                     /* tp_clear */
    0,                     /* tp_richcompare */
    0,                     /* tp_weaklistoffset */
    0,                     /* tp_iter */
    0,                     /* tp_iternext */
    BarycenterInput_methods,             /* tp_methods */
    BarycenterInput_members,             /* tp_members */
    BarycenterInput_getseters,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)BarycenterInput__init__,      /* tp_init */
    0,                         /* tp_alloc */
    BarycenterInput__new__,                 /* tp_new */
};

/*
 * ============================================================================
 *
 *                            LALIFOData
 *
 * ============================================================================
 */

/*Methods*/

 /* Destructor for LALIFOData */
static void LALIFOData_dealloc(li_LALIFOData *self)
{
    self->ob_type->tp_free((PyObject *)self);
}

static int LALIFOData__init__(li_LALIFOData *self, PyObject *args, PyObject *kwds)
{
    self->data=(LALIFOData*)malloc(sizeof(LALIFOData));
    memset((void*)self->data,0,sizeof(LALIFOData));
    return 0;
}

static PyObject* LALIFOData__new__(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    
    li_LALIFOData *obj = (li_LALIFOData*) PyType_GenericNew(type, args, kwds);
    if(!obj)
        return NULL;
    obj->data=(LALIFOData*)malloc(sizeof(LALIFOData));
    memset((void*)obj->data,0,sizeof(LALIFOData));
    obj->owner=NULL;
    obj->name=NULL;
    obj->modelParams=NULL;
    obj->dataParams=NULL;
    obj->modelDomain=NULL;
    obj->timeData=NULL;
    obj->timeModelhPlus=NULL;
    obj->timeModelhCross=NULL;
    obj->whiteTimeData=NULL;
    obj->windowedTimeData=NULL;
    obj->timeDomainNoiseWeights=NULL;
    obj->freqData=NULL;
    obj->freqModelhPlus=NULL;
    obj->freqModelhCross=NULL;
    obj->whiteFreqData=NULL;
    obj->compTimeData=NULL;
    obj->compModelData=NULL;
    obj->oneSidedNoisePowerSpectrum=NULL;
    obj->timeToFreqFFTPlan=NULL;
    obj->freqToTimeFFTPlan=NULL;
    obj->epoch=NULL;
    obj->window=NULL;
    obj->next=NULL;
    obj->bary=NULL;
    obj->ephem=NULL;
    
    return (PyObject*)obj;
}

/***********getsetters******************/

/*name*/

static PyObject* LALIFOData_getname(li_LALIFOData *self, void *closure)
{
    if(self->name){
        Py_INCREF(self->name);
        return self->name;
    }
    else{
        Py_INCREF(Py_None);
        return Py_None;
    }
}

static int LALIFOData_setname(li_LALIFOData *self, PyObject *value, void *closure)
{
    if(value==NULL){
        PyErr_SetString(PyExc_TypeError, "Cannot delete the attribute");
        return -1;
    }
    
    if (!PyString_Check(value)) {
        PyErr_SetString(PyExc_TypeError, 
                "The name attribute value must be a string");
        return -1;
    }

    if(self->name) Py_DECREF(self->name);
    Py_INCREF(value);
    self->name=value;
    strcpy(self->data->name,PyString_AsString(value));
    
    return 0;
}


PyObject* LALIFOData_getLALVariables(LALVariables* internallv,PyObject* owner){
    if(owner){
        return LALVariables_new(internallv,owner);
    }
    else {
        Py_INCREF(Py_None);
        return Py_None;
    }
}

int LALIFOData_setLALVariables(LALVariables* internallv,PyObject* oldlv,PyObject* newlv){
    if (!newlv) {
        PyErr_SetString(PyExc_TypeError, "Cannot delete the attribute");
        return -1;
    }
    if (!PyObject_TypeCheck(newlv, &li_LALVariables_Type)){
        PyErr_SetString(PyExc_TypeError, "The attribute must be a lalinference.LALVariables object.");
        return -1;
    }

    if(oldlv) Py_DECREF(oldlv);

    Py_INCREF(newlv);
    oldlv=newlv;

    li_LALVariables* valuelv=(li_LALVariables*)newlv;
    internallv=valuelv->vars;
    
    return 0;
}

/*modelParams*/
static PyObject* LALIFOData_getmodelParams(li_LALIFOData *self, void *closure) {return LALIFOData_getLALVariables(self->data->modelParams,self->modelParams);};
static int LALIFOData_setmodelParams(li_LALIFOData *self, PyObject *value, void *closure){return LALIFOData_setLALVariables(self->data->modelParams,(PyObject*)self->modelParams,value);}

/*dataParams*/
static PyObject* LALIFOData_getdataParams(li_LALIFOData *self, void *closure) {return LALIFOData_getLALVariables(self->data->dataParams,self->dataParams);};
static int LALIFOData_setdataParams(li_LALIFOData *self, PyObject *value, void *closure){return LALIFOData_setLALVariables(self->data->dataParams,(PyObject*)self->dataParams,value);}


static PyObject* LALIFOData_getmodelDomain(li_LALIFOData *self, void *closure)
{

    if(self->data->modelDomain==timeDomain){
        return PyString_FromString("timeDomain");
    }
    else if(self->data->modelDomain==frequencyDomain){
        return PyString_FromString("freqDomain");
    }
    else return NULL;
}

static int LALIFOData_setmodelDomain(li_LALIFOData *self, PyObject *value, void *closure)
{
    if (value == NULL) {
        PyErr_SetString(PyExc_TypeError, "Cannot delete the modelDomain attribute");
        return -1;
    }
    
    if (! PyString_Check(value)) {
        PyErr_SetString(PyExc_TypeError, 
                    "The modelDomain attribute value is a string.");
        return -1;
    }
    
    Py_INCREF(value);
    char* name=PyString_AsString(value);

    LALDomain var;

    if(!strcmp(name,"timeDomain")){
        var=timeDomain;
    }
    else if(!strcmp(name,"frequencyDomain")){
        var=frequencyDomain;
    }
    else{
        die("modelDomain must be one of the strings 'frequencyDomain' or 'timeDomain'");
    }

    return 0;
}

/* timeData */
static int LALIFOData_settimeData(li_LALIFOData *self, PyObject *value, void *closure) {return setREAL8TimeSeriesFromData(self->data->timeData,(PyObject*)self->timeData,value);}
static PyObject* LALIFOData_gettimeData(li_LALIFOData *self, void *closure) {return getREAL8TimeSeriesFromData(self->data->timeData,(PyObject*)self->timeData);}//{return getREAL8TimeSeriesFromData(self->data->timeData);}

/* timeModelhPlus */
static int LALIFOData_settimeModelhPlus(li_LALIFOData *self, PyObject *value, void *closure) {return setREAL8TimeSeriesFromData(self->data->timeModelhPlus,(PyObject*)self->timeModelhPlus,value);}
static PyObject* LALIFOData_gettimeModelhPlus(li_LALIFOData *self, void *closure) {return getREAL8TimeSeriesFromData(self->data->timeModelhPlus,(PyObject*)self->timeModelhPlus);}//{return getREAL8TimeSeriesFromData(self->data->timeModelhPlus);}

/* timeModelhCross */
static int LALIFOData_settimeModelhCross(li_LALIFOData *self, PyObject *value, void *closure) {return setREAL8TimeSeriesFromData(self->data->timeModelhCross,(PyObject*)self->timeModelhCross,value);}
static PyObject* LALIFOData_gettimeModelhCross(li_LALIFOData *self, void *closure) {return getREAL8TimeSeriesFromData(self->data->timeModelhCross,(PyObject*)self->timeModelhCross);}//{return getREAL8TimeSeriesFromData(self->data->timeModelhCross);}

/* whiteTimeData */
static int LALIFOData_setwhiteTimeData(li_LALIFOData *self, PyObject *value, void *closure) {return setREAL8TimeSeriesFromData(self->data->whiteTimeData,(PyObject*)self->whiteTimeData,value);}
static PyObject* LALIFOData_getwhiteTimeData(li_LALIFOData *self, void *closure) {return getREAL8TimeSeriesFromData(self->data->whiteTimeData,(PyObject*)self->whiteTimeData);}//{return getREAL8TimeSeriesFromData(self->data->whiteTimeData);}

/* windowedTimeData */
static int LALIFOData_setwindowedTimeData(li_LALIFOData *self, PyObject *value, void *closure) {return setREAL8TimeSeriesFromData(self->data->windowedTimeData,(PyObject*)self->windowedTimeData,value);}
static PyObject* LALIFOData_getwindowedTimeData(li_LALIFOData *self, void *closure) {return getREAL8TimeSeriesFromData(self->data->windowedTimeData,(PyObject*)self->windowedTimeData);}//getREAL8TimeSeriesFromData(self->data->windowedTimeData);}

/* timeDomainNoiseWeights */
static int LALIFOData_settimeDomainNoiseWeights(li_LALIFOData *self, PyObject *value, void *closure) {return setREAL8TimeSeriesFromData(self->data->timeDomainNoiseWeights,(PyObject*)self->timeDomainNoiseWeights,value);}
static PyObject* LALIFOData_gettimeDomainNoiseWeights(li_LALIFOData *self, void *closure) {return getREAL8TimeSeriesFromData(self->data->timeDomainNoiseWeights,(PyObject*)self->timeDomainNoiseWeights);}//{return getREAL8TimeSeriesFromData(self->data->timeDomainNoiseWeights);}

/* oneSidedNoisePowerSpectrum */
static int LALIFOData_setoneSidedNoisePowerSpectrum(li_LALIFOData *self, PyObject *value, void *closure) {return setREAL8FrequencySeriesFromData(self->data->oneSidedNoisePowerSpectrum,self->oneSidedNoisePowerSpectrum,value);}
static PyObject* LALIFOData_getoneSidedNoisePowerSpectrum(li_LALIFOData *self, void *closure) {return getREAL8FrequencySeriesFromData(self->data->oneSidedNoisePowerSpectrum,(PyObject*)self->oneSidedNoisePowerSpectrum);}

/*freqData*/
static int LALIFOData_setfreqData(li_LALIFOData *self, PyObject *value, void *closure) {return setCOMPLEX16FrequencySeriesFromData(self->data->freqData,self->freqData,value);}
static PyObject* LALIFOData_getfreqData(li_LALIFOData *self, void *closure) {return getCOMPLEX16FrequencySeriesFromData(self->data->freqData,(PyObject*)self->freqData);}

/*freqModelhPlus*/
static int LALIFOData_setfreqModelhPlus(li_LALIFOData *self, PyObject *value, void *closure) {return setCOMPLEX16FrequencySeriesFromData(self->data->freqModelhPlus,self->freqModelhPlus,value);}
static PyObject* LALIFOData_getfreqModelhPlus(li_LALIFOData *self, void *closure) {return getCOMPLEX16FrequencySeriesFromData(self->data->freqModelhPlus,(PyObject*)self->freqModelhPlus);}

/*freqModelhCross*/
static int LALIFOData_setfreqModelhCross(li_LALIFOData *self, PyObject *value, void *closure) {return setCOMPLEX16FrequencySeriesFromData(self->data->freqModelhCross,self->freqModelhCross,value);}
static PyObject* LALIFOData_getfreqModelhCross(li_LALIFOData *self, void *closure){return getCOMPLEX16FrequencySeriesFromData(self->data->freqModelhCross,(PyObject*)self->freqModelhCross);}

/*whiteFreqData*/
static int LALIFOData_setwhiteFreqData(li_LALIFOData *self, PyObject *value, void *closure) {return setCOMPLEX16FrequencySeriesFromData(self->data->whiteFreqData,self->whiteFreqData,value);}
static PyObject* LALIFOData_getwhiteFreqData(li_LALIFOData *self, void *closure){return getCOMPLEX16FrequencySeriesFromData(self->data->whiteFreqData,(PyObject*)self->whiteFreqData);}

/*compTimeData*/
static int LALIFOData_setcompTimeData(li_LALIFOData *self, PyObject *value, void *closure) {return setCOMPLEX16TimeSeriesFromData(self->data->compTimeData,self->compTimeData,value);}
static PyObject* LALIFOData_getcompTimeData(li_LALIFOData *self, void *closure) {return getCOMPLEX16TimeSeriesFromData(self->data->compTimeData,self->compTimeData);}

/*compModelData*/
static int LALIFOData_setcompModelData(li_LALIFOData *self, PyObject *value, void *closure) {return setCOMPLEX16TimeSeriesFromData(self->data->compModelData,self->compModelData,value);}
static PyObject* LALIFOData_getcompModelData(li_LALIFOData *self, void *closure) {return getCOMPLEX16TimeSeriesFromData(self->data->compModelData,(PyObject*)self->compModelData);}

/*timeToFreqFFTPlan*/
static int LALIFOData_settimeToFreqFFTPlan(li_LALIFOData *self, PyObject *value, void *closure) {return setREAL8FFTPlanFromData(self->data->timeToFreqFFTPlan,self->timeToFreqFFTPlan,value);}
static PyObject* LALIFOData_gettimeToFreqFFTPlan(li_LALIFOData *self, void *closure) {return getREAL8FFTPlanFromData(self->data->timeToFreqFFTPlan,(PyObject*)self->timeToFreqFFTPlan);}

/*freqToTimeFFTPlan*/
static int LALIFOData_setfreqToTimeFFTPlan(li_LALIFOData *self, PyObject *value, void *closure) {return setREAL8FFTPlanFromData(self->data->freqToTimeFFTPlan,self->freqToTimeFFTPlan,value);}
static PyObject* LALIFOData_getfreqToTimeFFTPlan(li_LALIFOData *self, void *closure) {return getREAL8FFTPlanFromData(self->data->freqToTimeFFTPlan,(PyObject*)self->freqToTimeFFTPlan);}

/* epoch */
static int LALIFOData_setepoch(li_LALIFOData *self, PyObject *value, void *closure) {return setLIGOTimeGPSFromData(self->data->epoch,value);}
static PyObject* LALIFOData_getepoch(li_LALIFOData *self, void *closure) {return pylal_LIGOTimeGPS_new(self->data->epoch);}

/* window */
static int LALIFOData_setwindow(li_LALIFOData *self, PyObject *value, void *closure) {return setREAL8WindowFromData(self->data->window,self->window,value);}
static PyObject* LALIFOData_getwindow(li_LALIFOData *self, void *closure) {return getREAL8WindowFromData(self->data->window,(PyObject*)self->window);}

/* next */
static int LALIFOData_setnext(li_LALIFOData *self, PyObject *value, void *closure) {return setLALIFODataFromData(self->data->next,self->next,value);}
static PyObject* LALIFOData_getnext(li_LALIFOData *self, void *closure) {return getLALIFODataFromData(self->data->next,(PyObject*)self->next);}

/* bary */
static int LALIFOData_setbary(li_LALIFOData *self, PyObject *value, void *closure) {return setBarycenterInputFromData(self->data->bary,self->bary,value);}
static PyObject* LALIFOData_getbary(li_LALIFOData *self, void *closure) {return getBarycenterInputFromData(self->data->bary,(PyObject*)self->bary);}

/* ephem */
static int LALIFOData_setephem(li_LALIFOData *self, PyObject *value, void *closure) {return setEphemerisDataFromData(self->data->ephem,self->ephem,value);}
static PyObject* LALIFOData_getephem(li_LALIFOData *self, void *closure) {return getEphemerisDataFromData(self->data->ephem,(PyObject*)self->ephem);}

/**getsetters registration struct**/

static PyGetSetDef LALIFOData_getseters[] = {
    {"name",(getter)LALIFOData_getname,(setter)LALIFOData_setname,"name",NULL},
    //LALVariables
    {"modelParams",(getter)LALIFOData_getmodelParams,(setter)LALIFOData_setmodelParams,"modelParams",NULL},
    {"dataParams",(getter)LALIFOData_getdataParams,(setter)LALIFOData_setdataParams,"dataParams",NULL},
    //LALDomain ...represented as string
    {"modelDomain",(getter)LALIFOData_getmodelDomain,(setter)LALIFOData_setmodelDomain,"modelDomain",NULL},
    //REAL8TimeSeries
    {"timeData",(getter)LALIFOData_gettimeData,(setter)LALIFOData_settimeData,"timeData",NULL},
    {"timeModelhPlus",(getter)LALIFOData_gettimeModelhPlus,(setter)LALIFOData_settimeModelhPlus,"timeModelhPlus",NULL},
    {"timeModelhCross",(getter)LALIFOData_gettimeModelhCross,(setter)LALIFOData_settimeModelhCross,"timeModelhCross",NULL},
    {"whiteTimeData",(getter)LALIFOData_getwhiteTimeData,(setter)LALIFOData_setwhiteTimeData,"whiteTimeData",NULL},
    {"windowedTimeData",(getter)LALIFOData_getwindowedTimeData,(setter)LALIFOData_setwindowedTimeData,"windowedTimeData",NULL},
    {"timeDomainNoiseWeights",(getter)LALIFOData_gettimeDomainNoiseWeights,(setter)LALIFOData_settimeDomainNoiseWeights,"timeDomainNoiseWeights",NULL},
    //COMPLEX16FrequencySeries
    {"freqData",(getter)LALIFOData_getfreqData,(setter)LALIFOData_setfreqData,"freqData",NULL},
    {"freqModelhPlus",(getter)LALIFOData_getfreqModelhPlus,(setter)LALIFOData_setfreqModelhPlus,"freqModelhPlus",NULL},
    {"freqModelhCross",(getter)LALIFOData_getfreqModelhCross,(setter)LALIFOData_setfreqModelhCross,"freqModelhCross",NULL},
    {"whiteFreqData",(getter)LALIFOData_getwhiteFreqData,(setter)LALIFOData_setwhiteFreqData,"whiteFreqData",NULL},
    //COMPLEX16TimeSeries
    {"compTimeData",(getter)LALIFOData_getcompTimeData,(setter)LALIFOData_setcompTimeData,"compTimeData",NULL},
    {"compModelData",(getter)LALIFOData_getcompModelData,(setter)LALIFOData_setcompModelData,"compModelData",NULL},
    //REAL8FrequencySeries
    {"oneSidedNoisePowerSpectrum",(getter)LALIFOData_getoneSidedNoisePowerSpectrum,(setter)LALIFOData_setoneSidedNoisePowerSpectrum,"oneSidedNoisePowerSpectrum",NULL},
    //REAL8FFTPlan
    {"timeToFreqFFTPlan",(getter)LALIFOData_gettimeToFreqFFTPlan,(setter)LALIFOData_settimeToFreqFFTPlan,"timeToFreqFFTPlan",NULL},
    {"freqToTimeFFTPlan",(getter)LALIFOData_getfreqToTimeFFTPlan,(setter)LALIFOData_setfreqToTimeFFTPlan,"freqToTimeFFTPlan",NULL},
    //LIGOTimeGPS
    {"epoch",(getter)LALIFOData_getepoch,(setter)LALIFOData_setepoch,"epoch",NULL},
    //REAL8Window
    {"window",(getter)LALIFOData_getwindow,(setter)LALIFOData_setwindow,"window",NULL},
    //LALIFOData
    {"next",(getter)LALIFOData_getnext,(setter)LALIFOData_setnext,"next",NULL},
    //BarycenterInput
    {"bary",(getter)LALIFOData_getbary,(setter)LALIFOData_setbary,"bary",NULL},
    //EphemerisData
    {"ephem",(getter)LALIFOData_getephem,(setter)LALIFOData_setephem,"ephem",NULL},
    {NULL}  /* Sentinel */
};

static struct PyMemberDef LALIFOData_members[] = {
    //REAL8's
    {"nullloglikelihood", T_DOUBLE, offsetof(li_LALIFOData, data)+offsetof(LALIFOData,nullloglikelihood), 0, "nullloglikelihood"},
    {"loglikelihood", T_DOUBLE, offsetof(li_LALIFOData, data)+offsetof(LALIFOData,loglikelihood), 0, "loglikelihood"},
    {"acceptedloglikelihood", T_DOUBLE, offsetof(li_LALIFOData, data)+offsetof(LALIFOData,acceptedloglikelihood), 0, "acceptedloglikelihood"},
    {"fPlus", T_DOUBLE, offsetof(li_LALIFOData, data)+offsetof(LALIFOData,fPlus), 0, "fPlus"},
    {"fCross", T_DOUBLE, offsetof(li_LALIFOData, data)+offsetof(LALIFOData,fCross), 0, "fCross"},
    {"timeshift", T_DOUBLE, offsetof(li_LALIFOData, data)+offsetof(LALIFOData,timeshift), 0, "timeshift"},
    {"fLow", T_DOUBLE, offsetof(li_LALIFOData, data)+offsetof(LALIFOData,fLow), 0, "fLow"},
    {"fHigh", T_DOUBLE,offsetof(li_LALIFOData, data)+offsetof(LALIFOData,fHigh), 0, "fHigh"},
    
    {NULL,}
};

static PyMethodDef LALIFOData_methods[]= {
    /* {"name", (PyCFunction)function, METH_NOARGS, "DESCR"}, */
    {NULL} /* Sentinel */
};

static PyTypeObject li_lalifodata_type = {
    PyObject_HEAD_INIT(NULL)
    0,              /* obj_size - unused (must be 0) */
    "lalinference.LALIFOData",    /* tp_name, name of type */
    sizeof(li_LALIFOData),  /* tp_basicsize */
    0,              /* tp_itemsize, need to check */
    (destructor)LALIFOData_dealloc,  /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,        /*tp_flags*/
    "LALInference LALIFOData object.", /* tp_doc */
    0,                     /* tp_traverse */
    0,                     /* tp_clear */
    0,                     /* tp_richcompare */
    0,                     /* tp_weaklistoffset */
    0,                     /* tp_iter */
    0,                     /* tp_iternext */
    LALIFOData_methods,             /* tp_methods */
    LALIFOData_members,             /* tp_members */
    LALIFOData_getseters,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)LALIFOData__init__,      /* tp_init */
    0,                         /* tp_alloc */
    LALIFOData__new__,                 /* tp_new */
    
};

/*
 * ============================================================================
 *
 *                            LALInferenceRunState
 *
 * ============================================================================
 */

/*Methods*/

 /* Destructor for LALInferenceRunState */
static void LALInferenceRunState_dealloc(li_LALInferenceRunState *self)
{
    self->ob_type->tp_free((PyObject *)self);
}

static int LALInferenceRunState__init__(li_LALInferenceRunState *self, PyObject *args, PyObject *kwds)
{
    self->data=(LALInferenceRunState*)malloc(sizeof(LALInferenceRunState));
    memset((void*)self->data,0,sizeof(LALInferenceRunState));
    return 0;
}

static PyObject* LALInferenceRunState__new__(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    
    li_LALInferenceRunState *obj = (li_LALInferenceRunState*) PyType_GenericNew(type, args, kwds);
    if(!obj)
        return NULL;
    obj->state=(LALInferenceRunState*)malloc(sizeof(LALInferenceRunState));
    memset((void*)obj->data,0,sizeof(LALInferenceRunState));
    obj->owner=NULL;
    //ProcessParams*
    obj->commandLine=NULL;
    //LALAlgorithm*
    obj->algorithm=NULL;
    //LALEvolveOneStepFunction*
    obj->evolve=NULL;
    //LALPriorFunction*
    obj->prior=NULL;
    //LALLikelihoodFunction*
    obj->likelihood=NULL;
    //LALProposalFunction*
    obj->proposal=NULL;
    //LALTemplateFunction*
    obj->template=NULL;
    //struct tagLALIFOData*
    obj->data=NULL;
    //LALVariables*
    obj->currentParams=NULL;
    obj->priorArgs=NULL;
    obj->proposalArgs=NULL;
    obj->algorithmParams=NULL; /* Parameters which control the running of the algorithm */
    //LALVariables**
    obj->livePoints=NULL; /* Array of live points for Nested Sampling */
    obj->differentialPoints=NULL;
    //gsl_rng*
    obj->GSLrandom=NULL; 
    return (PyObject*)obj;
}

/***********getsetters******************/
/*algorithm*/
static PyObject* LALInferenceRunState_getalgorithm(li_LALInferenceRunState *self, void *closure) {return getLALAlgorithmFromData(self->state->algorithm,self->algorithm);};
static int LALInferenceRunState_setalgorithm(li_LALInferenceRunState *self, PyObject *value, void *closure){return setLALAlgorithmFromData(self->state->algorithm,(PyObject*)self->algorithm,value);}
/*evolve*/
static PyObject* LALInferenceRunState_getevolve(li_LALInferenceRunState *self, void *closure) {return getLALEvolveOneStepFunctionFromData(self->state->evolve,self->evolve);};
static int LALInferenceRunState_setevolve(li_LALInferenceRunState *self, PyObject *value, void *closure){return setLALEvolveOneStepFunctionFromData(self->state->evolve,(PyObject*)self->evolve,value);}

/**getsetters registration struct**/

static PyGetSetDef LALInferenceRunState_getseters[] = {
    {"algorithm",(getter)LALInferenceRunState_getalgorithm,(setter)LALInferenceRunState_setalgorithm,"algorithm",NULL},
    {"evolve",(getter)LALInferenceRunState_getevolve,(setter)LALInferenceRunState_setevolve,"evolve",NULL},
    //{"prior",(getter)LALInferenceRunState_getprior,(setter)LALInferenceRunState_setprior,"prior",NULL},
    //{"likelihood",(getter)LALInferenceRunState_getlikelihood,(setter)LALInferenceRunState_setlikelihood,"likelihood",NULL},
    //{"proposal",(getter)LALInferenceRunState_getproposal,(setter)LALInferenceRunState_setproposal,"proposal",NULL},
    //{"template",(getter)LALInferenceRunState_gettemplate,(setter)LALInferenceRunState_settemplate,"template",NULL},
    ////LALVariables*
    //{"currentParams",(getter)LALInferenceRunState_getcurrentParams,(setter)LALInferenceRunState_setcurrentParams,"currentParams",NULL},
    //{"priorArgs",(getter)LALInferenceRunState_getpriorArgs,(setter)LALInferenceRunState_setpriorArgs,"priorArgs",NULL},
    //{"proposalArgs",(getter)LALInferenceRunState_getproposalArgs,(setter)LALInferenceRunState_setproposalArgs,"proposalArgs",NULL},
    //{"algorithmParams",(getter)LALInferenceRunState_getalgorithmParams,(setter)LALInferenceRunState_setalgorithmParams,"algorithmParams",NULL},
    ////LALVariables**
    //{"livePoints",(getter)LALInferenceRunState_getlivePoints,(setter)LALInferenceRunState_setlivePoints,"livePoints",NULL},
    //{"differentialPoints",(getter)LALInferenceRunState_getdifferentialPoints,(setter)LALInferenceRunState_setdifferentialPoints,"differentialPoints",NULL},
    ////LALIFOData*
    //{"data",(getter)LALInferenceRunState_getdata,(setter)LALInferenceRunState_setdata,"data",NULL},
    {NULL}  /* Sentinel */
};

static struct PyMemberDef LALInferenceRunState_members[] = {
    {"currentLikelihood", T_DOUBLE, offsetof(li_LALInferenceRunState,state)+offsetof(LALInferenceRunState,currentLikelihood), 0, "currentLikelihood"},
    {"currentPrior", T_DOUBLE, offsetof(li_LALInferenceRunState,state)+offsetof(LALInferenceRunState,currentPrior), 0, "currentPrior"},
    {"differentialPointsLength", T_UINT, offsetof(li_LALInferenceRunState,state)+offsetof(LALInferenceRunState,differentialPointsLength), 0, "differentialPointsLength"},
    {NULL}
};

static PyMethodDef LALInferenceRunState_methods[]= {
    /* {"name", (PyCFunction)function, METH_NOARGS, "DESCR"}, */
    {NULL} /* Sentinel */
};

static PyTypeObject li_lalinferencerunstate_type = {
    PyObject_HEAD_INIT(NULL)
    0,              /* obj_size - unused (must be 0) */
    "lalinference.LALInferenceRunState",    /* tp_name, name of type */
    sizeof(li_LALInferenceRunState),  /* tp_basicsize */
    0,              /* tp_itemsize, need to check */
    (destructor)LALInferenceRunState_dealloc,  /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,        /*tp_flags*/
    "LALInference LALInferenceRunState object.", /* tp_doc */
    0,                     /* tp_traverse */
    0,                     /* tp_clear */
    0,                     /* tp_richcompare */
    0,                     /* tp_weaklistoffset */
    0,                     /* tp_iter */
    0,                     /* tp_iternext */
    LALInferenceRunState_methods,             /* tp_methods */
    LALInferenceRunState_members,             /* tp_members */
    LALInferenceRunState_getseters,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)LALInferenceRunState__init__,      /* tp_init */
    0,                         /* tp_alloc */
    LALInferenceRunState__new__,                 /* tp_new */
    
};
/*
 * ============================================================================
 *
 *                            LALVariables
 *
 * ============================================================================
 */

/*Methods*/


/* Destructor for LALVariables */
static void LALVariables_dealloc(li_LALVariables *self)
{
    destroyVariables(self->vars);
    self->ob_type->tp_free((PyObject *)self);
}

VariableType LALVariables_infer_litype_from_pyvalue(PyObject* pyvalue){

    VariableType type;
    
    if(PyInt_Check(pyvalue)){
        type=INT8_t;
    }   
    else if(PyFloat_Check(pyvalue)){
        type=REAL8_t;
    }
    else{
        PyErr_SetObject(PyExc_TypeError, pyvalue);
        return -1;
    }

    return type;
}

VariableType LALVariables_convert_string_to_litype(char* typestring){

    VariableType type;
    if(!strcmp(typestring,"INT4")){
        type=INT4_t;
    }
    else if(!strcmp(typestring,"INT8")){
        type=INT8_t;
    }
    else if(!strcmp(typestring,"UINT4")){
        type=UINT4_t;
    }
    else if(!strcmp(typestring,"REAL4")){
        type=REAL4_t;
    }
    else if(!strcmp(typestring,"REAL8")){
        type=REAL8_t;
    }
    else{
        PyErr_SetString(PyExc_TypeError,"LALInference type not found!!");
        return -1;
    }
    return type;
}

void* LALVariables_convert_pyobj_to_livar_value(PyObject* pyvalue,VariableType type){
    void* value=(void *)malloc(typeSize[type]);
    
    if(type==INT4_t){
        INT4 cast_value=((INT4)PyInt_AsLong(pyvalue));
        INT4* cast_valuep=&cast_value;
        memcpy(value,(void*)cast_valuep,typeSize[type]);
    }
    else if(type==INT8_t){
        INT8 cast_value=((INT8)PyInt_AsLong(pyvalue));
        INT8* cast_valuep=&cast_value;
        memcpy(value,(void*)cast_valuep,typeSize[type]);
    }
    else if(type==UINT4_t){
        UINT4 cast_value=(UINT4)((unsigned long int)PyInt_AsLong(pyvalue));
        UINT4* cast_valuep=&cast_value;
        memcpy(value,(void*)cast_valuep,typeSize[type]);
    }   
    else if(type==REAL4_t){
    
        REAL4 cast_value=((REAL4)PyFloat_AsDouble(pyvalue));
        REAL4* cast_valuep=&cast_value;
        memcpy(value,(void*)cast_valuep,typeSize[type]);
        
    }
    else if(type==REAL8_t){
        REAL8 cast_value=((REAL8)PyFloat_AsDouble(pyvalue));
        REAL8* cast_valuep=&cast_value;
        memcpy(value,(void*)cast_valuep,typeSize[type]);
    }
    else{
        PyErr_SetObject(PyExc_TypeError, pyvalue);
        return NULL;
    }
    
    return value;
}

static PyObject* LALVariables_add_variable(li_LALVariables *self,PyObject* args,PyObject* kwds){
    PyObject *pyname=NULL,*pyvalue=NULL,*pyvarytype=NULL,*pytype=NULL;

    char *name;char* temp;void* value;VariableType type;ParamVaryType varytype;
    char *typestring=NULL;

    static char *kwlist[]={"name","value","varytype","type"};
    
    if (! PyArg_ParseTupleAndKeywords(args,kwds,"OOO|O",kwlist,&pyname,&pyvalue,&pyvarytype,&pytype)) return NULL;

    /*Extract name of variable*/
    if(PyString_Check(pyname)){
        name=PyString_AsString(pyname);
    }
    else{
        PyErr_SetObject(PyExc_TypeError, pyname);
        return NULL;
    }

    /*Extract and determine type of parameter value*/

    //If type given convert string to type...
    if(pytype) {
        typestring=PyString_AsString(pytype);
        type=LALVariables_convert_string_to_litype(typestring);
    }
    //...else infer type from python object type (this is more limited).
    else type=LALVariables_infer_litype_from_pyvalue(pyvalue);
    
    value=LALVariables_convert_pyobj_to_livar_value(pyvalue,type);
    
    /*Determine variable wrapping type from string*/
    if(PyString_Check(pyvarytype)){
        temp=PyString_AsString(pyvarytype);
        Py_INCREF(pyvarytype);

        if(!strcmp(temp,"linear")){
            varytype=PARAM_LINEAR;
        }
        else if(!strcmp(temp,"circular")){
            varytype=PARAM_CIRCULAR;
        }
        else if(!strcmp(temp,"fixed")){
            varytype=PARAM_FIXED;
        }
        else if(!strcmp(temp,"output")){
            varytype=PARAM_OUTPUT;
        }
        else {
            
            PyErr_SetObject(PyExc_ValueError,pyvarytype);
            return NULL;
        }
        Py_DECREF(pyvarytype);
    }
    else{
        PyErr_SetObject(PyExc_TypeError, pyvarytype);
        return NULL;
    }

    /* If we survived then call addVariable with self->vars ...*/
    
    addVariable((self->vars),name,value,type,varytype);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject* LALVariables_check_variable(li_LALVariables *self,PyObject* args)
/* Check for existence of name */
{
    PyObject* pyname;
    char* name;
    
    if (! PyArg_ParseTuple(args,"O",&pyname)) return NULL;

    name=PyString_AsString(pyname);

    if(getItem(self->vars,name)){
        Py_RETURN_TRUE;
    }
    else{
        Py_RETURN_FALSE;
    }
}


PyObject* LALVariables_convert_livar_value_to_pyobj(li_LALVariables *pyvars,char* name){
    VariableType type;
    
    PyObject* returnObj=NULL;
    
    type=getVariableType(pyvars->vars,name);
    void* uncastval=getVariable(pyvars->vars,name);
    if(type==INT4_t){
        
        long int cast_val=(long int)(*(INT4*)uncastval);
        returnObj=PyInt_FromLong(cast_val);
    }
    else if(type==INT8_t){
        long long cast_val=(long long)(*(INT8*)uncastval);
        returnObj=PyInt_FromLong(cast_val);
    }
    else if(type==UINT4_t){
        returnObj=PyInt_FromLong(*(UINT4*)uncastval);
    }
    else if(type==REAL4_t){
        float cast_val=(float)(*(REAL4*)uncastval);
        returnObj=PyFloat_FromDouble(cast_val);
    }
    else if(type==REAL8_t){
        returnObj=PyFloat_FromDouble(*(REAL8*)uncastval);
        
    }
    else {
        PyErr_SetString(PyExc_TypeError,"This type of variable cannot be converted to a python type!");
        return NULL;
    }
    
    return returnObj;
}

static int LALVariables_remove_variable(li_LALVariables *self,PyObject* args){
    PyObject* pyname;
    char* name;
    
    if (! PyArg_ParseTuple(args,"O",&pyname)) return -1;

    name=PyString_AsString(pyname);
    removeVariable(self->vars,name);
    
    return 0;
}

static PyObject* LALVariables_get_variable(li_LALVariables *self,PyObject* args){
    PyObject* pyname=NULL,*returnObj=NULL;
    char* name;
    
    if (! PyArg_ParseTuple(args,"O",&pyname)) return NULL;
    name=PyString_AsString(pyname);
    returnObj=LALVariables_convert_livar_value_to_pyobj(self,name);
    return returnObj;
}

static PyObject* LALVariables_get_variable_name(li_LALVariables *self,PyObject* args){
    char* name;
    long var_idx;
    
    if (! PyArg_ParseTuple(args,"i",&var_idx)) return NULL;
    
    name=getVariableName(self->vars,(int)var_idx);
    
    return PyString_FromString(name);
    
}

static int LALVariables_set_variable(li_LALVariables *self,PyObject* args){
    PyObject* pyname;char* name;
    PyObject* pynew_var_value;void* new_var_value;
    VariableType type;

    if (! PyArg_ParseTuple(args,"OO",&pyname,&pynew_var_value)) return -1;

    name=PyString_AsString(pyname);

    type=getVariableType(self->vars,name);

    new_var_value=LALVariables_convert_pyobj_to_livar_value(pynew_var_value,type);

    setVariable(self->vars,name,new_var_value);
    
    return 0;
    
}

char* LALVariables_get_type_string(VariableType type){
    char* type_name=NULL;
    
    if(type==INT4_t){
        type_name="INT4";
    }
    else if(type==UINT4_t){
        type_name="UINT4";
    }
    else if(type==INT8_t){
        type_name="INT8";
    }
    else if(type==REAL4_t){
        type_name="REAL4";
    }
    else if(type==REAL8_t){
        type_name="REAL8";
    }
    
    return type_name;
}

static PyObject* LALVariables_get_variable_type(li_LALVariables *self,PyObject* args){
    VariableType type;
    char* name;
    PyObject* pyname;
    char* type_name=NULL;
    if (! PyArg_ParseTuple(args,"O",&pyname)) return NULL;
    
    name=PyString_AsString(pyname);
    type=getVariableType(self->vars,name);
    
    type_name=LALVariables_get_type_string(type);
    
    if(type_name==NULL){
        Py_INCREF(Py_None);
        return Py_None;
    }    
    
    return PyString_FromString(type_name);
    
}
    
static PyObject* LALVariables_get_variable_type_by_index(li_LALVariables *self,PyObject* args){
    VariableType type;
    long var_idx;
    char* type_name=NULL;
    if (! PyArg_ParseTuple(args,"i",&var_idx)) return NULL;
    
    type=getVariableTypeByIndex(self->vars,(int)var_idx);
    type_name=LALVariables_get_type_string(type);
    
    return PyString_FromString(type_name);
    
}

char* LALVariables_get_varytype_string(ParamVaryType varytype){
    
    char* varytype_name=NULL;
    if(varytype==PARAM_LINEAR){
        varytype_name="linear";
    }
    else if(varytype==PARAM_CIRCULAR){
        varytype_name="circular";
    }
    else if(varytype==PARAM_FIXED){
        varytype_name="fixed";
    }
    else if(varytype==PARAM_OUTPUT){
        varytype_name="output";
    }
    
    return varytype_name;
    
}

static PyObject* LALVariables_get_variable_vary_type(li_LALVariables *self,PyObject* args){
    ParamVaryType varytype;
    char* name;
    PyObject* pyname;
    char* varytype_name=NULL;
    if (! PyArg_ParseTuple(args,"O",&pyname)) return NULL;
    
    name=PyString_AsString(pyname);
    varytype=getVariableVaryType(self->vars,name);
    varytype_name=LALVariables_get_varytype_string(varytype);
    
    if(varytype_name==NULL){
        Py_INCREF(Py_None);
        return Py_None;
    }    
    
    return PyString_FromString(varytype_name);
    
}

static PyObject* LALVariables_get_variable_dimension(li_LALVariables *self){
    long dimension;
    PyObject* pydimension=NULL;
    
    dimension=(long int)getVariableDimension(self->vars);

    pydimension=PyInt_FromLong(dimension);

    return pydimension;
    
}

static PyObject* LALVariables_get_variable_dimension_non_fixed(li_LALVariables *self){
    long dimension;
    PyObject* pydimension=NULL;
    
    dimension=(long int)getVariableDimensionNonFixed(self->vars);

    pydimension=PyInt_FromLong(dimension);

    return pydimension;
    
}

static int LALVariables_init(li_LALVariables *self, PyObject *args, PyObject *kwds)
{
    /* Should fill in the array using a dictionary as input */
    
    self->vars=(LALVariables*)malloc(sizeof(LALVariables));
    memset((void*)self->vars,0,sizeof(LALVariables));
    return 0;
}



static PyMethodDef LALVariables_methods[]= {
    /* {"name", (PyCFunction)function, METH_NOARGS, "DESCR"}, */
    {"addVariable",(PyCFunction)LALVariables_add_variable,METH_VARARGS,""},
    {"getVariable",(PyCFunction)LALVariables_get_variable,METH_VARARGS,""},
    {"checkVariable",(PyCFunction)LALVariables_check_variable,METH_VARARGS,""},
    {"removeVariable",(PyCFunction)LALVariables_remove_variable,METH_VARARGS,""},
    {"getVariableName",(PyCFunction)LALVariables_get_variable_name,METH_VARARGS,""},
    {"getVariableType",(PyCFunction)LALVariables_get_variable_type,METH_VARARGS,""},
    {"getVariableTypeByIndex",(PyCFunction)LALVariables_get_variable_type_by_index,METH_VARARGS,""},
    {"getVariableVaryType",(PyCFunction)LALVariables_get_variable_vary_type,METH_VARARGS,""},
    {"setVariable",(PyCFunction)LALVariables_set_variable,METH_VARARGS,""},
    {"getVariableDimension",(PyCFunction)LALVariables_get_variable_dimension,METH_NOARGS,""},
    {"getVariableDimensionNonFixed",(PyCFunction)LALVariables_get_variable_dimension_non_fixed,METH_NOARGS,""},
    
    {NULL} /* Sentinel */
};

static PyTypeObject li_lalvariables_type = {
    PyObject_HEAD_INIT(NULL)
    0,              /* obj_size - unused (must be 0) */
    "lalinference.BaseLALVariables",    /* tp_name, name of type */
    sizeof(li_LALVariables),  /* tp_basicsize */
    0,              /* tp_itemsize, need to check */
    (destructor)LALVariables_dealloc,  /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,        /*tp_flags*/
    "LALInference BaseLALVariables objects", /* tp_doc */
    0,                     /* tp_traverse */
    0,                     /* tp_clear */
    0,                     /* tp_richcompare */
    0,                     /* tp_weaklistoffset */
    0,                     /* tp_iter */
    0,                     /* tp_iternext */
    LALVariables_methods,             /* tp_methods */
    0,             /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)LALVariables_init,      /* tp_init */
    0,                         /* tp_alloc */
    PyType_GenericNew,                 /* tp_new */
};
/*
 * ============================================================================
 *
 *                            PosVelAcc
 *
 * ============================================================================
 */

/*Methods*/

 /* Destructor for LALIFOData */
static void PosVelAcc_dealloc(li_PosVelAcc *self)
{
    self->ob_type->tp_free((PyObject *)self);
}

static int PosVelAcc__init__(li_PosVelAcc *self, PyObject *args, PyObject *kwds)
{
    self->data=(PosVelAcc*)malloc(sizeof(PosVelAcc));
    memset((void*)self->data,0,sizeof(PosVelAcc));
    return 0;
}

static PyObject* PosVelAcc__new__(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    
    li_PosVelAcc *obj = (li_PosVelAcc*) PyType_GenericNew(type, args, kwds);
    if(!obj)
        return NULL;
    obj->data=(PosVelAcc*)malloc(sizeof(PosVelAcc));
    memset((void*)obj->data,0,sizeof(PosVelAcc));
    
    return (PyObject*)obj;
}

PyObject* getthreeTupleFromData(PyObject* tuple){
    if(tuple==NULL){
        Py_INCREF(Py_None);
        return Py_None;
    }else{
        Py_INCREF(tuple);
        return tuple;
    }
}

int setthreeTupleFromData(REAL8* internal,PyObject* old,PyObject* new){

    if(!PyTuple_Check(new)||PyTuple_Size(new)!=3){
        PyErr_SetObject(PyExc_TypeError, new);
        return -1;
    }
    Py_ssize_t i;
    for(i=0;i<3;i++){
        PyObject* item=PyTuple_GetItem(new,i);
        if(!PyFloat_Check(item)){
            PyErr_SetObject(PyExc_TypeError,item);
            return -1;
        }else{
            internal[i]=PyFloat_AsDouble(item);
        }
    }

    Py_DECREF(old);
    old=new;
    return 0;
}


/*pos*/
static PyObject* PosVelAcc_getpos(li_PosVelAcc *self, void *closure){return getthreeTupleFromData(self->pos);}
static int PosVelAcc_setpos(li_PosVelAcc *self, PyObject *value, void *closure){return setthreeTupleFromData(self->data->pos,self->pos,value);}

/*vel*/
static PyObject* PosVelAcc_getvel(li_PosVelAcc *self, void *closure){return getthreeTupleFromData(self->vel);}
static int PosVelAcc_setvel(li_PosVelAcc *self, PyObject *value, void *closure){return setthreeTupleFromData(self->data->vel,self->vel,value);}

/*acc*/
static PyObject* PosVelAcc_getacc(li_PosVelAcc *self, void *closure){return getthreeTupleFromData(self->acc);}
static int PosVelAcc_setacc(li_PosVelAcc *self, PyObject *value, void *closure){return setthreeTupleFromData(self->data->acc,self->acc,value);}

static PyMethodDef PosVelAcc_methods[]= {
    {NULL} /* Sentinel */
};

static PyGetSetDef PosVelAcc_getseters[] = {
    {"pos",(getter)PosVelAcc_getpos,(setter)PosVelAcc_setpos,"pos",NULL},
    {"vel",(getter)PosVelAcc_getvel,(setter)PosVelAcc_setvel,"vel",NULL},
    {"acc",(getter)PosVelAcc_getacc,(setter)PosVelAcc_setacc,"acc",NULL},
    {NULL}  /* Sentinel */
};

static struct PyMemberDef PosVelAcc_members[] = {
    {"gps", T_DOUBLE, offsetof(li_PosVelAcc, data)+offsetof(PosVelAcc,gps), 0, "gps"},
    {NULL,}
};

static PyTypeObject li_posvelacc_type = {
    PyObject_HEAD_INIT(NULL)
    0,              /* obj_size - unused (must be 0) */
    "lalinference.PosVelAcc",    /* tp_name, name of type */
    sizeof(li_PosVelAcc),  /* tp_basicsize */
    0,              /* tp_itemsize, need to check */
    (destructor)PosVelAcc_dealloc,  /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,        /*tp_flags*/
    "", /* tp_doc */
    0,                     /* tp_traverse */
    0,                     /* tp_clear */
    0,                     /* tp_richcompare */
    0,                     /* tp_weaklistoffset */
    0,                     /* tp_iter */
    0,                     /* tp_iternext */
    PosVelAcc_methods,             /* tp_methods */
    PosVelAcc_members,             /* tp_members */
    PosVelAcc_getseters,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)PosVelAcc__init__,      /* tp_init */
    0,                         /* tp_alloc */
    PosVelAcc__new__,                 /* tp_new */
};


/*
 * ============================================================================
 *
 *                            Module methods
 *
 * ============================================================================
 */
static PyMethodDef module_methods[] = {
    {NULL}  /* Sentinel */
};
/*
 * ============================================================================
 *
 *                            Module Registration
 *
 * ============================================================================
 */
 
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
    pylal_laldetector_import();
    
    _li_EphemerisData_Type = &li_ephemerisdata_type;
    if (PyType_Ready(&li_EphemerisData_Type) < 0)
        return;

    _li_BarycenterInput_Type = &li_barycenterinput_type;
    if (PyType_Ready(&li_BarycenterInput_Type) < 0)
        return;

    _li_LALVariables_Type = &li_lalvariables_type;
    if (PyType_Ready(&li_LALVariables_Type) < 0)
        return;  

    _li_PosVelAcc_Type = &li_posvelacc_type;
    if (PyType_Ready(&li_PosVelAcc_Type) < 0)
        return;
    
    _li_LALIFOData_Type = &li_lalifodata_type;
    if (PyType_Ready(&li_LALIFOData_Type) < 0)
        return;
    _li_LALInferenceRunState_Type = &li_lalinferencerunstate_type;
    if (PyType_Ready(&li_LALInferenceRunState_Type) < 0)
        return;

    m = Py_InitModule3(MODULE_NAME,module_methods,LIDocString);
    
        
    Py_INCREF(&li_EphemerisData_Type);
    PyModule_AddObject(m, "EphemerisData", (PyObject *)&li_EphemerisData_Type);
    
    Py_INCREF(&li_BarycenterInput_Type);
    PyModule_AddObject(m, "BarycenterInput", (PyObject *)&li_BarycenterInput_Type);
    
    Py_INCREF(&li_LALVariables_Type);
    PyModule_AddObject(m, "BaseLALVariables", (PyObject *)&li_LALVariables_Type);

    Py_INCREF(&li_PosVelAcc_Type);
    PyModule_AddObject(m, "PosVelAcc", (PyObject *)&li_PosVelAcc_Type);

    Py_INCREF(&li_LALIFOData_Type);
    PyModule_AddObject(m, "BaseLALIFOData", (PyObject *)&li_LALIFOData_Type);

    Py_INCREF(&li_LALIFOData_Type);
    PyModule_AddObject(m, "LALInferenceRunState", (PyObject *)&li_LALInferenceRunState_Type);
}
