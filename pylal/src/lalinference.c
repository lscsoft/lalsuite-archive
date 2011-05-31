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
#include <laldetector.h>

#include "lalinference/LALIFOData.h"
#include "lalinference/BarycenterInput.h"
#include "lalinference/EphemerisData.h"
#include "lalinference/LALVariables.h"
#include "lalinference/PosVelAcc.h"
#include "lalinference/LALInferenceRunState.h"
#include "lalinference/LALAlgorithm.h"
#include "lalinference/LALPriorFunction.h"
#include "lalinference/LALTemplateFunction.h"
#include "lalinference/LALEvolveOneStepFunction.h"
#include "lalinference/LALLikelihoodFunction.h"
#include "lalinference/LALProposalFunction.h"


#define MODULE_NAME "pylal._lalinference"

const char LIDocString[] =
"This module provides data types and function wrappers for"
"LALInference.";

int setLALInferenceIFODataFromData(LALInferenceIFOData* internal,PyObject* old,PyObject* new){
    /* require li_LALIFOData*/
    if (new == NULL) {
        PyErr_SetString(PyExc_TypeError, "Cannot delete the REAL8Window attribute");
        return -1;
    }
    if(!PyObject_TypeCheck(new, &li_LALInferenceIFOData_Type)){
        PyErr_SetObject(PyExc_TypeError, new);
        return -1;
    }

    if(old) Py_DECREF(old);

    Py_INCREF(new);
    old=new;

    li_LALInferenceIFOData* newdata=(li_LALInferenceIFOData*)new;
    internal=newdata->data;
    return 0;
}

PyObject* getLALInferenceIFODataFromData(LALInferenceIFOData* internal,PyObject* owner){
    if(!owner){
        Py_INCREF(Py_None);
        return Py_None;
    }else{
        
        return li_LALInferenceIFOData_new(internal,owner);
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

PyObject* getLALDetectorFromData(LALDetector* internal,PyObject* owner){
    if(!owner){
        Py_INCREF(Py_None);
        return Py_None;
    }else{
        
        return pylal_LALDetector_new(internal);
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
static PyObject* BarycenterInput_getsite(li_BarycenterInput *self, void *closure){return getLALDetectorFromData(&self->data->site,self->site);}
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
 *                            LALAlgorithm
 *
 * ============================================================================
 */

/*Methods*/

 /* Destructor for LALIFOData */
static void LALInferenceAlgorithm_dealloc(li_LALInferenceAlgorithm *self)
{
    self->ob_type->tp_free((PyObject *)self);
}

static int LALInferenceAlgorithm__init__(li_LALInferenceAlgorithm *self, PyObject *args, PyObject *kwds)
{    
    return 0;
}

static PyObject* LALInferenceAlgorithm__new__(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    
    li_LALInferenceAlgorithm *obj = (li_LALInferenceAlgorithm*) PyType_GenericNew(type, args, kwds);
    if(!obj)
        return NULL;
    obj->func=NULL;
    return (PyObject*)obj;
}

static PyObject* LALInferenceAlgorithm__call__(PyObject *self, PyObject *args, PyObject *kwds){
    return NULL;
}

static PyMethodDef LALInferenceAlgorithm_methods[]= {
    {NULL} /* Sentinel */
};

static PyGetSetDef LALInferenceAlgorithm_getseters[] = {
    {NULL}  /* Sentinel */
};

static struct PyMemberDef LALInferenceAlgorithm_members[] = {
    {NULL,}
};

static PyTypeObject li_lalalgorithm_type = {
    PyObject_HEAD_INIT(NULL)
    0,              /* obj_size - unused (must be 0) */
    "lalinference.LALAlgorithm",    /* tp_name, name of type */
    sizeof(li_LALInferenceAlgorithm),  /* tp_basicsize */
    0,              /* tp_itemsize, need to check */
    (destructor)LALInferenceAlgorithm_dealloc,  /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    LALInferenceAlgorithm__call__,                         /*tp_call*/
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
    LALInferenceAlgorithm_methods,             /* tp_methods */
    LALInferenceAlgorithm_members,             /* tp_members */
    LALInferenceAlgorithm_getseters,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)LALInferenceAlgorithm__init__,      /* tp_init */
    0,                         /* tp_alloc */
    LALInferenceAlgorithm__new__,                 /* tp_new */
};

/*
 * ============================================================================
 *
 *                            LALTemplateFunction
 *
 * ============================================================================
 */

/*Methods*/

 /* Destructor for LALInferenceIFOData */
static void LALInferenceTemplateFunction_dealloc(li_LALInferenceTemplateFunction *self)
{
    self->ob_type->tp_free((PyObject *)self);
}

static int LALInferenceTemplateFunction__init__(li_LALInferenceTemplateFunction *self, PyObject *args, PyObject *kwds)
{    
    return 0;
}

static PyObject* LALInferenceTemplateFunction__new__(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    
    li_LALInferenceTemplateFunction *obj = (li_LALInferenceTemplateFunction*) PyType_GenericNew(type, args, kwds);
    if(!obj)
        return NULL;
    obj->func=NULL;
    return (PyObject*)obj;
}

static PyObject* LALInferenceTemplateFunction__call__(PyObject *self, PyObject *args, PyObject *kwds){
    return NULL;
}

static PyMethodDef LALInferenceTemplateFunction_methods[]= {
    {NULL} /* Sentinel */
};

static PyGetSetDef LALInferenceTemplateFunction_getseters[] = {
    {NULL}  /* Sentinel */
};

static struct PyMemberDef LALInferenceTemplateFunction_members[] = {
    {NULL,}
};

static PyTypeObject li_laltemplatefunction_type = {
    PyObject_HEAD_INIT(NULL)
    0,              /* obj_size - unused (must be 0) */
    "lalinference.LALInferenceTemplateFunction",    /* tp_name, name of type */
    sizeof(li_LALInferenceTemplateFunction),  /* tp_basicsize */
    0,              /* tp_itemsize, need to check */
    (destructor)LALInferenceTemplateFunction_dealloc,  /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    LALInferenceTemplateFunction__call__,                         /*tp_call*/
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
    LALInferenceTemplateFunction_methods,             /* tp_methods */
    LALInferenceTemplateFunction_members,             /* tp_members */
    LALInferenceTemplateFunction_getseters,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)LALInferenceTemplateFunction__init__,      /* tp_init */
    0,                         /* tp_alloc */
    LALInferenceTemplateFunction__new__,                 /* tp_new */
};
/*
 * ============================================================================
 *
 *                            LALInferenceProposalFunction
 *
 * ============================================================================
 */

/*Methods*/

 /* Destructor for LALInferenceIFOData */
static void LALInferenceProposalFunction_dealloc(li_LALInferenceProposalFunction *self)
{
    self->ob_type->tp_free((PyObject *)self);
}

static int LALInferenceProposalFunction__init__(li_LALInferenceProposalFunction *self, PyObject *args, PyObject *kwds)
{    
    return 0;
}

static PyObject* LALInferenceProposalFunction__new__(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    
    li_LALInferenceProposalFunction *obj = (li_LALInferenceProposalFunction*) PyType_GenericNew(type, args, kwds);
    if(!obj)
        return NULL;
    obj->func=NULL;
    return (PyObject*)obj;
}

static PyObject* LALInferenceProposalFunction__call__(PyObject *self, PyObject *args, PyObject *kwds){
    return NULL;
}

static PyMethodDef LALInferenceProposalFunction_methods[]= {
    {NULL} /* Sentinel */
};

static PyGetSetDef LALInferenceProposalFunction_getseters[] = {
    {NULL}  /* Sentinel */
};

static struct PyMemberDef LALInferenceProposalFunction_members[] = {
    {NULL,}
};

static PyTypeObject li_lalproposalfunction_type = {
    PyObject_HEAD_INIT(NULL)
    0,              /* obj_size - unused (must be 0) */
    "lalinference.LALInferenceProposalFunction",    /* tp_name, name of type */
    sizeof(li_LALInferenceProposalFunction),  /* tp_basicsize */
    0,              /* tp_itemsize, need to check */
    (destructor)LALInferenceProposalFunction_dealloc,  /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    LALInferenceProposalFunction__call__,                         /*tp_call*/
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
    LALInferenceProposalFunction_methods,             /* tp_methods */
    LALInferenceProposalFunction_members,             /* tp_members */
    LALInferenceProposalFunction_getseters,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)LALInferenceProposalFunction__init__,      /* tp_init */
    0,                         /* tp_alloc */
    LALInferenceProposalFunction__new__,                 /* tp_new */
};

/*
 * ============================================================================
 *
 *                            LALInferenceEvolveOneStepFunction
 *
 * ============================================================================
 */

/*Methods*/

 /* Destructor for LALInferenceIFOData */
static void LALInferenceEvolveOneStepFunction_dealloc(li_LALInferenceEvolveOneStepFunction *self)
{
    self->ob_type->tp_free((PyObject *)self);
}

static int LALInferenceEvolveOneStepFunction__init__(li_LALInferenceEvolveOneStepFunction *self, PyObject *args, PyObject *kwds)
{    
    return 0;
}

static PyObject* LALInferenceEvolveOneStepFunction__new__(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    
    li_LALInferenceEvolveOneStepFunction *obj = (li_LALInferenceEvolveOneStepFunction*) PyType_GenericNew(type, args, kwds);
    if(!obj)
        return NULL;
    obj->func=NULL;
    return (PyObject*)obj;
}

static PyObject* LALInferenceEvolveOneStepFunction__call__(PyObject *self, PyObject *args, PyObject *kwds){
    return NULL;
}

static PyMethodDef LALInferenceEvolveOneStepFunction_methods[]= {
    {NULL} /* Sentinel */
};

static PyGetSetDef LALInferenceEvolveOneStepFunction_getseters[] = {
    {NULL}  /* Sentinel */
};

static struct PyMemberDef LALInferenceEvolveOneStepFunction_members[] = {
    {NULL,}
};

static PyTypeObject li_lalevolveonestepfunction_type = {
    PyObject_HEAD_INIT(NULL)
    0,              /* obj_size - unused (must be 0) */
    "lalinference.LALEvolveOneStepFunction",    /* tp_name, name of type */
    sizeof(li_LALInferenceEvolveOneStepFunction),  /* tp_basicsize */
    0,              /* tp_itemsize, need to check */
    (destructor)LALInferenceEvolveOneStepFunction_dealloc,  /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    LALInferenceEvolveOneStepFunction__call__,                         /*tp_call*/
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
    LALInferenceEvolveOneStepFunction_methods,             /* tp_methods */
    LALInferenceEvolveOneStepFunction_members,             /* tp_members */
    LALInferenceEvolveOneStepFunction_getseters,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)LALInferenceEvolveOneStepFunction__init__,      /* tp_init */
    0,                         /* tp_alloc */
    LALInferenceEvolveOneStepFunction__new__,                 /* tp_new */
};

/*
 * ============================================================================
 *
 *                            LALInferencePriorFunction
 *
 * ============================================================================
 */

/*Methods*/

 /* Destructor for LALInferenceIFOData */
static void LALInferencePriorFunction_dealloc(li_LALInferencePriorFunction *self)
{
    self->ob_type->tp_free((PyObject *)self);
}

static int LALInferencePriorFunction__init__(li_LALInferencePriorFunction *self, PyObject *args, PyObject *kwds)
{    
    return 0;
}

static PyObject* LALInferencePriorFunction__new__(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    
    li_LALInferencePriorFunction *obj = (li_LALInferencePriorFunction*) PyType_GenericNew(type, args, kwds);
    if(!obj)
        return NULL;
    obj->func=NULL;
    return (PyObject*)obj;
}

static PyObject* LALInferencePriorFunction__call__(PyObject *self, PyObject *args, PyObject *kwds){
    return NULL;
}

static PyMethodDef LALInferencePriorFunction_methods[]= {
    {NULL} /* Sentinel */
};

static PyGetSetDef LALInferencePriorFunction_getseters[] = {
    {NULL}  /* Sentinel */
};

static struct PyMemberDef LALInferencePriorFunction_members[] = {
    {NULL,}
};

static PyTypeObject li_lalpriorfunction_type = {
    PyObject_HEAD_INIT(NULL)
    0,              /* obj_size - unused (must be 0) */
    "lalinference.LALInferencePriorFunction",    /* tp_name, name of type */
    sizeof(li_LALInferencePriorFunction),  /* tp_basicsize */
    0,              /* tp_itemsize, need to check */
    (destructor)LALInferencePriorFunction_dealloc,  /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    LALInferencePriorFunction__call__,                         /*tp_call*/
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
    LALInferencePriorFunction_methods,             /* tp_methods */
    LALInferencePriorFunction_members,             /* tp_members */
    LALInferencePriorFunction_getseters,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)LALInferencePriorFunction__init__,      /* tp_init */
    0,                         /* tp_alloc */
    LALInferencePriorFunction__new__,                 /* tp_new */
};
/*
 * ============================================================================
 *
 *                            LALInferenceLikelihoodFunction
 *
 * ============================================================================
 */

/*Methods*/

 /* Destructor for LALInferenceIFOData */
static void LALInferenceLikelihoodFunction_dealloc(li_LALInferenceLikelihoodFunction *self)
{
    self->ob_type->tp_free((PyObject *)self);
}

static int LALInferenceLikelihoodFunction__init__(li_LALInferenceLikelihoodFunction *self, PyObject *args, PyObject *kwds)
{    
    return 0;
}

static PyObject* LALInferenceLikelihoodFunction__new__(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    
    li_LALInferenceLikelihoodFunction *obj = (li_LALInferenceLikelihoodFunction*) PyType_GenericNew(type, args, kwds);
    if(!obj)
        return NULL;
    obj->func=NULL;
    return (PyObject*)obj;
}

static PyObject* LALInferenceLikelihoodFunction__call__(PyObject *self, PyObject *args, PyObject *kwds){
    return NULL;
}

static PyMethodDef LALInferenceLikelihoodFunction_methods[]= {
    {NULL} /* Sentinel */
};

static PyGetSetDef LALInferenceLikelihoodFunction_getseters[] = {
    {NULL}  /* Sentinel */
};

static struct PyMemberDef LALInferenceLikelihoodFunction_members[] = {
    {NULL,}
};

static PyTypeObject li_lallikelihoodfunction_type = {
    PyObject_HEAD_INIT(NULL)
    0,              /* obj_size - unused (must be 0) */
    "lalinference.LALInferenceLikelihoodFunction",    /* tp_name, name of type */
    sizeof(li_LALInferenceLikelihoodFunction),  /* tp_basicsize */
    0,              /* tp_itemsize, need to check */
    (destructor)LALInferenceLikelihoodFunction_dealloc,  /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    LALInferenceLikelihoodFunction__call__,                         /*tp_call*/
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
    LALInferenceLikelihoodFunction_methods,             /* tp_methods */
    LALInferenceLikelihoodFunction_members,             /* tp_members */
    LALInferenceLikelihoodFunction_getseters,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)LALInferenceLikelihoodFunction__init__,      /* tp_init */
    0,                         /* tp_alloc */
    LALInferenceLikelihoodFunction__new__,                 /* tp_new */
};

/*
 * ============================================================================
 *
 *                            LALInferenceIFOData
 *
 * ============================================================================
 */

/*Methods*/

 /* Destructor for LALInferenceIFOData */
static void LALInferenceIFOData_dealloc(li_LALInferenceIFOData *self)
{
    self->ob_type->tp_free((PyObject *)self);
}

static int LALInferenceIFOData__init__(li_LALInferenceIFOData *self, PyObject *args, PyObject *kwds)
{
    self->data=(LALInferenceIFOData*)malloc(sizeof(LALInferenceIFOData));
    memset((void*)self->data,0,sizeof(LALInferenceIFOData));
    return 0;
}

static PyObject* LALInferenceIFOData__new__(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    
    li_LALInferenceIFOData *obj = (li_LALInferenceIFOData*) PyType_GenericNew(type, args, kwds);
    if(!obj)
        return NULL;
    obj->data=(LALInferenceIFOData*)malloc(sizeof(LALInferenceIFOData));
    memset((void*)obj->data,0,sizeof(LALInferenceIFOData));
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

static PyObject* LALInferenceIFOData_getname(li_LALInferenceIFOData *self, void *closure)
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

static int LALInferenceIFOData_setname(li_LALInferenceIFOData *self, PyObject *value, void *closure)
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


PyObject* getLALInferenceVariablesFromData(LALInferenceVariables* internallv,PyObject* owner){
    if(owner){
        return LALInferenceVariables_new(internallv,owner);
    }
    else {
        Py_INCREF(Py_None);
        return Py_None;
    }
}

int setLALInferenceVariablesFromData(LALInferenceVariables* internallv,PyObject* oldlv,PyObject* newlv){
    if (!newlv) {
        PyErr_SetString(PyExc_TypeError, "Cannot delete the attribute");
        return -1;
    }
    if (!PyObject_TypeCheck(newlv, &li_LALInferenceVariables_Type)){
        PyErr_SetString(PyExc_TypeError, "The attribute must be a lalinference.LALInferenceVariables object.");
        return -1;
    }

    if(oldlv) Py_DECREF(oldlv);

    Py_INCREF(newlv);
    oldlv=newlv;

    li_LALInferenceVariables* valuelv=(li_LALInferenceVariables*)newlv;
    internallv=valuelv->vars;
    
    return 0;
}

/*modelParams*/
static PyObject* LALInferenceIFOData_getmodelParams(li_LALInferenceIFOData *self, void *closure) {return getLALInferenceVariablesFromData(self->data->modelParams,self->modelParams);};
static int LALInferenceIFOData_setmodelParams(li_LALInferenceIFOData *self, PyObject *value, void *closure){return setLALInferenceVariablesFromData(self->data->modelParams,(PyObject*)self->modelParams,value);}

/*dataParams*/
static PyObject* LALInferenceIFOData_getdataParams(li_LALInferenceIFOData *self, void *closure) {return getLALInferenceVariablesFromData(self->data->dataParams,self->dataParams);};
static int LALInferenceIFOData_setdataParams(li_LALInferenceIFOData *self, PyObject *value, void *closure){return setLALInferenceVariablesFromData(self->data->dataParams,(PyObject*)self->dataParams,value);}


static PyObject* LALInferenceIFOData_getmodelDomain(li_LALInferenceIFOData *self, void *closure)
{

    if(self->data->modelDomain==timeDomain){
        return PyString_FromString("timeDomain");
    }
    else if(self->data->modelDomain==frequencyDomain){
        return PyString_FromString("freqDomain");
    }
    else return NULL;
}

static int LALInferenceIFOData_setmodelDomain(li_LALInferenceIFOData *self, PyObject *value, void *closure)
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

    LALInferenceDomain var;

    if(!strcmp(name,"timeDomain")){
        var=timeDomain;
    }
    else if(!strcmp(name,"frequencyDomain")){
        var=frequencyDomain;
    }
    else{
        XLALPrintError("modelDomain must be one of the strings 'frequencyDomain' or 'timeDomain'");
	XLAL_ERROR_INT("LALInference_IFOData_setmodelDomain",XLAL_EINVAL);
    }

    return 0;
}

/* timeData */
static int LALInferenceIFOData_settimeData(li_LALInferenceIFOData *self, PyObject *value, void *closure) {return setREAL8TimeSeriesFromData(self->data->timeData,(PyObject*)self->timeData,value);}
static PyObject* LALInferenceIFOData_gettimeData(li_LALInferenceIFOData *self, void *closure) {return getREAL8TimeSeriesFromData(self->data->timeData,(PyObject*)self->timeData);}//{return getREAL8TimeSeriesFromData(self->data->timeData);}

/* timeModelhPlus */
static int LALInferenceIFOData_settimeModelhPlus(li_LALInferenceIFOData *self, PyObject *value, void *closure) {return setREAL8TimeSeriesFromData(self->data->timeModelhPlus,(PyObject*)self->timeModelhPlus,value);}
static PyObject* LALInferenceIFOData_gettimeModelhPlus(li_LALInferenceIFOData *self, void *closure) {return getREAL8TimeSeriesFromData(self->data->timeModelhPlus,(PyObject*)self->timeModelhPlus);}//{return getREAL8TimeSeriesFromData(self->data->timeModelhPlus);}

/* timeModelhCross */
static int LALInferenceIFOData_settimeModelhCross(li_LALInferenceIFOData *self, PyObject *value, void *closure) {return setREAL8TimeSeriesFromData(self->data->timeModelhCross,(PyObject*)self->timeModelhCross,value);}
static PyObject* LALInferenceIFOData_gettimeModelhCross(li_LALInferenceIFOData *self, void *closure) {return getREAL8TimeSeriesFromData(self->data->timeModelhCross,(PyObject*)self->timeModelhCross);}//{return getREAL8TimeSeriesFromData(self->data->timeModelhCross);}

/* whiteTimeData */
static int LALInferenceIFOData_setwhiteTimeData(li_LALInferenceIFOData *self, PyObject *value, void *closure) {return setREAL8TimeSeriesFromData(self->data->whiteTimeData,(PyObject*)self->whiteTimeData,value);}
static PyObject* LALInferenceIFOData_getwhiteTimeData(li_LALInferenceIFOData *self, void *closure) {return getREAL8TimeSeriesFromData(self->data->whiteTimeData,(PyObject*)self->whiteTimeData);}//{return getREAL8TimeSeriesFromData(self->data->whiteTimeData);}

/* windowedTimeData */
static int LALInferenceIFOData_setwindowedTimeData(li_LALInferenceIFOData *self, PyObject *value, void *closure) {return setREAL8TimeSeriesFromData(self->data->windowedTimeData,(PyObject*)self->windowedTimeData,value);}
static PyObject* LALInferenceIFOData_getwindowedTimeData(li_LALInferenceIFOData *self, void *closure) {return getREAL8TimeSeriesFromData(self->data->windowedTimeData,(PyObject*)self->windowedTimeData);}//getREAL8TimeSeriesFromData(self->data->windowedTimeData);}

/* timeDomainNoiseWeights */
static int LALInferenceIFOData_settimeDomainNoiseWeights(li_LALInferenceIFOData *self, PyObject *value, void *closure) {return setREAL8TimeSeriesFromData(self->data->timeDomainNoiseWeights,(PyObject*)self->timeDomainNoiseWeights,value);}
static PyObject* LALInferenceIFOData_gettimeDomainNoiseWeights(li_LALInferenceIFOData *self, void *closure) {return getREAL8TimeSeriesFromData(self->data->timeDomainNoiseWeights,(PyObject*)self->timeDomainNoiseWeights);}//{return getREAL8TimeSeriesFromData(self->data->timeDomainNoiseWeights);}

/* oneSidedNoisePowerSpectrum */
static int LALInferenceIFOData_setoneSidedNoisePowerSpectrum(li_LALInferenceIFOData *self, PyObject *value, void *closure) {return setREAL8FrequencySeriesFromData(self->data->oneSidedNoisePowerSpectrum,self->oneSidedNoisePowerSpectrum,value);}
static PyObject* LALInferenceIFOData_getoneSidedNoisePowerSpectrum(li_LALInferenceIFOData *self, void *closure) {return getREAL8FrequencySeriesFromData(self->data->oneSidedNoisePowerSpectrum,(PyObject*)self->oneSidedNoisePowerSpectrum);}

/*freqData*/
static int LALInferenceIFOData_setfreqData(li_LALInferenceIFOData *self, PyObject *value, void *closure) {return setCOMPLEX16FrequencySeriesFromData(self->data->freqData,self->freqData,value);}
static PyObject* LALInferenceIFOData_getfreqData(li_LALInferenceIFOData *self, void *closure) {return getCOMPLEX16FrequencySeriesFromData(self->data->freqData,(PyObject*)self->freqData);}

/*freqModelhPlus*/
static int LALInferenceIFOData_setfreqModelhPlus(li_LALInferenceIFOData *self, PyObject *value, void *closure) {return setCOMPLEX16FrequencySeriesFromData(self->data->freqModelhPlus,self->freqModelhPlus,value);}
static PyObject* LALInferenceIFOData_getfreqModelhPlus(li_LALInferenceIFOData *self, void *closure) {return getCOMPLEX16FrequencySeriesFromData(self->data->freqModelhPlus,(PyObject*)self->freqModelhPlus);}

/*freqModelhCross*/
static int LALInferenceIFOData_setfreqModelhCross(li_LALInferenceIFOData *self, PyObject *value, void *closure) {return setCOMPLEX16FrequencySeriesFromData(self->data->freqModelhCross,self->freqModelhCross,value);}
static PyObject* LALInferenceIFOData_getfreqModelhCross(li_LALInferenceIFOData *self, void *closure){return getCOMPLEX16FrequencySeriesFromData(self->data->freqModelhCross,(PyObject*)self->freqModelhCross);}

/*whiteFreqData*/
static int LALInferenceIFOData_setwhiteFreqData(li_LALInferenceIFOData *self, PyObject *value, void *closure) {return setCOMPLEX16FrequencySeriesFromData(self->data->whiteFreqData,self->whiteFreqData,value);}
static PyObject* LALInferenceIFOData_getwhiteFreqData(li_LALInferenceIFOData *self, void *closure){return getCOMPLEX16FrequencySeriesFromData(self->data->whiteFreqData,(PyObject*)self->whiteFreqData);}

/*compTimeData*/
static int LALInferenceIFOData_setcompTimeData(li_LALInferenceIFOData *self, PyObject *value, void *closure) {return setCOMPLEX16TimeSeriesFromData(self->data->compTimeData,self->compTimeData,value);}
static PyObject* LALInferenceIFOData_getcompTimeData(li_LALInferenceIFOData *self, void *closure) {return getCOMPLEX16TimeSeriesFromData(self->data->compTimeData,self->compTimeData);}

/*compModelData*/
static int LALInferenceIFOData_setcompModelData(li_LALInferenceIFOData *self, PyObject *value, void *closure) {return setCOMPLEX16TimeSeriesFromData(self->data->compModelData,self->compModelData,value);}
static PyObject* LALInferenceIFOData_getcompModelData(li_LALInferenceIFOData *self, void *closure) {return getCOMPLEX16TimeSeriesFromData(self->data->compModelData,(PyObject*)self->compModelData);}

/*timeToFreqFFTPlan*/
static int LALInferenceIFOData_settimeToFreqFFTPlan(li_LALInferenceIFOData *self, PyObject *value, void *closure) {return setREAL8FFTPlanFromData(self->data->timeToFreqFFTPlan,self->timeToFreqFFTPlan,value);}
static PyObject* LALInferenceIFOData_gettimeToFreqFFTPlan(li_LALInferenceIFOData *self, void *closure) {return getREAL8FFTPlanFromData(self->data->timeToFreqFFTPlan,(PyObject*)self->timeToFreqFFTPlan);}

/*freqToTimeFFTPlan*/
static int LALInferenceIFOData_setfreqToTimeFFTPlan(li_LALInferenceIFOData *self, PyObject *value, void *closure) {return setREAL8FFTPlanFromData(self->data->freqToTimeFFTPlan,self->freqToTimeFFTPlan,value);}
static PyObject* LALInferenceIFOData_getfreqToTimeFFTPlan(li_LALInferenceIFOData *self, void *closure) {return getREAL8FFTPlanFromData(self->data->freqToTimeFFTPlan,(PyObject*)self->freqToTimeFFTPlan);}

/* epoch */
static int LALInferenceIFOData_setepoch(li_LALInferenceIFOData *self, PyObject *value, void *closure) {return setLIGOTimeGPSFromData(self->data->epoch,value);}
static PyObject* LALInferenceIFOData_getepoch(li_LALInferenceIFOData *self, void *closure) {return pylal_LIGOTimeGPS_new(self->data->epoch);}

/* window */
static int LALInferenceIFOData_setwindow(li_LALInferenceIFOData *self, PyObject *value, void *closure) {return setREAL8WindowFromData(self->data->window,self->window,value);}
static PyObject* LALInferenceIFOData_getwindow(li_LALInferenceIFOData *self, void *closure) {return getREAL8WindowFromData(self->data->window,(PyObject*)self->window);}

/* next */
static int LALInferenceIFOData_setnext(li_LALInferenceIFOData *self, PyObject *value, void *closure) {return setLALInferenceIFODataFromData(self->data->next,self->next,value);}
static PyObject* LALInferenceIFOData_getnext(li_LALInferenceIFOData *self, void *closure) {return getLALInferenceIFODataFromData(self->data->next,(PyObject*)self->next);}

/* bary */
static int LALInferenceIFOData_setbary(li_LALInferenceIFOData *self, PyObject *value, void *closure) {return setBarycenterInputFromData(self->data->bary,self->bary,value);}
static PyObject* LALInferenceIFOData_getbary(li_LALInferenceIFOData *self, void *closure) {return getBarycenterInputFromData(self->data->bary,(PyObject*)self->bary);}

/* ephem */
static int LALInferenceIFOData_setephem(li_LALInferenceIFOData *self, PyObject *value, void *closure) {return setEphemerisDataFromData(self->data->ephem,self->ephem,value);}
static PyObject* LALInferenceIFOData_getephem(li_LALInferenceIFOData *self, void *closure) {return getEphemerisDataFromData(self->data->ephem,(PyObject*)self->ephem);}

/**getsetters registration struct**/

static PyGetSetDef LALInferenceIFOData_getseters[] = {
    {"name",(getter)LALInferenceIFOData_getname,(setter)LALInferenceIFOData_setname,"name",NULL},
    //LALInferenceVariables
    {"modelParams",(getter)LALInferenceIFOData_getmodelParams,(setter)LALInferenceIFOData_setmodelParams,"modelParams",NULL},
    {"dataParams",(getter)LALInferenceIFOData_getdataParams,(setter)LALInferenceIFOData_setdataParams,"dataParams",NULL},
    //LALInferenceDomain ...represented as string
    {"modelDomain",(getter)LALInferenceIFOData_getmodelDomain,(setter)LALInferenceIFOData_setmodelDomain,"modelDomain",NULL},
    //REAL8TimeSeries
    {"timeData",(getter)LALInferenceIFOData_gettimeData,(setter)LALInferenceIFOData_settimeData,"timeData",NULL},
    {"timeModelhPlus",(getter)LALInferenceIFOData_gettimeModelhPlus,(setter)LALInferenceIFOData_settimeModelhPlus,"timeModelhPlus",NULL},
    {"timeModelhCross",(getter)LALInferenceIFOData_gettimeModelhCross,(setter)LALInferenceIFOData_settimeModelhCross,"timeModelhCross",NULL},
    {"whiteTimeData",(getter)LALInferenceIFOData_getwhiteTimeData,(setter)LALInferenceIFOData_setwhiteTimeData,"whiteTimeData",NULL},
    {"windowedTimeData",(getter)LALInferenceIFOData_getwindowedTimeData,(setter)LALInferenceIFOData_setwindowedTimeData,"windowedTimeData",NULL},
    {"timeDomainNoiseWeights",(getter)LALInferenceIFOData_gettimeDomainNoiseWeights,(setter)LALInferenceIFOData_settimeDomainNoiseWeights,"timeDomainNoiseWeights",NULL},
    //COMPLEX16FrequencySeries
    {"freqData",(getter)LALInferenceIFOData_getfreqData,(setter)LALInferenceIFOData_setfreqData,"freqData",NULL},
    {"freqModelhPlus",(getter)LALInferenceIFOData_getfreqModelhPlus,(setter)LALInferenceIFOData_setfreqModelhPlus,"freqModelhPlus",NULL},
    {"freqModelhCross",(getter)LALInferenceIFOData_getfreqModelhCross,(setter)LALInferenceIFOData_setfreqModelhCross,"freqModelhCross",NULL},
    {"whiteFreqData",(getter)LALInferenceIFOData_getwhiteFreqData,(setter)LALInferenceIFOData_setwhiteFreqData,"whiteFreqData",NULL},
    //COMPLEX16TimeSeries
    {"compTimeData",(getter)LALInferenceIFOData_getcompTimeData,(setter)LALInferenceIFOData_setcompTimeData,"compTimeData",NULL},
    {"compModelData",(getter)LALInferenceIFOData_getcompModelData,(setter)LALInferenceIFOData_setcompModelData,"compModelData",NULL},
    //REAL8FrequencySeries
    {"oneSidedNoisePowerSpectrum",(getter)LALInferenceIFOData_getoneSidedNoisePowerSpectrum,(setter)LALInferenceIFOData_setoneSidedNoisePowerSpectrum,"oneSidedNoisePowerSpectrum",NULL},
    //REAL8FFTPlan
    {"timeToFreqFFTPlan",(getter)LALInferenceIFOData_gettimeToFreqFFTPlan,(setter)LALInferenceIFOData_settimeToFreqFFTPlan,"timeToFreqFFTPlan",NULL},
    {"freqToTimeFFTPlan",(getter)LALInferenceIFOData_getfreqToTimeFFTPlan,(setter)LALInferenceIFOData_setfreqToTimeFFTPlan,"freqToTimeFFTPlan",NULL},
    //LIGOTimeGPS
    {"epoch",(getter)LALInferenceIFOData_getepoch,(setter)LALInferenceIFOData_setepoch,"epoch",NULL},
    //REAL8Window
    {"window",(getter)LALInferenceIFOData_getwindow,(setter)LALInferenceIFOData_setwindow,"window",NULL},
    //LALInferenceIFOData
    {"next",(getter)LALInferenceIFOData_getnext,(setter)LALInferenceIFOData_setnext,"next",NULL},
    //BarycenterInput
    {"bary",(getter)LALInferenceIFOData_getbary,(setter)LALInferenceIFOData_setbary,"bary",NULL},
    //EphemerisData
    {"ephem",(getter)LALInferenceIFOData_getephem,(setter)LALInferenceIFOData_setephem,"ephem",NULL},
    {NULL}  /* Sentinel */
};

static struct PyMemberDef LALInferenceIFOData_members[] = {
    //REAL8's
    {"nullloglikelihood", T_DOUBLE, offsetof(li_LALInferenceIFOData, data)+offsetof(LALInferenceIFOData,nullloglikelihood), 0, "nullloglikelihood"},
    {"loglikelihood", T_DOUBLE, offsetof(li_LALInferenceIFOData, data)+offsetof(LALInferenceIFOData,loglikelihood), 0, "loglikelihood"},
    {"acceptedloglikelihood", T_DOUBLE, offsetof(li_LALInferenceIFOData, data)+offsetof(LALInferenceIFOData,acceptedloglikelihood), 0, "acceptedloglikelihood"},
    {"fPlus", T_DOUBLE, offsetof(li_LALInferenceIFOData, data)+offsetof(LALInferenceIFOData,fPlus), 0, "fPlus"},
    {"fCross", T_DOUBLE, offsetof(li_LALInferenceIFOData, data)+offsetof(LALInferenceIFOData,fCross), 0, "fCross"},
    {"timeshift", T_DOUBLE, offsetof(li_LALInferenceIFOData, data)+offsetof(LALInferenceIFOData,timeshift), 0, "timeshift"},
    {"fLow", T_DOUBLE, offsetof(li_LALInferenceIFOData, data)+offsetof(LALInferenceIFOData,fLow), 0, "fLow"},
    {"fHigh", T_DOUBLE,offsetof(li_LALInferenceIFOData, data)+offsetof(LALInferenceIFOData,fHigh), 0, "fHigh"},
    
    {NULL,}
};

static PyMethodDef LALInferenceIFOData_methods[]= {
    /* {"name", (PyCFunction)function, METH_NOARGS, "DESCR"}, */
    {NULL} /* Sentinel */
};

static PyTypeObject li_lalifodata_type = {
    PyObject_HEAD_INIT(NULL)
    0,              /* obj_size - unused (must be 0) */
    "lalinference.LALInferenceIFOData",    /* tp_name, name of type */
    sizeof(li_LALInferenceIFOData),  /* tp_basicsize */
    0,              /* tp_itemsize, need to check */
    (destructor)LALInferenceIFOData_dealloc,  /*tp_dealloc*/
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
    "LALInference LALInferenceIFOData object.", /* tp_doc */
    0,                     /* tp_traverse */
    0,                     /* tp_clear */
    0,                     /* tp_richcompare */
    0,                     /* tp_weaklistoffset */
    0,                     /* tp_iter */
    0,                     /* tp_iternext */
    LALInferenceIFOData_methods,             /* tp_methods */
    LALInferenceIFOData_members,             /* tp_members */
    LALInferenceIFOData_getseters,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)LALInferenceIFOData__init__,      /* tp_init */
    0,                         /* tp_alloc */
    LALInferenceIFOData__new__,                 /* tp_new */
    
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
    self->state=(LALInferenceRunState*)malloc(sizeof(LALInferenceRunState));
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
    //struct tagLALInferenceIFOData*
    obj->data=NULL;
    //LALInferenceVariables*
    obj->currentParams=NULL;
    obj->priorArgs=NULL;
    obj->proposalArgs=NULL;
    obj->algorithmParams=NULL; /* Parameters which control the running of the algorithm */
    //LALInferenceVariables**
    obj->livePoints=NULL; /* Array of live points for Nested Sampling */
    obj->differentialPoints=NULL;
    //gsl_rng*
    obj->GSLrandom=NULL; 
    return (PyObject*)obj;
}

/***********getsetters******************/

int setLALInferenceAlgorithmFromData(LALInferenceAlgorithm* internal,PyObject* old,PyObject* new){
    /* require pylal_LALInferenceAlgorithm*/
    if (new == NULL) {
        PyErr_SetString(PyExc_TypeError, "Cannot delete the LALInferenceAlgorithm attribute");
        return -1;
    }
    if(!PyObject_TypeCheck(new, &li_LALInferenceAlgorithm_Type)){
        PyErr_SetObject(PyExc_TypeError, new);
        return -1;
    }

    if(old) Py_DECREF(old);

    Py_INCREF(new);
    old=new;

    li_LALInferenceAlgorithm* newalgo=(li_LALInferenceAlgorithm*)new;
    internal=newalgo->func;
    return 0;
}

PyObject* getLALInferenceAlgorithmFromData(LALInferenceAlgorithm* internal,PyObject* owner){
    if(!owner){
        Py_INCREF(Py_None);
        return Py_None;
    }else{
        
        return li_LALInferenceAlgorithm_new(internal,owner);
    }
}

int setLALInferenceEvolveOneStepFunctionFromData(LALInferenceEvolveOneStepFunction* internal,PyObject* old,PyObject* new){
    /* require pylal_LALInferenceEvolveOneStepFunction*/
    if (new == NULL) {
        PyErr_SetString(PyExc_TypeError, "Cannot delete the LALInferenceEvolveOneStepFunction attribute");
        return -1;
    }
    if(!PyObject_TypeCheck(new, &li_LALInferenceEvolveOneStepFunction_Type)){
        PyErr_SetObject(PyExc_TypeError, new);
        return -1;
    }

    if(old) Py_DECREF(old);

    Py_INCREF(new);
    old=new;

    li_LALInferenceEvolveOneStepFunction* newosf=(li_LALInferenceEvolveOneStepFunction*)new;
    internal=newosf->func;
    return 0;
}

PyObject* getLALInferenceEvolveOneStepFunctionFromData(LALInferenceEvolveOneStepFunction* internal,PyObject* owner){
    if(!owner){
        Py_INCREF(Py_None);
        return Py_None;
    }else{
        
        return li_LALInferenceEvolveOneStepFunction_new(internal,owner);
    }
}

int setLALInferencePriorFunctionFromData(LALInferencePriorFunction* internal,PyObject* old,PyObject* new){
    /* require pylal_LALInferencePriorFunction*/
    if (new == NULL) {
        PyErr_SetString(PyExc_TypeError, "Cannot delete the LALInferencePriorFunction attribute");
        return -1;
    }
    if(!PyObject_TypeCheck(new, &li_LALInferencePriorFunction_Type)){
        PyErr_SetObject(PyExc_TypeError, new);
        return -1;
    }

    if(old) Py_DECREF(old);

    Py_INCREF(new);
    old=new;

    li_LALInferencePriorFunction* newosf=(li_LALInferencePriorFunction*)new;
    internal=newosf->func;
    return 0;
}

PyObject* getLALInferencePriorFunctionFromData(LALInferencePriorFunction* internal,PyObject* owner){
    if(!owner){
        Py_INCREF(Py_None);
        return Py_None;
    }else{
        
        return li_LALInferencePriorFunction_new(internal,owner);
    }
}

int setLALInferenceLikelihoodFunctionFromData(LALInferenceLikelihoodFunction* internal,PyObject* old,PyObject* new){
    /* require pylal_LALInferenceLikelihoodFunction*/
    if (new == NULL) {
        PyErr_SetString(PyExc_TypeError, "Cannot delete the LALInferenceLikelihoodFunction attribute");
        return -1;
    }
    if(!PyObject_TypeCheck(new, &li_LALInferenceLikelihoodFunction_Type)){
        PyErr_SetObject(PyExc_TypeError, new);
        return -1;
    }

    if(old) Py_DECREF(old);

    Py_INCREF(new);
    old=new;

    li_LALInferenceLikelihoodFunction* newosf=(li_LALInferenceLikelihoodFunction*)new;
    internal=newosf->func;
    return 0;
}

PyObject* getLALInferenceLikelihoodFunctionFromData(LALInferenceLikelihoodFunction* internal,PyObject* owner){
    if(!owner){
        Py_INCREF(Py_None);
        return Py_None;
    }else{
        
        return li_LALInferenceLikelihoodFunction_new(internal,owner);
    }
}

int setLALInferenceProposalFunctionFromData(LALInferenceProposalFunction* internal,PyObject* old,PyObject* new){
    /* require pylal_LALInferenceProposalFunction*/
    if (new == NULL) {
        PyErr_SetString(PyExc_TypeError, "Cannot delete the LALInferenceProposalFunction attribute");
        return -1;
    }
    if(!PyObject_TypeCheck(new, &li_LALInferenceProposalFunction_Type)){
        PyErr_SetObject(PyExc_TypeError, new);
        return -1;
    }

    if(old) Py_DECREF(old);

    Py_INCREF(new);
    old=new;

    li_LALInferenceProposalFunction* newosf=(li_LALInferenceProposalFunction*)new;
    internal=newosf->func;
    return 0;
}

PyObject* getLALInferenceProposalFunctionFromData(LALInferenceProposalFunction* internal,PyObject* owner){
    if(!owner){
        Py_INCREF(Py_None);
        return Py_None;
    }else{
        
        return li_LALInferenceProposalFunction_new(internal,owner);
    }
}

int setLALInferenceTemplateFunctionFromData(LALInferenceTemplateFunction* internal,PyObject* old,PyObject* new){
    /* require pylal_LALInferenceTemplateFunction*/
    if (new == NULL) {
        PyErr_SetString(PyExc_TypeError, "Cannot delete the LALInferenceTemplateFunction attribute");
        return -1;
    }
    if(!PyObject_TypeCheck(new, &li_LALInferenceTemplateFunction_Type)){
        PyErr_SetObject(PyExc_TypeError, new);
        return -1;
    }

    if(old) Py_DECREF(old);

    Py_INCREF(new);
    old=new;

    li_LALInferenceTemplateFunction* newosf=(li_LALInferenceTemplateFunction*)new;
    internal=newosf->func;
    return 0;
}

PyObject* getLALInferenceTemplateFunctionFromData(LALInferenceTemplateFunction* internal,PyObject* owner){
    if(!owner){
        Py_INCREF(Py_None);
        return Py_None;
    }else{
        
        return li_LALInferenceTemplateFunction_new(internal,owner);
    }
}

int setLALInferenceVariablesListFromData(LALInferenceVariables** internal,PyObject* old,PyObject* new){
    /* require list of li_LALInferenceVariables*/
    if (new == NULL) {
        PyErr_SetString(PyExc_TypeError, "Cannot delete the LALInferenceTemplateFunction attribute");
        return -1;
    }
    if(!PyList_CheckExact(new)){
        PyErr_SetObject(PyExc_TypeError, new);
        return -1;
    }
    Py_ssize_t i;
    Py_ssize_t lst_n=PyList_Size(new);

    LALInferenceVariables** LVlist=(LALInferenceVariables**)malloc((int)lst_n*sizeof(LALInferenceVariables*));
    
    for(i=0;i<lst_n;i++){
        PyObject* item=PyList_GetItem(new,i);
        if(!PyObject_TypeCheck(new, &li_LALInferenceVariables_Type)){
            PyErr_SetObject(PyExc_TypeError,item);
            Py_ssize_t j;
            for(j=0;j<i;j++){
                item=PyList_GetItem(new,i);
                Py_DECREF(item);
                Py_DECREF(item);
            }
            return -1;
            
        }else{
            
            *(LVlist+(int)i)=((li_LALInferenceVariables*)item)->vars;
        }
    }
    if(old) {
        free(internal);
        Py_DECREF(old);
    }

    Py_INCREF(new);
    old=new;
    
    internal=LVlist;
    return 0;
}

PyObject* getLALInferenceVariablesListFromData(LALInferenceVariables** internal,PyObject* owner){
    if(!owner){
        Py_INCREF(Py_None);
        return Py_None;
    }else{
         
        return NULL;
    }
}

/*algorithm*/
static PyObject* LALInferenceRunState_getalgorithm(li_LALInferenceRunState *self, void *closure) {return getLALInferenceAlgorithmFromData(self->state->algorithm,self->algorithm);};
static int LALInferenceRunState_setalgorithm(li_LALInferenceRunState *self, PyObject *value, void *closure){return setLALInferenceAlgorithmFromData(self->state->algorithm,(PyObject*)self->algorithm,value);}
/*evolve*/
static PyObject* LALInferenceRunState_getevolve(li_LALInferenceRunState *self, void *closure) {return getLALInferenceEvolveOneStepFunctionFromData(self->state->evolve,self->evolve);};
static int LALInferenceRunState_setevolve(li_LALInferenceRunState *self, PyObject *value, void *closure){return setLALInferenceEvolveOneStepFunctionFromData(self->state->evolve,(PyObject*)self->evolve,value);}
/*prior*/
static PyObject* LALInferenceRunState_getprior(li_LALInferenceRunState *self, void *closure) {return getLALInferencePriorFunctionFromData(self->state->prior,self->prior);};
static int LALInferenceRunState_setprior(li_LALInferenceRunState *self, PyObject *value, void *closure){return setLALInferencePriorFunctionFromData(self->state->prior,(PyObject*)self->prior,value);}
/*likelihood*/
static PyObject* LALInferenceRunState_getlikelihood(li_LALInferenceRunState *self, void *closure) {return getLALInferenceLikelihoodFunctionFromData(self->state->likelihood,self->likelihood);};
static int LALInferenceRunState_setlikelihood(li_LALInferenceRunState *self, PyObject *value, void *closure){return setLALInferenceLikelihoodFunctionFromData(self->state->likelihood,(PyObject*)self->likelihood,value);}
/*proposal*/
static PyObject* LALInferenceRunState_getproposal(li_LALInferenceRunState *self, void *closure) {return getLALInferenceProposalFunctionFromData(self->state->proposal,self->proposal);};
static int LALInferenceRunState_setproposal(li_LALInferenceRunState *self, PyObject *value, void *closure){return setLALInferenceProposalFunctionFromData(self->state->proposal,(PyObject*)self->proposal,value);}
/*template*/
static PyObject* LALInferenceRunState_gettemplate(li_LALInferenceRunState *self, void *closure) {return getLALInferenceTemplateFunctionFromData(self->state->template,self->template);};
static int LALInferenceRunState_settemplate(li_LALInferenceRunState *self, PyObject *value, void *closure){return setLALInferenceTemplateFunctionFromData(self->state->template,(PyObject*)self->template,value);}
/*currentParams*/
static PyObject* LALInferenceRunState_getcurrentParams(li_LALInferenceRunState *self, void *closure) {return getLALInferenceVariablesFromData(self->state->currentParams,self->currentParams);};
static int LALInferenceRunState_setcurrentParams(li_LALInferenceRunState *self, PyObject *value, void *closure){return setLALInferenceVariablesFromData(self->state->currentParams,(PyObject*)self->currentParams,value);}
/*priorArgs*/
static PyObject* LALInferenceRunState_getpriorArgs(li_LALInferenceRunState *self, void *closure) {return getLALInferenceVariablesFromData(self->state->priorArgs,self->priorArgs);};
static int LALInferenceRunState_setpriorArgs(li_LALInferenceRunState *self, PyObject *value, void *closure){return setLALInferenceVariablesFromData(self->state->priorArgs,(PyObject*)self->priorArgs,value);}
/*proposalArgs*/
static PyObject* LALInferenceRunState_getproposalArgs(li_LALInferenceRunState *self, void *closure) {return getLALInferenceVariablesFromData(self->state->proposalArgs,self->proposalArgs);};
static int LALInferenceRunState_setproposalArgs(li_LALInferenceRunState *self, PyObject *value, void *closure){return setLALInferenceVariablesFromData(self->state->proposalArgs,(PyObject*)self->proposalArgs,value);}
/*algorithmParams*/
static PyObject* LALInferenceRunState_getalgorithmParams(li_LALInferenceRunState *self, void *closure) {return getLALInferenceVariablesFromData(self->state->algorithmParams,self->algorithmParams);};
static int LALInferenceRunState_setalgorithmParams(li_LALInferenceRunState *self, PyObject *value, void *closure){return setLALInferenceVariablesFromData(self->state->algorithmParams,(PyObject*)self->algorithmParams,value);}
/*data*/
static PyObject* LALInferenceRunState_getdata(li_LALInferenceRunState *self, void *closure) {return getLALInferenceIFODataFromData(self->state->data,self->data);};
static int LALInferenceRunState_setdata(li_LALInferenceRunState *self, PyObject *value, void *closure){return setLALInferenceIFODataFromData(self->state->data,(PyObject*)self->data,value);}
/*livePoints*/
static PyObject* LALInferenceRunState_getlivePoints(li_LALInferenceRunState *self, void *closure) {return getLALInferenceVariablesListFromData(self->state->livePoints,self->livePoints);};
static int LALInferenceRunState_setlivePoints(li_LALInferenceRunState *self, PyObject *value, void *closure){return setLALInferenceVariablesListFromData(self->state->livePoints,(PyObject*)self->livePoints,value);}
/*differentialPoints*/
static PyObject* LALInferenceRunState_getdifferentialPoints(li_LALInferenceRunState *self, void *closure) {return getLALInferenceVariablesListFromData(self->state->differentialPoints,self->differentialPoints);};
static int LALInferenceRunState_setdifferentialPoints(li_LALInferenceRunState *self, PyObject *value, void *closure){return setLALInferenceVariablesListFromData(self->state->differentialPoints,(PyObject*)self->differentialPoints,value);}
/**getsetters registration struct**/

static PyGetSetDef LALInferenceRunState_getseters[] = {
    {"algorithm",(getter)LALInferenceRunState_getalgorithm,(setter)LALInferenceRunState_setalgorithm,"algorithm",NULL},
    {"evolve",(getter)LALInferenceRunState_getevolve,(setter)LALInferenceRunState_setevolve,"evolve",NULL},
    {"prior",(getter)LALInferenceRunState_getprior,(setter)LALInferenceRunState_setprior,"prior",NULL},
    {"likelihood",(getter)LALInferenceRunState_getlikelihood,(setter)LALInferenceRunState_setlikelihood,"likelihood",NULL},
    {"proposal",(getter)LALInferenceRunState_getproposal,(setter)LALInferenceRunState_setproposal,"proposal",NULL},
    {"template",(getter)LALInferenceRunState_gettemplate,(setter)LALInferenceRunState_settemplate,"template",NULL},
    //LALVariables*
    {"currentParams",(getter)LALInferenceRunState_getcurrentParams,(setter)LALInferenceRunState_setcurrentParams,"currentParams",NULL},
    {"priorArgs",(getter)LALInferenceRunState_getpriorArgs,(setter)LALInferenceRunState_setpriorArgs,"priorArgs",NULL},
    {"proposalArgs",(getter)LALInferenceRunState_getproposalArgs,(setter)LALInferenceRunState_setproposalArgs,"proposalArgs",NULL},
    {"algorithmParams",(getter)LALInferenceRunState_getalgorithmParams,(setter)LALInferenceRunState_setalgorithmParams,"algorithmParams",NULL},
    //LALInferenceVariables**
    {"livePoints",(getter)LALInferenceRunState_getlivePoints,(setter)LALInferenceRunState_setlivePoints,"livePoints",NULL},
    {"differentialPoints",(getter)LALInferenceRunState_getdifferentialPoints,(setter)LALInferenceRunState_setdifferentialPoints,"differentialPoints",NULL},
    //LALInferenceIFOData*
    {"data",(getter)LALInferenceRunState_getdata,(setter)LALInferenceRunState_setdata,"data",NULL},
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
 *                            LALInferenceVariables
 *
 * ============================================================================
 */

/*Methods*/


/* Destructor for LALInferenceVariables */
static void LALInferenceVariables_dealloc(li_LALInferenceVariables *self)
{
    LALInferenceDestroyVariables(self->vars);
    self->ob_type->tp_free((PyObject *)self);
}

LALInferenceVariableType LALInferenceVariables_infer_litype_from_pyvalue(PyObject* pyvalue){

    LALInferenceVariableType type;
    
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

LALInferenceVariableType LALInferenceVariables_convert_string_to_litype(char* typestring){

    LALInferenceVariableType type;
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
        PyErr_SetString(PyExc_TypeError,"LALInferenceInference type not found!!");
        return -1;
    }
    return type;
}

void* LALInferenceVariables_convert_pyobj_to_livar_value(PyObject* pyvalue,LALInferenceVariableType type){
    void* value=(void *)malloc(LALInferenceTypeSize[type]);
    
    if(type==INT4_t){
        INT4 cast_value=((INT4)PyInt_AsLong(pyvalue));
        INT4* cast_valuep=&cast_value;
        memcpy(value,(void*)cast_valuep,LALInferenceTypeSize[type]);
    }
    else if(type==INT8_t){
        INT8 cast_value=((INT8)PyInt_AsLong(pyvalue));
        INT8* cast_valuep=&cast_value;
        memcpy(value,(void*)cast_valuep,LALInferenceTypeSize[type]);
    }
    else if(type==UINT4_t){
        UINT4 cast_value=(UINT4)((unsigned long int)PyInt_AsLong(pyvalue));
        UINT4* cast_valuep=&cast_value;
        memcpy(value,(void*)cast_valuep,LALInferenceTypeSize[type]);
    }   
    else if(type==REAL4_t){
    
        REAL4 cast_value=((REAL4)PyFloat_AsDouble(pyvalue));
        REAL4* cast_valuep=&cast_value;
        memcpy(value,(void*)cast_valuep,LALInferenceTypeSize[type]);
        
    }
    else if(type==REAL8_t){
        REAL8 cast_value=((REAL8)PyFloat_AsDouble(pyvalue));
        REAL8* cast_valuep=&cast_value;
        memcpy(value,(void*)cast_valuep,LALInferenceTypeSize[type]);
    }
    else{
        PyErr_SetObject(PyExc_TypeError, pyvalue);
        return NULL;
    }
    
    return value;
}

static PyObject* LALInferenceVariables_add_variable(li_LALInferenceVariables *self,PyObject* args,PyObject* kwds){
    PyObject *pyname=NULL,*pyvalue=NULL,*pyvarytype=NULL,*pytype=NULL;

    char *name;char* temp;void* value;LALInferenceVariableType type;LALInferenceParamVaryType varytype;
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
        type=LALInferenceVariables_convert_string_to_litype(typestring);
    }
    //...else infer type from python object type (this is more limited).
    else type=LALInferenceVariables_infer_litype_from_pyvalue(pyvalue);
    
    value=LALInferenceVariables_convert_pyobj_to_livar_value(pyvalue,type);
    
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

    /* If we survived then call LALInferenceAddVariable with self->vars ...*/
    
    LALInferenceAddVariable((self->vars),name,value,type,varytype);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject* LALInferenceVariables_check_variable(li_LALInferenceVariables *self,PyObject* args)
/* Check for existence of name */
{
    PyObject* pyname;
    char* name;
    
    if (! PyArg_ParseTuple(args,"O",&pyname)) return NULL;

    name=PyString_AsString(pyname);

    if(LALInferenceGetItem(self->vars,name)){
        Py_RETURN_TRUE;
    }
    else{
        Py_RETURN_FALSE;
    }
}


PyObject* LALInferenceVariables_convert_livar_value_to_pyobj(li_LALInferenceVariables *pyvars,char* name){
    LALInferenceVariableType type;
    
    PyObject* returnObj=NULL;
    
    type=LALInferenceGetVariableType(pyvars->vars,name);
    void* uncastval=LALInferenceGetVariable(pyvars->vars,name);
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

static int LALInferenceVariables_remove_variable(li_LALInferenceVariables *self,PyObject* args){
    PyObject* pyname;
    char* name;
    
    if (! PyArg_ParseTuple(args,"O",&pyname)) return -1;

    name=PyString_AsString(pyname);
    LALInferenceRemoveVariable(self->vars,name);
    
    return 0;
}

static PyObject* LALInferenceVariables_get_variable(li_LALInferenceVariables *self,PyObject* args){
    PyObject* pyname=NULL,*returnObj=NULL;
    char* name;
    
    if (! PyArg_ParseTuple(args,"O",&pyname)) return NULL;
    name=PyString_AsString(pyname);
    returnObj=LALInferenceVariables_convert_livar_value_to_pyobj(self,name);
    return returnObj;
}

static PyObject* LALInferenceVariables_get_variable_name(li_LALInferenceVariables *self,PyObject* args){
    char* name;
    long var_idx;
    
    if (! PyArg_ParseTuple(args,"i",&var_idx)) return NULL;
    
    name=LALInferenceGetVariableName(self->vars,(int)var_idx);
    
    return PyString_FromString(name);
    
}

static int LALInferenceVariables_set_variable(li_LALInferenceVariables *self,PyObject* args){
    PyObject* pyname;char* name;
    PyObject* pynew_var_value;void* new_var_value;
    LALInferenceVariableType type;

    if (! PyArg_ParseTuple(args,"OO",&pyname,&pynew_var_value)) return -1;

    name=PyString_AsString(pyname);

    type=LALInferenceGetVariableType(self->vars,name);

    new_var_value=LALInferenceVariables_convert_pyobj_to_livar_value(pynew_var_value,type);

    LALInferenceSetVariable(self->vars,name,new_var_value);
    
    return 0;
    
}

char* LALInferenceVariables_get_type_string(LALInferenceVariableType type){
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

static PyObject* LALInferenceVariables_get_variable_type(li_LALInferenceVariables *self,PyObject* args){
    LALInferenceVariableType type;
    char* name;
    PyObject* pyname;
    char* type_name=NULL;
    if (! PyArg_ParseTuple(args,"O",&pyname)) return NULL;
    
    name=PyString_AsString(pyname);
    type=LALInferenceGetVariableType(self->vars,name);
    
    type_name=LALInferenceVariables_get_type_string(type);
    
    if(type_name==NULL){
        Py_INCREF(Py_None);
        return Py_None;
    }    
    
    return PyString_FromString(type_name);
    
}
    
static PyObject* LALInferenceVariables_get_variable_type_by_index(li_LALInferenceVariables *self,PyObject* args){
    LALInferenceVariableType type;
    long var_idx;
    char* type_name=NULL;
    if (! PyArg_ParseTuple(args,"i",&var_idx)) return NULL;
    
    type=LALInferenceGetVariableTypeByIndex(self->vars,(int)var_idx);
    type_name=LALInferenceVariables_get_type_string(type);
    
    return PyString_FromString(type_name);
    
}

char* LALInferenceVariables_get_varytype_string(LALInferenceParamVaryType varytype){
    
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

static PyObject* LALInferenceVariables_get_variable_vary_type(li_LALInferenceVariables *self,PyObject* args){
    LALInferenceParamVaryType varytype;
    char* name;
    PyObject* pyname;
    char* varytype_name=NULL;
    if (! PyArg_ParseTuple(args,"O",&pyname)) return NULL;
    
    name=PyString_AsString(pyname);
    varytype=LALInferenceGetVariableVaryType(self->vars,name);
    varytype_name=LALInferenceVariables_get_varytype_string(varytype);
    
    if(varytype_name==NULL){
        Py_INCREF(Py_None);
        return Py_None;
    }    
    
    return PyString_FromString(varytype_name);
    
}

static PyObject* LALInferenceVariables_get_variable_dimension(li_LALInferenceVariables *self){
    long dimension;
    PyObject* pydimension=NULL;
    
    dimension=(long int)LALInferenceGetVariableDimension(self->vars);

    pydimension=PyInt_FromLong(dimension);

    return pydimension;
    
}

static PyObject* LALInferenceVariables_get_variable_dimension_non_fixed(li_LALInferenceVariables *self){
    long dimension;
    PyObject* pydimension=NULL;
    
    dimension=(long int)LALInferenceGetVariableDimensionNonFixed(self->vars);

    pydimension=PyInt_FromLong(dimension);

    return pydimension;
    
}

static int LALInferenceVariables_init(li_LALInferenceVariables *self, PyObject *args, PyObject *kwds)
{
    /* Should fill in the array using a dictionary as input */
    
    self->vars=(LALInferenceVariables*)malloc(sizeof(LALInferenceVariables));
    memset((void*)self->vars,0,sizeof(LALInferenceVariables));
    return 0;
}



static PyMethodDef LALInferenceVariables_methods[]= {
    /* {"name", (PyCFunction)function, METH_NOARGS, "DESCR"}, */
    {"LALInferenceAddVariable",(PyCFunction)LALInferenceVariables_add_variable,METH_VARARGS,""},
    {"LALInferenceGetVariable",(PyCFunction)LALInferenceVariables_get_variable,METH_VARARGS,""},
    {"LALInferenceCheckVariable",(PyCFunction)LALInferenceVariables_check_variable,METH_VARARGS,""},
    {"LALInferenceRemoveVariable",(PyCFunction)LALInferenceVariables_remove_variable,METH_VARARGS,""},
    {"LALInferenceGetVariableName",(PyCFunction)LALInferenceVariables_get_variable_name,METH_VARARGS,""},
    {"LALInferenceGetVariableType",(PyCFunction)LALInferenceVariables_get_variable_type,METH_VARARGS,""},
    {"LALInferenceGetVariableTypeByIndex",(PyCFunction)LALInferenceVariables_get_variable_type_by_index,METH_VARARGS,""},
    {"LALInferenceGetVariableVaryType",(PyCFunction)LALInferenceVariables_get_variable_vary_type,METH_VARARGS,""},
    {"LALInferenceSetVariable",(PyCFunction)LALInferenceVariables_set_variable,METH_VARARGS,""},
    {"LALInferenceGetVariableDimension",(PyCFunction)LALInferenceVariables_get_variable_dimension,METH_NOARGS,""},
    {"LALInferenceGetVariableDimensionNonFixed",(PyCFunction)LALInferenceVariables_get_variable_dimension_non_fixed,METH_NOARGS,""},
    
    {NULL} /* Sentinel */
};

static PyTypeObject li_lalvariables_type = {
    PyObject_HEAD_INIT(NULL)
    0,              /* obj_size - unused (must be 0) */
    "lalinference.BaseLALInferenceVariables",    /* tp_name, name of type */
    sizeof(li_LALInferenceVariables),  /* tp_basicsize */
    0,              /* tp_itemsize, need to check */
    (destructor)LALInferenceVariables_dealloc,  /*tp_dealloc*/
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
    "LALInference BaseLALInferenceVariables objects", /* tp_doc */
    0,                     /* tp_traverse */
    0,                     /* tp_clear */
    0,                     /* tp_richcompare */
    0,                     /* tp_weaklistoffset */
    0,                     /* tp_iter */
    0,                     /* tp_iternext */
    LALInferenceVariables_methods,             /* tp_methods */
    0,             /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)LALInferenceVariables_init,      /* tp_init */
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

 /* Destructor for LALInferenceIFOData */
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

    _li_LALInferenceVariables_Type = &li_lalvariables_type;
    if (PyType_Ready(&li_LALInferenceVariables_Type) < 0)
        return;  

    _li_PosVelAcc_Type = &li_posvelacc_type;
    if (PyType_Ready(&li_PosVelAcc_Type) < 0)
        return;
    
    _li_LALInferenceIFOData_Type = &li_lalifodata_type;
    if (PyType_Ready(&li_LALInferenceIFOData_Type) < 0)
        return;
        
    _li_LALInferenceRunState_Type = &li_lalinferencerunstate_type;
    if (PyType_Ready(&li_LALInferenceRunState_Type) < 0)
        return;

    _li_LALInferenceAlgorithm_Type = &li_lalalgorithm_type;
    if (PyType_Ready(&li_LALInferenceAlgorithm_Type) < 0)
        return;

    _li_LALInferencePriorFunction_Type = &li_lalpriorfunction_type;
    if (PyType_Ready(&li_LALInferencePriorFunction_Type) < 0)
        return;

    _li_LALInferenceProposalFunction_Type = &li_lalproposalfunction_type;
    if (PyType_Ready(&li_LALInferenceProposalFunction_Type) < 0)
        return;

    _li_LALInferenceTemplateFunction_Type = &li_laltemplatefunction_type;
    if (PyType_Ready(&li_LALInferenceTemplateFunction_Type) < 0)
        return;

    _li_LALInferenceLikelihoodFunction_Type = &li_lallikelihoodfunction_type;
    if (PyType_Ready(&li_LALInferenceLikelihoodFunction_Type) < 0)
        return;

    _li_LALInferenceEvolveOneStepFunction_Type = &li_lalevolveonestepfunction_type;
    if (PyType_Ready(&li_LALInferenceEvolveOneStepFunction_Type) < 0)
        return;

    m = Py_InitModule3(MODULE_NAME,module_methods,LIDocString);
    
        
    Py_INCREF(&li_EphemerisData_Type);
    PyModule_AddObject(m, "EphemerisData", (PyObject *)&li_EphemerisData_Type);
    
    Py_INCREF(&li_BarycenterInput_Type);
    PyModule_AddObject(m, "BarycenterInput", (PyObject *)&li_BarycenterInput_Type);
    
    Py_INCREF(&li_LALInferenceVariables_Type);
    PyModule_AddObject(m, "BaseLALInferenceVariables", (PyObject *)&li_LALInferenceVariables_Type);

    Py_INCREF(&li_PosVelAcc_Type);
    PyModule_AddObject(m, "PosVelAcc", (PyObject *)&li_PosVelAcc_Type);

    Py_INCREF(&li_LALInferenceIFOData_Type);
    PyModule_AddObject(m, "BaseLALInferenceIFOData", (PyObject *)&li_LALInferenceIFOData_Type);

    Py_INCREF(&li_LALInferenceIFOData_Type);
    PyModule_AddObject(m, "LALInferenceRunState", (PyObject *)&li_LALInferenceRunState_Type);

    Py_INCREF(&li_LALInferenceAlgorithm_Type);
    PyModule_AddObject(m, "LALInferenceAlgorithm", (PyObject *)&li_LALInferenceAlgorithm_Type);

    Py_INCREF(&li_LALInferencePriorFunction_Type);
    PyModule_AddObject(m, "LALInferencePriorFunction", (PyObject *)&li_LALInferencePriorFunction_Type);

    Py_INCREF(&li_LALInferenceProposalFunction_Type);
    PyModule_AddObject(m, "LALInferenceProposalFunction", (PyObject *)&li_LALInferenceProposalFunction_Type);

    Py_INCREF(&li_LALInferenceTemplateFunction_Type);
    PyModule_AddObject(m, "LALInferenceTemplateFunction", (PyObject *)&li_LALInferenceTemplateFunction_Type);

    Py_INCREF(&li_LALInferenceLikelihoodFunction_Type);
    PyModule_AddObject(m, "LALInferenceLikelihoodFunction", (PyObject *)&li_LALInferenceLikelihoodFunction_Type);

    Py_INCREF(&li_LALInferenceEvolveOneStepFunction_Type);
    PyModule_AddObject(m, "LALInferenceEvolveOneStepFunction", (PyObject *)&li_LALInferenceEvolveOneStepFunction_Type);
}
