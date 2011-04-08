//LALIFOData.c
#include <Python.h>
#include <structmember.h>

#include <numpy/arrayobject.h>

#include <lal/LALInference.h>

#include <lal/Sequence.h>
#include <lal/LALBarycenter.h>
#include <lal/LALDetectors.h>

#include <ligotimegps.h>

#include "lalinference.h"

#include "LALIFOData.h"
#include "LALVariables.h"


/*
 * ============================================================================
 *
 *                            LALIFOData
 *
 * ============================================================================
 */

#define MODULE_NAME "pylal.lalinference.lalifodata"

const char LIDDocString[] =
"This module provides data types and function wrappers for"
"LALIFOData.";

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
    "lalinference.lalifodata.LALIFOData",    /* tp_name, name of type */
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
 *                            Module Registration
 *
 * ============================================================================
 */
 
static PyMethodDef module_methods[] = {
    {NULL}  /* Sentinel */
};

PyMODINIT_FUNC
init_lalifodata(void)
{
    PyObject *m;

    _li_LALIFOData_Type = &li_lalifodata_type;
    
    if (PyType_Ready(&li_LALIFOData_Type) < 0)
        return;
    
    m = Py_InitModule3(MODULE_NAME,module_methods,LIDDocString);

    import_array();
    pylal_ligotimegps_import();
    LALVariables_import();
    Py_INCREF(&li_LALIFOData_Type);
    PyModule_AddObject(m, "LALIFOData", (PyObject *)&li_LALIFOData_Type);
}
