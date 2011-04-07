//EphemerisData.c
#include <Python.h>
#include <structmember.h>

#include <numpy/arrayobject.h>

#include <lal/LALInference.h>
#include <lal/LALBarycenter.h>

#include <tools.h>


#include "EphemerisData.h"

/*
 * ============================================================================
 *
 *                            EphemerisData
 *
 * ============================================================================
 */

#define MODULE_NAME LALINFERENCE_EPHEMERISDATA_MODULE_NAME

const char EDDocString[] =
"This module provides data types and function wrappers for"
"EphemerisData.";

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

PyObject* getPosVelAccFromData(PosVelAcc* internal,PyObject* owner){
    PyObject *empty_tuple = PyTuple_New(0);
    pylal_PosVelAcc *obj = (pylal_PosVelAcc *) PyType_GenericNew(&pylal_PosVelAcc_Type,empty_tuple,NULL);
    Py_DECREF(empty_tuple);

    if(!obj) {
        return NULL;
    }
    
    memcpy((void*)&obj->detector,(void*)&internal,sizeof(PosVelAcc));
    return (PyObject *) obj;
}
/*site*/

static PyObject* EphemerisData_getephemE(li_EphemerisData *self, void *closure){return getLALDetectorFromData(self->data->site);}
static int EphemerisData_setephemE(li_EphemerisData *self, PyObject *value, void *closure){return setLALDetectorFromData(&self->data->site,self->site,value);}


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
    "lalinference.Ephemerisdata.EphemerisData",    /* tp_name, name of type */
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

PyMODINIT_FUNC
init_lalinference(void)
{
    PyObject *m;

    _li_EphemerisData_Type = &li_ephemerisdata_type;
    //li_LALVariables_Type.tp_new = PyType_GenericNew;
    if (PyType_Ready(&li_EphemerisData_Type) < 0)
        return;
    
    m = Py_InitModule3(MODULE_NAME,module_methods,EDDocString);

    import_array();
    PosVelAcc_import();
    
    Py_INCREF(&li_EphemerisData_Type);
    PyModule_AddObject(m, "EphemerisData", (PyObject *)&li_EphemerisData_Type);
}
