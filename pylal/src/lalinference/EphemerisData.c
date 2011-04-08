//EphemerisData.c
#include <Python.h>
#include <structmember.h>

#include <numpy/arrayobject.h>

#include <lal/LALInference.h>
#include <lal/LALBarycenter.h>

#include <tools.h>

#include "PosVelAcc.h"
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
init_ephemerisdata(void)
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
