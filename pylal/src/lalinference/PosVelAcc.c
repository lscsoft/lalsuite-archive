//PosVelAcc.c
#include <Python.h>
#include <structmember.h>

#include <numpy/arrayobject.h>

#include <lal/LALInference.h>
#include <lal/LALBarycenter.h>

#include <tools.h>

#include "PosVelAcc.h"

/*
 * ============================================================================
 *
 *                            PosVelAcc
 *
 * ============================================================================
 */

#define MODULE_NAME LALINFERENCE_POSVELACC_MODULE_NAME

const char EDDocString[] =
"This module provides data types and function wrappers for"
"PosVelAcc.";

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

    if(!PyTuple_Check(new)|PyTuple_Size(new)!=3){
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
    {"vel",(getter)PosVelAcc_getvelsetter)PosVelAcc_setvel,"vel",NULL},
    {"acc",(getter)PosVelAcc_getaccsetter)PosVelAcc_setacc,"acc",NULL},
    {NULL}  /* Sentinel */
};

static struct PyMemberDef PosVelAcc_members[] = {
    {"gps", T_DOUBLE, offsetof(li_PosVelAcc, data)+offsetof(PosVelAcc,gps), 0, "gps"},
    {NULL,}
};

static PyTypeObject li_posvelacc_type = {
    PyObject_HEAD_INIT(NULL)
    0,              /* obj_size - unused (must be 0) */
    "lalinference.PosVelAcc.PosVelAcc",    /* tp_name, name of type */
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

PyMODINIT_FUNC
init_posvelacc(void)
{
    PyObject *m;

    _li_PosVelAcc_Type = &li_posvelacc_type;
    if (PyType_Ready(&li_PosVelAcc_Type) < 0)
        return;
    
    m = Py_InitModule3(MODULE_NAME,module_methods,EDDocString);

    import_array();
    
    Py_INCREF(&li_PosVelAcc_Type);
    PyModule_AddObject(m, "PosVelAcc", (PyObject *)&li_PosVelAcc_Type);
}
