//PosVelAcc.h
#include <lal/LALBarycenter.h>
#define LALINFERENCE_POSVELACC_MODULE_NAME "pylal._lalinference"

static PyTypeObject *_li_PosVelAcc_Type = NULL;
#define li_PosVelAcc_Type (*_li_PosVelAcc_Type)

typedef struct {
    PyObject_HEAD
    PyObject* owner;
    /* Type-specific fields go here */
    PosVelAcc* data;
    PyObject* pos;
    PyObject* vel;
    PyObject* acc;
}li_PosVelAcc;

static PyObject* li_PosVelAcc_new(PosVelAcc *bci, PyObject *owner){
    PyObject *empty_tuple = PyTuple_New(0);
    li_PosVelAcc *obj = (li_PosVelAcc *) PyType_GenericNew(&li_PosVelAcc_Type,empty_tuple,NULL);
    Py_DECREF(empty_tuple);
    
    if(!obj) {
        if(!owner)
            return NULL;
    }
    if(owner)
        Py_INCREF(owner);
    obj->owner = owner;
    
    obj->data = bci;
    return (PyObject *) obj;
}

static PyObject *PosVelAcc_import(void)
{
    PyObject *name = PyString_FromString(LALINFERENCE_POSVELACC_MODULE_NAME);
    PyObject *module = PyImport_Import(name);
    Py_DECREF(name);

    name = PyString_FromString("PosVelAcc");
    _li_PosVelAcc_Type = (PyTypeObject *) PyDict_GetItem(PyModule_GetDict(module), name);
    Py_INCREF(&li_PosVelAcc_Type);
    Py_DECREF(name);

    return module;
}
