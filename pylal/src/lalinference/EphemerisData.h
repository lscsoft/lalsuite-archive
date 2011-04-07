//EmpherisData.h
#include <lal/LALBarycenter.h>
#define LALINFERENCE_EPHEMERISDATA_MODULE_NAME "pylal.lalinference.ephemerisdata"

static PyTypeObject *_li_EphemerisData_Type = NULL;
#define li_EphemerisData_Type (*_li_EphemerisData_Type)

typedef struct {
    PyObject_HEAD
    PyObject* owner;
    /* Type-specific fields go here */
    EphemerisData* data;
    PyObject* ephemE;
    PyObject* ephemS;
    PyObject* ephiles;
}li_EphemerisData;

static PyObject* li_EphemerisData_new(EphemerisData *bci, PyObject *owner){
    PyObject *empty_tuple = PyTuple_New(0);
    li_EphemerisData *obj = (li_EphemerisData *) PyType_GenericNew(&li_EphemerisData_Type,empty_tuple,NULL);
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

static PyObject *EphemerisData_import(void)
{
    PyObject *name = PyString_FromString(LALINFERENCE_EPHEMERISDATA_MODULE_NAME);
    PyObject *module = PyImport_Import(name);
    Py_DECREF(name);

    name = PyString_FromString("EphemerisData");
    _li_EphemerisData_Type = (PyTypeObject *) PyDict_GetItem(PyModule_GetDict(module), name);
    Py_INCREF(&li_EphemerisData_Type);
    Py_DECREF(name);

    return module;
}
