//BarycenterInput.h
#include <lal/LALBarycenter.h>
#define LALINFERENCE_BARYCENTERINPUT_MODULE_NAME "pylal.lalinference.barycenterinput"

static PyTypeObject *_li_BarycenterInput_Type = NULL;
#define li_BarycenterInput_Type (*_li_BarycenterInput_Type)

typedef struct {
    PyObject_HEAD
    PyObject* owner;
    /* Type-specific fields go here */
    BarycenterInput* data;
    PyObject* tgps;
    PyObject* site;
}li_BarycenterInput;

static PyObject* li_BarycenterInput_new(BarycenterInput *bci, PyObject *owner){
    PyObject *empty_tuple = PyTuple_New(0);
    li_BarycenterInput *obj = (li_BarycenterInput *) PyType_GenericNew(&li_BarycenterInput_Type,empty_tuple,NULL);
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

static PyObject *BarycenterInput_import(void)
{
    PyObject *name = PyString_FromString(LALINFERENCE_BARYCENTERINPUT_MODULE_NAME);
    PyObject *module = PyImport_Import(name);
    Py_DECREF(name);

    name = PyString_FromString("BarycenterInput");
    _li_BarycenterInput_Type = (PyTypeObject *) PyDict_GetItem(PyModule_GetDict(module), name);
    Py_INCREF(&li_BarycenterInput_Type);
    Py_DECREF(name);

    return module;
}
