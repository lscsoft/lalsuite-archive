//LALEvolveOneStepFunction.h

static PyTypeObject *_li_LALEvolveOneStepFunction_Type = NULL;
#define li_LALEvolveOneStepFunction_Type (*_li_LALEvolveOneStepFunction_Type)

typedef struct {
    PyObject_HEAD
    PyObject* owner;
    /* Type-specific fields go here */
    LALEvolveOneStepFunction* func;
} li_LALEvolveOneStepFunction;

static PyObject* li_LALEvolveOneStepFunction_new(LALEvolveOneStepFunction* func, PyObject *owner){
    PyObject *empty_tuple = PyTuple_New(0);
    li_LALEvolveOneStepFunction *obj = (li_LALEvolveOneStepFunction *) PyType_GenericNew(&li_LALEvolveOneStepFunction_Type,empty_tuple,NULL);
    Py_DECREF(empty_tuple);
    
    if(!obj) {
        if(!owner)
            return NULL;
    }
    if(owner)
        Py_INCREF(owner);
    obj->owner = owner;
    
    obj->func = func;
    return (PyObject *) obj;
}
