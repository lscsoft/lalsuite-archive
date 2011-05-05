//LALPriorFunction.h

static PyTypeObject *_li_LALInferencePriorFunction_Type = NULL;
#define li_LALInferencePriorFunction_Type (*_li_LALInferencePriorFunction_Type)

typedef struct {
    PyObject_HEAD
    PyObject* owner;
    /* Type-specific fields go here */
    LALInferencePriorFunction* func;
} li_LALInferencePriorFunction;

static PyObject* li_LALInferencePriorFunction_new(LALInferencePriorFunction* func, PyObject *owner){
    PyObject *empty_tuple = PyTuple_New(0);
    li_LALInferencePriorFunction *obj = (li_LALInferencePriorFunction *) PyType_GenericNew(&li_LALInferencePriorFunction_Type,empty_tuple,NULL);
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
