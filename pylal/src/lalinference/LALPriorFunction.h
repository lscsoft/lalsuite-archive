//LALPriorFunction.h

static PyTypeObject *_li_LALPriorFunction_Type = NULL;
#define li_LALPriorFunction_Type (*_li_LALPriorFunction_Type)

typedef struct {
    PyObject_HEAD
    PyObject* owner;
    /* Type-specific fields go here */
    LALPriorFunction* func;
} li_LALPriorFunction;

static PyObject* li_LALPriorFunction_new(LALPriorFunction* func, PyObject *owner){
    PyObject *empty_tuple = PyTuple_New(0);
    li_LALPriorFunction *obj = (li_LALPriorFunction *) PyType_GenericNew(&li_LALPriorFunction_Type,empty_tuple,NULL);
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
