//LALEvolveOneStepFunction.h

static PyTypeObject *_li_LALInferenceEvolveOneStepFunction_Type = NULL;
#define li_LALInferenceEvolveOneStepFunction_Type (*_li_LALInferenceEvolveOneStepFunction_Type)

typedef struct {
    PyObject_HEAD
    PyObject* owner;
    /* Type-specific fields go here */
    LALInferenceEvolveOneStepFunction* func;
} li_LALInferenceEvolveOneStepFunction;

static PyObject* li_LALInferenceEvolveOneStepFunction_new(LALInferenceEvolveOneStepFunction* func, PyObject *owner){
    PyObject *empty_tuple = PyTuple_New(0);
    li_LALInferenceEvolveOneStepFunction *obj = (li_LALInferenceEvolveOneStepFunction *) PyType_GenericNew(&li_LALInferenceEvolveOneStepFunction_Type,empty_tuple,NULL);
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
