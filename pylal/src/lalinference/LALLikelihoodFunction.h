//LALLikelihoodFunction.h

static PyTypeObject *_li_LALInferenceLikelihoodFunction_Type = NULL;
#define li_LALInferenceLikelihoodFunction_Type (*_li_LALInferenceLikelihoodFunction_Type)

typedef struct {
    PyObject_HEAD
    PyObject* owner;
    /* Type-specific fields go here */
    LALInferenceLikelihoodFunction* func;
} li_LALInferenceLikelihoodFunction;

static PyObject* li_LALInferenceLikelihoodFunction_new(LALInferenceLikelihoodFunction* func, PyObject *owner){
    PyObject *empty_tuple = PyTuple_New(0);
    li_LALInferenceLikelihoodFunction *obj = (li_LALInferenceLikelihoodFunction *) PyType_GenericNew(&li_LALInferenceLikelihoodFunction_Type,empty_tuple,NULL);
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
