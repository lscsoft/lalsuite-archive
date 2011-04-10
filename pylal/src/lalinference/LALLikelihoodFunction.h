//LALLikelihoodFunction.h

static PyTypeObject *_li_LALLikelihoodFunction_Type = NULL;
#define li_LALLikelihoodFunction_Type (*_li_LALLikelihoodFunction_Type)

typedef struct {
    PyObject_HEAD
    PyObject* owner;
    /* Type-specific fields go here */
    LALLikelihoodFunction* func;
} li_LALLikelihoodFunction;

static PyObject* li_LALLikelihoodFunction_new(LALLikelihoodFunction* func, PyObject *owner){
    PyObject *empty_tuple = PyTuple_New(0);
    li_LALLikelihoodFunction *obj = (li_LALLikelihoodFunction *) PyType_GenericNew(&li_LALLikelihoodFunction_Type,empty_tuple,NULL);
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
