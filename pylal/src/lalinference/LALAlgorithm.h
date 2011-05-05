//LALAlgorithm.h

static PyTypeObject *_li_LALInferenceAlgorithm_Type = NULL;
#define li_LALInferenceAlgorithm_Type (*_li_LALInferenceAlgorithm_Type)

typedef struct {
    PyObject_HEAD
    PyObject* owner;
    /* Type-specific fields go here */
    LALInferenceAlgorithm* func;
} li_LALInferenceAlgorithm;

static PyObject* li_LALInferenceAlgorithm_new(LALInferenceAlgorithm* func, PyObject *owner){
    PyObject *empty_tuple = PyTuple_New(0);
    li_LALInferenceAlgorithm *obj = (li_LALInferenceAlgorithm *) PyType_GenericNew(&li_LALInferenceAlgorithm_Type,empty_tuple,NULL);
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
