//LALAlgorithm.h

static PyTypeObject *_li_LALAlgorithm_Type = NULL;
#define li_LALAlgorithm_Type (*_li_LALAlgorithm_Type)

typedef struct {
    PyObject_HEAD
    PyObject* owner;
    /* Type-specific fields go here */
    LALAlgorithm* func;
} li_LALAlgorithm;

static PyObject* li_LALAlgorithm_new(LALAlgorithm* func, PyObject *owner){
    PyObject *empty_tuple = PyTuple_New(0);
    li_LALAlgorithm *obj = (li_LALAlgorithm *) PyType_GenericNew(&li_LALAlgorithm_Type,empty_tuple,NULL);
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
