//LALTemplateFunction.h

static PyTypeObject *_li_LALTemplateFunction_Type = NULL;
#define li_LALTemplateFunction_Type (*_li_LALTemplateFunction_Type)

typedef struct {
    PyObject_HEAD
    PyObject* owner;
    /* Type-specific fields go here */
    LALTemplateFunction* func;
} li_LALTemplateFunction;

static PyObject* li_LALTemplateFunction_new(LALTemplateFunction* func, PyObject *owner){
    PyObject *empty_tuple = PyTuple_New(0);
    li_LALTemplateFunction *obj = (li_LALTemplateFunction *) PyType_GenericNew(&li_LALTemplateFunction_Type,empty_tuple,NULL);
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
