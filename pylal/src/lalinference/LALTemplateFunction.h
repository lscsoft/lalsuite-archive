//LALTemplateFunction.h

static PyTypeObject *_li_LALInferenceTemplateFunction_Type = NULL;
#define li_LALInferenceTemplateFunction_Type (*_li_LALInferenceTemplateFunction_Type)

typedef struct {
    PyObject_HEAD
    PyObject* owner;
    /* Type-specific fields go here */
    LALInferenceTemplateFunction* func;
} li_LALInferenceTemplateFunction;

static PyObject* li_LALInferenceTemplateFunction_new(LALInferenceTemplateFunction* func, PyObject *owner){
    PyObject *empty_tuple = PyTuple_New(0);
    li_LALInferenceTemplateFunction *obj = (li_LALInferenceTemplateFunction *) PyType_GenericNew(&li_LALInferenceTemplateFunction_Type,empty_tuple,NULL);
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
