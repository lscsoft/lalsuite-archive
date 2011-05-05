//LALProposalFunction.h

static PyTypeObject *_li_LALInferenceProposalFunction_Type = NULL;
#define li_LALInferenceProposalFunction_Type (*_li_LALInferenceProposalFunction_Type)

typedef struct {
    PyObject_HEAD
    PyObject* owner;
    /* Type-specific fields go here */
    LALInferenceProposalFunction* func;
} li_LALInferenceProposalFunction;

static PyObject* li_LALInferenceProposalFunction_new(LALInferenceProposalFunction* func, PyObject *owner){
    PyObject *empty_tuple = PyTuple_New(0);
    li_LALInferenceProposalFunction *obj = (li_LALInferenceProposalFunction *) PyType_GenericNew(&li_LALInferenceProposalFunction_Type,empty_tuple,NULL);
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
