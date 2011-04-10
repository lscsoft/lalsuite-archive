//LALProposalFunction.h

static PyTypeObject *_li_LALProposalFunction_Type = NULL;
#define li_LALProposalFunction_Type (*_li_LALProposalFunction_Type)

typedef struct {
    PyObject_HEAD
    PyObject* owner;
    /* Type-specific fields go here */
    LALProposalFunction* func;
} li_LALProposalFunction;

static PyObject* li_LALProposalFunction_new(LALProposalFunction* func, PyObject *owner){
    PyObject *empty_tuple = PyTuple_New(0);
    li_LALProposalFunction *obj = (li_LALProposalFunction *) PyType_GenericNew(&li_LALProposalFunction_Type,empty_tuple,NULL);
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
