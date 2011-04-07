//LALVariables.h

static PyTypeObject *_li_LALVariables_Type = NULL;
#define li_LALVariables_Type (*_li_LALVariables_Type)

typedef struct {
    PyObject_HEAD
    PyObject* owner;
    /* Type-specific fields go here */
    LALVariables* vars;
} li_LALVariables;

static PyObject* LALVariables_new(LALVariables *vars, PyObject *owner){
    PyObject *empty_tuple = PyTuple_New(0);
    li_LALVariables *obj = (li_LALVariables *) PyType_GenericNew(&li_LALVariables_Type,empty_tuple,NULL);
    Py_DECREF(empty_tuple);
    
    if(!obj) {
        if(!owner)
            return NULL;
    }
    if(owner)
        Py_INCREF(owner);
    obj->owner = owner;
    destroyVariables(obj->vars);
    obj->vars = vars;
    return (PyObject *) obj;
}
