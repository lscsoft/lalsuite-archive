//LALVariables.h

#define LALINFERENCE_LALVARIABLES_MODULE_NAME "pylal.lalinference.lalvariables"

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


static PyObject *LALVariables_import(void)
{
    PyObject *name = PyString_FromString(LALINFERENCE_LALVARIABLES_MODULE_NAME);
    PyObject *module = PyImport_Import(name);
    Py_DECREF(name);

    name = PyString_FromString("LALVariables");
    _li_LALVariables_Type = (PyTypeObject *) PyDict_GetItem(PyModule_GetDict(module), name);
    Py_INCREF(&li_LALVariables_Type);
    Py_DECREF(name);

    return module;
}
