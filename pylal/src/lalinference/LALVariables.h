//LALVariables.h

#define LALINFERENCE_LALVARIABLES_MODULE_NAME "pylal._lalinference"

static PyTypeObject *_li_LALInferenceVariables_Type = NULL;
#define li_LALInferenceVariables_Type (*_li_LALInferenceVariables_Type)

typedef struct {
    PyObject_HEAD
    PyObject* owner;
    /* Type-specific fields go here */
    LALInferenceVariables* vars;
} li_LALInferenceVariables;

static PyObject* LALInferenceVariables_new(LALInferenceVariables *vars, PyObject *owner){
    PyObject *empty_tuple = PyTuple_New(0);
    li_LALInferenceVariables *obj = (li_LALInferenceVariables *) PyType_GenericNew(&li_LALInferenceVariables_Type,empty_tuple,NULL);
    Py_DECREF(empty_tuple);
    
    if(!obj) {
        if(!owner)
            return NULL;
    }
    if(owner)
        Py_INCREF(owner);
    obj->owner = owner;
    LALInferenceDestroyVariables(obj->vars);
    obj->vars = vars;
    return (PyObject *) obj;
}


static PyObject *LALInferenceVariables_import(void)
{
    PyObject *name = PyString_FromString(LALINFERENCE_LALVARIABLES_MODULE_NAME);
    PyObject *module = PyImport_Import(name);
    Py_DECREF(name);

    name = PyString_FromString("LALInferenceVariables");
    _li_LALInferenceVariables_Type = (PyTypeObject *) PyDict_GetItem(PyModule_GetDict(module), name);
    Py_INCREF(&li_LALInferenceVariables_Type);
    Py_DECREF(name);

    return module;
}
