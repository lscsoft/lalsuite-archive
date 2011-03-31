#include <Python.h>
#include <lal/LALInference.h>

#define MODULE_NAME "_lalinference"

const char LIDocString[] =
"This module provides data types and function wrappers for"
"LALInference.";

typedef struct {
    PyObject_HEAD
    /* Type-specific fields go here */
    LALVariables var;
} li_LALVariablesObject;

/* Destructor for LALVariables */
static void LALVariables_dealloc(li_LALVariablesObject *self)
{
    destroyVariables(&(self->var));
    self->ob_type->tp_free((PyObject *)self);
}

static PyMethodDef LALVariables_methods[]= {
    /* {"name", (PyCFunction)function, METH_NOARGS, "DESCR"}, */
    {NULL} /* Sentinel */
};

static int LALVariables_init(li_LALVariablesObject *self, PyObject *args, PyObject *kwds)
{
    /* Should fill in the array using a dictionary as input */
    
    return 0;
}

static PyTypeObject li_LALVariablesType = {
    PyObject_HEAD_INIT(NULL)
    0,              /* obj_size - unused (must be 0) */
    "lalinference.LALVariables",    /* tp_name, name of type */
    sizeof(li_LALVariablesObject),  /* tp_basicsize */
    0,              /* tp_itemsize, need to check */
    (destructor)LALVariables_dealloc,  /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT,        /*tp_flags*/
    "LALInference LALVariables objects", /* tp_doc */
    0,                     /* tp_traverse */
    0,                     /* tp_clear */
    0,                     /* tp_richcompare */
    0,                     /* tp_weaklistoffset */
    0,                     /* tp_iter */
    0,                     /* tp_iternext */
    LALVariables_methods,             /* tp_methods */
    0,             /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)LALVariables_init,      /* tp_init */
    0,                         /* tp_alloc */
    PyType_GenericNew,                 /* tp_new */
};

static PyMethodDef li_methods[]={
    {NULL} /* Sentinel */
};

PyMODINIT_FUNC
init_lalinference(void)
{
    PyObject *m;
    li_LALVariablesType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&li_LALVariablesType) < 0)
        return;
    
    m = Py_InitModule3(MODULE_NAME,li_methods,LIDocString);
    Py_INCREF(&li_LALVariablesType);
    PyModule_AddObject(m, "LALVariables", (PyObject *)&li_LALVariablesType);
}
