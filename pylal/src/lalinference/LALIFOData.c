#include <Python.h>
#include <structmember.h>

#include <lal/LALInference.h>

#define MODULE_NAME "pylal._lalinference.lalifodata"

const char LIDocString[] =
"This module provides data types and function wrappers for"
"LALIFOData.";

typedef struct {
    PyObject_HEAD
    /* Type-specific fields go here */
    LALIFOData data;
} li_LALIFODataObject;

/* Destructor for LALIFOData */
static void LALIFOData_dealloc(li_LALIFODataObject *self)
{
    self->ob_type->tp_free((PyObject *)self);
}

static int LALIFOData_init(li_LALIFODataObject *self, PyObject *args, PyObject *kwds)
{
    return 0;
}

static PyMethodDef LALIFOData_methods[]= {
    /* {"name", (PyCFunction)function, METH_NOARGS, "DESCR"}, */
    {NULL} /* Sentinel */
};

static PyTypeObject li_LALIFODataType = {
    PyObject_HEAD_INIT(NULL)
    0,              /* obj_size - unused (must be 0) */
    "lalinference.BaseLALIFOData",    /* tp_name, name of type */
    sizeof(li_LALIFODataObject),  /* tp_basicsize */
    0,              /* tp_itemsize, need to check */
    (destructor)LALIFOData_dealloc,  /*tp_dealloc*/
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
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,        /*tp_flags*/
    "LALInference BaseLALIFOData object.", /* tp_doc */
    0,                     /* tp_traverse */
    0,                     /* tp_clear */
    0,                     /* tp_richcompare */
    0,                     /* tp_weaklistoffset */
    0,                     /* tp_iter */
    0,                     /* tp_iternext */
    LALIFOData_methods,             /* tp_methods */
    0,             /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)LALIFOData_init,      /* tp_init */
    0,                         /* tp_alloc */
    PyType_GenericNew,                 /* tp_new */
};

static PyMethodDef module_methods[] = {
    {NULL}  /* Sentinel */
};

PyMODINIT_FUNC
init_lalifodata(void)
{
    PyObject *m;
    li_LALIFODataType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&li_LALIFODataType) < 0)
        return;
    
    m = Py_InitModule3(MODULE_NAME,module_methods,LIDocString);
    Py_INCREF(&li_LALIFODataType);
    PyModule_AddObject(m, "BaseLALIFOData", (PyObject *)&li_LALIFODataType);
}
