//BarycenterInput.c
#include <Python.h>
#include <structmember.h>

#include <numpy/arrayobject.h>

#include <lal/LALInference.h>

#include <lal/LALBarycenter.h>
#include <lal/LALDetectors.h>

#include <ligotimegps.h>
#include <tools.h>

#include "lalvariables.h"
#include "BarycenterInput.h"

/*
 * ============================================================================
 *
 *                            BarycenterInput
 *
 * ============================================================================
 */

#define MODULE_NAME LALINFERENCE_BARYCENTERINPUT_MODULE_NAME

/*Methods*/

 /* Destructor for LALIFOData */
static void BarycenterInput_dealloc(li_BarycenterInput *self)
{
    self->ob_type->tp_free((PyObject *)self);
}

static int BarycenterInput__init__(li_BarycenterInput *self, PyObject *args, PyObject *kwds)
{
    self->data=(BarycenterInput*)malloc(sizeof(BarycenterInput));
    memset((void*)self->data,0,sizeof(BarycenterInput));
    return 0;
}

static PyObject* BarycenterInput__new__(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    
    li_BarycenterInput *obj = (li_BarycenterInput*) PyType_GenericNew(type, args, kwds);
    if(!obj)
        return NULL;
    obj->data=(BarycenterInput*)malloc(sizeof(BarycenterInput));
    memset((void*)obj->data,0,sizeof(BarycenterInput));
    obj->tgps=NULL;
    obj->site=NULL;
    return (PyObject*)obj;
}

/*tgps*/
static PyObject* BarycenterInput_gettgps(li_BarycenterInput *self, void *closure){return pylal_LIGOTimeGPS_new(self->data->tgps);}
static int BarycenterInput_settgps(li_BarycenterInput *self, PyObject *value, void *closure){return setLIGOTimeGPSFromLALIFOData(self->data->tgps,value);}

/*site*/
static PyObject* BarycenterInput_getsite(li_BarycenterInput *self, void *closure){return getLALDetectorFromData(self->data->site);}
static int BarycenterInput_setsite(li_BarycenterInput *self, PyObject *value, void *closure){return setLALDetectorFromData(&self->data->site,self->site,value);}


static PyMethodDef BarycenterInput_methods[]= {
    {NULL} /* Sentinel */
};

static PyGetSetDef BarycenterInput_getseters[] = {
    {"tgps",(getter)BarycenterInput_gettgps,(setter)BarycenterInput_settgps,"name",NULL},
    {"site",(getter)BarycenterInput_getsite,(setter)BarycenterInput_setsite,"name",NULL},
    
    {NULL}  /* Sentinel */
};

static struct PyMemberDef BarycenterInput_members[] = {
    //REAL8's
    {"alpha", T_DOUBLE, offsetof(li_BarycenterInput, data)+offsetof(BarycenterInput,alpha), 0, "alpha"},
    {"delta", T_DOUBLE, offsetof(li_BarycenterInput, data)+offsetof(BarycenterInput,delta), 0, "delta"},
    {"dInv", T_DOUBLE, offsetof(li_BarycenterInput, data)+offsetof(BarycenterInput,dInv), 0, "dInv"},
    
    {NULL,}
};

static PyTypeObject li_barycenterinput_type = {
    PyObject_HEAD_INIT(NULL)
    0,              /* obj_size - unused (must be 0) */
    "lalinference.barycenterinput.BarycenterInput",    /* tp_name, name of type */
    sizeof(li_BarycenterInput),  /* tp_basicsize */
    0,              /* tp_itemsize, need to check */
    (destructor)BarycenterInput_dealloc,  /*tp_dealloc*/
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
    "", /* tp_doc */
    0,                     /* tp_traverse */
    0,                     /* tp_clear */
    0,                     /* tp_richcompare */
    0,                     /* tp_weaklistoffset */
    0,                     /* tp_iter */
    0,                     /* tp_iternext */
    BarycenterInput_methods,             /* tp_methods */
    BarycenterInput_members,             /* tp_members */
    BarycenterInput_getseters,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)BarycenterInput__init__,      /* tp_init */
    0,                         /* tp_alloc */
    BarycenterInput__new__,                 /* tp_new */
};

/*
 * ============================================================================
 *
 *                            Module Registration
 *
 * ============================================================================
 */
 
static PyMethodDef module_methods[] = {
    {NULL}  /* Sentinel */
};

PyMODINIT_FUNC
init_barycenterinput(void)
{
    PyObject *m;

    _li_BarycenterInput_Type = &li_barycenterinput_type;
    if (PyType_Ready(&li_BarycenterInput_Type) < 0)
        return;

    m = Py_InitModule3(MODULE_NAME,module_methods,BIDocString);

    import_array();
    
    pylal_ligotimegps_import();
    
    Py_INCREF(&li_LALIFOData_Type);
    PyModule_AddObject(m, "BarycenterInput", (PyObject *)&li_BarycenterInput_Type);
}
