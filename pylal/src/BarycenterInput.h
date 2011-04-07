//BarycenterInput.h

static PyTypeObject *_li_BarycenterInput_Type = NULL;
#define li_BarycenterInput_Type (*_li_BarycenterInput_Type)

typedef struct {
    PyObject_HEAD
    PyObject* owner;
    /* Type-specific fields go here */
    BarycenterInput* data;
    PyObject* tgps;
    PyObject* site;
}li_BarycenterInput;

static PyObject* BarycenterInput_new(BarycenterInput *bci, PyObject *owner){
    PyObject *empty_tuple = PyTuple_New(0);
    li_BarycenterInput *obj = (li_BarycenterInput *) PyType_GenericNew(&li_BarycenterInput_Type,empty_tuple,NULL);
    Py_DECREF(empty_tuple);
    
    if(!obj) {
        if(!owner)
            return NULL;
    }
    if(owner)
        Py_INCREF(owner);
    obj->owner = owner;
    
    obj->data = bci;
    return (PyObject *) obj;
}
