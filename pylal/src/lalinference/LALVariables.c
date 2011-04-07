//LALVariables.c
#include <Python.h>
#include <structmember.h>

#include <numpy/arrayobject.h>

#include <lal/LALInference.h>

#include "LALVariables.h"

/*
 * ============================================================================
 *
 *                            LALVariables
 *
 * ============================================================================
 */

#define MODULE_NAME LALINFERENCE_LALVARIABLES_MODULE_NAME
const char LVDocString[] =
"This module provides data types and function wrappers for"
"LALVariables.";
/*Methods*/


/* Destructor for LALVariables */
static void LALVariables_dealloc(li_LALVariables *self)
{
    destroyVariables(self->vars);
    self->ob_type->tp_free((PyObject *)self);
}

VariableType LALVariables_infer_litype_from_pyvalue(PyObject* pyvalue){

    VariableType type;
    
    if(PyInt_Check(pyvalue)){
        type=INT8_t;
    }   
    else if(PyFloat_Check(pyvalue)){
        type=REAL8_t;
    }
    else{
        PyErr_SetObject(PyExc_TypeError, pyvalue);
        return -1;
    }

    return type;
}

VariableType LALVariables_convert_string_to_litype(char* typestring){

    VariableType type;
    if(!strcmp(typestring,"INT4")){
        type=INT4_t;
    }
    else if(!strcmp(typestring,"INT8")){
        type=INT8_t;
    }
    else if(!strcmp(typestring,"UINT4")){
        type=UINT4_t;
    }
    else if(!strcmp(typestring,"REAL4")){
        type=REAL4_t;
    }
    else if(!strcmp(typestring,"REAL8")){
        type=REAL8_t;
    }
    else{
        PyErr_SetString(PyExc_TypeError,"LALInference type not found!!");
        return -1;
    }
    return type;
}

void* LALVariables_convert_pyobj_to_livar_value(PyObject* pyvalue,VariableType type){
    void* value=(void *)malloc(typeSize[type]);
    
    if(type==INT4_t){
        INT4 cast_value=((INT4)PyInt_AsLong(pyvalue));
        INT4* cast_valuep=&cast_value;
        memcpy(value,(void*)cast_valuep,typeSize[type]);
    }
    else if(type==INT8_t){
        INT8 cast_value=((INT8)PyInt_AsLong(pyvalue));
        INT8* cast_valuep=&cast_value;
        memcpy(value,(void*)cast_valuep,typeSize[type]);
    }
    else if(type==UINT4_t){
        UINT4 cast_value=(UINT4)((unsigned long int)PyInt_AsLong(pyvalue));
        UINT4* cast_valuep=&cast_value;
        memcpy(value,(void*)cast_valuep,typeSize[type]);
    }   
    else if(type==REAL4_t){
    
        REAL4 cast_value=((REAL4)PyFloat_AsDouble(pyvalue));
        REAL4* cast_valuep=&cast_value;
        memcpy(value,(void*)cast_valuep,typeSize[type]);
        
    }
    else if(type==REAL8_t){
        REAL8 cast_value=((REAL8)PyFloat_AsDouble(pyvalue));
        REAL8* cast_valuep=&cast_value;
        memcpy(value,(void*)cast_valuep,typeSize[type]);
    }
    else{
        PyErr_SetObject(PyExc_TypeError, pyvalue);
        return NULL;
    }
    
    return value;
}

static PyObject* LALVariables_add_variable(li_LALVariables *self,PyObject* args,PyObject* kwds){
    PyObject *pyname=NULL,*pyvalue=NULL,*pyvarytype=NULL,*pytype=NULL;

    char *name;char* temp;void* value;VariableType type;ParamVaryType varytype;
    char *typestring=NULL;

    static char *kwlist[]={"name","value","varytype","type"};
    
    if (! PyArg_ParseTupleAndKeywords(args,kwds,"OOO|O",kwlist,&pyname,&pyvalue,&pyvarytype,&pytype)) return NULL;

    /*Extract name of variable*/
    if(PyString_Check(pyname)){
        name=PyString_AsString(pyname);
    }
    else{
        PyErr_SetObject(PyExc_TypeError, pyname);
        return NULL;
    }

    /*Extract and determine type of parameter value*/

    //If type given convert string to type...
    if(pytype) {
        typestring=PyString_AsString(pytype);
        type=LALVariables_convert_string_to_litype(typestring);
    }
    //...else infer type from python object type (this is more limited).
    else type=LALVariables_infer_litype_from_pyvalue(pyvalue);
    
    value=LALVariables_convert_pyobj_to_livar_value(pyvalue,type);
    
    /*Determine variable wrapping type from string*/
    if(PyString_Check(pyvarytype)){
        temp=PyString_AsString(pyvarytype);
        Py_INCREF(pyvarytype);

        if(!strcmp(temp,"linear")){
            varytype=PARAM_LINEAR;
        }
        else if(!strcmp(temp,"circular")){
            varytype=PARAM_CIRCULAR;
        }
        else if(!strcmp(temp,"fixed")){
            varytype=PARAM_FIXED;
        }
        else if(!strcmp(temp,"output")){
            varytype=PARAM_OUTPUT;
        }
        else {
            
            PyErr_SetObject(PyExc_ValueError,pyvarytype);
            return NULL;
        }
        Py_DECREF(pyvarytype);
    }
    else{
        PyErr_SetObject(PyExc_TypeError, pyvarytype);
        return NULL;
    }

    /* If we survived then call addVariable with self->vars ...*/
    
    addVariable((self->vars),name,value,type,varytype);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject* LALVariables_check_variable(li_LALVariables *self,PyObject* args)
/* Check for existence of name */
{
    PyObject* pyname;
    char* name;
    
    if (! PyArg_ParseTuple(args,"O",&pyname)) return NULL;

    name=PyString_AsString(pyname);

    if(getItem(self->vars,name)){
        Py_RETURN_TRUE;
    }
    else{
        Py_RETURN_FALSE;
    }
}


PyObject* LALVariables_convert_livar_value_to_pyobj(li_LALVariables *pyvars,char* name){
    VariableType type;
    
    PyObject* returnObj=NULL;
    
    type=getVariableType(pyvars->vars,name);
    void* uncastval=getVariable(pyvars->vars,name);
    if(type==INT4_t){
        
        long int cast_val=(long int)(*(INT4*)uncastval);
        returnObj=PyInt_FromLong(cast_val);
    }
    else if(type==INT8_t){
        long long cast_val=(long long)(*(INT8*)uncastval);
        returnObj=PyInt_FromLong(cast_val);
    }
    else if(type==UINT4_t){
        returnObj=PyInt_FromLong(*(UINT4*)uncastval);
    }
    else if(type==REAL4_t){
        float cast_val=(float)(*(REAL4*)uncastval);
        returnObj=PyFloat_FromDouble(cast_val);
    }
    else if(type==REAL8_t){
        returnObj=PyFloat_FromDouble(*(REAL8*)uncastval);
        
    }
    else {
        PyErr_SetString(PyExc_TypeError,"This type of variable cannot be converted to a python type!");
        return NULL;
    }
    
    return returnObj;
}

static int LALVariables_remove_variable(li_LALVariables *self,PyObject* args){
    PyObject* pyname;
    char* name;
    
    if (! PyArg_ParseTuple(args,"O",&pyname)) return -1;

    name=PyString_AsString(pyname);
    removeVariable(self->vars,name);
    
    return 0;
}

static PyObject* LALVariables_get_variable(li_LALVariables *self,PyObject* args){
    PyObject* pyname=NULL,*returnObj=NULL;
    char* name;
    
    if (! PyArg_ParseTuple(args,"O",&pyname)) return NULL;
    name=PyString_AsString(pyname);
    returnObj=LALVariables_convert_livar_value_to_pyobj(self,name);
    return returnObj;
}

static PyObject* LALVariables_get_variable_name(li_LALVariables *self,PyObject* args){
    char* name;
    long var_idx;
    
    if (! PyArg_ParseTuple(args,"i",&var_idx)) return NULL;
    
    name=getVariableName(self->vars,(int)var_idx);
    
    return PyString_FromString(name);
    
}

static int LALVariables_set_variable(li_LALVariables *self,PyObject* args){
    PyObject* pyname;char* name;
    PyObject* pynew_var_value;void* new_var_value;
    VariableType type;

    if (! PyArg_ParseTuple(args,"OO",&pyname,&pynew_var_value)) return -1;

    name=PyString_AsString(pyname);

    type=getVariableType(self->vars,name);

    new_var_value=LALVariables_convert_pyobj_to_livar_value(pynew_var_value,type);

    setVariable(self->vars,name,new_var_value);
    
    return 0;
    
}

char* LALVariables_get_type_string(VariableType type){
    char* type_name=NULL;
    
    if(type==INT4_t){
        type_name="INT4";
    }
    else if(type==UINT4_t){
        type_name="UINT4";
    }
    else if(type==INT8_t){
        type_name="INT8";
    }
    else if(type==REAL4_t){
        type_name="REAL4";
    }
    else if(type==REAL8_t){
        type_name="REAL8";
    }
    
    return type_name;
}

static PyObject* LALVariables_get_variable_type(li_LALVariables *self,PyObject* args){
    VariableType type;
    char* name;
    PyObject* pyname;
    char* type_name=NULL;
    if (! PyArg_ParseTuple(args,"O",&pyname)) return NULL;
    
    name=PyString_AsString(pyname);
    type=getVariableType(self->vars,name);
    
    type_name=LALVariables_get_type_string(type);
    
    if(type_name==NULL){
        Py_INCREF(Py_None);
        return Py_None;
    }    
    
    return PyString_FromString(type_name);
    
}
    
static PyObject* LALVariables_get_variable_type_by_index(li_LALVariables *self,PyObject* args){
    VariableType type;
    long var_idx;
    char* type_name=NULL;
    if (! PyArg_ParseTuple(args,"i",&var_idx)) return NULL;
    
    type=getVariableTypeByIndex(self->vars,(int)var_idx);
    type_name=LALVariables_get_type_string(type);
    
    return PyString_FromString(type_name);
    
}

char* LALVariables_get_varytype_string(ParamVaryType varytype){
    
    char* varytype_name=NULL;
    if(varytype==PARAM_LINEAR){
        varytype_name="linear";
    }
    else if(varytype==PARAM_CIRCULAR){
        varytype_name="circular";
    }
    else if(varytype==PARAM_FIXED){
        varytype_name="fixed";
    }
    else if(varytype==PARAM_OUTPUT){
        varytype_name="output";
    }
    
    return varytype_name;
    
}

static PyObject* LALVariables_get_variable_vary_type(li_LALVariables *self,PyObject* args){
    ParamVaryType varytype;
    char* name;
    PyObject* pyname;
    char* varytype_name=NULL;
    if (! PyArg_ParseTuple(args,"O",&pyname)) return NULL;
    
    name=PyString_AsString(pyname);
    varytype=getVariableVaryType(self->vars,name);
    varytype_name=LALVariables_get_varytype_string(varytype);
    
    if(varytype_name==NULL){
        Py_INCREF(Py_None);
        return Py_None;
    }    
    
    return PyString_FromString(varytype_name);
    
}

static PyObject* LALVariables_get_variable_dimension(li_LALVariables *self){
    long dimension;
    PyObject* pydimension=NULL;
    
    dimension=(long int)getVariableDimension(self->vars);

    pydimension=PyInt_FromLong(dimension);

    return pydimension;
    
}

static PyObject* LALVariables_get_variable_dimension_non_fixed(li_LALVariables *self){
    long dimension;
    PyObject* pydimension=NULL;
    
    dimension=(long int)getVariableDimensionNonFixed(self->vars);

    pydimension=PyInt_FromLong(dimension);

    return pydimension;
    
}

static int LALVariables_init(li_LALVariables *self, PyObject *args, PyObject *kwds)
{
    /* Should fill in the array using a dictionary as input */
    
    self->vars=(LALVariables*)malloc(sizeof(LALVariables));
    memset((void*)self->vars,0,sizeof(LALVariables));
    return 0;
}



static PyMethodDef LALVariables_methods[]= {
    /* {"name", (PyCFunction)function, METH_NOARGS, "DESCR"}, */
    {"addVariable",(PyCFunction)LALVariables_add_variable,METH_VARARGS,""},
    {"getVariable",(PyCFunction)LALVariables_get_variable,METH_VARARGS,""},
    {"checkVariable",(PyCFunction)LALVariables_check_variable,METH_VARARGS,""},
    {"removeVariable",(PyCFunction)LALVariables_remove_variable,METH_VARARGS,""},
    {"getVariableName",(PyCFunction)LALVariables_get_variable_name,METH_VARARGS,""},
    {"getVariableType",(PyCFunction)LALVariables_get_variable_type,METH_VARARGS,""},
    {"getVariableTypeByIndex",(PyCFunction)LALVariables_get_variable_type_by_index,METH_VARARGS,""},
    {"getVariableVaryType",(PyCFunction)LALVariables_get_variable_vary_type,METH_VARARGS,""},
    {"setVariable",(PyCFunction)LALVariables_set_variable,METH_VARARGS,""},
    {"getVariableDimension",(PyCFunction)LALVariables_get_variable_dimension,METH_NOARGS,""},
    {"getVariableDimensionNonFixed",(PyCFunction)LALVariables_get_variable_dimension_non_fixed,METH_NOARGS,""},
    
    {NULL} /* Sentinel */
};

static PyTypeObject li_lalvariables_type = {
    PyObject_HEAD_INIT(NULL)
    0,              /* obj_size - unused (must be 0) */
    "lalinference.lalvariables.BaseLALVariables",    /* tp_name, name of type */
    sizeof(li_LALVariables),  /* tp_basicsize */
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
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,        /*tp_flags*/
    "LALInference BaseLALVariables objects", /* tp_doc */
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
init_lalvariables(void)
{
    PyObject *m;

    _li_LALVariables_Type = &li_lalvariables_type;
    //li_LALVariables_Type.tp_new = PyType_GenericNew;
    if (PyType_Ready(&li_LALVariables_Type) < 0)
        return;    
    
    m = Py_InitModule3(MODULE_NAME,module_methods,LVDocString);

    import_array();
    
    Py_INCREF(&li_LALVariables_Type);
    PyModule_AddObject(m, "BaseLALVariables", (PyObject *)&li_LALVariables_Type);
}
