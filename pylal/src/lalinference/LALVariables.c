#include <Python.h>
#include <structmember.h>

#include <lal/LALInference.h>

#define MODULE_NAME "pylal._lalinference.lalvariables"

const char LIDocString[] =
"This module provides data types and function wrappers for"
"LALVariables.";

typedef struct {
    PyObject_HEAD
    /* Type-specific fields go here */
    LALVariables vars;
} li_LALVariablesObject;

/* Destructor for LALVariables */
static void LALVariables_dealloc(li_LALVariablesObject *self)
{
    destroyVariables(&self->vars);
    self->ob_type->tp_free((PyObject *)self);
}

VariableType infer_litype_from_pyvalue(PyObject* pyvalue){

    VariableType type;
    
    if(PyInt_Check(pyvalue)){
        type=INT8_t;
    }   
    else if(PyFloat_Check(pyvalue)){
        type=REAL8_t;
    }
    else{
        PyErr_SetObject(PyExc_TypeError, pyvalue);
        return NULL;
    }

    return type;
}

VariableType convert_string_to_litype(char* typestring){

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
        return NULL;
    }
    return type;
}


static PyObject* add_variable(li_LALVariablesObject *self,PyObject* args){
    PyObject *pyname=NULL,*pyvalue=NULL,*pyvarytype=NULL;

    char *name;char* temp;void* value;VariableType type;ParamVaryType varytype;
    char *typestring=NULL;
    
    if (! PyArg_ParseTuple(args,"OOO",&pyname,&pyvalue,&pyvarytype)) {
        PyErr_SetString(PyExc_TypeError, "Input");
        return NULL;
    }

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
    //if(typestring) type=convert_string_to_litype(typestring);
    //...else infer type from python object type (this is more limited).
    //else{
        type=infer_litype_from_pyvalue(pyvalue);
    //}
    
    if(PyInt_Check(pyvalue)){
        value=(void*)((INT8*)PyInt_AsLong(pyvalue));
        type=INT8_t;
    }   
    else if(PyFloat_Check(pyvalue)){
        REAL8 cast_value=((REAL8)PyFloat_AsDouble(pyvalue));
        REAL8* cast_valuep=&cast_value;
        value=(void*)cast_valuep;
        type=REAL8_t;
    }
    else{
        PyErr_SetObject(PyExc_TypeError, pyvalue);
        return NULL;
    }

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
    //DEBUG//printf("%s %lf %d %d\n",name,*(double*)value,type,varytype);
    addVariable(&(self->vars),name,value,type,varytype);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject* check_variable(li_LALVariablesObject *self,PyObject* args)
/* Check for existence of name */
{
    PyObject* pyname;
    char* name;
    
    if (! PyArg_ParseTuple(args,"O",&pyname)) return NULL;

    name=PyString_AsString(pyname);

    if(getItem(&self->vars,name)){
        Py_RETURN_TRUE;
    }
    else{
        Py_RETURN_FALSE;
    }
}


PyObject* convert_livar_value_to_pyobj(li_LALVariablesObject *pyvars,char* name){
    VariableType type;
    
    PyObject* returnObj;
    
    
    type=getVariableType(&pyvars->vars,name);
    
    if(type==INT4_t){
        returnObj=PyInt_FromLong(*(long int*)getVariable(&pyvars->vars,name));
    }
    else if(type==INT8_t){
        returnObj=PyInt_FromLong(*(long int*)getVariable(&pyvars->vars,name));
    }
    else if(type==UINT4_t){
        returnObj=PyInt_FromLong(*(long int*)getVariable(&pyvars->vars,name));
    }
    else if(type==REAL4_t){
        returnObj=PyFloat_FromDouble(*(double*)getVariable(&pyvars->vars,name));
    }
    else if(type==REAL8_t){
        returnObj=PyFloat_FromDouble(*(double*)getVariable(&pyvars->vars,name));
    }
    else {
        PyErr_SetString(PyExc_TypeError,"This type of variable cannot be converted to a python type!");
        return NULL;
    }
    
    return Py_BuildValue("O",returnObj);
}

void* convert_pyobj_to_livar_value(PyObject* inObj){
    void* value=NULL;
    if(PyInt_Check(inObj)){
        INT8 cast_value=((INT8)PyInt_AsLong(inObj));
        value=(void*)&cast_value;
        
    }
    else if(PyFloat_Check(inObj)){
        REAL8 cast_value=((REAL8)PyFloat_AsDouble(inObj));
        value=(void*)&cast_value;
    }
    else{
    }
    return value;
}

static int remove_variable(li_LALVariablesObject *self,PyObject* args){
    PyObject* pyname;
    char* name;
    
    if (! PyArg_ParseTuple(args,"O",&pyname)) return -1;

    name=PyString_AsString(pyname);
    removeVariable(&self->vars,name);
    
    return 0;
}

static PyObject* get_variable(li_LALVariablesObject *self,PyObject* args){
    PyObject* pyname;
    char* name;
    
    if (! PyArg_ParseTuple(args,"O",&pyname)) return NULL;
    name=PyString_AsString(pyname);
    return convert_livar_value_to_pyobj(self,name);
    
}

static PyObject* get_variable_name(li_LALVariablesObject *self,PyObject* args){
    char* name;
    long var_idx;
    
    if (! PyArg_ParseTuple(args,"i",&var_idx)) return NULL;
    
    name=getVariableName(&self->vars,(int)var_idx);
    
    return PyString_FromString(name);
    
}

char* get_type_string(VariableType type){
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

static PyObject* get_variable_type(li_LALVariablesObject *self,PyObject* args){
    VariableType type;
    char* name;
    PyObject* pyname;
    char* type_name=NULL;
    if (! PyArg_ParseTuple(args,"O",&pyname)) return NULL;
    
    name=PyString_FromString(pyname);
    type=getVariableType(&self->vars,name);
    
    type_name=get_type_string(type);
    
    if(type_name==NULL){
        Py_INCREF(Py_None);
        return Py_None;
    }    
    
    return PyString_FromString(type_name);
    
}
    
static PyObject* get_variable_type_by_index(li_LALVariablesObject *self,PyObject* args){
    VariableType type;
    long var_idx;
    char* type_name=NULL;
    if (! PyArg_ParseTuple(args,"i",&var_idx)) return NULL;
    
    type=getVariableTypeByIndex(&self->vars,(int)var_idx);
    type_name=get_type_string(type);
    
    return PyString_FromString(type_name);
    
}

char* get_varytype_string(ParamVaryType varytype){
    
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

static PyObject* get_variable_vary_type(li_LALVariablesObject *self,PyObject* args){
    ParamVaryType varytype;
    char* name;
    PyObject* pyname;
    char* varytype_name=NULL;
    if (! PyArg_ParseTuple(args,"O",&pyname)) return NULL;
    
    name=PyString_FromString(pyname);
    varytype=getVariableVaryType(&self->vars,name);
    varytype_name=get_varytype_string(varytype);
    
    if(varytype_name==NULL){
        Py_INCREF(Py_None);
        return Py_None;
    }    
    
    return PyString_FromString(varytype_name);
    
}

static int LALVariables_init(li_LALVariablesObject *self, PyObject *args, PyObject *kwds)
{
    /* Should fill in the array using a dictionary as input */
    
    //memset(self->vars,0,sizeof(LALVariables));

    return 0;
}

static PyMethodDef LALVariables_methods[]= {
    /* {"name", (PyCFunction)function, METH_NOARGS, "DESCR"}, */
    {"addVariable",(PyCFunction)add_variable,METH_VARARGS,""},
    {"getVariable",(PyCFunction)get_variable,METH_VARARGS,""},
    {"checkVariable",(PyCFunction)check_variable,METH_VARARGS,""},
    {"removeVariable",(PyCFunction)remove_variable,METH_VARARGS,""},
    {"getVariableName",(PyCFunction)get_variable_name,METH_VARARGS,""},
    {"getVariableType",(PyCFunction)get_variable_type,METH_VARARGS,""},
    {"getVariableTypeByIndex",(PyCFunction)get_variable_type_by_index,METH_VARARGS,""},
    {"getVariableVaryType",(PyCFunction)get_variable_vary_type,METH_VARARGS,""},
    
    {NULL} /* Sentinel */
};

static PyTypeObject li_LALVariablesType = {
    PyObject_HEAD_INIT(NULL)
    0,              /* obj_size - unused (must be 0) */
    "lalinference.BaseLALVariables",    /* tp_name, name of type */
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

static PyMethodDef module_methods[] = {
    {NULL}  /* Sentinel */
};

PyMODINIT_FUNC
init_lalvariables(void)
{
    PyObject *m;
    li_LALVariablesType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&li_LALVariablesType) < 0)
        return;
    
    m = Py_InitModule3(MODULE_NAME,module_methods,LIDocString);
    Py_INCREF(&li_LALVariablesType);
    PyModule_AddObject(m, "BaseLALVariables", (PyObject *)&li_LALVariablesType);
}
