#include <Python.h>
#include <structmember.h>

#include <lal/LALInference.h>

#define MODULE_NAME "pylal._lalinference"

const char LIDocString[] =
"This module provides data types and function wrappers for"
"LALInference.";

/*
 * ============================================================================
 *
 *                            LALVariables
 *
 * ============================================================================
 */

typedef struct {
    PyObject_HEAD
    /* Type-specific fields go here */
    LALVariables* vars;
} li_LALVariablesObject;

/*Methods*/


/* Destructor for LALVariables */
static void LALVariables_dealloc(li_LALVariablesObject *self)
{
    destroyVariables(self->vars);
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
        return -1;
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
        return -1;
    }
    return type;
}

void* convert_pyobj_to_livar_value(PyObject* pyvalue,VariableType type){
    void* value=NULL;
    
    if(type==INT4_t){
        INT4 cast_value=((INT4)PyInt_AsLong(pyvalue));
        INT4* cast_valuep=&cast_value;
        value=(void*)cast_valuep;
    }
    else if(type==INT8_t){
        INT8 cast_value=((INT8)PyInt_AsLong(pyvalue));
        //printf("%li\n",(long int)cast_value);
        INT8* cast_valuep=&cast_value;
        value=(void*)cast_valuep;
    }
    else if(type==UINT4_t){
        UINT4 cast_value=(UINT4)((unsigned long int)PyInt_AsLong(pyvalue));
        //printf("%lu\n",(unsigned long int)cast_value);
        UINT4* cast_valuep=&cast_value;
        
        value=(void*)cast_valuep;
    }   
    else if(type==REAL4_t){
    
        REAL4 cast_value=((REAL4)PyFloat_AsDouble(pyvalue));
        REAL4* cast_valuep=&cast_value;
        value=(void*)cast_valuep;
        
    }
    else if(type==REAL8_t){
        REAL8 cast_value=((REAL8)PyFloat_AsDouble(pyvalue));
        
        REAL8* cast_valuep=&cast_value;
        printf("%f\n",*cast_valuep);
        value=(void*)cast_valuep;
    }
    else{
        PyErr_SetObject(PyExc_TypeError, pyvalue);
        return NULL;
    }
    return value;
}

static PyObject* add_variable(li_LALVariablesObject *self,PyObject* args,PyObject* kwds){
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
        type=convert_string_to_litype(typestring);
    }
    //...else infer type from python object type (this is more limited).
    else type=infer_litype_from_pyvalue(pyvalue);
    
    value=convert_pyobj_to_livar_value(pyvalue,type);

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
    printf("%s %f %d %d\n",name,*((REAL8*)value),type,varytype);//DEBUG//
    addVariable((self->vars),name,value,type,varytype);
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

    if(getItem(self->vars,name)){
        Py_RETURN_TRUE;
    }
    else{
        Py_RETURN_FALSE;
    }
}


PyObject* convert_livar_value_to_pyobj(li_LALVariablesObject *pyvars,char* name){
    VariableType type;
    
    PyObject* returnObj;
    
    type=getVariableType(pyvars->vars,name);
    
    if(type==INT4_t){
        returnObj=PyInt_FromLong((long int)(*(INT4*)getVariable(pyvars->vars,name)));
    }
    else if(type==INT8_t){
        returnObj=PyInt_FromLong(*(INT8*)getVariable(pyvars->vars,name));
    }
    else if(type==UINT4_t){
        returnObj=PyInt_FromLong((long int)(*(UINT4*)getVariable(pyvars->vars,name)));
    }
    else if(type==REAL4_t){
        returnObj=PyFloat_FromDouble((double)(*(REAL4*)getVariable(pyvars->vars,name)));
    }
    else if(type==REAL8_t){
        REAL8 val=*((REAL8*)getVariable(pyvars->vars,name));
        printf("%f",val);
        
        returnObj=PyFloat_FromDouble(val);
        
    }
    else {
        PyErr_SetString(PyExc_TypeError,"This type of variable cannot be converted to a python type!");
        return NULL;
    }
    
    return returnObj;
}

static int remove_variable(li_LALVariablesObject *self,PyObject* args){
    PyObject* pyname;
    char* name;
    
    if (! PyArg_ParseTuple(args,"O",&pyname)) return -1;

    name=PyString_AsString(pyname);
    removeVariable(self->vars,name);
    
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
    
    name=getVariableName(self->vars,(int)var_idx);
    
    return PyString_FromString(name);
    
}

static int set_variable(li_LALVariablesObject *self,PyObject* args){
    PyObject* pyname;char* name;
    PyObject* pynew_var_value;void* new_var_value;
    VariableType type;

    if (! PyArg_ParseTuple(args,"OO",&pyname,&pynew_var_value)) return -1;

    name=PyString_AsString(pyname);

    type=getVariableType(self->vars,name);

    new_var_value=convert_pyobj_to_livar_value(pynew_var_value,type);

    setVariable(self->vars,name,new_var_value);
    
    return 0;
    
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
    
    name=PyString_AsString(pyname);
    type=getVariableType(self->vars,name);
    
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
    
    type=getVariableTypeByIndex(self->vars,(int)var_idx);
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
    
    name=PyString_AsString(pyname);
    varytype=getVariableVaryType(self->vars,name);
    varytype_name=get_varytype_string(varytype);
    
    if(varytype_name==NULL){
        Py_INCREF(Py_None);
        return Py_None;
    }    
    
    return PyString_FromString(varytype_name);
    
}

static PyObject* get_variable_dimension(li_LALVariablesObject *self){
    long dimension;
    PyObject* pydimension=NULL;
    
    dimension=(long int)getVariableDimension(self->vars);

    pydimension=PyInt_FromLong(dimension);

    return pydimension;
    
}

static PyObject* get_variable_dimension_non_fixed(li_LALVariablesObject *self){
    long dimension;
    PyObject* pydimension=NULL;
    
    dimension=(long int)getVariableDimensionNonFixed(self->vars);

    pydimension=PyInt_FromLong(dimension);

    return pydimension;
    
}

static int LALVariables_init(li_LALVariablesObject *self, PyObject *args, PyObject *kwds)
{
    /* Should fill in the array using a dictionary as input */
    
    self->vars=(LALVariables*)malloc(sizeof(LALVariables));

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
    {"setVariable",(PyCFunction)set_variable,METH_VARARGS,""},
    {"getVariableDimension",(PyCFunction)get_variable_dimension,METH_NOARGS,""},
    {"getVariableDimensionNonFixed",(PyCFunction)get_variable_dimension_non_fixed,METH_NOARGS,""},
    
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

/*
 * ============================================================================
 *
 *                            LALIFOData
 *
 * ============================================================================
 */

typedef struct {
    PyObject_HEAD
    /* Type-specific fields go here */
    LALIFOData* data;
} li_LALIFODataObject;


/*Methods*/

 /* Destructor for LALIFOData */
static void LALIFOData_dealloc(li_LALIFODataObject *self)
{
    self->ob_type->tp_free((PyObject *)self);
}

static int LALIFOData_init(li_LALIFODataObject *self, PyObject *args, PyObject *kwds)
{
    self->data=(LALIFOData*)malloc(sizeof(LALIFOData));
    return 0;
}

/***********getsetters******************/

/*name*/

static PyObject* LALIFOData_getname(li_LALIFODataObject *self, void *closure)
{
    return PyString_FromString(self->data->name);
}

static int LALIFOData_setname(li_LALIFODataObject *self, PyObject *value, void *closure)
{
  if (value == NULL) {
    PyErr_SetString(PyExc_TypeError, "Cannot delete the name attribute");
    return -1;
  }
  
  if (! PyString_Check(value)) {
    PyErr_SetString(PyExc_TypeError, 
                    "The name attribute value must be a string");
    return -1;
  }
  Py_INCREF(value);
  
  strcpy(self->data->name,PyString_AsString(value));

  return 0;
}

/*nullloglikelihood*/

static PyObject* LALIFOData_getnullloglikelihood(li_LALIFODataObject *self, void *closure)
{
    return PyFloat_FromDouble((double)(self->data->nullloglikelihood));
}

static int LALIFOData_setnullloglikelihood(li_LALIFODataObject *self, PyObject *value, void *closure)
{
    if (value == NULL) {
    PyErr_SetString(PyExc_TypeError, "Cannot delete the nullloglikelihood attribute");
    return -1;
  }
  
  if (! PyFloat_Check(value)) {
    PyErr_SetString(PyExc_TypeError, 
                    "The nullloglikelihood attribute value must be a float");
    return -1;
  }
      
  self->data->nullloglikelihood=(REAL8)PyFloat_AsDouble(value);

  return 0;
}

/*loglikelihood*/

static PyObject* LALIFOData_getloglikelihood(li_LALIFODataObject *self, void *closure)
{
    return PyFloat_FromDouble((double)(self->data->loglikelihood));
}

static int LALIFOData_setloglikelihood(li_LALIFODataObject *self, PyObject *value, void *closure)
{
    if (value == NULL) {
    PyErr_SetString(PyExc_TypeError, "Cannot delete the loglikelihood attribute");
    return -1;
  }
  
  if (! PyFloat_Check(value)) {
    PyErr_SetString(PyExc_TypeError, 
                    "The loglikelihood attribute value must be a float");
    return -1;
  }
      
  self->data->loglikelihood=(REAL8)PyFloat_AsDouble(value);

  return 0;
}

/*acceptedloglikelihood*/


static PyObject* LALIFOData_getacceptedloglikelihood(li_LALIFODataObject *self, void *closure)
{
    return PyFloat_FromDouble((double)(self->data->acceptedloglikelihood));
}

static int LALIFOData_setacceptedloglikelihood(li_LALIFODataObject *self, PyObject *value, void *closure)
{
    if (value == NULL) {
    PyErr_SetString(PyExc_TypeError, "Cannot delete the acceptedloglikelihood attribute");
    return -1;
  }
  
  if (! PyFloat_Check(value)) {
    PyErr_SetString(PyExc_TypeError, 
                    "The acceptedloglikelihood attribute value must be a float");
    return -1;
  }
      
  self->data->acceptedloglikelihood=(REAL8)PyFloat_AsDouble(value);

  return 0;
}

/*fPlus*/

static PyObject* LALIFOData_getfPlus(li_LALIFODataObject *self, void *closure)
{
    return PyFloat_FromDouble((double)(self->data->fPlus));
}

static int LALIFOData_setfPlus(li_LALIFODataObject *self, PyObject *value, void *closure)
{
    if (value == NULL) {
    PyErr_SetString(PyExc_TypeError, "Cannot delete the fPlus attribute");
    return -1;
  }
  
  if (! PyFloat_Check(value)) {
    PyErr_SetString(PyExc_TypeError, 
                    "The fPlus attribute value must be a float");
    return -1;
  }
      
  self->data->fPlus=(REAL8)PyFloat_AsDouble(value);

  return 0;
}

/*fCross*/

static PyObject* LALIFOData_getfCross(li_LALIFODataObject *self, void *closure)
{
    return PyFloat_FromDouble((double)(self->data->fCross));
}

static int LALIFOData_setfCross(li_LALIFODataObject *self, PyObject *value, void *closure)
{
    if (value == NULL) {
    PyErr_SetString(PyExc_TypeError, "Cannot delete the fCross attribute");
    return -1;
  }
  
  if (! PyFloat_Check(value)) {
    PyErr_SetString(PyExc_TypeError, 
                    "The fCross attribute value must be a float");
    return -1;
  }
      
  self->data->fCross=(REAL8)PyFloat_AsDouble(value);

  return 0;
}

/*timeshift*/

static PyObject* LALIFOData_gettimeshift(li_LALIFODataObject *self, void *closure)
{
    return PyFloat_FromDouble((double)(self->data->timeshift));
}

static int LALIFOData_settimeshift(li_LALIFODataObject *self, PyObject *value, void *closure)
{
    if (value == NULL) {
    PyErr_SetString(PyExc_TypeError, "Cannot delete the timeshift attribute");
    return -1;
  }
  
  if (! PyFloat_Check(value)) {
    PyErr_SetString(PyExc_TypeError, 
                    "The timeshift attribute value must be a float");
    return -1;
  }
      
  self->data->timeshift=(REAL8)PyFloat_AsDouble(value);

  return 0;
}

/*fLow*/
static PyObject* LALIFOData_getfLow(li_LALIFODataObject *self, void *closure)
{
    return PyFloat_FromDouble((double)(self->data->fLow));
}

static int LALIFOData_setfLow(li_LALIFODataObject *self, PyObject *value, void *closure)
{
    if (value == NULL) {
    PyErr_SetString(PyExc_TypeError, "Cannot delete the fLow attribute");
    return -1;
  }
  
  if (! PyFloat_Check(value)) {
    PyErr_SetString(PyExc_TypeError, 
                    "The fLow attribute value must be a float");
    return -1;
  }
      
  self->data->fLow=(REAL8)PyFloat_AsDouble(value);

  return 0;
}

/*fHigh*/
static PyObject* LALIFOData_getfHigh(li_LALIFODataObject *self, void *closure)
{
    return PyFloat_FromDouble((double)(self->data->fHigh));
}

static int LALIFOData_setfHigh(li_LALIFODataObject *self, PyObject *value, void *closure)
{
    if (value == NULL) {
    PyErr_SetString(PyExc_TypeError, "Cannot delete the fHigh attribute");
    return -1;
  }
  
  if (! PyFloat_Check(value)) {
    PyErr_SetString(PyExc_TypeError, 
                    "The fHigh attribute value must be a float");
    return -1;
  }
      
  self->data->fHigh=(REAL8)PyFloat_AsDouble(value);

  return 0;
}

/*modelParams*/
static PyObject* LALIFOData_getmodelParams(li_LALIFODataObject *self, void *closure)
{
    LALVariables* val=(LALVariables*)malloc(sizeof(LALVariables));
    copyVariables(val,self->data->modelParams);
    return NULL;
}

static int LALIFOData_setmodelParams(li_LALIFODataObject *self, PyObject *value, void *closure)
{
    if (value == NULL) {
        PyErr_SetString(PyExc_TypeError, "Cannot delete the modelParams attribute");
        return -1;
    }

    
    
  return 0;
}

/**getsetters registration struct**/

static PyGetSetDef LALIFOData_getseters[] = {
    {"name", 
     (getter)LALIFOData_getname, (setter)LALIFOData_setname,
     "name",
     NULL},
     {"nullloglikelihood", 
     (getter)LALIFOData_getnullloglikelihood, (setter)LALIFOData_setnullloglikelihood,
     "nullloglikelihood",
     NULL},
     {"loglikelihood", 
     (getter)LALIFOData_getloglikelihood, (setter)LALIFOData_setloglikelihood,
     "loglikelihood",
     NULL},
     {"acceptedloglikelihood", 
     (getter)LALIFOData_getacceptedloglikelihood, (setter)LALIFOData_setacceptedloglikelihood,
     "acceptedloglikelihood",
     NULL},
     {"fPlus", 
     (getter)LALIFOData_getfPlus, (setter)LALIFOData_setfPlus,
     "fPlus",
     NULL},
     {"fCross", 
     (getter)LALIFOData_getfCross, (setter)LALIFOData_setfCross,
     "fCross",
     NULL},
     {"timeshift", 
     (getter)LALIFOData_gettimeshift, (setter)LALIFOData_settimeshift,
     "timeshift",
     NULL},
     {"fLow", 
     (getter)LALIFOData_getfLow, (setter)LALIFOData_setfLow,
     "fLow",
     NULL},
     {"fHigh", 
     (getter)LALIFOData_getfHigh, (setter)LALIFOData_setfHigh,
     "fHigh",
     NULL},
    {"modelParams", 
     (getter)LALIFOData_getmodelParams, (setter)LALIFOData_setmodelParams,
     "modelParams",
     NULL},
    {NULL}  /* Sentinel */
};

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
    LALIFOData_getseters,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)LALIFOData_init,      /* tp_init */
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
init_lalinference(void)
{
    PyObject *m;
    li_LALVariablesType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&li_LALVariablesType) < 0)
        return;
    
    li_LALIFODataType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&li_LALIFODataType) < 0)
        return;
    
    m = Py_InitModule3(MODULE_NAME,module_methods,LIDocString);
    
    Py_INCREF(&li_LALVariablesType);
    PyModule_AddObject(m, "BaseLALVariables", (PyObject *)&li_LALVariablesType);
    Py_INCREF(&li_LALIFODataType);
    PyModule_AddObject(m, "BaseLALIFOData", (PyObject *)&li_LALIFODataType);
}
