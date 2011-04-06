#include <Python.h>
#include <structmember.h>

#include <numpy/arrayobject.h>

#include <lal/LALInference.h>
#include <lal/Sequence.h>

#include <complex16frequencyseries.h>
#include <complex16timeseries.h>
#include <real8fftplan.h>
#include <real8frequencyseries.h>
#include <real8timeseries.h>

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

static PyObject* LALVariables_add_variable(li_LALVariablesObject *self,PyObject* args,PyObject* kwds){
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

static PyObject* LALVariables_check_variable(li_LALVariablesObject *self,PyObject* args)
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


PyObject* LALVariables_convert_livar_value_to_pyobj(li_LALVariablesObject *pyvars,char* name){
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

static int LALVariables_remove_variable(li_LALVariablesObject *self,PyObject* args){
    PyObject* pyname;
    char* name;
    
    if (! PyArg_ParseTuple(args,"O",&pyname)) return -1;

    name=PyString_AsString(pyname);
    removeVariable(self->vars,name);
    
    return 0;
}

static PyObject* LALVariables_get_variable(li_LALVariablesObject *self,PyObject* args){
    PyObject* pyname=NULL,*returnObj=NULL;
    char* name;
    
    if (! PyArg_ParseTuple(args,"O",&pyname)) return NULL;
    name=PyString_AsString(pyname);
    returnObj=LALVariables_convert_livar_value_to_pyobj(self,name);
    return returnObj;
}

static PyObject* LALVariables_get_variable_name(li_LALVariablesObject *self,PyObject* args){
    char* name;
    long var_idx;
    
    if (! PyArg_ParseTuple(args,"i",&var_idx)) return NULL;
    
    name=getVariableName(self->vars,(int)var_idx);
    
    return PyString_FromString(name);
    
}

static int LALVariables_set_variable(li_LALVariablesObject *self,PyObject* args){
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

static PyObject* LALVariables_get_variable_type(li_LALVariablesObject *self,PyObject* args){
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
    
static PyObject* LALVariables_get_variable_type_by_index(li_LALVariablesObject *self,PyObject* args){
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

static PyObject* LALVariables_get_variable_vary_type(li_LALVariablesObject *self,PyObject* args){
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

static PyObject* LALVariables_get_variable_dimension(li_LALVariablesObject *self){
    long dimension;
    PyObject* pydimension=NULL;
    
    dimension=(long int)getVariableDimension(self->vars);

    pydimension=PyInt_FromLong(dimension);

    return pydimension;
    
}

static PyObject* LALVariables_get_variable_dimension_non_fixed(li_LALVariablesObject *self){
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
    LALIFOData data;
    PyObject* window;
     
} li_LALIFODataObject;


/*Methods*/

 /* Destructor for LALIFOData */
static void LALIFOData_dealloc(li_LALIFODataObject *self)
{
    self->ob_type->tp_free((PyObject *)self);
}

static int LALIFOData_init(li_LALIFODataObject *self, PyObject *args, PyObject *kwds)
{
    //self->data=(LALIFOData*)malloc(sizeof(LALIFOData));
    memset((void*)&self->data,0,sizeof(LALIFOData));
    return 0;
}

/***********getsetters******************/

/*name*/

static PyObject* LALIFOData_getname(li_LALIFODataObject *self, void *closure)
{
    return PyString_FromString(self->data.name);
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
  
  strcpy(self->data.name,PyString_AsString(value));

  return 0;
}

/*modelParams*/
static PyObject* LALIFOData_getmodelParams(li_LALIFODataObject *self, void *closure)
{
    if(self->data.modelParams){
        
        li_LALVariablesObject* newLV=(li_LALVariablesObject*)PyObject_CallFunction((PyObject *)&li_LALVariablesType,NULL);
        copyVariables(self->data.modelParams,newLV->vars);
        
        return (PyObject*)newLV;
    }
    else {
        Py_INCREF(Py_None);
        return Py_None;
    }
}

static int LALIFOData_setmodelParams(li_LALIFODataObject *self, PyObject *value, void *closure)
{   
    if (!value) {
        PyErr_SetString(PyExc_TypeError, "Cannot delete the modelParams attribute");
        return -1;
    }
    
    li_LALVariablesObject* existingLV=(li_LALVariablesObject*)value;
    
    if(!self->data.modelParams){
        self->data.modelParams=(LALVariables*)malloc(sizeof(LALVariables));
    }
    
    copyVariables(existingLV->vars,self->data.modelParams);
    
  return 0;
}

/*dataParams*/
static PyObject* LALIFOData_getdataParams(li_LALIFODataObject *self, void *closure)
{
    if(self->data.dataParams){
        
        li_LALVariablesObject* newLV=(li_LALVariablesObject*)PyObject_CallFunction((PyObject *)&li_LALVariablesType,NULL);
        copyVariables(self->data.dataParams,newLV->vars);
        
        return (PyObject*)newLV;
    }
    else {
        Py_INCREF(Py_None);
        return Py_None;
    }
}

static int LALIFOData_setdataParams(li_LALIFODataObject *self, PyObject *value, void *closure)
{   
    if (!value) {
        PyErr_SetString(PyExc_TypeError, "Cannot delete the dataParams attribute");
        return -1;
    }
    
    li_LALVariablesObject* existingLV=(li_LALVariablesObject*)value;
    
    if(!self->data.dataParams){
        self->data.dataParams=(LALVariables*)malloc(sizeof(LALVariables));
    }
    
    copyVariables(existingLV->vars,self->data.dataParams);
    
  return 0;
}

static PyObject* LALIFOData_getmodelDomain(li_LALIFODataObject *self, void *closure)
{

    if(self->data.modelDomain==timeDomain){
        return PyString_FromString("timeDomain");
    }
    else if(self->data.modelDomain==frequencyDomain){
        return PyString_FromString("freqDomain");
    }
    else return NULL;
}

static int LALIFOData_setmodelDomain(li_LALIFODataObject *self, PyObject *value, void *closure)
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

    char* name=PyString_AsString(value);
    LALDomain var;

    if(!strcmp(name,"timeDomain")){
        var=timeDomain;
    }
    else if(!strcmp(name,"freqDomain")){
        var=frequencyDomain;
    }
    else{
        die("modelDomain must be one of 'freqDomain' or 'timeDomain'");
    }

    return 0;
}

int setREAL8TimeSeriesFromLALIFOData(REAL8TimeSeries* target,PyObject* origin){
    /* require pylal_REAL8TimeSeries */
    if(!PyObject_TypeCheck(origin, &pylal_REAL8TimeSeries_Type)){
        PyErr_SetObject(PyExc_TypeError, origin);
        return -1;
    }

    pylal_REAL8TimeSeries* originvec=(pylal_REAL8TimeSeries*)origin;

    int n=originvec->series->data->length;
    if(n != target->data->length)
        XLALResizeREAL8Sequence(target->data, 0, n);
        
    memcpy(target->data->data,originvec->series->data->data,n * sizeof(*originvec->series->data->data));
    return 0;
}

//PyObject* getREAL8TimeSeriesFromLALIFOData(REAL8TimeSeries* target){

    //pylal_REAL8TimeSeries* new=(pylal_REAL8TimeSeries*)PyObject_CallFunction((PyObject *)&pylal_REAL8TimeSeries_Type,NULL);

    //int n=target->data->length;
    
    //memcpy(new->series,target,sizeof(REAL8TimeSeries));
    //XLALResizeREAL8Sequence(new->series->data, 0, n);
    //memcpy(new->series->data->data,target->data->data,n * sizeof(*target->data->data));
    //return (PyObject*)new;

    
//}

/* timeData */
static int LALIFOData_settimeData(li_LALIFODataObject *self, PyObject *value, void *closure) {return setREAL8TimeSeriesFromLALIFOData(self->data.timeData,value);}
static PyObject* LALIFOData_gettimeData(li_LALIFODataObject *self, void *closure) {return pylal_REAL8TimeSeries_new(self->data.timeData,(PyObject*)self);}//{return getREAL8TimeSeriesFromLALIFOData(self->data.timeData);}
/* timeModelhPlus */
static int LALIFOData_settimeModelhPlus(li_LALIFODataObject *self, PyObject *value, void *closure) {return setREAL8TimeSeriesFromLALIFOData(self->data.timeModelhPlus,value);}
static PyObject* LALIFOData_gettimeModelhPlus(li_LALIFODataObject *self, void *closure) {return pylal_REAL8TimeSeries_new(self->data.timeModelhPlus,(PyObject*)self);}//{return getREAL8TimeSeriesFromLALIFOData(self->data.timeModelhPlus);}
/* timeModelhCross */
static int LALIFOData_settimeModelhCross(li_LALIFODataObject *self, PyObject *value, void *closure) {return setREAL8TimeSeriesFromLALIFOData(self->data.timeModelhCross,value);}
static PyObject* LALIFOData_gettimeModelhCross(li_LALIFODataObject *self, void *closure) {return pylal_REAL8TimeSeries_new(self->data.timeModelhCross,(PyObject*)self);}//{return getREAL8TimeSeriesFromLALIFOData(self->data.timeModelhCross);}
/* whiteTimeData */
static int LALIFOData_setwhiteTimeData(li_LALIFODataObject *self, PyObject *value, void *closure) {return setREAL8TimeSeriesFromLALIFOData(self->data.whiteTimeData,value);}
static PyObject* LALIFOData_getwhiteTimeData(li_LALIFODataObject *self, void *closure) {return pylal_REAL8TimeSeries_new(self->data.whiteTimeData,(PyObject*)self);}//{return getREAL8TimeSeriesFromLALIFOData(self->data.whiteTimeData);}
/* windowedTimeData */
static int LALIFOData_setwindowedTimeData(li_LALIFODataObject *self, PyObject *value, void *closure) {return setREAL8TimeSeriesFromLALIFOData(self->data.windowedTimeData,value);}
static PyObject* LALIFOData_getwindowedTimeData(li_LALIFODataObject *self, void *closure) {return pylal_REAL8TimeSeries_new(self->data.windowedTimeData,(PyObject*)self);}//getREAL8TimeSeriesFromLALIFOData(self->data.windowedTimeData);}
/* timeDomainNoiseWeights */
static int LALIFOData_settimeDomainNoiseWeights(li_LALIFODataObject *self, PyObject *value, void *closure) {return setREAL8TimeSeriesFromLALIFOData(self->data.timeDomainNoiseWeights,value);}
static PyObject* LALIFOData_gettimeDomainNoiseWeights(li_LALIFODataObject *self, void *closure) {return pylal_REAL8TimeSeries_new(self->data.timeDomainNoiseWeights,(PyObject*)self);}//{return getREAL8TimeSeriesFromLALIFOData(self->data.timeDomainNoiseWeights);}

int setREAL8FrequencySeriesFromLALIFOData(REAL8FrequencySeries* target,PyObject* origin){
    /* require pylal_REAL8TimeSeries */
    if(!PyObject_TypeCheck(origin, &pylal_REAL8FrequencySeries_Type)){
        PyErr_SetObject(PyExc_TypeError, origin);
        return -1;
    }

    pylal_REAL8FrequencySeries* originvec=(pylal_REAL8FrequencySeries*)origin;

    int n=originvec->series->data->length;
    if(n != target->data->length)
        XLALResizeREAL8Sequence(target->data, 0, n);
        
    memcpy(target->data->data,originvec->series->data->data,n * sizeof(*originvec->series->data->data));
    return 0;
}

//PyObject* getREAL8FrequencySeriesFromLALIFOData(REAL8FrequencySeries* target){

    //pylal_REAL8FrequencySeries* new=(pylal_REAL8FrequencySeries*)PyObject_CallFunction((PyObject *)&pylal_REAL8FrequencySeries_Type,NULL);

    //int n=target->data->length;
    
    //memcpy(new->series,target,sizeof(REAL8FrequencySeries));
    //XLALResizeREAL8Sequence(new->series->data, 0, n);
    //memcpy(new->series->data->data,target->data->data,n * sizeof(*target->data->data));
    //return (PyObject*)new;

    
//}

/* oneSidedNoisePowerSpectrum */
static int LALIFOData_setoneSidedNoisePowerSpectrum(li_LALIFODataObject *self, PyObject *value, void *closure) {return setREAL8FrequencySeriesFromLALIFOData(self->data.oneSidedNoisePowerSpectrum,value);}
static PyObject* LALIFOData_getoneSidedNoisePowerSpectrum(li_LALIFODataObject *self, void *closure) {return pylal_REAL8FrequencySeries_new(self->data.oneSidedNoisePowerSpectrum,(PyObject*)self);}//{return getREAL8FrequencySeriesFromLALIFOData(self->data.oneSidedNoisePowerSpectrum);}


int setCOMPLEX16FrequencySeriesFromLALIFOData(COMPLEX16FrequencySeries* target,PyObject* origin){
	/* require pylal_COMPLEX16TimeSeries */
    if(!PyObject_TypeCheck(origin, &pylal_COMPLEX16FrequencySeries_Type)){
        PyErr_SetObject(PyExc_TypeError, origin);
        return -1;
    }

    pylal_COMPLEX16FrequencySeries* originvec=(pylal_COMPLEX16FrequencySeries*)origin;

    int n=originvec->series->data->length;
    if(n != target->data->length)
        XLALResizeCOMPLEX16Sequence(target->data, 0, n);
        
    memcpy(target->data->data,originvec->series->data->data,n * sizeof(*originvec->series->data->data));
    return 0;
}
//PyObject* getCOMPLEX16FrequencySeriesFromLALIFOData(COMPLEX16FrequencySeries* target){
	//pylal_COMPLEX16FrequencySeries* new=(pylal_COMPLEX16FrequencySeries*)PyObject_CallFunction((PyObject *)&pylal_COMPLEX16FrequencySeries_Type,NULL);

    //int n=target->data->length;
    
    //memcpy(new->series,target,sizeof(COMPLEX16FrequencySeries));
    //XLALResizeCOMPLEX16Sequence(new->series->data, 0, n);
    //memcpy(new->series->data->data,target->data->data,n * sizeof(*target->data->data));
    //return (PyObject*)new;
//}

/*freqData*/
static int LALIFOData_setfreqData(li_LALIFODataObject *self, PyObject *value, void *closure) {return setCOMPLEX16FrequencySeriesFromLALIFOData(self->data.freqData,value);}
static PyObject* LALIFOData_getfreqData(li_LALIFODataObject *self, void *closure) {return pylal_COMPLEX16FrequencySeries_new(self->data.freqData,(PyObject*)self);}//{return getCOMPLEX16FrequencySeriesFromLALIFOData(self->data.freqData);}
/*freqModelhPlus*/
static int LALIFOData_setfreqModelhPlus(li_LALIFODataObject *self, PyObject *value, void *closure) {return setCOMPLEX16FrequencySeriesFromLALIFOData(self->data.freqModelhPlus,value);}
static PyObject* LALIFOData_getfreqModelhPlus(li_LALIFODataObject *self, void *closure) {return pylal_COMPLEX16FrequencySeries_new(self->data.freqModelhPlus,(PyObject*)self);}//{return getCOMPLEX16FrequencySeriesFromLALIFOData(self->data.freqModelhPlus);}
/*freqModelhCross*/
static int LALIFOData_setfreqModelhCross(li_LALIFODataObject *self, PyObject *value, void *closure) {return setCOMPLEX16FrequencySeriesFromLALIFOData(self->data.freqModelhCross,value);}
static PyObject* LALIFOData_getfreqModelhCross(li_LALIFODataObject *self, void *closure) {return pylal_COMPLEX16FrequencySeries_new(self->data.freqModelhCross,(PyObject*)self);}//{return getCOMPLEX16FrequencySeriesFromLALIFOData(self->data.freqModelhCross);}
/*whiteFreqData*/
static int LALIFOData_setwhiteFreqData(li_LALIFODataObject *self, PyObject *value, void *closure) {return setCOMPLEX16FrequencySeriesFromLALIFOData(self->data.whiteFreqData,value);}
static PyObject* LALIFOData_getwhiteFreqData(li_LALIFODataObject *self, void *closure) {return pylal_COMPLEX16FrequencySeries_new(self->data.whiteFreqData,(PyObject*)self);}//{return getCOMPLEX16FrequencySeriesFromLALIFOData(self->data.whiteFreqData);}

int setCOMPLEX16TimeSeriesFromLALIFOData(COMPLEX16TimeSeries* target,PyObject* origin){
	/* require pylal_COMPLEX16TimeSeries */
    if(!PyObject_TypeCheck(origin, &pylal_COMPLEX16TimeSeries_Type)){
        PyErr_SetObject(PyExc_TypeError, origin);
        return -1;
    }

    pylal_COMPLEX16TimeSeries* originvec=(pylal_COMPLEX16TimeSeries*)origin;

    int n=originvec->series->data->length;
    if(n != target->data->length)
        XLALResizeCOMPLEX16Sequence(target->data, 0, n);
        
    memcpy(target->data->data,originvec->series->data->data,n * sizeof(*originvec->series->data->data));
    return 0;
}
//PyObject* getCOMPLEX16TimeSeriesFromLALIFOData(COMPLEX16TimeSeries* target){
	//pylal_COMPLEX16TimeSeries* new=(pylal_COMPLEX16TimeSeries*)PyObject_CallFunction((PyObject *)&pylal_COMPLEX16TimeSeries_Type,NULL);

    //int n=target->data->length;
    
    //memcpy(new->series,target,sizeof(COMPLEX16TimeSeries));
    //XLALResizeCOMPLEX16Sequence(new->series->data, 0, n);
    //memcpy(new->series->data->data,target->data->data,n * sizeof(*target->data->data));
    //return (PyObject*)new;
//}

/*compTimeData*/
static int LALIFOData_setcompTimeData(li_LALIFODataObject *self, PyObject *value, void *closure) {return setCOMPLEX16TimeSeriesFromLALIFOData(self->data.compTimeData,value);}
static PyObject* LALIFOData_getcompTimeData(li_LALIFODataObject *self, void *closure) {return pylal_COMPLEX16TimeSeries_new(self->data.compTimeData,(PyObject*)self);}

/*compModelData*/
static int LALIFOData_setcompModelData(li_LALIFODataObject *self, PyObject *value, void *closure) {return setCOMPLEX16TimeSeriesFromLALIFOData(self->data.compModelData,value);}
static PyObject* LALIFOData_getcompModelData(li_LALIFODataObject *self, void *closure) {return pylal_COMPLEX16TimeSeries_new(self->data.compModelData,(PyObject*)self);}

int setREAL8FFTPlanFromLALIFOData(REAL8FFTPlan* target,PyObject* origin){
    
    return 0;
}

/*timeToFreqFFTPlan*/
static int LALIFOData_settimeToFreqFFTPlan(li_LALIFODataObject *self, PyObject *value, void *closure) {return setREAL8FFTPlanFromLALIFOData(self->data.timeToFreqFFTPlan,value);}
static PyObject* LALIFOData_gettimeToFreqFFTPlan(li_LALIFODataObject *self, void *closure) {return pylal_REAL8FFTPlan_new(self->data.timeToFreqFFTPlan,(PyObject*)self);}

/*freqToTimeFFTPlan*/
static int LALIFOData_setfreqToTimeFFTPlan(li_LALIFODataObject *self, PyObject *value, void *closure) {return setREAL8FFTPlanFromLALIFOData(self->data.freqToTimeFFTPlan,value);}
static PyObject* LALIFOData_getfreqToTimeFFTPlan(li_LALIFODataObject *self, void *closure) {return pylal_REAL8FFTPlan_new(self->data.freqToTimeFFTPlan,(PyObject*)self);}


/**getsetters registration struct**/

static PyGetSetDef LALIFOData_getseters[] = {
    {"name",(getter)LALIFOData_getname,(setter)LALIFOData_setname,"name",NULL},
    //LALVariables
    {"modelParams",(getter)LALIFOData_getmodelParams,(setter)LALIFOData_setmodelParams,"modelParams",NULL},
    {"dataParams",(getter)LALIFOData_getdataParams,(setter)LALIFOData_setdataParams,"dataParams",NULL},
    //LALDomain ...represented as string
    {"modelDomain",(getter)LALIFOData_getmodelDomain,(setter)LALIFOData_setmodelDomain,"modelDomain",NULL},
    //REAL8TimeSeries
    {"timeData",(getter)LALIFOData_gettimeData,(setter)LALIFOData_settimeData,"timeData",NULL},
    {"timeModelhPlus",(getter)LALIFOData_gettimeModelhPlus,(setter)LALIFOData_settimeModelhPlus,"timeModelhPlus",NULL},
    {"timeModelhCross",(getter)LALIFOData_gettimeModelhCross,(setter)LALIFOData_settimeModelhCross,"timeModelhCross",NULL},
    {"whiteTimeData",(getter)LALIFOData_getwhiteTimeData,(setter)LALIFOData_setwhiteTimeData,"whiteTimeData",NULL},
    {"windowedTimeData",(getter)LALIFOData_getwindowedTimeData,(setter)LALIFOData_setwindowedTimeData,"windowedTimeData",NULL},
    {"timeDomainNoiseWeights",(getter)LALIFOData_gettimeDomainNoiseWeights,(setter)LALIFOData_settimeDomainNoiseWeights,"timeDomainNoiseWeights",NULL},
    //COMPLEX16FrequencySeries
    {"freqData",(getter)LALIFOData_getfreqData,(setter)LALIFOData_setfreqData,"freqData",NULL},
    {"freqModelhPlus",(getter)LALIFOData_getfreqModelhPlus,(setter)LALIFOData_setfreqModelhPlus,"freqModelhPlus",NULL},
    {"freqModelhCross",(getter)LALIFOData_getfreqModelhCross,(setter)LALIFOData_setfreqModelhCross,"freqModelhCross",NULL},
    {"whiteFreqData",(getter)LALIFOData_getwhiteFreqData,(setter)LALIFOData_setwhiteFreqData,"whiteFreqData",NULL},
    //COMPLEX16TimeSeries
    {"compTimeData",(getter)LALIFOData_getcompTimeData,(setter)LALIFOData_setcompTimeData,"compTimeData",NULL},
    {"compModelData",(getter)LALIFOData_getcompModelData,(setter)LALIFOData_setcompModelData,"compModelData",NULL},
    //REAL8FrequencySeries
    {"oneSidedNoisePowerSpectrum",(getter)LALIFOData_getoneSidedNoisePowerSpectrum,(setter)LALIFOData_setoneSidedNoisePowerSpectrum,"oneSidedNoisePowerSpectrum",NULL},
    //REAL8FFTPlan
    {"timeToFreqFFTPlan",(getter)LALIFOData_gettimeToFreqFFTPlan,(setter)LALIFOData_settimeToFreqFFTPlan,"timeToFreqFFTPlan",NULL},
    {"freqToTimeFFTPlan",(getter)LALIFOData_getfreqToTimeFFTPlan,(setter)LALIFOData_setfreqToTimeFFTPlan,"freqToTimeFFTPlan",NULL},
    {NULL}  /* Sentinel */
};

static struct PyMemberDef LALIFOData_members[] = {
    //REAL8's
    {"nullloglikelihood", T_DOUBLE, offsetof(li_LALIFODataObject, data.nullloglikelihood), 0, "nullloglikelihood"},
    {"loglikelihood", T_DOUBLE, offsetof(li_LALIFODataObject, data.loglikelihood), 0, "loglikelihood"},
    {"acceptedloglikelihood", T_DOUBLE, offsetof(li_LALIFODataObject, data.acceptedloglikelihood), 0, "acceptedloglikelihood"},
    {"fPlus", T_DOUBLE, offsetof(li_LALIFODataObject, data.fPlus), 0, "fPlus"},
    {"fCross", T_DOUBLE, offsetof(li_LALIFODataObject, data.fCross), 0, "fCross"},
    {"timeshift", T_DOUBLE, offsetof(li_LALIFODataObject, data.timeshift), 0, "timeshift"},
    {"fLow", T_DOUBLE, offsetof(li_LALIFODataObject, data.fLow), 0, "fLow"},
    {"fHigh", T_DOUBLE, offsetof(li_LALIFODataObject, data.acceptedloglikelihood), 0, "fHigh"},
    {"acceptedloglikelihood", T_DOUBLE, offsetof(li_LALIFODataObject, data.acceptedloglikelihood), 0, "acceptedloglikelihood"},
	
    {"window", T_OBJECT, offsetof(li_LALIFODataObject,window), 0, "window"},
    
    {NULL,}
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
    LALIFOData_members,             /* tp_members */
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

    import_array();
    pylal_complex16frequencyseries_import();
    pylal_complex16timeseries_import();
    pylal_real8fftplan_import();
    pylal_real8frequencyseries_import();
    pylal_real8timeseries_import();
    
    Py_INCREF(&li_LALVariablesType);
    PyModule_AddObject(m, "BaseLALVariables", (PyObject *)&li_LALVariablesType);
    Py_INCREF(&li_LALIFODataType);
    PyModule_AddObject(m, "BaseLALIFOData", (PyObject *)&li_LALIFODataType);
}
