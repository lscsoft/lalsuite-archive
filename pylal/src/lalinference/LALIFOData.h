//LALIFOData.h

static PyTypeObject *_li_LALIFOData_Type = NULL;
#define li_LALIFOData_Type (*_li_LALIFOData_Type)

typedef struct {
    PyObject_HEAD
    PyObject* owner;
    /* Type-specific fields go here */
    LALIFOData *    data;
    PyObject* name;
    //LALVariables
    PyObject* modelParams;
    PyObject* dataParams;
    //LALDomain ...represented as string
    PyObject* modelDomain;
    //REAL8TimeSeries
    PyObject* timeData;
    PyObject* timeModelhPlus;
    PyObject* timeModelhCross;
    PyObject* whiteTimeData;
    PyObject* windowedTimeData;
    PyObject* timeDomainNoiseWeights;
    //COMPLEX16FrequencySeries
    PyObject* freqData;
    PyObject* freqModelhPlus;
    PyObject* freqModelhCross;
    PyObject* whiteFreqData;
    //COMPLEX16TimeSeries
    PyObject* compTimeData;
    PyObject* compModelData;
    //REAL8FrequencySeries
    PyObject* oneSidedNoisePowerSpectrum;
    //REAL8FFTPlan
    PyObject* timeToFreqFFTPlan;
    PyObject* freqToTimeFFTPlan;
    PyObject* epoch;
    //REAL8Window
    PyObject* window;
    //LALIFOData
    PyObject* next;
    //BarycenterInput
    PyObject* bary;
    //EphemerisData
    PyObject* ephem;
    
     
} li_LALIFOData;

static PyObject* li_LALIFOData_new(LALIFOData *data, PyObject *owner){
    PyObject *empty_tuple = PyTuple_New(0);
    li_LALIFOData *obj = (li_LALIFOData *) PyType_GenericNew(&li_LALIFOData_Type,empty_tuple,NULL);
    Py_DECREF(empty_tuple);
    
    if(!obj) {
        if(!owner)
            return NULL;
    }
    if(owner)
        Py_INCREF(owner);
    obj->owner = owner;
    
    obj->data = data;
    return (PyObject *) obj;
}
