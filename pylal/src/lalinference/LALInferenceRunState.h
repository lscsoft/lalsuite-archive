//LALInferenceRunState.h

static PyTypeObject *_li_LALInferenceRunState_Type = NULL;
#define li_LALInferenceRunState_Type (*_li_LALInferenceRunState_Type)

typedef struct {
    PyObject_HEAD
    PyObject* owner;
    /* Type-specific fields go here */
    LALInferenceRunState *    state;
    //ProcessParams*
    PyObject* commandLine;
    //LALAlgorithm*
    PyObject* algorithm;
    //LALEvolveOneStepFunction*
    PyObject* evolve;
    //LALPriorFunction*
    PyObject* prior;
    //LALLikelihoodFunction*
    PyObject* likelihood;
    //LALProposalFunction*
    PyObject* proposal;
    //LALTemplateFunction*
    PyObject* template;
    //struct tagLALIFOData*
    PyObject* data;
    //LALVariables*
    PyObject* currentParams;
    PyObject* priorArgs;
    PyObject* proposalArgs;
    PyObject* algorithmParams; /* Parameters which control the running of the algorithm */
    //LALVariables**
    PyObject* livePoints; /* Array of live points for Nested Sampling */
    PyObject* differentialPoints;
    //size_t
    PyObject* differentialPointsLength;
    //gsl_rng*
    PyObject* GSLrandom; 
} li_LALInferenceRunState;

static PyObject* li_LALInferenceRunState_new(LALInferenceRunState *state, PyObject *owner){
    PyObject *empty_tuple = PyTuple_New(0);
    li_LALInferenceRunState *obj = (li_LALInferenceRunState *) PyType_GenericNew(&li_LALInferenceRunState_Type,empty_tuple,NULL);
    Py_DECREF(empty_tuple);
    
    if(!obj) {
        if(!owner)
            return NULL;
    }
    if(owner)
        Py_INCREF(owner);
    obj->owner = owner;
    
    obj->state = state;
    return (PyObject *) obj;
}
