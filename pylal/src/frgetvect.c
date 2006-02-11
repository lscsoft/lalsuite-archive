/*
    Nick Fotopoulos (nvf@mit.edu)
    Created: 2005-05-24
    Python port of frgetvect (based on Matlab frgetvect)

    Requires: Numarray or Numeric (change #include line appropriately)
    To do:  Check inputs, raise appropriate exceptions, test complex-valued
            frames, use keyword arguments, allow reading with time offset
            and without knowing t0 (perhaps useful if you know the first 5
            seconds of every run are noisy or something)
*/

#include "Python.h"
#include "numarray/arrayobject.h"
#include "FrameL.h"

static PyObject *frgetvect(PyObject *self, PyObject *args)
{
    import_libnumeric();

    int nrhs = 5;  // I don't actually use this at present.
    PyArrayObject *out1, *out2, *out3;
    double out4;
    char *out5, *out6, *out7;

    struct FrFile *iFile;
    struct FrVect *vect;	
    long  debugLvl, i, nData, utc;
    double *data, dt, df, t0, duration;
    char  *fileName, *vectName;

    int shape[1];
    int ok;

    /*--------------- unpack arguments --------------------*/
    /* Give bad default arguments as a way of counting how many
       arguments were passed to this function. */
    t0 = -1.;
    duration = -1.;
    debugLvl = -1;

    /* The | in the format string indicates the next arguments are
       optional.  They are simply not assigned anything. */
    ok = PyArg_ParseTuple(args, "ss|ddi", &fileName, &vectName,
                           &t0, &duration, &debugLvl);

    // Handle defaults
    if (debugLvl < 0)
    {
        debugLvl = 0;
        nrhs = 4;
    }
    if (duration < 0.)
    {
        duration = 1.;
        nrhs = 3;
    }
    if (t0 < 0.)
    {
        t0 = 0.;
        nrhs = 2;
    }

    FrLibSetLvl(debugLvl);

    /*-------------- open file --------------------------*/

    if (debugLvl > 0)
    {
        const char *msg = "Opening %s for channel %s (t0=%.2f, duration=%.2f).\n";
        printf(msg, fileName, vectName, t0, duration);
    }

    iFile = FrFileINew(fileName);
    if (iFile == NULL)
    {
        printf("%s\n", FrErrorGetHistory());
        return;
    }

    if (debugLvl > 0)
    {
        const char *msg = "Opened %s for channel %s!\n";
        printf(msg, fileName, vectName);
    }

    /*-------------- get vector --------------------------*/

    vect = FrFileIGetV(iFile, vectName, t0, duration);

    if(debugLvl > 0) FrVectDump(vect, stdout, debugLvl);
    if(vect == NULL)
    {
        printf("In file %s, vector not found: %s\n",fileName,vectName);
        FrFileIEnd(iFile);
        return;
    }

    if (debugLvl > 0)
    {
        const char *msg = "Extracted channel %s successfully!\n";
        printf(msg, vectName);
    }

    nData = vect->nData;
    dt = vect->dx[0];
    if (0.0 == dt)
    {
        printf("dt==0\n");
        df = 0.0;
    }
    else
    {
        df = 1.0/(dt*nData);
    }

    /*-------- copy data ------*/

    /* Numarray does not support chars. */

    /* Assign data first.  Copy to numarray afterwards. */

    if (vect->type == FR_VECT_8C || vect->type == FR_VECT_16C) {
        data = (double *)malloc(2*nData*sizeof(double));
    } else {
        data = (double *)malloc(nData*sizeof(double));
    }
    if (data==NULL) {
        printf("Unable to allocate space for data.\n");
        return;
    }

    if(vect->type == FR_VECT_2S)
    {for(i=0; i<nData; i++) {data[i] = vect->dataS[i];}}
    else if(vect->type == FR_VECT_4S)
    {for(i=0; i<nData; i++) {data[i] = vect->dataI[i];}}
    else if(vect->type == FR_VECT_8S)
    {for(i=0; i<nData; i++) {data[i] = vect->dataL[i];}}
    else if(vect->type == FR_VECT_1U)
    {for(i=0; i<nData; i++) {data[i] = vect->dataU[i];}}
    else if(vect->type == FR_VECT_2U)
    {for(i=0; i<nData; i++) {data[i] = vect->dataUS[i];}}
    else if(vect->type == FR_VECT_4U)
    {for(i=0; i<nData; i++) {data[i] = vect->dataUI[i];}}
    else if(vect->type == FR_VECT_8U)
    {for(i=0; i<nData; i++) {data[i] = vect->dataUL[i];}}
    else if(vect->type == FR_VECT_4R)
    {for(i=0; i<nData; i++) {data[i] = vect->dataF[i];}}
    else if(vect->type == FR_VECT_8R)
    {for(i=0; i<nData; i++) {data[i] = vect->dataD[i];}}
    // Note the 2*nData in the for loop for complex types
    else if(vect->type == FR_VECT_8C)
    {for(i=0; i<2*nData; i++) {data[i] = vect->dataF[i];}}
    else if(vect->type == FR_VECT_16C)
    {for(i=0; i<2*nData; i++) {data[i] = vect->dataD[i];}}
    else
    {
        printf("No numarray type for this channel");
        FrVectFree(vect);
        FrFileIEnd(iFile);
        return;
    }

    shape[0] = nData;

    if (vect->type == FR_VECT_8C || vect->type == FR_VECT_16C) {
        out1 = (PyArrayObject *)
                PyArray_FromDimsAndData(1,shape,PyArray_CDOUBLE,(char *)data);
    } else {
        out1 = (PyArrayObject *)
                PyArray_FromDimsAndData(1,shape,PyArray_DOUBLE,(char *)data);
    }

    /*------------- fill time and frequency array --------*/

    // output2 = x-axis values relative to first data point.
    // Usually time, but could be frequency in the case of a frequency
    // series
    data = malloc(nData*sizeof(double));
    for(i=0; i<nData; i++) {
      data[i] = vect->startX[0]+(double)i*dt;
    }
    shape[0] = nData;
    out2 = (PyArrayObject *)
            PyArray_FromDimsAndData(1,shape,PyArray_DOUBLE,(char *)data);

    // output3 = frequency values in the case of time series
    data = malloc(nData/2*sizeof(double));
    for (i=0;i<nData/2;i++) {
        data[i] = (double)i*df;
    }
    shape[0] = nData/2;
    out3 = (PyArrayObject *)
            PyArray_FromDimsAndData(1,shape,PyArray_DOUBLE,(char *)data);

    // output4 = gps start time
    out4 = (double)vect->GTime;

    // output5 = gps start time as a string
    utc = vect->GTime - vect->ULeapS + FRGPSTAI;
    out5 = malloc(200*sizeof(char));
    sprintf(out5,"Starting GPS time:%.1f UTC=%s",vect->GTime,FrStrGTime(utc));

    // output6 = unitX as a string
    // I don't really understand why vect->unitX is a char**. (see FrVect.h)
    out6 = malloc((strlen(vect->unitX[0])+1)*sizeof(char));
    strcpy(out6, vect->unitX[0]);

    // output7 = unitY as a string
    out7 = malloc((strlen(vect->unitY)+1)*sizeof(char));
    strcpy(out7, vect->unitY);

    /*------------- clean up -----------------------------*/
    /* One should not free(data).  That memory is in use! */

    FrVectFree(vect);
    FrFileIEnd(iFile);

    return Py_BuildValue("(OOOdsss)",out1,out2,out3,out4,out5,out6,out7);
}


static PyMethodDef frgetvectMethods[] = {
    {"frgetvect", frgetvect, METH_VARARGS, "Docstring goes here."},
        {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC initfrgetvect(void){
    (void) Py_InitModule("frgetvect", frgetvectMethods);
}
