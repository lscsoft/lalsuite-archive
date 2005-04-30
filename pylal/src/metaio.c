/******************************************************************** 
 * A simple C extension module for Python, called "hello"; compile
 * this into a ".so" on python path, import and call hello.message;
 ********************************************************************/

#include <Python.h>
#include <string.h>
#include <numarray/numarray.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOLwXMLHeaders.h>
#include <lal/LIGOLwXMLRead.h>

/* module functions */
static PyObject *                                 /* returns object */
read_sngl_burst(PyObject *self, PyObject *args)   /* self unused in modules */
{                                                 /* args from python call */
    char *fromPython;
    SnglBurstTable *eventHead=NULL;
    SnglBurstTable *event=NULL;
    int j, m=0, n=0;
    double *result;
    
    if (! PyArg_Parse(args, "(s)", &fromPython))  /* convert Python -> C */
        return NULL;                              /* null=raise exception */
    else {
        eventHead = XLALSnglBurstTableFromLIGOLw ( fromPython );
        m = XLALCountSnglBurst ( eventHead );
        n=9;
        result = (double *) malloc( m*n*sizeof(double) );
        for ( j=0, event = eventHead; event ; j++, event = event->next )
        {
          result[j*n] = (double)event->start_time.gpsSeconds;
          result[j*n+1] = (double)event->start_time.gpsNanoSeconds;
          result[j*n+2] = (double)event->peak_time.gpsSeconds;
          result[j*n+3] = (double)event->peak_time.gpsNanoSeconds;
          result[j*n+4] = (double)event->duration;
          result[j*n+5] = (double)event->central_freq;
          result[j*n+6] = (double)event->bandwidth;
          result[j*n+7] = (double)event->snr;
          result[j*n+8] = (double)event->confidence;
        }
    }
    return NA_NewArray((void *)result, tFloat64, 2, m, n);
}

/* module functions */
static PyObject *                                 /* returns object */
read_sngl_inspiral(PyObject *self, PyObject *args)           /* self unused in modules */
{                                                 /* args from python call */
    char *fromPython;
    SnglInspiralTable *eventHead=NULL;
    SnglInspiralTable *event=NULL;
    int j, m=0, n=0;
    int startEvent = 0, stopEvent = -1;
    double *result;
    
    if (! PyArg_Parse(args, "(s)", &fromPython))  /* convert Python -> C */
        return NULL;                              /* null=raise exception */
    else {
        m = LALSnglInspiralTableFromLIGOLw ( &eventHead, fromPython, 
              startEvent, stopEvent);
        n=12;
        result = (double *) malloc( m*n*sizeof(double) );
        for ( j=0, event = eventHead; event ; j++, event = event->next )
        {
          result[j*n+0] = (double)event->end_time.gpsSeconds;
          result[j*n+1] = (double)event->end_time.gpsNanoSeconds;
          result[j*n+2] = (double)event->eff_distance;
          result[j*n+3] = (double)event->coa_phase;
          result[j*n+4] = (double)event->mass1;
          result[j*n+5] = (double)event->mass2;
          result[j*n+6] = (double)event->mchirp;
          result[j*n+7] = (double)event->eta;
          result[j*n+8] = (double)event->snr;
          result[j*n+9] = (double)event->chisq;
          result[j*n+10] = (double)event->chisq_dof;
          result[j*n+11] = (double)event->sigmasq;
        }
        return NA_NewArray((void *)result, tFloat64, 2, m, n);
    }
}

/* module functions */
static PyObject *                                 /* returns object */
read_sim_inspiral(PyObject *self, PyObject *args)           /* self unused in modules */
{                                                 /* args from python call */
    char *fromPython;
    SimInspiralTable *eventHead=NULL;
    SimInspiralTable *event=NULL;
    int j, m=0, n=0;
    int startEvent = 0, stopEvent = 0;
    double *result;
    
    if (! PyArg_Parse(args, "(s)", &fromPython))  /* convert Python -> C */
        return NULL;                              /* null=raise exception */
    else {
        m = SimInspiralTableFromLIGOLw ( &eventHead, fromPython, 
              startEvent, stopEvent);
        n=12;
        result = (double *) malloc( m*n*sizeof(double) );
        for ( j=0, event = eventHead; event ; j++, event = event->next )
        {
          result[j*n+0] = (double)event->geocent_end_time.gpsSeconds;
          result[j*n+1] = (double)event->geocent_end_time.gpsNanoSeconds;
          result[j*n+2] = (double)event->distance;
          result[j*n+3] = (double)event->f_final;
          result[j*n+4] = (double)event->mass1;
          result[j*n+5] = (double)event->mass2;
          result[j*n+6] = (double)event->mchirp;
          result[j*n+7] = (double)event->eta;
          result[j*n+8] = (double)event->eff_dist_h;
          result[j*n+9] = (double)event->eff_dist_l;
          result[j*n+10] = (double)event->eff_dist_g;
          result[j*n+11] = (double)event->eff_dist_v;
        }
        return NA_NewArray((void *)result, tFloat64, 2, m, n);
    }
}

/* registration table  */
static struct PyMethodDef metaio_methods[] = {
    {"read_sngl_burst", read_sngl_burst, 1},  
    {"read_sngl_inspiral", read_sngl_inspiral, 1}, 
    {"read_sim_inspiral", read_sim_inspiral, 1}, 
    {NULL, NULL}                   /* end of table marker */
};

/* module initializer */
void initmetaio()                       /* called on first import */
{                                      /* name matters if loaded dynamically */
    (void) Py_InitModule("lgen.metaio", metaio_methods);   /* mod name, table ptr */
      import_libnumarray();
}
