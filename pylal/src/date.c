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
gps_to_utc(PyObject *self, PyObject *args)   /* self unused in modules */
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

/* registration table  */
static struct PyMethodDef date_methods[] = {
    {"gps_to_utc", gps_to_utc, 1},  
    {NULL, NULL}                   /* end of table marker */
};

/* module initializer */
void initdate()                       /* called on first import */
{                                      /* name matters if loaded dynamically */
    (void) Py_InitModule("lgen.date", date_methods);   /* mod name, table ptr */
      import_libnumarray();
}
