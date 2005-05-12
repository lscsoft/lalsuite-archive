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

#if defined(Py_DEBUG) || defined(DEBUG)
extern void _Py_CountReferences(FILE *);
#define CURIOUS(x) { fprintf(stderr,__FILE__":%d ",__LINE__); x; }
#else
#define CURIOUS(x)
#endif
#define MARKER()     CURIOUS(fprintf(stderr,"\n"))
#define DESCRIBE(x)     CURIOUS(fprintf(stderr,"  " #x "=%d\n",x))
#define DESCRIBE_HEX(x)     CURIOUS(fprintf(stderr,"  " #x "=%08x\n",x))
#define COUNTREFS()     CURIOUS(_Py_CountReferences(stderr))


/* module functions */
static PyObject *                                 /* returns object */
read_search_summary(PyObject *self, PyObject *args)           /* self unused in modules */
{ 
  int j, m=0, n=0, nelement=0, len;
  SearchSummaryTable *eventHead=NULL;
  SearchSummaryTable **addevent=&eventHead;
  SearchSummaryTable *event=NULL;
  PyObject *fromPython;
  PyObject *outlist;
  PyObject *tmpvalue;

  if (! PyArg_ParseTuple(args, "O", &fromPython))  /* convert Python -> C */
    return NULL;                              /* null=raise exception */

  if (! PyList_Check(fromPython) )
    return NULL;

  len = PyList_Size(fromPython);
  for(n=0;n<len;n++)
  {
    m = SearchSummaryTableFromLIGOLw ( addevent,
        PyString_AsString(PyList_GetItem(fromPython, n)));

    while(*addevent)
      addevent = &(*addevent)->next;

    nelement += m;
  }

  outlist = PyList_New(nelement);
  for ( j=0, event = eventHead; event ; j++, event = event->next )
  {
    tmpvalue = Py_BuildValue(
        "{s:s, s:d, s:d, s:d, s:d, s:i, s:i, s:s}",
        "comment", event->comment,
        "in_start_time", (double)event->in_start_time.gpsSeconds,
        "in_end_time", (double)event->in_end_time.gpsSeconds,
        "out_start_time", (double)event->out_start_time.gpsSeconds,
        "out_end_time", (double)event->out_end_time.gpsSeconds,
        "nevents", event->nevents,
        "nnodes", event->nnodes,
        "ifos", event->ifos);
    PyList_SetItem(outlist, j, tmpvalue);
  }

  while(eventHead) {
    event = eventHead;
    eventHead = eventHead->next;
    LALFree(event);
  }

  return outlist;
}

/* module functions */
static PyObject *                                 /* returns object */
read_sngl_burst(PyObject *self, PyObject *args)           /* self unused in modules */
{ 
  SnglBurstTable *eventHead=NULL;
  SnglBurstTable **addevent=&eventHead;
  SnglBurstTable *event=NULL;
  PyObject *fromPython;
  int j, m=0, n=0, nelement=0, len;
  int startEvent = 0, stopEvent = -1;
  PyObject *outlist;
  PyObject *tmpvalue;

  if (! PyArg_ParseTuple(args, "O", &fromPython))  /* convert Python -> C */
    return NULL;                              /* null=raise exception */

  if (! PyList_Check(fromPython) )
    return NULL;

  len = PyList_Size(fromPython);
  for(n=0;n<len;n++)
  {
    *addevent = XLALSnglBurstTableFromLIGOLw (
        PyString_AsString(PyList_GetItem(fromPython, n)) );

    while(*addevent)
      addevent = &(*addevent)->next;

  }
  nelement = XLALCountSnglBurst ( eventHead );

  outlist = PyList_New(nelement);
  for ( j=0, event = eventHead; event ; j++, event = event->next )
  {
    tmpvalue = Py_BuildValue(
        "{s:s, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d}",
        "ifo", event->ifo,
        "start_time", (double)event->start_time.gpsSeconds,
        "peak_time", (double)event->peak_time.gpsSeconds,
        "duration", event->duration,
        "central_freq", event->central_freq,
        "bandwidth", event->bandwidth,
        "amplitude", event->amplitude,
        "snr", event->snr,
        "confidence", event->confidence);
    PyList_SetItem(outlist, j, tmpvalue);
  }

  while(eventHead) {
    event = eventHead;
    eventHead = eventHead->next;
    LALFree(event);
  }

  return outlist;
}

/* module functions */
static PyObject *                                 /* returns object */
read_sngl_inspiral(PyObject *self, PyObject *args)           /* self unused in modules */
{ 
  SnglInspiralTable *eventHead=NULL;
  SnglInspiralTable **addevent=&eventHead;
  SnglInspiralTable *event=NULL;
  PyObject *fromPython;
  int j, m=0, n=0, nelement=0, len;
  int startEvent = 0, stopEvent = -1;
  PyObject *outlist;
  PyObject *tmpvalue;

  if (! PyArg_ParseTuple(args, "O", &fromPython))  /* convert Python -> C */
    return NULL;                              /* null=raise exception */

  if (! PyList_Check(fromPython) )
    return NULL;

  len = PyList_Size(fromPython);
  for(n=0;n<len;n++)
  {
    m = LALSnglInspiralTableFromLIGOLw ( addevent,
        PyString_AsString(PyList_GetItem(fromPython, n)), 
        startEvent, stopEvent);

    while(*addevent)
      addevent = &(*addevent)->next;

    nelement += m;
  }

  outlist = PyList_New(nelement);
  for ( j=0, event = eventHead; event ; j++, event = event->next )
  {
    tmpvalue = Py_BuildValue(
        "{s:s, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:i, s:d}",
        "ifo", event->ifo,
        "eff_distance", event->eff_distance,
        "coa_phase", event->coa_phase,
        "mass1", event->mass1,
        "mass2", event->mass2,
        "mchirp", event->mchirp,
        "mtotal", event->mtotal,
        "eta", event->eta,
        "snr", event->snr,
        "chisq", event->chisq,
        "chisq_dof", event->chisq_dof,
        "sigmasq", event->sigmasq);
    PyList_SetItem(outlist, j, tmpvalue);
  }

  while(eventHead) {
    event = eventHead;
    eventHead = eventHead->next;
    LALFree(event);
  }

  return outlist;
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
    {"read_search_summary", read_search_summary, 1},  
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
