/******************************************************************** 
 * A simple C extension module for Python, called "hello"; compile
 * this into a ".so" on python path, import and call hello.message;
 ********************************************************************/

#include <Python.h>
#include <string.h>
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
        "{s:s, s:i, s:i, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:i, s:d}",
        "ifo", event->ifo,
        "end_time", event->end_time.gpsSeconds,
        "end_time_ns", event->end_time.gpsNanoSeconds,
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

/* registration table  */
static struct PyMethodDef metaio_methods[] = {
    {"read_search_summary", read_search_summary, 1},  
    {"read_sngl_burst", read_sngl_burst, 1},  
    {"read_sngl_inspiral", read_sngl_inspiral, 1}, 
    {NULL, NULL}                   /* end of table marker */
};

/* module initializer */
void initmetaio()                       /* called on first import */
{                                      /* name matters if loaded dynamically */
    (void) Py_InitModule("lgen.metaio", metaio_methods);   /* mod name, table ptr */
      /* import_libnumarray(); */
}
