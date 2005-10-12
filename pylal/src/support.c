/******************************************************************** 
 * A C extension module for Python, called "support"
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

/******************************************************************** 
 * Process Params Table Reading Function
 ********************************************************************/
static PyObject *   
read_process_params(PyObject *self, PyObject *args)
{ 
  ProcessParamsTable *eventHead=NULL;
  ProcessParamsTable **addevent=&eventHead;
  ProcessParamsTable *event=NULL;
  PyObject *fromPython;
  int j, m=0, n=0, nelement=0, len;
  int startEvent = 0, stopEvent = -1;
  PyObject *outlist;
  PyObject *tmpvalue;

  if (! PyArg_ParseTuple(args, "O", &fromPython))
    return NULL;                             

  if (! PyList_Check(fromPython) )
    return NULL;

  len = PyList_Size(fromPython);
  for(n=0;n<len;n++)
  {
    *addevent = XLALProcessParamsTableFromLIGOLw (
        PyString_AsString(PyList_GetItem(fromPython, n)) );

    while(*addevent)
      addevent = &(*addevent)->next;

  }
  nelement = XLALCountProcessParamsTable ( eventHead );

  outlist = PyList_New(nelement);
  for ( j=0, event = eventHead; event ; j++, event = event->next )
  {
    tmpvalue = Py_BuildValue(
        "{s:s, s:s, s:s, s:s}",
        "program", event->program,
        "param", event->param,
        "type", event->type,
        "value", event->value
        );
    PyList_SetItem(outlist, j, tmpvalue);
  }

  while(eventHead) {
    event = eventHead;
    eventHead = eventHead->next;
    LALFree(event);
  }

  return outlist;
}

/******************************************************************** 
 * Process Table Reading Function
 ********************************************************************/
static PyObject *   
read_process(PyObject *self, PyObject *args)
{ 
  ProcessTable *eventHead=NULL;
  ProcessTable **addevent=&eventHead;
  ProcessTable *event=NULL;
  PyObject *fromPython;
  int j, m=0, n=0, nelement=0, len;
  int startEvent = 0, stopEvent = -1;
  PyObject *outlist;
  PyObject *tmpvalue;

  if (! PyArg_ParseTuple(args, "O", &fromPython))
    return NULL;                             

  if (! PyList_Check(fromPython) )
    return NULL;

  len = PyList_Size(fromPython);
  for(n=0;n<len;n++)
  {
    *addevent = XLALProcessTableFromLIGOLw (
        PyString_AsString(PyList_GetItem(fromPython, n)) );

    while(*addevent)
      addevent = &(*addevent)->next;

  }
  nelement = XLALCountProcessTable ( eventHead );

  outlist = PyList_New(nelement);
  for ( j=0, event = eventHead; event ; j++, event = event->next )
  {
    tmpvalue = Py_BuildValue(
        "{s:s, s:s, s:s, s:d, s:s, s:d, s:s, s:s, s:d, s:d, s:d, s:s, s:s}",
        "program", event->program,
        "version", event->version,
        "cvs_repository", event->version,
        "cvs_entry_time", (double)event->cvs_entry_time.gpsSeconds,
        "comment", event->comment,
        "is_online", (double)event->is_online,
        "node", event->node,
        "username", event->username,
        "start_time", (double)event->start_time.gpsSeconds,
        "end_time", (double)event->end_time.gpsSeconds,
        "jobid", (double)event->jobid,
        "domain", event->domain,
        "ifos", event->ifos
        );
    PyList_SetItem(outlist, j, tmpvalue);
  }

  while(eventHead) {
    event = eventHead;
    eventHead = eventHead->next;
    LALFree(event);
  }

  return outlist;
}

/******************************************************************** 
 * Search Summary Reading Function
 ********************************************************************/
static PyObject *                   
read_search_summary(PyObject *self, PyObject *args)
{ 
  int j, m=0, n=0, nelement=0, len;
  SearchSummaryTable *eventHead=NULL;
  SearchSummaryTable **addevent=&eventHead;
  SearchSummaryTable *event=NULL;
  PyObject *fromPython;
  PyObject *outlist;
  PyObject *tmpvalue;

  if (! PyArg_ParseTuple(args, "O", &fromPython)) 
    return NULL;                             

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

/******************************************************************** 
 * SummValue Reading Function
 ********************************************************************/
static PyObject *     
read_summ_value(PyObject *self, PyObject *args)
{ 
  int j, m=0, n=0, nelement=0, len;
  SummValueTable *eventHead=NULL;
  SummValueTable **addevent=&eventHead;
  SummValueTable *event=NULL;
  PyObject *fromPython;
  PyObject *outlist;
  PyObject *tmpvalue;

  if (! PyArg_ParseTuple(args, "O", &fromPython))
    return NULL;                           

  if (! PyList_Check(fromPython) )
    return NULL;

  len = PyList_Size(fromPython);
  for(n=0;n<len;n++)
  {
    m = SummValueTableFromLIGOLw ( addevent,
        PyString_AsString(PyList_GetItem(fromPython, n)));

    while(*addevent)
      addevent = &(*addevent)->next;

    nelement += m;
  }

  outlist = PyList_New(nelement);
  for ( j=0, event = eventHead; event ; j++, event = event->next )
  {
    tmpvalue = Py_BuildValue(
        "{s:s, s:d, s:d, s:s, s:s, s:d, s:s}",
        "program", event->program,
        "start_time", (double)event->start_time.gpsSeconds,
        "end_time", (double)event->end_time.gpsSeconds,
        "ifo", event->ifo,
        "name", event->name,
        "value", event->value,
        "comment", event->comment);
    PyList_SetItem(outlist, j, tmpvalue);
  }

  while(eventHead) {
    event = eventHead;
    eventHead = eventHead->next;
    LALFree(event);
  }

  return outlist;
}

/******************************************************************** 
 * Single Burst Reading Function
 ********************************************************************/
static PyObject *   
read_sngl_burst(PyObject *self, PyObject *args)
{ 
  SnglBurstTable *eventHead;
  SnglBurstTable *event;
  PyObject *fromPython;
  PyObject *outlist;
  int n, len;

  if(!PyArg_ParseTuple(args, "O", &fromPython))
    return NULL;

  if(!PyList_Check(fromPython)) {
    PyErr_SetString(PyExc_TypeError, "argument is not a list");
    return NULL;
  }

  len = PyList_Size(fromPython);
  outlist = PyList_New(0);
  for(n = 0; n < len; n++)
  {
    eventHead = XLALSnglBurstTableFromLIGOLw(PyString_AsString(PyList_GetItem(fromPython, n)));

    while(eventHead)
    {
      event = eventHead;
      PyList_Append(outlist, Py_BuildValue(
          "{s:s, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d}",
          "ifo", event->ifo,
          "start_time", event->start_time.gpsSeconds + event->start_time.gpsNanoSeconds * (double) 1e-9,
          "peak_time", event->peak_time.gpsSeconds + event->peak_time.gpsNanoSeconds * (double) 1e-9,
          "duration", event->duration,
          "central_freq", event->central_freq,
          "bandwidth", event->bandwidth,
          "amplitude", event->amplitude,
          "snr", event->snr,
          "confidence", event->confidence));
      eventHead = eventHead->next;
      LALFree(event);
    }
  }

  return outlist;
}

/******************************************************************** 
 * Single Inspiral Reading Function
 ********************************************************************/
static PyObject *   
read_sngl_inspiral(PyObject *self, PyObject *args) 
{ 
  SnglInspiralTable *eventHead=NULL;
  SnglInspiralTable **addevent=&eventHead;
  SnglInspiralTable *event=NULL;
  PyObject *fromPython;
  int j, m=0, n=0, nelement=0, len;
  int startEvent = 0, stopEvent = -1;
  PyObject *outlist;
  PyObject *tmpvalue;

  if (! PyArg_ParseTuple(args, "O", &fromPython)) 
    return NULL;                             

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

    long long tmpid = 0;
    
    if (  event->event_id )
      tmpid =  event->event_id->id;  
  
    tmpvalue = Py_BuildValue(
        "{s:s, s:i, s:i, s:d, s:d, s:d, s:d, s:d, s:d, s:d,\
        s:d, s:d, s:i, s:d, s:L, s:s, s:s, s:d, s:i, s:i,\
        s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d,\
        s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d}",
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
        "sigmasq", event->sigmasq,
        "event_id", tmpid,
        "search", event->search,
        "channel",event->channel,
        "end_time_gmst", event->end_time_gmst,
        "impulse_time", event->impulse_time.gpsSeconds,
        "impulse_time_ns", event->impulse_time.gpsNanoSeconds,
        
        "template_duration", event->template_duration,
        "event_duration", event->event_duration,
        "amplitude", event->amplitude,
        "tau0", event->tau0,
        "tau2", event->tau2,
        "tau3", event->tau3,
        "tau4", event->tau4,
        "tau5", event->tau5,
        "ttotal", event->ttotal,
        "psi0", event->psi0,
        
        "psi3", event->psi3,
        "alpha", event->alpha,
        "alpha1", event->alpha1,
        "alpha2", event->alpha2,
        "alpha3", event->alpha3,
        "alpha4", event->alpha4,
        "alpha5", event->alpha5,
        "alpha6", event->alpha6,
        "beta", event->beta,
        "f_final", event->f_final,
        "rsqveto_duration", event->rsqveto_duration);

    PyList_SetItem(outlist, j, tmpvalue);

  }

  while(eventHead) {
    event = eventHead;
    eventHead = eventHead->next;
    XLALFreeSnglInspiral ( &event );
  }

  return outlist;
}

/******************************************************************** 
 * Simulated Inspiral Reading Function
 ********************************************************************/
static PyObject *   
read_sim_inspiral(PyObject *self, PyObject *args) 
{ 
  SimInspiralTable *eventHead=NULL;
  SimInspiralTable **addevent=&eventHead;
  SimInspiralTable *event=NULL;
  PyObject *fromPython;
  int j, m=0, n=0, nelement=0, len;
  int startEvent = 0, stopEvent = 0;
  PyObject *outlist;
  PyObject *tmpvalue;

  if (! PyArg_ParseTuple(args, "O", &fromPython)) 
    return NULL;                             

  if (! PyList_Check(fromPython) )
    return NULL;

  len = PyList_Size(fromPython);
  for(n=0;n<len;n++)
  {
    m = SimInspiralTableFromLIGOLw ( addevent,
        PyString_AsString(PyList_GetItem(fromPython, n)), 
        startEvent, stopEvent);

    while(*addevent)
      addevent = &(*addevent)->next;

    nelement += m;
  }

  /* quit proceeding if no trigger to precess */
  if (nelement<=0) {
    return NULL;
  }

  outlist = PyList_New(nelement);
  for ( j=0, event = eventHead; event ; j++, event = event->next )
  {
    tmpvalue = Py_BuildValue(
        "{s:s, s:s, s:i, s:i, s:i, s:i, s:i, s:i, s:i, s:i,\
        s:i, s:i, s:i, s:i, s:d, s:d, s:d, s:d, s:d, s:d,\
        s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d,\
        s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d,\
        s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d}",
        "waveform", event->waveform,
        "source", event->source,
        "geocent_end_time", event->geocent_end_time.gpsSeconds,
        "geocent_end_time_ns", event->geocent_end_time.gpsNanoSeconds,
        "h_end_time", event->h_end_time.gpsSeconds,
        "h_end_time_ns", event->h_end_time.gpsNanoSeconds,
        "l_end_time", event->l_end_time.gpsSeconds,
        "l_end_time_ns", event->l_end_time.gpsNanoSeconds,
        "g_end_time", event->g_end_time.gpsSeconds,
        "g_end_time_ns", event->g_end_time.gpsNanoSeconds,

        "t_end_time", event->t_end_time.gpsSeconds,
        "t_end_time_ns", event->t_end_time.gpsNanoSeconds,
        "v_end_time", event->v_end_time.gpsSeconds,
        "v_end_time_ns", event->v_end_time.gpsNanoSeconds,
        "end_time_gmst", event->end_time_gmst,
        "mass1", event->mass1, 
        "mass2", event->mass2, 
        "eta", event->eta, 
        "distance", event->distance, 
        "longitude", event->longitude, 

        "latitude", event->latitude, 
        "inclination", event->inclination, 
        "coa_phase", event->coa_phase, 
        "polarization", event->polarization,
        "psi0", event->psi0,
        "psi3", event->psi3,
        "alpha", event->alpha,
        "alpha1", event->alpha1,
        "alpha2", event->alpha2,
        "alpha3", event->alpha3,

        "alpha4", event->alpha4,
        "alpha5", event->alpha5,
        "alpha6", event->alpha6,
        "beta", event->beta,
        "spin1x", event->spin1x,
        "spin1y", event->spin1y,
        "spin1z", event->spin1z,
        "spin2x", event->spin2x,
        "spin2y", event->spin2y,
        "spin2z", event->spin2z,

        "theta0", event->theta0,
        "phi0", event->phi0,
        "f_lower", event->f_lower,
        "f_final", event->f_final,
        "mchirp", event->mchirp,
        "eff_dist_h", event->eff_dist_h,
        "eff_dist_l", event->eff_dist_l,
        "eff_dist_g", event->eff_dist_g,
        "eff_dist_t", event->eff_dist_t,
        "eff_dist_v", event->eff_dist_v);
    PyList_SetItem(outlist, j, tmpvalue);
  }

  while(eventHead) {
    event = eventHead;
    eventHead = eventHead->next;
    XLALFreeSimInspiral ( &event );
  }

  return outlist;
}


/******************************************************************** 
 * Multi Inspiral Table Reading Function
 ********************************************************************/
static PyObject *   
read_multi_inspiral(PyObject *self, PyObject *args)
{ 
  MultiInspiralTable *eventHead=NULL;
  MultiInspiralTable **addevent=&eventHead;
  MultiInspiralTable *event=NULL;
  PyObject *fromPython;
  int j, m=0, n=0, nelement=0, len;
  int startEvent = 0, stopEvent = -1;
  PyObject *outlist;
  PyObject *tmpvalue;

  if (! PyArg_ParseTuple(args, "O", &fromPython))
    return NULL;                             

  if (! PyList_Check(fromPython) )
    return NULL;

  len = PyList_Size(fromPython);
  
  for(n=0;n<len;n++)
  {
    *addevent = XLALMultiInspiralTableFromLIGOLw (
        PyString_AsString(PyList_GetItem(fromPython, n)) );

    while(*addevent)
      addevent = &(*addevent)->next;

  }
   nelement = XLALCountMultiInspiralTable ( eventHead );

  outlist = PyList_New(nelement);
  for ( j=0, event = eventHead; event ; j++, event = event->next )
  {
    tmpvalue = Py_BuildValue(
        "{s:s, s:i, s:i, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:i, s:L}",
        "ifos", event->ifos,
        "end_time", event->end_time.gpsSeconds,
        "end_time_ns", event->end_time.gpsNanoSeconds,
        "coa_phase", event->coa_phase,
        "mass1", event->mass1,
        "mass2", event->mass2,
        "mchirp", event->mchirp,
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


/******************************************************************** 
 * XML opening Function
 ********************************************************************/
static PyObject*
open_xml(PyObject *self, PyObject *args)
{
  char* filename;
  LIGOLwXMLStream* xmlStream;
  static LALStatus status;
  PyObject* fileObj;
  
  /* extract arguments */
  if (! PyArg_ParseTuple(args, "s",  &filename))
    return NULL; 

  /* setting xlal error to zero */
  XLALClearErrno();
  
  /* open output file */ 
  xmlStream=(LIGOLwXMLStream *) LALCalloc( 1, sizeof(LIGOLwXMLStream) );
  LALOpenLIGOLwXMLFile ( &status, xmlStream, filename );

  /* convert to python object */
  fileObj=PyCObject_FromVoidPtr((void*)(xmlStream),NULL);
  return fileObj;
}

/******************************************************************** 
 * XML closing Function
 ********************************************************************/
static PyObject* 
close_xml(PyObject *self, PyObject *args) 
{  
  PyObject* fileObj;
  LIGOLwXMLStream* xmlStream;
  static LALStatus status;
  
  /* extract arguments */
  if (! PyArg_ParseTuple(args, "O",  &fileObj))
    return NULL; 

  /* get the file stream */
  xmlStream=(LIGOLwXMLStream*)PyCObject_AsVoidPtr(fileObj);

  /* and close the file */
  LALCloseLIGOLwXMLFile( &status, xmlStream );
  LALFree(xmlStream);

  return PyInt_FromLong(1);
}


/******************************************************************** 
 * Process Writing Function
 ********************************************************************/
static PyObject* 
write_process(PyObject *self, PyObject *args) 
{  
  ProcessTable *eventHead=NULL;
  ProcessTable *event=NULL;
  int n=0, len;
  PyObject* fileObj;
  PyObject* inlist;
  PyObject* tmpvalue;
  LIGOLwXMLStream* xmlStream;
  MetadataTable outputTable;
  static LALStatus status;

  /* extract arguments */
  if (! PyArg_ParseTuple(args, "OO",  &fileObj, &inlist))
    return NULL; 
  
  /* check input file */
  if ( !PyList_Check(inlist) ) 
    return NULL;

  /* get the file stream */
  xmlStream=(LIGOLwXMLStream*)PyCObject_AsVoidPtr(fileObj);

  /* fill SnglInspiralTable */
  len = PyList_Size(inlist);
  for(n=0;n<len;n++)  {
    
    /* allocate new memory for next table */
    if ( ! eventHead ) {
      event = eventHead = (ProcessTable *) LALCalloc( 1, sizeof(ProcessTable) );
    } else {
      event = event->next = (ProcessTable *) LALCalloc( 1, sizeof(ProcessTable) );
    }

    /* get a 'row' from the python list */
    tmpvalue=PyList_GetItem(inlist, n);
    LALSnprintf( event->program, LIGOMETA_PROGRAM_MAX  * sizeof(CHAR),
                 "%s", PyString_AsString(PyDict_GetItem(tmpvalue,  PyString_FromString("program") )) );
    LALSnprintf( event->version, LIGOMETA_VERSION_MAX  * sizeof(CHAR),
                 "%s", PyString_AsString(PyDict_GetItem(tmpvalue,  PyString_FromString("version") )) );
    LALSnprintf( event->cvs_repository, LIGOMETA_CVS_REPOSITORY_MAX  * sizeof(CHAR),
                 "%s", PyString_AsString(PyDict_GetItem(tmpvalue,  PyString_FromString("cvs_repository") )) );
    event->cvs_entry_time.gpsSeconds = (int)PyFloat_AsDouble(PyDict_GetItem(tmpvalue, PyString_FromString("cvs_entry_time") ));
    LALSnprintf( event->comment, LIGOMETA_COMMENT_MAX  * sizeof(CHAR),
                 "%s", PyString_AsString(PyDict_GetItem(tmpvalue,  PyString_FromString("comment") )) );
    event->is_online = (int)PyFloat_AsDouble(PyDict_GetItem(tmpvalue, PyString_FromString("is_online") ));
    LALSnprintf( event->node, LIGOMETA_NODE_MAX  * sizeof(CHAR),
                 "%s", PyString_AsString(PyDict_GetItem(tmpvalue,  PyString_FromString("node") )) );
    LALSnprintf( event->username, LIGOMETA_USERNAME_MAX  * sizeof(CHAR),
                 "%s", PyString_AsString(PyDict_GetItem(tmpvalue,  PyString_FromString("username") )) );
    event->start_time.gpsSeconds = (int)PyFloat_AsDouble( PyDict_GetItem( tmpvalue, PyString_FromString("start_time") ));
    event->end_time.gpsSeconds   = (int)PyFloat_AsDouble( PyDict_GetItem( tmpvalue, PyString_FromString("end_time") ));
    event->jobid = (int)PyFloat_AsDouble( PyDict_GetItem( tmpvalue, PyString_FromString("jobid") ));
    LALSnprintf( event->domain, LIGOMETA_DOMAIN_MAX  * sizeof(CHAR),
                 "%s", PyString_AsString(PyDict_GetItem(tmpvalue,  PyString_FromString("domain") )) );
    LALSnprintf( event->ifos, LIGOMETA_IFOS_MAX  * sizeof(CHAR),
                 "%s", PyString_AsString(PyDict_GetItem(tmpvalue,  PyString_FromString("ifos") )) );

  }

  /* write data to XML file */
  outputTable.processTable = eventHead;
  LALBeginLIGOLwXMLTable( &status, xmlStream, process_table );
  LALWriteLIGOLwXMLTable( &status, xmlStream, outputTable, process_table );
  LALEndLIGOLwXMLTable(   &status, xmlStream ); 
  
  /* clearing memory */
  while(eventHead) {
    event = eventHead;
    eventHead = eventHead->next;
    LALFree(event);
  }

  return PyInt_FromLong(1);
 }


/******************************************************************** 
 * Process Params Writing Function
 ********************************************************************/
static PyObject* 
write_process_params(PyObject *self, PyObject *args) 
{  
  ProcessParamsTable *eventHead=NULL;
  ProcessParamsTable *event=NULL;
  int n=0, len;
  PyObject* fileObj;
  PyObject* inlist;
  PyObject* tmpvalue;
  LIGOLwXMLStream* xmlStream;
  MetadataTable outputTable;
  static LALStatus status;

  /* extract arguments */
  if (! PyArg_ParseTuple(args, "OO",  &fileObj, &inlist))
    return NULL; 
  
  /* check input file */
  if ( !PyList_Check(inlist) ) 
    return NULL;

  /* get the file stream */
  xmlStream=(LIGOLwXMLStream*)PyCObject_AsVoidPtr(fileObj);

  /* fill SnglInspiralTable */
  len = PyList_Size(inlist);
  for(n=0;n<len;n++)  {
    
    /* allocate new memory for next table */
    if ( ! eventHead ) {
      event = eventHead = (ProcessParamsTable *) LALCalloc( 1, sizeof(ProcessParamsTable) );
    } else {
      event = event->next = (ProcessParamsTable *) LALCalloc( 1, sizeof(ProcessParamsTable) );
    }

    /* get a 'row' from the python list */
    tmpvalue=PyList_GetItem(inlist, n);
    LALSnprintf( event->program, LIGOMETA_PROGRAM_MAX  * sizeof(CHAR),
                 "%s", PyString_AsString(PyDict_GetItem(tmpvalue,  PyString_FromString("program") )) );
    LALSnprintf( event->param, LIGOMETA_PARAM_MAX  * sizeof(CHAR),
                 "%s", PyString_AsString(PyDict_GetItem(tmpvalue,  PyString_FromString("param") )) );
    LALSnprintf( event->type, LIGOMETA_TYPE_MAX  * sizeof(CHAR),
                 "%s", PyString_AsString(PyDict_GetItem(tmpvalue,  PyString_FromString("type") )) );
    LALSnprintf( event->value, LIGOMETA_VALUE_MAX  * sizeof(CHAR),
                 "%s", PyString_AsString(PyDict_GetItem(tmpvalue,  PyString_FromString("value") )) );
  }

  /* write data to XML file */
  outputTable.processParamsTable = eventHead;
  LALBeginLIGOLwXMLTable( &status, xmlStream, process_params_table );
  LALWriteLIGOLwXMLTable( &status, xmlStream, outputTable, process_params_table );
  LALEndLIGOLwXMLTable( &status, xmlStream );
  
  
  /* clearing memory */
  while(eventHead) {
    event = eventHead;
    eventHead = eventHead->next;
    LALFree(event);
  }

  return PyInt_FromLong(1);
 }


/******************************************************************** 
 * Search Summary Writing Function
 ********************************************************************/
static PyObject* 
write_search_summary(PyObject *self, PyObject *args) 
{  
  SearchSummaryTable *eventHead=NULL;
  SearchSummaryTable *event=NULL;
  int n=0, len;
  PyObject* fileObj;
  PyObject* inlist;
  PyObject* tmpvalue;
  LIGOLwXMLStream* xmlStream;
  MetadataTable outputTable;
  static LALStatus status;

  /* extract arguments */
  if (! PyArg_ParseTuple(args, "OO",  &fileObj, &inlist))
    return NULL; 
  
  /* check input file */
  if ( !PyList_Check(inlist) ) 
    return NULL;

  /* get the file stream */
  xmlStream=(LIGOLwXMLStream*)PyCObject_AsVoidPtr(fileObj);

  /* fill SearchSummaryTable */
  len = PyList_Size(inlist);
  for(n=0;n<len;n++)  {
    
    /* allocate new memory for next table */
    if ( ! eventHead ) {
      event = eventHead = (SearchSummaryTable *) LALCalloc( 1, sizeof(SearchSummaryTable) );
    } else {
      event = event->next = (SearchSummaryTable *) LALCalloc( 1, sizeof(SearchSummaryTable) );
    }

    /* get a 'row' from the python list */
    tmpvalue=PyList_GetItem(inlist, n);
    LALSnprintf( event->comment, LIGOMETA_COMMENT_MAX  * sizeof(CHAR),
                 "%s", PyString_AsString(PyDict_GetItem(tmpvalue,  PyString_FromString("comment") )) );
    LALSnprintf( event->ifos, LIGOMETA_IFOS_MAX  * sizeof(CHAR),
                 "%s", PyString_AsString(PyDict_GetItem(tmpvalue,  PyString_FromString("ifos") )) );
    event->in_start_time.gpsSeconds    = PyInt_AsLong(PyDict_GetItem(tmpvalue, PyString_FromString("in_start_time") ));
    event->in_start_time.gpsNanoSeconds= 0;
    event->in_end_time.gpsSeconds    = PyInt_AsLong(PyDict_GetItem(tmpvalue, PyString_FromString("in_end_time") ));
    event->in_end_time.gpsNanoSeconds= 0;
    event->out_start_time.gpsSeconds    = PyInt_AsLong(PyDict_GetItem(tmpvalue, PyString_FromString("out_start_time") ));
    event->out_start_time.gpsNanoSeconds= 0;
    event->out_end_time.gpsSeconds    = PyInt_AsLong(PyDict_GetItem(tmpvalue, PyString_FromString("out_end_time") ));
    event->out_end_time.gpsNanoSeconds= 0;
    event->nevents= PyInt_AsLong(PyDict_GetItem(tmpvalue, PyString_FromString("nevents") ));
    event->nnodes= PyInt_AsLong(PyDict_GetItem(tmpvalue, PyString_FromString("nnodes") ));
  }

  /* write data to XML file */
  outputTable.searchSummaryTable = eventHead;
  LALBeginLIGOLwXMLTable( &status, xmlStream, search_summary_table );
  LALWriteLIGOLwXMLTable( &status, xmlStream, outputTable, search_summary_table );
  LALEndLIGOLwXMLTable( &status, xmlStream );
  
  
  /* clearing memory */
  while(eventHead) {
    event = eventHead;
    eventHead = eventHead->next;
    LALFree(event);
  }

  return PyInt_FromLong(1);
 }


/******************************************************************** 
 * Sim Inspiral Writing Function
 ********************************************************************/
static PyObject* 
write_sim_inspiral(PyObject *self, PyObject *args) 
{  
  SimInspiralTable *eventHead=NULL;
  SimInspiralTable *event=NULL;
  int n=0, len;
  PyObject* fileObj;
  PyObject* inlist;
  PyObject* tmpvalue;
  LIGOLwXMLStream* xmlStream;
  MetadataTable outputTable;
  static LALStatus status;

  /* extract arguments */
  if (! PyArg_ParseTuple(args, "OO",  &fileObj, &inlist))
    return NULL; 
  
  /* check input file */
  if ( !PyList_Check(inlist) ) 
    return NULL;

  /* get the file stream */
  xmlStream=(LIGOLwXMLStream*)PyCObject_AsVoidPtr(fileObj);

  /* fill SnglInspiralTable */
  len = PyList_Size(inlist);
  if (!len) return PyInt_FromLong(1);
  for(n=0;n<len;n++)  {
    
    /* allocate new memory for next table */
    if ( ! eventHead ) {
      event = eventHead = (SimInspiralTable *) LALCalloc( 1, sizeof(SimInspiralTable) );
    } else {
      event = event->next = (SimInspiralTable *) LALCalloc( 1, sizeof(SimInspiralTable) );
    }

    /* get a 'row' from the python list */
    tmpvalue=PyList_GetItem(inlist, n);

    LALSnprintf( event->waveform, LIGOMETA_WAVEFORM_MAX * sizeof(CHAR),
                 "%s", PyString_AsString(PyDict_GetItem(tmpvalue,  PyString_FromString("waveform"))) );
    LALSnprintf( event->source, LIGOMETA_SOURCE_MAX * sizeof(CHAR),
                 "%s", PyString_AsString(PyDict_GetItem(tmpvalue,  PyString_FromString("source"))) );
    event->geocent_end_time.gpsSeconds     = PyInt_AsLong(PyDict_GetItem(tmpvalue,  PyString_FromString("geocent_end_time")));
    event->geocent_end_time.gpsNanoSeconds = PyInt_AsLong(PyDict_GetItem(tmpvalue,  PyString_FromString("geocent_end_time_ns")));
    event->h_end_time.gpsSeconds= PyInt_AsLong(PyDict_GetItem(tmpvalue,  PyString_FromString("h_end_time")));
    event->h_end_time.gpsNanoSeconds= PyInt_AsLong(PyDict_GetItem(tmpvalue,  PyString_FromString("h_end_time_ns")));
    event->l_end_time.gpsSeconds    = PyInt_AsLong(PyDict_GetItem(tmpvalue,  PyString_FromString("l_end_time")));
    event->l_end_time.gpsNanoSeconds= PyInt_AsLong(PyDict_GetItem(tmpvalue, PyString_FromString("l_end_time_ns") ));
    event->g_end_time.gpsSeconds    = PyInt_AsLong(PyDict_GetItem(tmpvalue, PyString_FromString("g_end_time") ));
    event->g_end_time.gpsNanoSeconds= PyInt_AsLong(PyDict_GetItem(tmpvalue,  PyString_FromString("g_end_time_ns")));    

    event->t_end_time.gpsSeconds    = PyInt_AsLong(PyDict_GetItem(tmpvalue, PyString_FromString("t_end_time") ));
    event->t_end_time.gpsNanoSeconds= PyInt_AsLong(PyDict_GetItem(tmpvalue, PyString_FromString("t_end_time_ns") ));
    event->v_end_time.gpsSeconds    = PyInt_AsLong(PyDict_GetItem(tmpvalue, PyString_FromString("v_end_time") ));
    event->v_end_time.gpsNanoSeconds= PyInt_AsLong(PyDict_GetItem(tmpvalue,  PyString_FromString("v_end_time_ns")));
    event->end_time_gmst            = PyFloat_AsDouble(PyDict_GetItem(tmpvalue, PyString_FromString("end_time_gmst") ));
    event->mass1= PyFloat_AsDouble(PyDict_GetItem(tmpvalue,  PyString_FromString("mass1")));
    event->mass2=  PyFloat_AsDouble(PyDict_GetItem(tmpvalue,  PyString_FromString("mass2")));
    event->eta = PyFloat_AsDouble(PyDict_GetItem(tmpvalue, PyString_FromString("eta") ));
    event->distance = PyFloat_AsDouble(PyDict_GetItem(tmpvalue, PyString_FromString("distance") ));
    event->longitude = PyFloat_AsDouble(PyDict_GetItem(tmpvalue, PyString_FromString("longitude") ));

    event->latitude = PyFloat_AsDouble(PyDict_GetItem(tmpvalue,  PyString_FromString("latitude") ));
    event->inclination = PyFloat_AsDouble(PyDict_GetItem(tmpvalue, PyString_FromString("inclination") ));
    event->coa_phase = PyFloat_AsDouble(PyDict_GetItem(tmpvalue, PyString_FromString("coa_phase") ));
    event->polarization = PyFloat_AsDouble(PyDict_GetItem(tmpvalue, PyString_FromString("polarization") ));
    event->psi0 = PyFloat_AsDouble(PyDict_GetItem(tmpvalue,  PyString_FromString("psi0") ));
    event->psi3 = PyFloat_AsDouble(PyDict_GetItem(tmpvalue,  PyString_FromString("psi3") ));
    event->alpha = PyFloat_AsDouble(PyDict_GetItem(tmpvalue,  PyString_FromString("alpha") ));
    event->alpha1 = PyFloat_AsDouble(PyDict_GetItem(tmpvalue,  PyString_FromString("alpha1") ));
    event->alpha2 = PyFloat_AsDouble(PyDict_GetItem(tmpvalue,  PyString_FromString("alpha2") ));
    event->alpha3 = PyFloat_AsDouble(PyDict_GetItem(tmpvalue,  PyString_FromString("alpha3") ));

    event->alpha4 = PyFloat_AsDouble(PyDict_GetItem(tmpvalue,  PyString_FromString("alpha4") ));
    event->alpha5 = PyFloat_AsDouble(PyDict_GetItem(tmpvalue,  PyString_FromString("alpha5") ));
    event->alpha6 = PyFloat_AsDouble(PyDict_GetItem(tmpvalue,  PyString_FromString("alpha6") ));
    event->beta = PyFloat_AsDouble(PyDict_GetItem(tmpvalue,  PyString_FromString("beta") ));
    event->spin1x= PyFloat_AsDouble(PyDict_GetItem(tmpvalue,  PyString_FromString("spin1x") ));
    event->spin1y= PyFloat_AsDouble(PyDict_GetItem(tmpvalue,  PyString_FromString("spin1y") ));
    event->spin1z= PyFloat_AsDouble(PyDict_GetItem(tmpvalue,  PyString_FromString("spin1z") ));
    event->spin2x= PyFloat_AsDouble(PyDict_GetItem(tmpvalue,  PyString_FromString("spin2x") ));
    event->spin2y= PyFloat_AsDouble(PyDict_GetItem(tmpvalue,  PyString_FromString("spin2y") ));
    event->spin2z= PyFloat_AsDouble(PyDict_GetItem(tmpvalue,  PyString_FromString("spin2z") ));

    event->theta0 = PyFloat_AsDouble(PyDict_GetItem(tmpvalue,  PyString_FromString("theta0") ));
    event->phi0= PyFloat_AsDouble(PyDict_GetItem(tmpvalue, PyString_FromString("phi0") ));
    event->f_lower= PyFloat_AsDouble(PyDict_GetItem(tmpvalue, PyString_FromString("f_lower") ));
    event->f_final= PyFloat_AsDouble(PyDict_GetItem(tmpvalue, PyString_FromString("f_final") ));
    event->mchirp= PyFloat_AsDouble(PyDict_GetItem(tmpvalue, PyString_FromString("mchirp") ));
    event->eff_dist_h= PyFloat_AsDouble(PyDict_GetItem(tmpvalue, PyString_FromString("eff_dist_h") ));
    event->eff_dist_l= PyFloat_AsDouble(PyDict_GetItem(tmpvalue, PyString_FromString("eff_dist_l") ));  
    event->eff_dist_g= PyFloat_AsDouble(PyDict_GetItem(tmpvalue, PyString_FromString("eff_dist_g") ));     
    event->eff_dist_t= PyFloat_AsDouble(PyDict_GetItem(tmpvalue, PyString_FromString("eff_dist_t") ));    
    event->eff_dist_v= PyFloat_AsDouble(PyDict_GetItem(tmpvalue, PyString_FromString("eff_dist_v") ));
  }

  /* write data to XML file */
  outputTable.simInspiralTable = eventHead;
  LALBeginLIGOLwXMLTable( &status, xmlStream, sim_inspiral_table );
  LALWriteLIGOLwXMLTable( &status, xmlStream, outputTable, sim_inspiral_table );
  LALEndLIGOLwXMLTable(   &status, xmlStream );

  /* deleting entries in the SimInspiralTable */
  while(eventHead) {
    event = eventHead;
    eventHead = eventHead->next;
    XLALFreeSimInspiral ( &event );
  }

  LALCheckMemoryLeaks();
 
  return PyInt_FromLong(1);
}

/******************************************************************** 
 * Sngl Inspiral Writing Function (BEGIN)
 ********************************************************************/
static PyObject* 
write_sngl_inspiral_begin(PyObject *self, PyObject *args) 
{  
  SnglInspiralTable *eventHead=NULL;
  SnglInspiralTable *event=NULL;
  int n=0, len;
  PyObject* fileObj;
  PyObject* inlist;
  PyObject* tmpvalue;
  LIGOLwXMLStream* xmlStream;
  MetadataTable outputTable;
  static LALStatus status;

  /* extract arguments */
  if (! PyArg_ParseTuple(args, "O",  &fileObj))
    return NULL; 
  
  /* open the table */
  xmlStream=(LIGOLwXMLStream*)PyCObject_AsVoidPtr(fileObj);
  LALBeginLIGOLwXMLTable( &status, xmlStream, sngl_inspiral_table );
  
  return PyInt_FromLong(1);
}

/******************************************************************** 
 * Sngl Inspiral Writing Function(WRITE)
 ********************************************************************/
static PyObject* 
write_sngl_inspiral_write(PyObject *self, PyObject *args) 
{  
  SnglInspiralTable *eventHead=NULL;
  SnglInspiralTable *event=NULL;
  int n=0, len;
  PyObject* fileObj;
  PyObject* inlist;
  PyObject* tmpvalue;
  LIGOLwXMLStream* xmlStream;
  MetadataTable outputTable;
  static LALStatus status;

  /* extract arguments */
  if (! PyArg_ParseTuple(args, "OO",  &fileObj, &inlist))
    return NULL; 
  
  /* check input file */
  if ( !PyList_Check(inlist) ) 
    return NULL;

  /* get the file stream */
  xmlStream=(LIGOLwXMLStream*)PyCObject_AsVoidPtr(fileObj);

  /* fill SnglInspiralTable */
  len = PyList_Size(inlist);
  for(n=0;n<len;n++)  {
    
    /* allocate new memory for next table */
    if ( ! eventHead ) {
      event = eventHead = (SnglInspiralTable *) LALCalloc( 1, sizeof(SnglInspiralTable) );
    } else {
      event = event->next = (SnglInspiralTable *) LALCalloc( 1, sizeof(SnglInspiralTable) );
    }

    /* get a 'row' from the python list */
    tmpvalue=PyList_GetItem(inlist, n);

    /* fill eht values into the event structure */  
    LALSnprintf( event->ifo, LIGOMETA_IFO_MAX  * sizeof(CHAR),
                 "%s", PyString_AsString(PyDict_GetItem(tmpvalue,  PyString_FromString("ifo") )) );
    LALSnprintf( event->search, LIGOMETA_SEARCH_MAX  * sizeof(CHAR),
                 "%s", PyString_AsString(PyDict_GetItem(tmpvalue,  PyString_FromString("search") )) );
    LALSnprintf( event->channel, LIGOMETA_CHANNEL_MAX  * sizeof(CHAR),
                 "%s", PyString_AsString(PyDict_GetItem(tmpvalue,  PyString_FromString("channel") )) );
    event->end_time.gpsSeconds    = PyInt_AsLong(PyDict_GetItem(tmpvalue, PyString_FromString("end_time") ));
    event->end_time.gpsNanoSeconds = PyInt_AsLong(PyDict_GetItem(tmpvalue, PyString_FromString("end_time_ns") ));
    event->end_time_gmst = PyFloat_AsDouble(PyDict_GetItem(tmpvalue, PyString_FromString("end_time_gmst") ));
    event->impulse_time.gpsSeconds    = PyInt_AsLong(PyDict_GetItem(tmpvalue, PyString_FromString("impulse_time") ));
    event->impulse_time.gpsNanoSeconds = PyInt_AsLong(PyDict_GetItem(tmpvalue, PyString_FromString("impulse_time_ns") ));
    event->template_duration = PyFloat_AsDouble(PyDict_GetItem(tmpvalue, PyString_FromString("template_duration") ));
    event->event_duration = PyFloat_AsDouble(PyDict_GetItem(tmpvalue, PyString_FromString("event_duration") ));

    event->amplitude = PyFloat_AsDouble(PyDict_GetItem(tmpvalue, PyString_FromString("amplitude") ));
    event->eff_distance= PyFloat_AsDouble(PyDict_GetItem(tmpvalue, PyString_FromString("eff_distance") ));
    event->coa_phase = PyFloat_AsDouble( PyDict_GetItem( tmpvalue, PyString_FromString("coa_phase") ));
    event->mass1= PyFloat_AsDouble(PyDict_GetItem(tmpvalue, PyString_FromString("mass1")));
    event->mass2=  PyFloat_AsDouble(PyDict_GetItem(tmpvalue, PyString_FromString("mass2")));
    event->mchirp=  PyFloat_AsDouble(PyDict_GetItem(tmpvalue, PyString_FromString("mchirp")));
    event->mtotal=  PyFloat_AsDouble(PyDict_GetItem(tmpvalue, PyString_FromString("mtotal")));
    event->eta = PyFloat_AsDouble(PyDict_GetItem(tmpvalue, PyString_FromString("eta")));
    event->tau0 = PyFloat_AsDouble(PyDict_GetItem(tmpvalue, PyString_FromString("tau0") ));
    event->tau2 = PyFloat_AsDouble(PyDict_GetItem(tmpvalue, PyString_FromString("tau2") ));

    event->tau3   = PyFloat_AsDouble(PyDict_GetItem(tmpvalue, PyString_FromString("tau3") ));
    event->tau4   = PyFloat_AsDouble(PyDict_GetItem(tmpvalue, PyString_FromString("tau4") ));
    event->tau5   = PyFloat_AsDouble(PyDict_GetItem(tmpvalue, PyString_FromString("tau5") ));
    event->ttotal = PyFloat_AsDouble(PyDict_GetItem(tmpvalue, PyString_FromString("ttotal") ));
    event->psi0   = PyFloat_AsDouble(PyDict_GetItem(tmpvalue, PyString_FromString("psi0") ));
    event->psi3   = PyFloat_AsDouble(PyDict_GetItem(tmpvalue, PyString_FromString("psi3") ));
    event->alpha  = PyFloat_AsDouble(PyDict_GetItem(tmpvalue, PyString_FromString("alpha") ));
    event->alpha1 = PyFloat_AsDouble(PyDict_GetItem(tmpvalue, PyString_FromString("alpha1") ));
    event->alpha2 = PyFloat_AsDouble(PyDict_GetItem(tmpvalue, PyString_FromString("alpha2") ));
    event->alpha3 = PyFloat_AsDouble(PyDict_GetItem(tmpvalue, PyString_FromString("alpha3") ));

    event->alpha4 = PyFloat_AsDouble(PyDict_GetItem(tmpvalue, PyString_FromString("alpha4") ));
    event->alpha5 = PyFloat_AsDouble(PyDict_GetItem(tmpvalue, PyString_FromString("alpha5") ));
    event->alpha6 = PyFloat_AsDouble(PyDict_GetItem(tmpvalue, PyString_FromString("alpha6") ));
    event->beta   = PyFloat_AsDouble(PyDict_GetItem(tmpvalue, PyString_FromString("beta")));
    event->f_final= PyFloat_AsDouble(PyDict_GetItem(tmpvalue, PyString_FromString("f_final")));
    event->snr    = PyFloat_AsDouble(PyDict_GetItem(tmpvalue, PyString_FromString("snr")));
    event->chisq  = PyFloat_AsDouble(PyDict_GetItem(tmpvalue, PyString_FromString("chisq")));
    event->chisq_dof = PyInt_AsLong(PyDict_GetItem(tmpvalue, PyString_FromString("chisq_dof")));
    event->sigmasq   = PyFloat_AsDouble(PyDict_GetItem(tmpvalue, PyString_FromString("sigmasq")));

    event->rsqveto_duration = PyFloat_AsDouble(PyDict_GetItem(tmpvalue, PyString_FromString("rsqveto_duration")));
    event->event_id = (EventIDColumn *) LALCalloc( 1, sizeof(EventIDColumn) );
    event->event_id->id = PyLong_AsLongLong(PyDict_GetItem(tmpvalue, PyString_FromString("event_id")));
    event->event_id->snglInspiralTable = event;
  }

  /* write data to XML file */
  outputTable.snglInspiralTable = eventHead;
  LALWriteLIGOLwXMLTable( &status, xmlStream, outputTable, sngl_inspiral_table );
  
  /* deleting entries in the SnglInspiralTable */
  while(eventHead) {
    event = eventHead;
    eventHead = eventHead->next;
    XLALFreeSnglInspiral ( &event );
  }
  LALCheckMemoryLeaks();
 
  return PyInt_FromLong(1);
}

/******************************************************************** 
 * Sngl Inspiral Writing Function (END)
 ********************************************************************/
static PyObject* 
write_sngl_inspiral_end(PyObject *self, PyObject *args) 
{  
  SnglInspiralTable *eventHead=NULL;
  SnglInspiralTable *event=NULL;
  int n=0, len;
  PyObject* fileObj;
  PyObject* inlist;
  PyObject* tmpvalue;
  LIGOLwXMLStream* xmlStream;
  MetadataTable outputTable;
  static LALStatus status;

  /* extract arguments */
  if (! PyArg_ParseTuple(args, "O",  &fileObj))
    return NULL; 
  
  /* close the table */
  xmlStream=(LIGOLwXMLStream*)PyCObject_AsVoidPtr(fileObj);
  LALEndLIGOLwXMLTable( &status, xmlStream);
  
  return PyInt_FromLong(1);
}


/******************************************************************** 
 * supporting snippets
 ********************************************************************/
/* registration table  */
static struct PyMethodDef support_methods[] = {
    {"read_process_params", read_process_params, 1},  
    {"read_process", read_process, 1},  
    {"read_search_summary", read_search_summary, 1},  
    {"read_summ_value", read_summ_value, 1},  
    {"read_sngl_burst", read_sngl_burst, 1},  
    {"read_sngl_inspiral", read_sngl_inspiral, 1}, 
    {"read_sim_inspiral", read_sim_inspiral, 1}, 
    {"read_multi_inspiral", read_multi_inspiral, 1},
    {"open_xml", open_xml, 1}, 
    {"close_xml", close_xml, 1}, 
    {"write_process", write_process, 1}, 
    {"write_process_params", write_process_params, 1},
    {"write_search_summary", write_search_summary, 1}, 
    {"write_sim_inspiral", write_sim_inspiral, 1}, 
    {"write_sngl_inspiral_begin", write_sngl_inspiral_begin, 1},
    {"write_sngl_inspiral_write", write_sngl_inspiral_write, 1},
    {"write_sngl_inspiral_end", write_sngl_inspiral_end, 1},
    {NULL, NULL} 
};

/* module initializer */
void initsupport()               /* called on first import */
{                               /* name matters if loaded dynamically */
    /* mod name, table ptr */
    (void) Py_InitModule("pylal.support", support_methods);
}
