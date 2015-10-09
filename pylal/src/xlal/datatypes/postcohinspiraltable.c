/*
 * Copyright (C) 2010  Kipp Cannon
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */


/*
 * ============================================================================
 *
 *                                  Preamble
 *
 * ============================================================================
 */


#include <Python.h>
#include <structmember.h>
#include <lal/LIGOMetadataTables.h>
#include <ligotimegps.h>
#include <misc.h>
#include <postcohinspiraltable.h>


#define MODULE_NAME PYLAL_POSTCOHINSPIRALTABLE_MODULE_NAME


/*
 * ============================================================================
 *
 *                                    Type
 *
 * ============================================================================
 */


/*
 * Cached ID types
 */


//static PyObject *postcoh_inspiral_event_id_type = NULL;
static PyObject *process_id_type = NULL;


/*
 * Member access
 */


static struct PyMemberDef members[] = {
	{"end_time", T_INT, offsetof(pylal_PostcohInspiralTable, postcoh_inspiral.end_time.gpsSeconds), 0, "end_time"},
	{"end_time_ns", T_INT, offsetof(pylal_PostcohInspiralTable, postcoh_inspiral.end_time.gpsNanoSeconds), 0, "end_time_ns"},
	{"chisq", T_FLOAT, offsetof(pylal_PostcohInspiralTable, postcoh_inspiral.chisq), 0, "chisq"},
	{"tmplt_idx", T_INT, offsetof(pylal_PostcohInspiralTable, postcoh_inspiral.tmplt_idx), 0, "tmplt_idx"},
	{"pix_idx", T_INT, offsetof(pylal_PostcohInspiralTable, postcoh_inspiral.pix_idx), 0, "pix_idx"},
	{"maxsnglsnr", T_FLOAT, offsetof(pylal_PostcohInspiralTable, postcoh_inspiral.maxsnglsnr), 0, "maxsnglsnr"},
	{"cohsnr", T_FLOAT, offsetof(pylal_PostcohInspiralTable, postcoh_inspiral.cohsnr), 0, "cohsnr"},
	{"nullsnr", T_FLOAT, offsetof(pylal_PostcohInspiralTable, postcoh_inspiral.nullsnr), 0, "nullsnr"},
	{NULL,}
};


static PyObject *end_get(PyObject *obj, void *data)
{
	return pylal_LIGOTimeGPS_new(((pylal_PostcohInspiralTable*)obj)->postcoh_inspiral.end_time);
}


static int end_set(PyObject *obj, PyObject *val, void *data)
{
	int seconds = 0;
	int nanoseconds = 0;

	if(val != Py_None) {
		PyObject *attr = PyObject_GetAttrString(val, "gpsSeconds");
		if(!attr)
			return -1;
		seconds = PyInt_AsLong(attr);
		Py_DECREF(attr);
		if(PyErr_Occurred())
			return -1;
		attr = PyObject_GetAttrString(val, "gpsNanoSeconds");
		if(!attr)
			return -1;
		nanoseconds = PyInt_AsLong(attr);
		Py_DECREF(attr);
		if(PyErr_Occurred())
			return -1;
	}

	((pylal_PostcohInspiralTable*)obj)->postcoh_inspiral.end_time.gpsSeconds = seconds;
	((pylal_PostcohInspiralTable*)obj)->postcoh_inspiral.end_time.gpsNanoSeconds = nanoseconds;

	return 0;
}


static struct PyGetSetDef getset[] = {
	{"ifos", pylal_inline_string_get, pylal_inline_string_set, "ifos", &(struct pylal_inline_string_description) {offsetof(pylal_PostcohInspiralTable, postcoh_inspiral.ifos), MAX_ALLIFO_LEN}},
	{"pivotal_ifo", pylal_inline_string_get, pylal_inline_string_set, "pivotal_ifo", &(struct pylal_inline_string_description) {offsetof(pylal_PostcohInspiralTable, postcoh_inspiral.pivotal_ifo), MAX_IFO_LEN}},
	{"is_background", pylal_inline_string_get, pylal_inline_string_set, "is_background", &(struct pylal_inline_string_description) {offsetof(pylal_PostcohInspiralTable, postcoh_inspiral.is_background), 1}},
	{"skymap_fname", pylal_inline_string_get, pylal_inline_string_set, "skymap_fname", &(struct pylal_inline_string_description) {offsetof(pylal_PostcohInspiralTable, postcoh_inspiral.skymap_fname), MAX_SKYMAP_FNAME_LEN}},
	{"end", end_get, end_set, "end", NULL},
	{NULL,}
};


static Py_ssize_t getreadbuffer(PyObject *self, Py_ssize_t segment, void **ptrptr)
{
	if(segment) {
		PyErr_SetString(PyExc_SystemError, "bad segment");
		return -1;
	}
	*ptrptr = &((pylal_PostcohInspiralTable*)self)->postcoh_inspiral;
	return sizeof(((pylal_PostcohInspiralTable*)self)->postcoh_inspiral);
}


static Py_ssize_t getsegcount(PyObject *self, Py_ssize_t *lenp)
{
	if(lenp)
		*lenp = sizeof(((pylal_PostcohInspiralTable*)self)->postcoh_inspiral);
	return 1;
}


static PyBufferProcs as_buffer = {
	.bf_getreadbuffer = getreadbuffer,
	.bf_getsegcount = getsegcount,
	.bf_getwritebuffer = NULL,
	.bf_getcharbuffer = NULL
};


/*
 * Methods
 */


static PyObject *__new__(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
	pylal_PostcohInspiralTable *new = (pylal_PostcohInspiralTable *) PyType_GenericNew(type, args, kwds);

	if(!new)
		return NULL;

	/* link the event_id pointer in the postcoh_inspiral table structure
	 * to the event_id structure */
	//new->postcoh_inspiral.event_id = &new->event_id;

	//new->process_id_i = 0;
	//new->event_id.id = 0;

	/* done */
	return (PyObject *) new;
}


static PyObject *from_buffer(PyObject *cls, PyObject *args)
{
	const PostcohInspiralTable *data;
	Py_ssize_t length;
	unsigned i;
	PyObject *result;

	if(!PyArg_ParseTuple(args, "s#", (const char **) &data, &length))
		return NULL;

	if(length % sizeof(PostcohInspiralTable)) {
		PyErr_SetString(PyExc_ValueError, "buffer size is not an integer multiple of PostcohInspiralTable struct size");
		return NULL;
	}
	length /= sizeof(PostcohInspiralTable);

	result = PyTuple_New(length);
	if(!result)
		return NULL;
	for(i = 0; i < length; i++) {
		PyObject *item = PyType_GenericNew((PyTypeObject *) cls, NULL, NULL);
		if(!item) {
			Py_DECREF(result);
			return NULL;
		}
		/* memcpy postcoh_inspiral row */
		((pylal_PostcohInspiralTable*)item)->postcoh_inspiral = *data++;
		/* repoint event_id to event_id structure */
		//((pylal_PostcohInspiralTable*)item)->postcoh_inspiral.event_id = &((pylal_PostcohInspiralTable*)item)->event_id;

		PyTuple_SET_ITEM(result, i, item);
	}

	return result;
}


static struct PyMethodDef methods[] = {
	{"from_buffer", from_buffer, METH_VARARGS | METH_CLASS, "Construct a tuple of PostcohInspiralTable objects from a buffer object.  The buffer is interpreted as a C array of PostcohInspiralTable structures."},
	{NULL,}
};


/*
 * Type
 */


static PyTypeObject pylal_postcohinspiraltable_type = {
	PyObject_HEAD_INIT(NULL)
	.tp_basicsize = sizeof(pylal_PostcohInspiralTable),
	.tp_doc = "LAL's PostcohInspiralTable structure",
	.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_CHECKTYPES,
	.tp_members = members,
	.tp_methods = methods,
	.tp_getset = getset,
	.tp_as_buffer = &as_buffer,
	.tp_name = MODULE_NAME ".PostcohInspiralTable",
	.tp_new = __new__,
};


/*
 * ============================================================================
 *
 *                            Module Registration
 *
 * ============================================================================
 */


PyMODINIT_FUNC initpostcohinspiraltable(void)
{
	PyObject *module = Py_InitModule3(MODULE_NAME, NULL, "Wrapper for LAL's PostcohInspiralTable type.");

	pylal_ligotimegps_import();

	/* Cached ID types */
	process_id_type = pylal_get_ilwdchar_class("process", "process_id");
	//postcoh_inspiral_event_id_type = pylal_get_ilwdchar_class("postcoh_inspiral", "event_id");

	/* PostcohInspiralTable */
	_pylal_PostcohInspiralTable_Type = &pylal_postcohinspiraltable_type;
	if(PyType_Ready(&pylal_PostcohInspiralTable_Type) < 0)
		return;
	Py_INCREF(&pylal_PostcohInspiralTable_Type);
	PyModule_AddObject(module, "PostcohInspiralTable", (PyObject *) &pylal_PostcohInspiralTable_Type);
}
