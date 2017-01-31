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
	{"end_time_L", T_INT, offsetof(pylal_PostcohInspiralTable, postcoh_inspiral.end_time_L.gpsSeconds), 0, "end_time_L"},
	{"end_time_ns_L", T_INT, offsetof(pylal_PostcohInspiralTable, postcoh_inspiral.end_time_L.gpsNanoSeconds), 0, "end_time_ns_L"},
	{"end_time_H", T_INT, offsetof(pylal_PostcohInspiralTable, postcoh_inspiral.end_time_H.gpsSeconds), 0, "end_time_H"},
	{"end_time_ns_H", T_INT, offsetof(pylal_PostcohInspiralTable, postcoh_inspiral.end_time_H.gpsNanoSeconds), 0, "end_time_ns_H"},
	{"end_time_V", T_INT, offsetof(pylal_PostcohInspiralTable, postcoh_inspiral.end_time_V.gpsSeconds), 0, "end_time_V"},
	{"end_time_ns_V", T_INT, offsetof(pylal_PostcohInspiralTable, postcoh_inspiral.end_time_V.gpsNanoSeconds), 0, "end_time_ns_V"},
	{"snglsnr_L", T_FLOAT, offsetof(pylal_PostcohInspiralTable, postcoh_inspiral.snglsnr_L), 0, "snglsnr_L"},
	{"snglsnr_H", T_FLOAT, offsetof(pylal_PostcohInspiralTable, postcoh_inspiral.snglsnr_H), 0, "snglsnr_H"},
	{"snglsnr_V", T_FLOAT, offsetof(pylal_PostcohInspiralTable, postcoh_inspiral.snglsnr_V), 0, "snglsnr_V"},
	{"coaphase_L", T_FLOAT, offsetof(pylal_PostcohInspiralTable, postcoh_inspiral.coaphase_L), 0, "coaphase_L"},
	{"coaphase_H", T_FLOAT, offsetof(pylal_PostcohInspiralTable, postcoh_inspiral.coaphase_H), 0, "coaphase_H"},
	{"coaphase_V", T_FLOAT, offsetof(pylal_PostcohInspiralTable, postcoh_inspiral.coaphase_V), 0, "coaphase_V"},
	{"chisq_L", T_FLOAT, offsetof(pylal_PostcohInspiralTable, postcoh_inspiral.chisq_L), 0, "chisq_L"},
	{"chisq_H", T_FLOAT, offsetof(pylal_PostcohInspiralTable, postcoh_inspiral.chisq_H), 0, "chisq_H"},
	{"chisq_V", T_FLOAT, offsetof(pylal_PostcohInspiralTable, postcoh_inspiral.chisq_V), 0, "chisq_V"},
	{"is_background", T_INT, offsetof(pylal_PostcohInspiralTable, postcoh_inspiral.is_background), 0, "is_background"},
	{"livetime", T_INT, offsetof(pylal_PostcohInspiralTable, postcoh_inspiral.livetime), 0, "livetime"},
	{"tmplt_idx", T_INT, offsetof(pylal_PostcohInspiralTable, postcoh_inspiral.tmplt_idx), 0, "tmplt_idx"},
	{"pix_idx", T_INT, offsetof(pylal_PostcohInspiralTable, postcoh_inspiral.pix_idx), 0, "pix_idx"},
	{"cohsnr", T_FLOAT, offsetof(pylal_PostcohInspiralTable, postcoh_inspiral.cohsnr), 0, "cohsnr"},
	{"nullsnr", T_FLOAT, offsetof(pylal_PostcohInspiralTable, postcoh_inspiral.nullsnr), 0, "nullsnr"},
	{"cmbchisq", T_FLOAT, offsetof(pylal_PostcohInspiralTable, postcoh_inspiral.cmbchisq), 0, "cmbchisq"},
	{"spearman_pval", T_FLOAT, offsetof(pylal_PostcohInspiralTable, postcoh_inspiral.spearman_pval), 0, "spearman_pval"},
	{"fap", T_FLOAT, offsetof(pylal_PostcohInspiralTable, postcoh_inspiral.fap), 0, "fap"},
	{"fap_h", T_FLOAT, offsetof(pylal_PostcohInspiralTable, postcoh_inspiral.fap_h), 0, "fap_h"},
	{"fap_l", T_FLOAT, offsetof(pylal_PostcohInspiralTable, postcoh_inspiral.fap_l), 0, "fap_l"},
	{"fap_v", T_FLOAT, offsetof(pylal_PostcohInspiralTable, postcoh_inspiral.fap_v), 0, "fap_v"},
	{"far", T_FLOAT, offsetof(pylal_PostcohInspiralTable, postcoh_inspiral.far), 0, "far"},
	{"far_h", T_FLOAT, offsetof(pylal_PostcohInspiralTable, postcoh_inspiral.far_h), 0, "far_h"},
	{"far_l", T_FLOAT, offsetof(pylal_PostcohInspiralTable, postcoh_inspiral.far_l), 0, "far_l"},
	{"far_v", T_FLOAT, offsetof(pylal_PostcohInspiralTable, postcoh_inspiral.far_v), 0, "far_v"},
	{"template_duration", T_DOUBLE, offsetof(pylal_PostcohInspiralTable, postcoh_inspiral.template_duration), 0, "template_duration"},
	{"mass1", T_FLOAT, offsetof(pylal_PostcohInspiralTable, postcoh_inspiral.mass1), 0, "mass1"},
	{"mass2", T_FLOAT, offsetof(pylal_PostcohInspiralTable, postcoh_inspiral.mass2), 0, "mass2"},
	{"mchirp", T_FLOAT, offsetof(pylal_PostcohInspiralTable, postcoh_inspiral.mchirp), 0, "mchirp"},
	{"mtotal", T_FLOAT, offsetof(pylal_PostcohInspiralTable, postcoh_inspiral.mtotal), 0, "mtotal"},
	{"eta", T_FLOAT, offsetof(pylal_PostcohInspiralTable, postcoh_inspiral.eta), 0, "eta"},
	{"spin1x", T_FLOAT, offsetof(pylal_PostcohInspiralTable, postcoh_inspiral.spin1x), 0, "spin1x"},
	{"spin1y", T_FLOAT, offsetof(pylal_PostcohInspiralTable, postcoh_inspiral.spin1y), 0, "spin1y"},
	{"spin1z", T_FLOAT, offsetof(pylal_PostcohInspiralTable, postcoh_inspiral.spin1z), 0, "spin1z"},
	{"spin2x", T_FLOAT, offsetof(pylal_PostcohInspiralTable, postcoh_inspiral.spin2x), 0, "spin2x"},
	{"spin2y", T_FLOAT, offsetof(pylal_PostcohInspiralTable, postcoh_inspiral.spin2y), 0, "spin2y"},
	{"spin2z", T_FLOAT, offsetof(pylal_PostcohInspiralTable, postcoh_inspiral.spin2z), 0, "spin2z"},
	{"ra", T_DOUBLE, offsetof(pylal_PostcohInspiralTable, postcoh_inspiral.ra), 0, "ra"},
	{"dec", T_DOUBLE, offsetof(pylal_PostcohInspiralTable, postcoh_inspiral.dec), 0, "dec"},
	{"deff_L", T_DOUBLE, offsetof(pylal_PostcohInspiralTable, postcoh_inspiral.deff_L), 0, "deff_L"},
	{"deff_H", T_DOUBLE, offsetof(pylal_PostcohInspiralTable, postcoh_inspiral.deff_H), 0, "deff_H"},
	{"deff_V", T_DOUBLE, offsetof(pylal_PostcohInspiralTable, postcoh_inspiral.deff_V), 0, "deff_V"},
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
