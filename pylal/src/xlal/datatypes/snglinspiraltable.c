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
#include <misc.h>
#include <snglinspiraltable.h>


#define MODULE_NAME PYLAL_SNGLINSPIRALTABLE_MODULE_NAME


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


static PyObject *sngl_inspiral_event_id_type = NULL;
static PyObject *process_id_type = NULL;


/*
 * Member access
 */


static struct PyMemberDef members[] = {
	{"end_time", T_INT, offsetof(pylal_SnglInspiralTable, sngl_inspiral.end_time.gpsSeconds), 0, "end_time"},
	{"end_time_ns", T_INT, offsetof(pylal_SnglInspiralTable, sngl_inspiral.end_time.gpsNanoSeconds), 0, "end_time_ns"},
	{"end_time_gmst", T_DOUBLE, offsetof(pylal_SnglInspiralTable, sngl_inspiral.end_time_gmst), 0, "end_time_gmst"},
	{"impulse_time", T_INT, offsetof(pylal_SnglInspiralTable, sngl_inspiral.impulse_time.gpsSeconds), 0, "impulse_time"},
	{"impulse_time_ns", T_INT, offsetof(pylal_SnglInspiralTable, sngl_inspiral.impulse_time.gpsNanoSeconds), 0, "impulse_time_ns"},
	{"template_duration", T_DOUBLE, offsetof(pylal_SnglInspiralTable, sngl_inspiral.template_duration), 0, "template_duration"},
	{"event_duration", T_DOUBLE, offsetof(pylal_SnglInspiralTable, sngl_inspiral.event_duration), 0, "event_duration"},
	{"amplitude", T_FLOAT, offsetof(pylal_SnglInspiralTable, sngl_inspiral.amplitude), 0, "amplitude"},
	{"eff_distance", T_FLOAT, offsetof(pylal_SnglInspiralTable, sngl_inspiral.eff_distance), 0, "eff_distance"},
	{"coa_phase", T_FLOAT, offsetof(pylal_SnglInspiralTable, sngl_inspiral.coa_phase), 0, "coa_phase"},
	{"mass1", T_FLOAT, offsetof(pylal_SnglInspiralTable, sngl_inspiral.mass1), 0, "mass1"},
	{"mass2", T_FLOAT, offsetof(pylal_SnglInspiralTable, sngl_inspiral.mass2), 0, "mass2"},
	{"mchirp", T_FLOAT, offsetof(pylal_SnglInspiralTable, sngl_inspiral.mchirp), 0, "mchirp"},
	{"mtotal", T_FLOAT, offsetof(pylal_SnglInspiralTable, sngl_inspiral.mtotal), 0, "mtotal"},
	{"eta", T_FLOAT, offsetof(pylal_SnglInspiralTable, sngl_inspiral.eta), 0, "eta"},
	{"kappa", T_FLOAT, offsetof(pylal_SnglInspiralTable, sngl_inspiral.kappa), 0, "kappa"},
	{"chi", T_FLOAT, offsetof(pylal_SnglInspiralTable, sngl_inspiral.chi), 0, "chi"},
	{"tau0", T_FLOAT, offsetof(pylal_SnglInspiralTable, sngl_inspiral.tau0), 0, "tau0"},
	{"tau2", T_FLOAT, offsetof(pylal_SnglInspiralTable, sngl_inspiral.tau2), 0, "tau2"},
	{"tau3", T_FLOAT, offsetof(pylal_SnglInspiralTable, sngl_inspiral.tau3), 0, "tau3"},
	{"tau4", T_FLOAT, offsetof(pylal_SnglInspiralTable, sngl_inspiral.tau4), 0, "tau4"},
	{"tau5", T_FLOAT, offsetof(pylal_SnglInspiralTable, sngl_inspiral.tau5), 0, "tau5"},
	{"ttotal", T_FLOAT, offsetof(pylal_SnglInspiralTable, sngl_inspiral.ttotal), 0, "ttotal"},
	{"psi0", T_FLOAT, offsetof(pylal_SnglInspiralTable, sngl_inspiral.psi0), 0, "psi0"},
	{"psi3", T_FLOAT, offsetof(pylal_SnglInspiralTable, sngl_inspiral.psi3), 0, "psi3"},
	{"alpha", T_FLOAT, offsetof(pylal_SnglInspiralTable, sngl_inspiral.alpha), 0, "alpha"},
	{"alpha1", T_FLOAT, offsetof(pylal_SnglInspiralTable, sngl_inspiral.alpha1), 0, "alpha1"},
	{"alpha2", T_FLOAT, offsetof(pylal_SnglInspiralTable, sngl_inspiral.alpha2), 0, "alpha2"},
	{"alpha3", T_FLOAT, offsetof(pylal_SnglInspiralTable, sngl_inspiral.alpha3), 0, "alpha3"},
	{"alpha4", T_FLOAT, offsetof(pylal_SnglInspiralTable, sngl_inspiral.alpha4), 0, "alpha4"},
	{"alpha5", T_FLOAT, offsetof(pylal_SnglInspiralTable, sngl_inspiral.alpha5), 0, "alpha5"},
	{"alpha6", T_FLOAT, offsetof(pylal_SnglInspiralTable, sngl_inspiral.alpha6), 0, "alpha6"},
	{"beta", T_FLOAT, offsetof(pylal_SnglInspiralTable, sngl_inspiral.beta), 0, "beta"},
	{"f_final", T_FLOAT, offsetof(pylal_SnglInspiralTable, sngl_inspiral.f_final), 0, "f_final"},
	{"snr", T_FLOAT, offsetof(pylal_SnglInspiralTable, sngl_inspiral.snr), 0, "snr"},
	{"chisq", T_FLOAT, offsetof(pylal_SnglInspiralTable, sngl_inspiral.chisq), 0, "chisq"},
	{"chisq_dof", T_INT, offsetof(pylal_SnglInspiralTable, sngl_inspiral.chisq_dof), 0, "chisq_dof"},
	{"bank_chisq", T_FLOAT, offsetof(pylal_SnglInspiralTable, sngl_inspiral.bank_chisq), 0, "bank_chisq"},
	{"bank_chisq_dof", T_INT, offsetof(pylal_SnglInspiralTable, sngl_inspiral.bank_chisq_dof), 0, "bank_chisq_dof"},
	{"cont_chisq", T_FLOAT, offsetof(pylal_SnglInspiralTable, sngl_inspiral.cont_chisq), 0, "cont_chisq"},
	{"cont_chisq_dof", T_INT, offsetof(pylal_SnglInspiralTable, sngl_inspiral.cont_chisq_dof), 0, "cont_chisq_dof"},
	{"sigmasq", T_DOUBLE, offsetof(pylal_SnglInspiralTable, sngl_inspiral.sigmasq), 0, "sigmasq"},
	{"rsqveto_duration", T_FLOAT, offsetof(pylal_SnglInspiralTable, sngl_inspiral.rsqveto_duration), 0, "rsqveto_duration"},
	{"Gamma0", T_FLOAT, offsetof(pylal_SnglInspiralTable, sngl_inspiral.Gamma[0]), 0, "Gamma0"},
	{"Gamma1", T_FLOAT, offsetof(pylal_SnglInspiralTable, sngl_inspiral.Gamma[1]), 0, "Gamma1"},
	{"Gamma2", T_FLOAT, offsetof(pylal_SnglInspiralTable, sngl_inspiral.Gamma[2]), 0, "Gamma2"},
	{"Gamma3", T_FLOAT, offsetof(pylal_SnglInspiralTable, sngl_inspiral.Gamma[3]), 0, "Gamma3"},
	{"Gamma4", T_FLOAT, offsetof(pylal_SnglInspiralTable, sngl_inspiral.Gamma[4]), 0, "Gamma4"},
	{"Gamma5", T_FLOAT, offsetof(pylal_SnglInspiralTable, sngl_inspiral.Gamma[5]), 0, "Gamma5"},
	{"Gamma6", T_FLOAT, offsetof(pylal_SnglInspiralTable, sngl_inspiral.Gamma[6]), 0, "Gamma6"},
	{"Gamma7", T_FLOAT, offsetof(pylal_SnglInspiralTable, sngl_inspiral.Gamma[7]), 0, "Gamma7"},
	{"Gamma8", T_FLOAT, offsetof(pylal_SnglInspiralTable, sngl_inspiral.Gamma[8]), 0, "Gamma8"},
	{"Gamma9", T_FLOAT, offsetof(pylal_SnglInspiralTable, sngl_inspiral.Gamma[9]), 0, "Gamma9"},
        {"spin1x", T_FLOAT, offsetof(pylal_SnglInspiralTable, sngl_inspiral.spin1x), 0, "spin1x"},
	{"spin1y", T_FLOAT, offsetof(pylal_SnglInspiralTable, sngl_inspiral.spin1y), 0, "spin1y"},
	{"spin1z", T_FLOAT, offsetof(pylal_SnglInspiralTable, sngl_inspiral.spin1z), 0, "spin1z"},
	{"spin2x", T_FLOAT, offsetof(pylal_SnglInspiralTable, sngl_inspiral.spin2x), 0, "spin2x"},
	{"spin2y", T_FLOAT, offsetof(pylal_SnglInspiralTable, sngl_inspiral.spin2y), 0, "spin2y"},
	{"spin2z", T_FLOAT, offsetof(pylal_SnglInspiralTable, sngl_inspiral.spin2z), 0, "spin2z"},
	{NULL,}
};


static struct PyGetSetDef getset[] = {
	{"ifo", pylal_inline_string_get, pylal_inline_string_set, "ifo", &(struct pylal_inline_string_description) {offsetof(pylal_SnglInspiralTable, sngl_inspiral.ifo), LIGOMETA_IFO_MAX}},
	{"search", pylal_inline_string_get, pylal_inline_string_set, "search", &(struct pylal_inline_string_description) {offsetof(pylal_SnglInspiralTable, sngl_inspiral.search), LIGOMETA_SEARCH_MAX}},
	{"channel", pylal_inline_string_get, pylal_inline_string_set, "channel", &(struct pylal_inline_string_description) {offsetof(pylal_SnglInspiralTable, sngl_inspiral.channel), LIGOMETA_CHANNEL_MAX}},
	{"process_id", pylal_ilwdchar_id_get, pylal_ilwdchar_id_set, "process_id", &(struct pylal_ilwdchar_id_description) {offsetof(pylal_SnglInspiralTable, process_id_i), &process_id_type}},
	{"event_id", pylal_ilwdchar_id_get, pylal_ilwdchar_id_set, "event_id", &(struct pylal_ilwdchar_id_description) {offsetof(pylal_SnglInspiralTable, event_id.id), &sngl_inspiral_event_id_type}},
	{NULL,}
};


#if PY_VERSION_HEX < 0x02050000
static int getreadbuffer(PyObject *self, int segment, void **ptrptr)
#else
static Py_ssize_t getreadbuffer(PyObject *self, Py_ssize_t segment, void **ptrptr)
#endif
{
	*ptrptr = &((pylal_SnglInspiralTable*)self)->sngl_inspiral;
	return sizeof(((pylal_SnglInspiralTable*)self)->sngl_inspiral);
}


#if PY_VERSION_HEX < 0x02050000
static int getsegcount(PyObject *self, int *lenp)
#else
static Py_ssize_t getsegcount(PyObject *self, Py_ssize_t *lenp)
#endif
{
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
	pylal_SnglInspiralTable *new = (pylal_SnglInspiralTable *) PyType_GenericNew(type, args, kwds);

	if(!new)
		return NULL;

	/* link the event_id pointer in the sngl_inspiral table structure
	 * to the event_id structure */
	new->sngl_inspiral.event_id = &new->event_id;

	new->process_id_i = 0;
	new->event_id.id = 0;

	/* done */
	return (PyObject *) new;
}


/*
 * Type
 */


static PyTypeObject pylal_snglinspiraltable_type = {
	PyObject_HEAD_INIT(NULL)
	.tp_basicsize = sizeof(pylal_SnglInspiralTable),
	.tp_doc = "LAL's SnglInspiralTable structure",
	.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_CHECKTYPES,
	.tp_members = members,
	.tp_getset = getset,
	.tp_as_buffer = &as_buffer,
	.tp_name = MODULE_NAME ".SnglInspiralTable",
	.tp_new = __new__,
};


/*
 * ============================================================================
 *
 *                                 Functions
 *
 * ============================================================================
 */


static PyObject *from_buffer(PyObject *self, PyObject *args)
{
	const SnglInspiralTable *data;
#if PY_VERSION_HEX < 0x02050000
	int length;
#else
	Py_ssize_t length;
#endif
	unsigned i;
	PyObject *result;

	if(!PyArg_ParseTuple(args, "s#", (const char **) &data, &length))
		return NULL;

	if(length % sizeof(SnglInspiralTable)) {
		PyErr_SetString(PyExc_ValueError, "buffer size is not an integer multiple of SnglInspiralTable struct size");
		return NULL;
	}
	length /= sizeof(SnglInspiralTable);

	result = PyTuple_New(length);
	if(!result)
		return NULL;
	for(i = 0; i < length; i++) {
		PyObject *item = pylal_SnglInspiralTable_new(data++);
		if(!item) {
			Py_DECREF(result);
			return NULL;
		}
		PyTuple_SET_ITEM(result, i, item);
	}

	return result;
}


/*
 * ============================================================================
 *
 *                            Module Registration
 *
 * ============================================================================
 */


static struct PyMethodDef functions[] = {
	{"from_buffer", from_buffer, METH_VARARGS, "Construct a tuple of SnglInspiralTable objects from a buffer object.  The buffer is interpreted as a C array of SnglInspiralTable structures."},
	{NULL, }
};


PyMODINIT_FUNC initsnglinspiraltable(void)
{
	PyObject *module = Py_InitModule3(MODULE_NAME, functions, "Wrapper for LAL's SnglInspiralTable type.");

	/* Cached ID types */
	process_id_type = pylal_get_ilwdchar_class("process", "process_id");
	sngl_inspiral_event_id_type = pylal_get_ilwdchar_class("sngl_inspiral", "event_id");

	/* SnglInspiralTable */
	_pylal_SnglInspiralTable_Type = &pylal_snglinspiraltable_type;
	if(PyType_Ready(&pylal_SnglInspiralTable_Type) < 0)
		return;
	Py_INCREF(&pylal_SnglInspiralTable_Type);
	PyModule_AddObject(module, "SnglInspiralTable", (PyObject *) &pylal_SnglInspiralTable_Type);
}
