/*
 * $Id$
 *
 * Copyright (C) 2006  Kipp C. Cannon
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
 *                   Python Wrapper For LAL's Tools Package
 *
 * ============================================================================
 */


#include <Python.h>
#include <structmember.h>
#include <string.h>
#include <numpy/arrayobject.h>
#include <lal/DetectorSite.h>
#include <misc.h>
#include <tools.h>

/* for inspiral-related things */
#include <lal/LIGOMetadataUtils.h>
#include <lal/CoincInspiralEllipsoid.h>
#include <lal/TrigScanEThincaCommon.h>


#define MODULE_NAME "pylal.xlal.tools"


/*
 * ============================================================================
 *
 *                         Attribute Get/Set Helpers
 *
 * ============================================================================
 */


struct inline_string_description {
	Py_ssize_t offset;
	Py_ssize_t length;
};


static PyObject *pylal_inline_string_get(PyObject *obj, void *data)
{
	struct inline_string_description *desc = data;
	char *s = (void *) obj + desc->offset;

	if(strlen(s) >= desc->length) {
		/* something's wrong, obj probably isn't a valid address */
	}

	return PyString_FromString(s);
}


static int pylal_inline_string_set(PyObject *obj, PyObject *val, void *data)
{
	struct inline_string_description *desc = data;
	char *v = PyString_AsString(val);
	char *s = (void *) obj + desc->offset;

	if(!v)
		return -1;
	if(strlen(v) >= desc->length) {
		PyErr_Format(PyExc_ValueError, "string too long \'%s\'", v);
		return -1;
	}

	strncpy(s, v, desc->length - 1);
	s[desc->length - 1] = '\0';

	return 0;
}


struct ilwdchar_id_description {
	Py_ssize_t offset;
	PyObject **id_type;
};


static int pylal_ilwdchar_id_set(PyObject *obj, PyObject *val, void *data)
{
	const struct ilwdchar_id_description *desc = (struct ilwdchar_id_description *) data;
	void *location = (void *) obj + desc->offset;
	long i = PyInt_AsLong(val);

	if(PyErr_Occurred())
		return -1;

	if((PyObject *) val->ob_type != *desc->id_type) {
		PyErr_SetObject(PyExc_TypeError, val);
		return -1;
	}

	*(long *) location = i;

	return 0;
}


static PyObject *pylal_ilwdchar_id_get(PyObject *obj, void *data)
{
	const struct ilwdchar_id_description *desc = (struct ilwdchar_id_description *) data;
	void *location = (void *) obj + desc->offset;

	return PyObject_CallFunction(*desc->id_type, "l", *(long *) location);
}


/*
 * ============================================================================
 *
 *                              LALDetector Type
 *
 * ============================================================================
 */


/*
 * Member access
 */


static struct PyMemberDef pylal_LALDetector_members[] = {
	{"name", T_STRING_INPLACE, offsetof(pylal_LALDetector, detector.frDetector.name), READONLY, "name"},
	{"prefix", T_STRING_INPLACE, offsetof(pylal_LALDetector, detector.frDetector.prefix), READONLY, "prefix"},
	{"vertexLongitudeRadians", T_DOUBLE, offsetof(pylal_LALDetector, detector.frDetector.vertexLongitudeRadians), READONLY, "vertexLongitudeRadians"},
	{"vertexLatitudeRadians", T_DOUBLE, offsetof(pylal_LALDetector, detector.frDetector.vertexLatitudeRadians), READONLY, "vertexLatitudeRadians"},
	{"vertexElevation", T_FLOAT, offsetof(pylal_LALDetector, detector.frDetector.vertexElevation), READONLY, "vertexElevation"},
	{"xArmAltitudeRadians", T_FLOAT, offsetof(pylal_LALDetector, detector.frDetector.xArmAltitudeRadians), READONLY, "xArmAltitudeRadians"},
	{"xArmAzimuthRadians", T_FLOAT, offsetof(pylal_LALDetector, detector.frDetector.xArmAzimuthRadians), READONLY, "xArmAzimuthRadians"},
	{"yArmAltitudeRadians", T_FLOAT, offsetof(pylal_LALDetector, detector.frDetector.yArmAltitudeRadians), READONLY, "yArmAltitudeRadians"},
	{"yArmAzimuthRadians", T_FLOAT, offsetof(pylal_LALDetector, detector.frDetector.yArmAzimuthRadians), READONLY, "yArmAzimuthRadians"},
	{"xArmMidpoint", T_FLOAT, offsetof(pylal_LALDetector, detector.frDetector.xArmMidpoint), READONLY, "xArmMidpoint"},
	{"yArmMidpoint", T_FLOAT, offsetof(pylal_LALDetector, detector.frDetector.yArmMidpoint), READONLY, "yArmMidpoint"},
	{"location", T_OBJECT, offsetof(pylal_LALDetector, location), READONLY, "location"},
	{"response", T_OBJECT, offsetof(pylal_LALDetector, response), READONLY, "response"},
	{NULL,}
};


/*
 * Type
 */


PyTypeObject pylal_LALDetector_Type = {
	PyObject_HEAD_INIT(NULL)
	.tp_basicsize = sizeof(pylal_LALDetector),
	.tp_doc = "LALDetector structure",
	.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_CHECKTYPES,
	.tp_members = pylal_LALDetector_members,
	.tp_name = MODULE_NAME ".LALDetector",
	.tp_new = PyType_GenericNew,
};


/*
 * ============================================================================
 *
 *                           SnglInspiralTable Type
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


static struct PyMemberDef pylal_SnglInspiralTable_members[] = {
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
	{NULL,}
};


static int pylal_SnglInspiralTable_event_id_set(PyObject *obj, PyObject *val, void *data)
{
	pylal_SnglInspiralTable *row = (pylal_SnglInspiralTable *) obj;
	long i = PyInt_AsLong(val);

	if(PyErr_Occurred())
		return -1;

	if((PyObject *) val->ob_type != sngl_inspiral_event_id_type) {
		PyErr_SetObject(PyExc_TypeError, val);
		return -1;
	}

	row->event_id.id = i;

	return 0;
}


static PyObject *pylal_SnglInspiralTable_event_id_get(PyObject *obj, void *data)
{
	pylal_SnglInspiralTable *row = (pylal_SnglInspiralTable *) obj;

	return PyObject_CallFunction(sngl_inspiral_event_id_type, "l", (long) row->event_id.id);
}


static struct PyGetSetDef pylal_SnglInspiralTable_getset[] = {
	{"ifo", pylal_inline_string_get, pylal_inline_string_set, "ifo", &(struct inline_string_description) {offsetof(pylal_SnglInspiralTable, sngl_inspiral.ifo), LIGOMETA_IFO_MAX}},
	{"search", pylal_inline_string_get, pylal_inline_string_set, "search", &(struct inline_string_description) {offsetof(pylal_SnglInspiralTable, sngl_inspiral.search), LIGOMETA_SEARCH_MAX}},
	{"channel", pylal_inline_string_get, pylal_inline_string_set, "channel", &(struct inline_string_description) {offsetof(pylal_SnglInspiralTable, sngl_inspiral.channel), LIGOMETA_CHANNEL_MAX}},
	{"process_id", pylal_ilwdchar_id_get, pylal_ilwdchar_id_set, "process_id", &(struct ilwdchar_id_description) {offsetof(pylal_SnglInspiralTable, process_id_i), &process_id_type}},
	{"event_id", pylal_SnglInspiralTable_event_id_get, pylal_SnglInspiralTable_event_id_set, "event_id", NULL},
	{NULL,}
};


/*
 * Methods
 */


static PyObject *pylal_SnglInspiralTable___new__(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
	/* call the generic __new__() */
	pylal_SnglInspiralTable *new = (pylal_SnglInspiralTable *) PyType_GenericNew(type, args, kwds);

	if(!new)
		return NULL;

	/* link the event_id pointer in the sngl_inspiral table structure
	 * to the event_id structure */
	new->sngl_inspiral.event_id = &new->event_id;

	/* done */
	return (PyObject *) new;
}


/*
 * Type
 */


PyTypeObject pylal_SnglInspiralTable_Type = {
	PyObject_HEAD_INIT(NULL)
	.tp_basicsize = sizeof(pylal_SnglInspiralTable),
	.tp_doc = "LAL's SnglInspiralTable structure",
	.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_CHECKTYPES,
	.tp_members = pylal_SnglInspiralTable_members,
	.tp_getset = pylal_SnglInspiralTable_getset,
	.tp_name = MODULE_NAME ".SnglInspiralTable",
	.tp_new = pylal_SnglInspiralTable___new__,
};


/*
 * ============================================================================
 *
 *                           SimInspiralTable Type
 *
 * ============================================================================
 */


/*
 * Cached ID types
 */


static PyObject *sim_inspiral_simulation_id_type = NULL;


/*
 * Member access
 */


static struct PyMemberDef pylal_SimInspiralTable_members[] = {
	{"geocent_end_time", T_INT, offsetof(pylal_SimInspiralTable, sim_inspiral.geocent_end_time.gpsSeconds), 0, "geocent_end_time"},
	{"geocent_end_time_ns", T_INT, offsetof(pylal_SimInspiralTable, sim_inspiral.geocent_end_time.gpsNanoSeconds), 0, "geocent_end_time_ns"},
	{"h_end_time", T_INT, offsetof(pylal_SimInspiralTable, sim_inspiral.h_end_time.gpsSeconds), 0, "h_end_time"},
	{"h_end_time_ns", T_INT, offsetof(pylal_SimInspiralTable, sim_inspiral.h_end_time.gpsNanoSeconds), 0, "h_end_time_ns"},
	{"l_end_time", T_INT, offsetof(pylal_SimInspiralTable, sim_inspiral.l_end_time.gpsSeconds), 0, "l_end_time"},
	{"l_end_time_ns", T_INT, offsetof(pylal_SimInspiralTable, sim_inspiral.l_end_time.gpsNanoSeconds), 0, "l_end_time_ns"},
	{"g_end_time", T_INT, offsetof(pylal_SimInspiralTable, sim_inspiral.g_end_time.gpsSeconds), 0, "g_end_time"},
	{"g_end_time_ns", T_INT, offsetof(pylal_SimInspiralTable, sim_inspiral.g_end_time.gpsNanoSeconds), 0, "g_end_time_ns"},
	{"t_end_time", T_INT, offsetof(pylal_SimInspiralTable, sim_inspiral.t_end_time.gpsSeconds), 0, "t_end_time"},
	{"t_end_time_ns", T_INT, offsetof(pylal_SimInspiralTable, sim_inspiral.t_end_time.gpsNanoSeconds), 0, "t_end_time_ns"},
	{"v_end_time", T_INT, offsetof(pylal_SimInspiralTable, sim_inspiral.v_end_time.gpsSeconds), 0, "v_end_time"},
	{"v_end_time_ns", T_INT, offsetof(pylal_SimInspiralTable, sim_inspiral.v_end_time.gpsNanoSeconds), 0, "v_end_time_ns"},
	{"end_time_gmst", T_DOUBLE, offsetof(pylal_SimInspiralTable, sim_inspiral.end_time_gmst), 0, "end_time_gmst"},
	{"mass1", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.mass1), 0, "mass1"},
	{"mass2", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.mass2), 0, "mass2"},
	{"mchirp", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.mchirp), 0, "mchirp"},
	{"eta", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.eta), 0, "eta"},
	{"distance", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.distance), 0, "distance"},
	{"longitude", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.longitude), 0, "longitude"},
	{"latitude", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.latitude), 0, "latitude"},
	{"inclination", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.inclination), 0, "inclination"},
	{"coa_phase", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.coa_phase), 0, "coa_phase"},
	{"polarization", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.polarization), 0, "polarization"},
	{"psi0", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.psi0), 0, "psi0"},
	{"psi3", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.psi3), 0, "psi3"},
	{"alpha", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.alpha), 0, "alpha"},
	{"alpha1", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.alpha1), 0, "alpha1"},
	{"alpha2", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.alpha2), 0, "alpha2"},
	{"alpha3", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.alpha3), 0, "alpha3"},
	{"alpha4", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.alpha4), 0, "alpha4"},
	{"alpha5", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.alpha5), 0, "alpha5"},
	{"alpha6", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.alpha6), 0, "alpha6"},
	{"beta", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.beta), 0, "beta"},
	{"spin1x", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.spin1x), 0, "spin1x"},
	{"spin1y", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.spin1y), 0, "spin1y"},
	{"spin1z", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.spin1z), 0, "spin1z"},
	{"spin2x", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.spin2x), 0, "spin2x"},
	{"spin2y", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.spin2y), 0, "spin2y"},
	{"spin2z", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.spin2z), 0, "spin2z"},
	{"theta0", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.theta0), 0, "theta0"},
	{"phi0", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.phi0), 0, "phi0"},
	{"f_lower", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.f_lower), 0, "f_lower"},
	{"f_final", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.f_final), 0, "f_final"},
	{"eff_dist_h", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.eff_dist_h), 0, "eff_dist_h"},
	{"eff_dist_l", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.eff_dist_l), 0, "eff_dist_l"},
	{"eff_dist_g", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.eff_dist_t), 0, "eff_dist_t"},
	{"eff_dist_t", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.eff_dist_g), 0, "eff_dist_g"},
	{"eff_dist_v", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.eff_dist_v), 0, "eff_dist_v"},
	{"numrel_mode_min", T_INT, offsetof(pylal_SimInspiralTable, sim_inspiral.numrel_mode_min), 0, "numrel_mode_min"},
	{"numrel_mode_max", T_INT, offsetof(pylal_SimInspiralTable, sim_inspiral.numrel_mode_max), 0, "numrel_mode_max"},
	{"amp_order", T_INT, offsetof(pylal_SimInspiralTable, sim_inspiral.amp_order), 0, "amp_order"},
	{"bandpass", T_INT, offsetof(pylal_SimInspiralTable, sim_inspiral.bandpass), 0, "bandpass"},
	{NULL,}
};


static int pylal_SimInspiralTable_simulation_id_set(PyObject *obj, PyObject *val, void *data)
{
	pylal_SimInspiralTable *row = (pylal_SimInspiralTable *) obj;
	long i = PyInt_AsLong(val);

	if(PyErr_Occurred())
		return -1;

	if((PyObject *) val->ob_type != sim_inspiral_simulation_id_type) {
		PyErr_SetObject(PyExc_TypeError, val);
		return -1;
	}

	row->event_id.id = i;

	return 0;
}


static PyObject *pylal_SimInspiralTable_simulation_id_get(PyObject *obj, void *data)
{
	pylal_SimInspiralTable *row = (pylal_SimInspiralTable *) obj;

	return PyObject_CallFunction(sim_inspiral_simulation_id_type, "l", (long) row->event_id.id);
}


static struct PyGetSetDef pylal_SimInspiralTable_getset[] = {
	{"waveform", pylal_inline_string_get, pylal_inline_string_set, "waveform", &(struct inline_string_description) {offsetof(pylal_SimInspiralTable, sim_inspiral.waveform), LIGOMETA_WAVEFORM_MAX}},
	{"source", pylal_inline_string_get, pylal_inline_string_set, "source", &(struct inline_string_description) {offsetof(pylal_SimInspiralTable, sim_inspiral.source), LIGOMETA_SOURCE_MAX}},
	{"numrel_data", pylal_inline_string_get, pylal_inline_string_set, "numrel_data", &(struct inline_string_description) {offsetof(pylal_SimInspiralTable, sim_inspiral.numrel_data), LIGOMETA_STRING_MAX}},
	{"taper", pylal_inline_string_get, pylal_inline_string_set, "taper", &(struct inline_string_description) {offsetof(pylal_SimInspiralTable, sim_inspiral.taper), LIGOMETA_INSPIRALTAPER_MAX}},
	{"process_id", pylal_ilwdchar_id_get, pylal_ilwdchar_id_set, "process_id", &(struct ilwdchar_id_description) {offsetof(pylal_SimInspiralTable, process_id_i), &process_id_type}},
	{"simulation_id", pylal_SimInspiralTable_simulation_id_get, pylal_SimInspiralTable_simulation_id_set, "simulation_id", NULL},
	{NULL,}
};


/*
 * Methods
 */


static PyObject *pylal_SimInspiralTable___new__(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
	/* call the generic __new__() */
	pylal_SimInspiralTable *new = (pylal_SimInspiralTable *) PyType_GenericNew(type, args, kwds);

	if(!new)
		return NULL;

	/* link the event_id pointer in the sngl_inspiral table structure
	 * to the event_id structure */
	new->sim_inspiral.event_id = &new->event_id;

	/* done */
	return (PyObject *) new;
}


/*
 * Type
 */


PyTypeObject pylal_SimInspiralTable_Type = {
	PyObject_HEAD_INIT(NULL)
	.tp_basicsize = sizeof(pylal_SimInspiralTable),
	.tp_doc = "LAL's SimInspiralTable structure",
	.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_CHECKTYPES,
	.tp_members = pylal_SimInspiralTable_members,
	.tp_getset = pylal_SimInspiralTable_getset,
	.tp_name = MODULE_NAME ".SimInspiralTable",
	.tp_new = pylal_SimInspiralTable___new__,
};


/*
 * ============================================================================
 *
 *                               CoincMap Type
 *
 * ============================================================================
 */


/*
 * Structure
 */


static PyObject *coinc_event_id_type = NULL;


typedef struct {
	PyObject_HEAD
	PyObject *event_id_type;
	long event_id_i;
	long coinc_event_id_i;
} ligolw_CoincMap;


/*
 * Attributes
 */


static int ligolw_CoincMap_event_id_set(PyObject *obj, PyObject *val, void *data)
{
	ligolw_CoincMap *coinc_map = (ligolw_CoincMap *) obj;
	long i = PyInt_AsLong(val);

	if(PyErr_Occurred())
		return -1;

	Py_XDECREF(coinc_map->event_id_type);
	coinc_map->event_id_type = (PyObject *) val->ob_type;
	Py_INCREF(coinc_map->event_id_type);
	coinc_map->event_id_i = i;

	return 0;
}


static PyObject *ligolw_CoincMap_event_id_get(PyObject *obj, void *data)
{
	ligolw_CoincMap *coinc_map = (ligolw_CoincMap *) obj;

	if(!coinc_map->event_id_type) {
		Py_INCREF(Py_None);
		return Py_None;
	}

	return PyObject_CallFunction(coinc_map->event_id_type, "l", coinc_map->event_id_i);
}


static int ligolw_CoincMap_table_name_set(PyObject *obj, PyObject *val, void *data)
{
	/* ignored */
	return 0;
}


static PyObject *ligolw_CoincMap_table_name_get(PyObject *obj, void *data)
{
	ligolw_CoincMap *coinc_map = (ligolw_CoincMap *) obj;

	if(!coinc_map->event_id_type) {
		Py_INCREF(Py_None);
		return Py_None;
	}

	return PyObject_GetAttrString(coinc_map->event_id_type, "table_name");
}


static struct PyGetSetDef ligolw_CoincMap_getset[] = {
	{"event_id", ligolw_CoincMap_event_id_get, ligolw_CoincMap_event_id_set, "event_id", NULL},
	{"table_name", ligolw_CoincMap_table_name_get, ligolw_CoincMap_table_name_set, "table_name", NULL},
	{"coinc_event_id", pylal_ilwdchar_id_get, pylal_ilwdchar_id_set, "coinc_event_id", &(struct ilwdchar_id_description) {offsetof(ligolw_CoincMap, coinc_event_id_i), &coinc_event_id_type}},
	{NULL,}
};


/*
 * Methods
 */


static void ligolw_CoincMap___del__(PyObject *self)
{
	ligolw_CoincMap *coinc_map = (ligolw_CoincMap *) self;

	Py_XDECREF(coinc_map->event_id_type);

	self->ob_type->tp_free(self);
}


/*
 * Type
 */


PyTypeObject ligolw_CoincMap_Type = {
	PyObject_HEAD_INIT(NULL)
	.tp_basicsize = sizeof(ligolw_CoincMap),
	.tp_dealloc = ligolw_CoincMap___del__,
	.tp_flags = Py_TPFLAGS_DEFAULT,
	.tp_name = MODULE_NAME ".CoincMap",
	.tp_new = PyType_GenericNew,
	.tp_getset = ligolw_CoincMap_getset,
};


/*
 * ============================================================================
 *
 *                                 Functions
 *
 * ============================================================================
 */


static PyObject *pylal_XLALSnglInspiralTimeError(PyObject *self, PyObject *args)
{
	pylal_SnglInspiralTable *row;
	double e_thinca_threshold;
	double delta_t;

	if(!PyArg_ParseTuple(args, "O!d", &pylal_SnglInspiralTable_Type, &row, &e_thinca_threshold))
		return NULL;

	delta_t = XLALSnglInspiralTimeError(&row->sngl_inspiral, e_thinca_threshold);
	if(XLAL_IS_REAL8_FAIL_NAN(delta_t)) {
		pylal_set_exception_from_xlalerrno();
		return NULL;
	}

	return PyFloat_FromDouble(delta_t);
}


static PyObject *pylal_XLALCalculateEThincaParameter(PyObject *self, PyObject *args)
{
	static InspiralAccuracyList accuracyparams;
	static int accuracyparams_set = 0;
	pylal_SnglInspiralTable *row1, *row2;
	double result;

	if(!PyArg_ParseTuple(args, "O!O!", &pylal_SnglInspiralTable_Type, &row1, &pylal_SnglInspiralTable_Type, &row2))
		return NULL;

	if(!accuracyparams_set) {
		memset(&accuracyparams, 0, sizeof(accuracyparams));
		XLALPopulateAccuracyParams(&accuracyparams);
		accuracyparams_set = 1;
	}

	result = XLALCalculateEThincaParameter(&row1->sngl_inspiral, &row2->sngl_inspiral, &accuracyparams);

	if(XLAL_IS_REAL8_FAIL_NAN(result)) {
		XLALClearErrno();
		PyErr_SetString(PyExc_ValueError, "not coincident");
		return NULL;
	}

	return PyFloat_FromDouble(result);
}


/*
 * ============================================================================
 *
 *                            Module Registration
 *
 * ============================================================================
 */


static PyObject *make_cached_detectors(void)
{
	PyObject *cached_detector = PyDict_New();
	int i;

	for(i = 0; i < LALNumCachedDetectors; i++) {
		pylal_LALDetector *new = (pylal_LALDetector *) _PyObject_New(&pylal_LALDetector_Type);
		memcpy(&new->detector, &lalCachedDetectors[i], sizeof(new->detector));
		{
		npy_intp dims[] = {3};
		new->location = PyArray_SimpleNewFromData(1, dims, NPY_FLOAT64, new->detector.location);
		}
		{
		npy_intp dims[] = {3, 3};
		new->response = PyArray_SimpleNewFromData(2, dims, NPY_FLOAT32, new->detector.response);
		}

		PyDict_SetItemString(cached_detector, new->detector.frDetector.name, (PyObject *) new);
	}

	return cached_detector;
}


static PyObject *get_ilwdchar_class(char *table_name, char *column_name)
{
	PyObject *module_name;
	PyObject *module;
	PyObject *func;
	PyObject *class;

	module_name = PyString_FromString("glue.ligolw.ilwd");
	if(!module_name)
		return NULL;

	module = PyImport_Import(module_name);
	Py_DECREF(module_name);
	if(!module)
		return NULL;

	func = PyMapping_GetItemString(PyModule_GetDict(module), "get_ilwdchar_class");
	Py_DECREF(module);
	if(!func)
		return NULL;

	class = PyObject_CallFunction(func, "ss", table_name, column_name);
	Py_DECREF(func);

	return class;
}


static struct PyMethodDef methods[] = {
	{"XLALSnglInspiralTimeError", pylal_XLALSnglInspiralTimeError, METH_VARARGS, "XLALSnglInspiralTimeError(row, threshold)\n\nFrom a sngl_inspiral event and compute the \\Delta t interval corresponding to the given e-thinca threshold."},
	{"XLALCalculateEThincaParameter", pylal_XLALCalculateEThincaParameter, METH_VARARGS, "XLALCalculateEThincaParameter(row1, row2)\n\nTakes two SnglInspiralTable objects and\ncalculates the overlap factor between them."},
	{NULL,}
};


void inittools(void)
{
	PyObject *module = Py_InitModule3(MODULE_NAME, methods, "Wrapper for LAL's tools package.");

	import_array();

	/* LALDetector */
	if(PyType_Ready(&pylal_LALDetector_Type) < 0)
		return;
	Py_INCREF(&pylal_LALDetector_Type);
	PyModule_AddObject(module, "LALDetector", (PyObject *) &pylal_LALDetector_Type);
	PyModule_AddObject(module, "cached_detector", make_cached_detectors());

	/* SnglInspiralTable */
	if(PyType_Ready(&pylal_SnglInspiralTable_Type) < 0)
		return;
	Py_INCREF(&pylal_SnglInspiralTable_Type);
	PyModule_AddObject(module, "SnglInspiralTable", (PyObject *) &pylal_SnglInspiralTable_Type);
	process_id_type = get_ilwdchar_class("process", "process_id");
	sngl_inspiral_event_id_type = get_ilwdchar_class("sngl_inspiral", "event_id");

	/* SimInspiralTable */
	if(PyType_Ready(&pylal_SimInspiralTable_Type) < 0)
		return;
	Py_INCREF(&pylal_SimInspiralTable_Type);
	PyModule_AddObject(module, "SimInspiralTable", (PyObject *) &pylal_SimInspiralTable_Type);
	sim_inspiral_simulation_id_type = get_ilwdchar_class("sim_inspiral", "simulation_id");

	/* CoincMap */
	if(PyType_Ready(&ligolw_CoincMap_Type) < 0)
		return;
	Py_INCREF(&ligolw_CoincMap_Type);
	PyModule_AddObject(module, "CoincMap", (PyObject *) &ligolw_CoincMap_Type);
	coinc_event_id_type = get_ilwdchar_class("coinc_event", "coinc_event_id");
}
