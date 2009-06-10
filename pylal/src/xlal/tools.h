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
#include <lal/LALDetectors.h>
#include <lal/LIGOMetadataTables.h>


/*
 * ============================================================================
 *
 *                              LALDetector Type
 *
 * ============================================================================
 */


/*
 * Type
 */


extern PyTypeObject pylal_LALDetector_Type;


/*
 * Structure
 */


typedef struct {
	PyObject_HEAD
	LALDetector detector;
	PyObject *location;
	PyObject *response;
} pylal_LALDetector;


/*
 * ============================================================================
 *
 *                           SnglInspiralTable Type
 *
 * ============================================================================
 */


/*
 * Type
 */


extern PyTypeObject pylal_SnglInspiralTable_Type;


/*
 * Structure
 */


typedef struct {
	PyObject_HEAD
	SnglInspiralTable sngl_inspiral;
	/* FIXME:  these should be incorporated into the LAL structure */
	long process_id_i;
	EventIDColumn event_id;
} pylal_SnglInspiralTable;


/*
 * ============================================================================
 *
 *                           SimInspiralTable Type
 *
 * ============================================================================
 */


/*
 * Type
 */


extern PyTypeObject pylal_SimInspiralTable_Type;


/*
 * Structure
 */


typedef struct {
	PyObject_HEAD
	SimInspiralTable sim_inspiral;
	/* FIXME:  these should be incorporated into the LAL structure */
	long process_id_i;
	EventIDColumn event_id;
} pylal_SimInspiralTable;
