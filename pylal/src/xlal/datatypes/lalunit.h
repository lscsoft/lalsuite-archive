/*
 * Copyright (C) 2009  Kipp Cannon
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
 *                              Helper Utilities
 *
 * ============================================================================
 */


#include <Python.h>
#include <lal/LALDatatypes.h>


/*
 * ============================================================================
 *
 *                                  LALUnit
 *
 * ============================================================================
 */


extern PyTypeObject pylal_LALUnit_Type;


typedef struct {
	PyObject_HEAD
	LALUnit unit;
} pylal_LALUnit;


PyObject *pylal_LALUnit_new(int power_of_ten, LALUnit unit);


extern PyObject *pylal_LALUnitMeter;
extern PyObject *pylal_LALUnitKiloGram;
extern PyObject *pylal_LALUnitSecond;
extern PyObject *pylal_LALUnitAmpere;
extern PyObject *pylal_LALUnitKelvin;
extern PyObject *pylal_LALUnitStrain;
extern PyObject *pylal_LALUnitADCCount;
