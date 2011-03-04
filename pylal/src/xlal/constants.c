/*
 * Copyright (C) 2010 Nickolas Fotopoulos
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
 *                   Python Wrapper For LAL's constants header
 *
 * ============================================================================
 */

#include <Python.h>
#include <lal/LALConstants.h>

#define MODULE_NAME "pylal.xlal.constants"

/*
 *  Function to add a C double as a Python float to a module
 */
static int add_const_double(PyObject *m, const char *name, double val) {
    PyObject *val_obj = PyFloat_FromDouble(val);
    if (val_obj && !PyModule_AddObject(m, name, val_obj)) return 0;
    Py_DECREF(val_obj);
    return -1;
}
#define add_double_macro(m, c) add_const_double(m, #c, c)

PyMODINIT_FUNC initconstants(void) {
    PyObject *m;

    m = Py_InitModule3("constants", NULL, "Mathematical and physical "\
        "constants from LAL");
    if (m == NULL)
        return;

    /* The following lines were generated from:
    grep -e "^#define LAL_" LALConstants.h \
      | awk '{ print "add_double_macro(m, " $2 ");" }'
    Then a few were corrected to integers as necessary.
    */
    PyModule_AddIntMacro(m, LAL_REAL4_MANT);
    add_double_macro(m, LAL_REAL4_MAX);
    add_double_macro(m, LAL_REAL4_MIN);
    add_double_macro(m, LAL_REAL4_EPS);
    PyModule_AddIntMacro(m, LAL_REAL8_MANT);
    add_double_macro(m, LAL_REAL8_MAX);
    add_double_macro(m, LAL_REAL8_MIN);
    add_double_macro(m, LAL_REAL8_EPS);
    add_double_macro(m, LAL_E);
    add_double_macro(m, LAL_LOG2E);
    add_double_macro(m, LAL_LOG10E);
    add_double_macro(m, LAL_LN2);
    add_double_macro(m, LAL_LN10);
    add_double_macro(m, LAL_SQRT2);
    add_double_macro(m, LAL_SQRT1_2);
    add_double_macro(m, LAL_GAMMA);
    add_double_macro(m, LAL_EXPGAMMA);
    add_double_macro(m, LAL_PI);
    add_double_macro(m, LAL_TWOPI);
    add_double_macro(m, LAL_PI_2);
    add_double_macro(m, LAL_PI_4);
    add_double_macro(m, LAL_1_PI);
    add_double_macro(m, LAL_2_PI);
    add_double_macro(m, LAL_2_SQRTPI);
    add_double_macro(m, LAL_PI_180);
    add_double_macro(m, LAL_180_PI);
    PyModule_AddIntMacro(m, LAL_C_SI);
    add_double_macro(m, LAL_EPSILON0_SI);
    add_double_macro(m, LAL_MU0_SI);
    add_double_macro(m, LAL_GEARTH_SI);
    PyModule_AddIntMacro(m, LAL_PATM_SI);
    PyModule_AddIntMacro(m, LAL_YRJUL_SI);
    add_double_macro(m, LAL_LYR_SI);
    add_double_macro(m, LAL_G_SI);
    add_double_macro(m, LAL_H_SI);
    add_double_macro(m, LAL_HBAR_SI);
    add_double_macro(m, LAL_MPL_SI);
    add_double_macro(m, LAL_LPL_SI);
    add_double_macro(m, LAL_TPL_SI);
    add_double_macro(m, LAL_K_SI);
    add_double_macro(m, LAL_R_SI);
    add_double_macro(m, LAL_MOL);
    add_double_macro(m, LAL_BWIEN_SI);
    add_double_macro(m, LAL_SIGMA_SI);
    add_double_macro(m, LAL_AMU_SI);
    add_double_macro(m, LAL_MP_SI);
    add_double_macro(m, LAL_ME_SI);
    add_double_macro(m, LAL_QE_SI);
    add_double_macro(m, LAL_ALPHA);
    add_double_macro(m, LAL_RE_SI);
    add_double_macro(m, LAL_LAMBDAE_SI);
    add_double_macro(m, LAL_AB_SI);
    add_double_macro(m, LAL_MUB_SI);
    add_double_macro(m, LAL_MUN_SI);
    add_double_macro(m, LAL_REARTH_SI);
    add_double_macro(m, LAL_AWGS84_SI);
    add_double_macro(m, LAL_BWGS84_SI);
    add_double_macro(m, LAL_MEARTH_SI);
    add_double_macro(m, LAL_IEARTH);
    add_double_macro(m, LAL_EEARTH);
    add_double_macro(m, LAL_RSUN_SI);
    add_double_macro(m, LAL_MSUN_SI);
    add_double_macro(m, LAL_MRSUN_SI);
    add_double_macro(m, LAL_MTSUN_SI);
    add_double_macro(m, LAL_LSUN_SI);
    add_double_macro(m, LAL_AU_SI);
    add_double_macro(m, LAL_PC_SI);
    add_double_macro(m, LAL_YRTROP_SI);
    add_double_macro(m, LAL_YRSID_SI);
    add_double_macro(m, LAL_DAYSID_SI);
    add_double_macro(m, LAL_H0FAC_SI);
    add_double_macro(m, LAL_H0_SI);
    add_double_macro(m, LAL_RHOCFAC_SI);
    add_double_macro(m, LAL_RHOC_SI);
    add_double_macro(m, LAL_TCBR_SI);
    add_double_macro(m, LAL_VCBR_SI);
    add_double_macro(m, LAL_RHOCBR_SI);
    add_double_macro(m, LAL_NCBR_SI);
    add_double_macro(m, LAL_SCBR_SI);
}
