/* standard includes */
#include <stdio.h>
#include <math.h>
#include <complex.h>

/* LAL Includes */

#include <lal/LALDatatypes.h>
#include <lal/LALStdlib.h>
#include <lal/Units.h>
#include <lal/LALInspiral.h>


#include <Python.h>
#include <numpy/arrayobject.h>


static int SPAWaveform (double mass1, double mass2, int order, double deltaF, double deltaT, double fLower, double fFinal, int numPoints, complex double *expPsi);
static double chirp_time (double m1, double m2, double fLower, int order);
static double chirp_time_between_f1_and_f2(double m1, double m2, double fLower, double fUpper, int order);
static double schwarz_isco(double m1, double m2);
static double bkl_isco(double m1, double m2);
static double light_ring(double m1, double m2);


const char SPADocstring[] =
"This module wraps SPA inspiral waveform generation and some useful functions "
"related to them.\n\n"
"EXAMPLE CODE: (you could cut and paste this into the interpreter)\n"
"\n"
"from pylal import spawaveform\n"
"import numpy\n"
"import pylab\n"
"import scipy\n"
"import time\n"
"\n"
"# some parameters\n"
"m1 = 25.0\n"
"m2 = 25.0\n"
"order = 7\n"
"fLower = 30.0\n"
"fFinal = spawaveform.ffinal(m1,m2,'light_ring')\n"
"# find next highest powers of 2\n"
"sr = 2**numpy.ceil(numpy.log2(fFinal*2))\n"
"dur = 2**numpy.ceil(numpy.log2(spawaveform.chirptime(m1,m2,order,fLower)))\n"
"deltaF = 1.0 / dur\n"
"deltaT = 1.0 / sr\n"
"\n"
"pylab.figure(1)\n"
"for j, endfreq in enumerate(['schwarz_isco','bkl_isco','light_ring']):\n"
"        fFinal = spawaveform.ffinal(m1,m2,endfreq)\n"
"        print endfreq, fFinal\n"
"        # time stamp vector\n"
"        timestamp = numpy.arange(dur * sr) / sr\n"
"        # initialise data for the template\n"
"        z = numpy.empty(sr * dur, 'complex128')\n"
"        # make a spawaveform\n"
"        spawaveform.waveform(m1, m2, order, deltaF, deltaT, fLower, fFinal, z)\n"
"        z = scipy.ifft(z)\n"
"        # plot it\n"
"        pylab.subplot(3,1,j+1)\n"
"        pylab.plot(timestamp, numpy.real(z))\n"
"        pylab.hold(1)\n"
"        pylab.plot(timestamp, numpy.imag(z),'r')\n"
"        pylab.legend([endfreq + ' hc(t)', endfreq + ' hs(t)'],loc='lower left')\n"
"        pylab.xlim([dur - spawaveform.chirptime(m1,m2,order,40.0,fFinal), dur])\n"
"pylab.hold(0)\n"
"pylab.show()\n";


/* static PyObject* SPAWaveformError; */

static PyObject *PyFFinal(PyObject *self, PyObject *args)
	{
	double mass1, mass2;
	const char *s = NULL;
	if(!PyArg_ParseTuple(args, "dd|s", &mass1, &mass2, &s)) return NULL;
	/* Default is schwarz isco */
	if ( !strcmp(s, "schwarz_isco") || !s) return Py_BuildValue("d", schwarz_isco(mass1,mass2));
	if ( !strcmp(s, "bkl_isco") ) return Py_BuildValue("d", bkl_isco(mass1,mass2));
	if ( !strcmp(s, "light_ring") ) return Py_BuildValue("d", light_ring(mass1,mass2));
	PyErr_SetString(PyExc_ValueError, "Unrecognized ending frequency, must be schwarz_isco | bkl_isco | light_ring");
	return NULL;
	}

static PyObject *PySPAWaveform(PyObject *self, PyObject *args) 
	{
	/* Generate a SPA (frequency domain) waveform at a given PN order */
	PyObject *arg9, *py_spa_array;
	double mass1, mass2, deltaF, deltaT, fLower, fFinal;
	int order;
	npy_intp *dims = NULL;
	complex double *data = NULL;
	/* FIXME properly handle references */

	if(!PyArg_ParseTuple(args, "ddiddddO", &mass1, &mass2, &order, &deltaF, &deltaT, &fLower, &fFinal, &arg9)) return NULL;

	/* this gets a contiguous memory numpy array */
        py_spa_array = PyArray_FROM_OTF(arg9, NPY_CDOUBLE, NPY_IN_ARRAY);
	if (py_spa_array == NULL) return NULL;

	/* Actually call the SPA waveform C function */
	/* FIXME no checking of the array dimensions, this could be done in a python wrapper */

	dims = PyArray_DIMS(py_spa_array);
	data = PyArray_DATA(py_spa_array);
	
	SPAWaveform(mass1, mass2, order, deltaF, deltaT, fLower, fFinal, dims[0], data);

	Py_DECREF(py_spa_array);
        Py_INCREF(Py_None);
        return Py_None;
	}

static PyObject *PyChirpTime(PyObject *self, PyObject *args) 
	{

	/* Generate a SPA (frequency domain) waveform at a given PN order */
	double mass1, mass2, fLower, fFinal, time;
	int order;
	fFinal = 0;
	/* FIXME properly handle references */

	if (!PyArg_ParseTuple(args, "ddid|d", &mass1, &mass2, &order, &fLower, &fFinal)) return NULL;
	if (fFinal)
		time = chirp_time_between_f1_and_f2(mass1, mass2, fLower, fFinal, order);
	else
		time = chirp_time(mass1, mass2, fLower, order);
	return Py_BuildValue("d", time);
	}


static struct PyMethodDef methods[] = {
	{"waveform", PySPAWaveform, METH_VARARGS, 
         "This function produces a frequency domain waveform at a "
         "specified mass1, mass2 and PN order.\n\n"
         "spawaveform(m1, m2, order, deltaF, deltaT, fLower, fFinal, signalArray)\n\n"
        },
	{"chirptime", PyChirpTime, METH_VARARGS, 
         "This function calculates the SPA chirptime at a specified mass1, mass2 " 
	 "and PN order between two frequencies.  If the second frequency is omitted "
	 "it is assumed to be infinite.\n\n"
         "chirptime(m1, m2, order, fLower, [fFinal])\n\n"
        },
	{"ffinal", PyFFinal, METH_VARARGS, 
         "This function calculates the Schwarzschild ISCO frequency specified by " 
	 "mass1 and mass2.\n\n"  
         "schwarzisco(m1, m2, ['schwarz_isco'|'bkl_isco'|'light_ring'])\n\n"
        },
	{NULL, NULL, 0, NULL}	
	};

/* FIXME make this function exist in LAL and have the LAL SPA waveform generator call it? */
static int SPAWaveform (double mass1, double mass2, int order, double deltaF, double deltaT, double fLower, double fFinal, int numPoints,  complex double *expPsi)
	{
	double m = mass1 + mass2;
	double eta = mass1 * mass2 / m / m;
	double mu = mass1 * mass2 / m;

	double x1 = pow (LAL_PI * m * LAL_MTSUN_SI * deltaF, -1.0 / 3.0);
	double psi = 0.0;
	double psi0 = 0.0;
	int k = 0;
	int kmin = fLower / deltaF > 1 ? fLower / deltaF : 1;
	int kmax = fFinal / deltaF < numPoints / 2 ? fFinal / deltaF : numPoints / 2;

	const double cannonDist = 1.0; /* Mpc */
	double distNorm = 2.0 * LAL_MRSUN_SI / (cannonDist * 1.0e6 * LAL_PC_SI);
	double tNorm = sqrt ((5.0 * mu) / 96.0) * pow (m / (LAL_PI * LAL_PI), 1.0 / 3.0) * pow (LAL_MTSUN_SI / (double) deltaT, -1.0 / 6.0);

	/* pn constants */
	double c0, c10, c15, c20, c25, c25Log, c30, c30Log, c35, c40P;
	double x;

	/* chebychev coefficents for expansion of sin and cos */
	const double s2 = -0.16605;
	const double s4 = 0.00761;
	const double c2 = -0.49670;
	const double c4 = 0.03705;

	complex double value;

	/* template norm */
	tNorm *= tNorm;
	tNorm *= distNorm * distNorm;

	/* zero output */
	memset (expPsi, 0, numPoints * sizeof (complex double));

	/* Initialize all PN phase coeffs to zero. */
	c0 = c10 = c15 = c20 = c25 = c25Log = c30 = c30Log = c35 = c40P = 0.;

	/* Switch on PN order, set the appropriate phase coeffs for that order */
	switch (order)
		{
		case 8:
			c40P = 3923.0;
		case 7:
			c35 = LAL_PI * (77096675.0 / 254016.0 + eta * 378515.0 / 1512.0 - eta * eta * 74045.0 / 756.0);
		case 6: 
			c30 = 11583231236531.0 / 4694215680.0 - LAL_GAMMA * 6848.0 / 21.0 - LAL_PI * LAL_PI * 640.0 / 3.0 + eta * (LAL_PI * LAL_PI * 2255.0 / 12.0 - 15737765635.0 / 3048192.0) + eta * eta * 76055.0 / 1728.0 - eta * eta * eta * 127825.0 / 1296.0 - 6848.0 * log (4.0) / 21.0;
			c30Log = -6848.0 / 21.0;
		case 5:
			c25 = LAL_PI * 38645.0 / 756.0 - LAL_PI * eta * 65.0 / 9.0;
			c25Log = 3 * c25;
		case 4:
			c20 = 15293365.0 / 508032.0 + eta * (27145.0 / 504.0 + eta * 3085.0 / 72.0);
			c15 = -16 * LAL_PI;
			c10 = 3715.0 / 756.0 + eta * 55.0 / 9.0;
			c0 = 3.0 / (eta * 128.0);
			break;
		default:
			fprintf (stderr, "unknown PN order (use 8 for 4PN 7 for 3.5 PN ... 4 for 2PN ) nothing above 4PN or below 2PN is recognized");
			break;
		}

	/* x1 */
	x1 = pow (LAL_PI * m * LAL_MTSUN_SI * deltaF, -1.0 / 3.0);
	x = x1 * pow ((double) kmin, -1.0 / 3.0);

	psi = c0 * (x * (c20 + x * (c15 + x * (c10 + x * x))) + c25 - c25Log * log (x) + (1.0 / x) * (c30 - c30Log * log (x) + (1.0 / x) * (c35 - (1.0 / x) * c40P * log (x))));
	psi0 = -2 * LAL_PI * (floor (0.5 * psi / LAL_PI));
	
	/* Chirp Time */
	/* This formula works for any PN order, because */
	/* higher order coeffs will be set to zero.     */
	for (k = kmin; k < kmax; ++k)
		{
		double x = x1 * pow ((double) k, -1.0 / 3.0);
		double psi = c0 * (x * (c20 + x * (c15 + x * (c10 + x * x))) + c25 - c25Log * log (x) + (1.0 / x) * (c30 - c30Log * log (x) + (1.0 / x) * (c35 - (1.0 / x) * c40P * log (x))));

		double psi1 = psi + psi0;
		double psi2;

		/* range reduction of psi1 */
		while (psi1 < -LAL_PI)
			{
			psi1 += 2 * LAL_PI;
			psi0 += 2 * LAL_PI;
			}
		while (psi1 > LAL_PI)
			{
			psi1 -= 2 * LAL_PI;
			psi0 -= 2 * LAL_PI;
			}

		/* compute approximate sine and cosine of psi1 */
		if (psi1 < -LAL_PI / 2)
			{
			psi1 = -LAL_PI - psi1;
			psi2 = psi1 * psi1;
			/* XXX minus sign added because of new sign convention for fft */
			/* FIXME minus sign put back because it makes a reverse chirp with scipy's ifft */
			value = psi1 * (1 + psi2 * (s2 + psi2 * s4)) +  I * (0. - 1. - psi2 * (c2 + psi2 * c4));
			expPsi[k]  = value;
			}
		else if (psi1 > LAL_PI / 2)
			{
			psi1 = LAL_PI - psi1;
			psi2 = psi1 * psi1;
			/* XXX minus sign added because of new sign convention for fft */
			/* FIXME minus sign put back because it makes a reverse chirp with scipy's ifft */
			value = psi1 * (1 + psi2 * (s2 + psi2 * s4)) + I * (0. - 1. - psi2 * (c2 + psi2 * c4));
			expPsi[k] = value;
			}
		else
			{
			psi2 = psi1 * psi1;
			/* XXX minus sign added because of new sign convention for fft */
			/* FIXME minus sign put back because it makes a reverse chirp with scipy's ifft */
			value = psi1 * (1 + psi2 * (s2 + psi2 * s4)) + I * (1. + psi2 * (c2 + psi2 * c4));
			expPsi[k] = value;
			}
		/* put in the first order amplitude factor */
		expPsi[k] *= pow (k, -7.0 / 6.0);
		}
	return 0;
	}

/* A function to compute the time of a spa waveform starting at a given low
 * frequency and going to infinity
 */

static double chirp_time (double m1, double m2, double fLower, int order)
	{

	/* variables used to compute chirp time */
	double c0T, c2T, c3T, c4T, c5T, c6T, c6LogT, c7T;
	double xT, x2T, x3T, x4T, x5T, x6T, x7T, x8T;
	double m = m1 + m2;
	double eta = m1 * m2 / m / m;

	c0T = c2T = c3T = c4T = c5T = c6T = c6LogT = c7T = 0.;

	/* Switch on PN order, set the chirp time coeffs for that order */
	switch (order)
		{
		case 8:
		case 7:
			c7T = LAL_PI * (14809.0 * eta * eta - 75703.0 * eta / 756.0 - 15419335.0 / 127008.0);
		case 6:
			c6T = LAL_GAMMA * 6848.0 / 105.0 - 10052469856691.0 / 23471078400.0 + LAL_PI * LAL_PI * 128.0 / 3.0 + eta * (3147553127.0 / 3048192.0 - LAL_PI * LAL_PI * 451.0 / 12.0) - eta * eta * 15211.0 / 1728.0 + eta * eta * eta * 25565.0 / 1296.0 + log (4.0) * 6848.0 / 105.0;
     			c6LogT = 6848.0 / 105.0;
		case 5:
			c5T = 13.0 * LAL_PI * eta / 3.0 - 7729.0 / 252.0;
		case 4:
			c4T = 3058673.0 / 508032.0 + eta * (5429.0 / 504.0 + eta * 617.0 / 72.0);
			c3T = -32.0 * LAL_PI / 5.0;
			c2T = 743.0 / 252.0 + eta * 11.0 / 3.0;
			c0T = 5.0 * m * LAL_MTSUN_SI / (256.0 * eta);	
			break;
		default:
			fprintf (stderr, "ERROR!!!\n");
			break;
		}

	/* This is the PN parameter v evaluated at the lower freq. cutoff */
	xT = pow (LAL_PI * m * LAL_MTSUN_SI * fLower, 1.0 / 3.0);
	x2T = xT * xT;
	x3T = xT * x2T;
	x4T = x2T * x2T;
	x5T = x2T * x3T;
	x6T = x3T * x3T;
	x7T = x3T * x4T;
	x8T = x4T * x4T;

	/* Computes the chirp time as tC = t(v_low)    */
	/* tC = t(v_low) - t(v_upper) would be more    */
	/* correct, but the difference is negligble.   */

	/* This formula works for any PN order, because */
	/* higher order coeffs will be set to zero.     */

	return c0T * (1 + c2T * x2T + c3T * x3T + c4T * x4T + c5T * x5T + (c6T + c6LogT * log (xT)) * x6T + c7T * x7T) / x8T;
}

/* A convenience function to compute the time between two frequencies */
static double chirp_time_between_f1_and_f2(double m1, double m2, double fLower, double fUpper, int order)
	{
	return chirp_time(m1,m2,fLower,order) - chirp_time(m1,m2,fUpper,order);
	}


static double schwarz_isco(double m1, double m2)
	{
	double m = LAL_MTSUN_SI * (m1 + m2);
	return 1.0 / pow(6.0, 1.5) / m / LAL_PI; 
	}

static double light_ring(double m1, double m2)
	{
	double m = LAL_MTSUN_SI * (m1 + m2);
	return 1.0 / pow(3.0, 1.5) / m / LAL_PI; 
	}

static double bkl_isco(double m1, double m2)
	{
	double q;
	q = (m1 < m2) ? (m1/m2) : (m2/m1);
	return (0.8 * q * q * q - 2.6 * q * q + 2.8 * q + 1.0 ) * schwarz_isco(m1,m2);
	}
	
/* The init function for this module */
void initspawaveform(void) 
	{
	(void) Py_InitModule3("pylal.spawaveform", methods, SPADocstring);
	import_array();
	/* FIXME someday handle errors
	 * SVMError = PyErr_NewException("spawaveform.SPAWaveformError", NULL, NULL);
	 * Py_INCREF(SPAWaveformError);
	 * PyModule_AddObject(m, "SPAWaveformError", SPAWaveformError);
         */
	}
