/*                                                                      
 *  Copyright (C) 2012 Nickolas Fotopoulos, Leo Singer, Alexander Dietz
 *                                                                     
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or   
 *  (at your option) any later version.                                 
 *                                                                      
 *  This program is distributed in the hope that it will be useful,     
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of      
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       
 *  GNU General Public License for more details.                        
 *                                                                      
 *  You should have received a copy of the GNU General Public License   
 *  along with with program; see the file COPYING. If not, write to the 
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, 
 *  MA  02111-1307  USA                                                 
 */


#include <Python.h>
#include <numpy/arrayobject.h>
#include <omp.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <lal/DetResponse.h>
#include <lal/LALConstants.h>
#include <lal/LALDetectors.h>
#include <time.h>
#include <ctype.h>
void str2network(LALDetector network[], size_t net_size, char *str);
static double rchisq_2(gsl_rng *rng, double lambda);
void Ssq(double *S2, gsl_rng *rng, double beam_fac, double *response, size_t network_size);

void simulate(double *efficiency, LALDetector *network, size_t network_size, double* BNS_horizon_dist, double jet_semiangle_deg, double m1_min, double m1_max, double m2_min, double m2_max, double rho_thresh, double *distances, size_t dist_bins, size_t samples){

  const static double beam_fac=1-cos(jet_semiangle_deg)*M_PI/180.0;
  memset(efficiency, 0, dist_bins*sizeof(double));
  const static double twopi=2 * M_PI;
  const static double BNSchirp=pow( pow(1.4*1.4,3.0)/(1.4+1.4), 0.2);
      

#pragma omp parallel
  {
    /*Create and seed random number generator*/
    gsl_rng *rng=gsl_rng_alloc(gsl_rng_mt19937);
#ifdef _OPENMP
    gsl_rng_set(rng, omp_get_thread_num()*time(NULL));//Seed with thread number.
#else
    gsl_rng_set(rng, time(NULL));
#endif

    double *response=malloc(2 * network_size * sizeof(double));

    /*Total detections for each distance step for the given thread.*/
    long *threadTotals=calloc(dist_bins, sizeof(long));

    /*Preallocate variables needed inside loops*/
    size_t j, k, l;
    /*size_t objects are inherently the size of a pointer on the system.*/

#pragma omp for schedule(static)
    for(k=1; k<=samples; k++){
      /* draw orientations and determine network response */
      const double lon=twopi*gsl_rng_uniform(rng);
      const double lat=M_PI_2 - acos(2*gsl_rng_uniform(rng) - 1);
      const double zeta=twopi*gsl_rng_uniform(rng);
      const double m1=m1_min+gsl_rng_uniform(rng)*(m1_max-m1_min);
      const double m2=m2_min+gsl_rng_uniform(rng)*(m2_max-m2_min);
      const double mass_correction=pow( pow(m1*m2,3.0)/(m1+m2), 0.2)/BNSchirp;
      const double rho_thresh_2=rho_thresh*rho_thresh;

      for(l=0; l<network_size; l++)
	XLALComputeDetAMResponse(response+2*l, response+2*l+1, network[l].response, lon, lat, zeta, 0.0);

      double S2[network_size];
      Ssq(S2, rng, beam_fac, response, network_size);
      double lambda_const[network_size];
      for(l=network_size; l--;)
	/* Some of the values in lambda are constant for a given mass pair/detector
	   regardless of distance. These have been pre-evaluated for speed */
	lambda_const[l]=(rho_thresh*rho_thresh-2) * BNS_horizon_dist[l]*BNS_horizon_dist[l] * S2[l] * mass_correction;

      for(j=dist_bins; j--;){
	/* Require two or more detectors to detect a signal
	   for a valid detection to be added to the total. */
	int successes=0;
	for(l=network_size; l--;){
	  const double lambda=(1/distances[j]) * (1/distances[j]) * lambda_const[l];
	  double rand=rchisq_2(rng,lambda);
	  successes+=rand>rho_thresh*rho_thresh;
	}
	threadTotals[j]+=(successes>=2);
      }
    }
    
#pragma omp critical
    for(j=dist_bins; j--;)
      efficiency[j]+=threadTotals[j]/((double)samples);

    gsl_rng_free(rng);
    free(threadTotals);
    free(response);
  }
}

static double rchisq_2(gsl_rng *rng, double lambda) {
  /*Generates a random value from a degree 2
    chi_squared with non-centrality parameter of lambda*/
  double a=sqrt(0.5*lambda);
  const double temp=gsl_ran_gaussian(rng, 1.0)+a;
  const double temp2=gsl_ran_gaussian(rng, 1.0)+a;
  return (temp*temp)+(temp2*temp2);
}

void Ssq(double *S2, gsl_rng *rng, double beam_fac, double *response, size_t network_size){
  /* Calculates the antenna factor for each detector for a random source orientation*/
  const double cosiota=1-beam_fac*gsl_rng_uniform(rng);//beam_fac determines max iota.
  const double cosiotasq=cosiota*cosiota;
  const double iotafac=0.25*(1+cosiotasq)*(1+cosiotasq);

  for(size_t l=0; l<network_size; l++) {
    double fplus=response[2*l], fcross=response[2*l+1];
    S2[l]=(fplus*fplus)*iotafac+(fcross*fcross)*cosiotasq;
  }
}

void str2network(LALDetector network[LAL_NUM_DETECTORS], size_t net_size, char *str){
  /*Convert string like "HLVK" to an array of LALDetectors.
    Return the size of the network.*/
  size_t k=0;
  while (k < net_size) {
    /* fprintf(stderr, "Detector '%c' horizon distance: %g\n", str[k], reach[k]); */
    if (str[k]=='H') {
      network[k++]=lalCachedDetectors[LAL_LHO_4K_DETECTOR];
    }
    else if (str[k]=='L') {
      network[k++]=lalCachedDetectors[LAL_LLO_4K_DETECTOR];
    }
    else if (str[k]=='V') {
      network[k++]=lalCachedDetectors[LAL_VIRGO_DETECTOR];
    }
    else if (str[k]=='K') {
      /* numbers from private communication with Koji Arai */
      LALDetector *detector=network+k;
      LALFrDetector *frDetector=&(detector->frDetector);
      strncpy(frDetector->name, "KAGRA", LALNameLength);
      strncpy(frDetector->prefix, "K1", 3);
      frDetector->vertexLatitudeRadians=2.396511595913414;
      frDetector->vertexLongitudeRadians=0.6354743806511354;
      frDetector->vertexElevation=372.0;
      frDetector->xArmAltitudeRadians=0.0;
      frDetector->xArmAzimuthRadians=1.076693615555302;
      frDetector->yArmAltitudeRadians=0.0;
      frDetector->yArmAzimuthRadians=5.789082595939991;
      detector=XLALCreateDetector(network+k, frDetector, LALDETECTORTYPE_IFODIFF);
      if (!detector) {
	fprintf(stderr, "Failed to create KAGRA detector\n");
      }
      k++;
    }
    else if (str[k] == 'I') {
      /* numbers from Schutz 2011 network FOMs */
      LALDetector *detector=network+k;
      LALFrDetector *frDetector=&(detector->frDetector);
      strncpy(frDetector->name, "Indigo", LALNameLength);
      strncpy(frDetector->prefix, "I1", 3);
      frDetector->vertexLatitudeRadians=1.3098647554849334;
      frDetector->vertexLongitudeRadians=0.33329486135237268;
      frDetector->vertexElevation=0.0;
      frDetector->xArmAltitudeRadians=0.0;
      frDetector->xArmAzimuthRadians=3.9269908169872414;
      frDetector->yArmAltitudeRadians=0.0;
      frDetector->yArmAzimuthRadians=5.497787143782138;
      detector=XLALCreateDetector(network+k, frDetector, LALDETECTORTYPE_IFODIFF);
      if (!detector) {
	fprintf(stderr, "Failed to create Indigo detector\n");
      }
      k++;
    }
    else {
      fprintf(stderr, "Unrecognized site: %c\n", str[k]);
    }
  }
}

static PyObject *pylalsimulate(PyObject *self, PyObject *args){
  PyObject *py_network, *py_distances;
  double m1_min, m1_max;
  double m2_min, m2_max;
  double rho_thresh, jet_semiangle_deg;
  size_t samples;

  if(!PyArg_ParseTuple(args, "OdddddOdi", &py_network, &jet_semiangle_deg, &m1_min, &m1_max, &m2_min, &m2_max, &py_distances, &rho_thresh, &samples)) return NULL;

  Py_ssize_t distances_length=PySequence_Length(py_distances);
  double distances[distances_length];
  for(Py_ssize_t i=0; i<distances_length; i++){
    distances[i]=PyFloat_AsDouble(PySequence_GetItem(py_distances, i));
  }

  Py_ssize_t network_length=PySequence_Length(py_network);
  char net_names[network_length];
  double net_hor_dists[network_length];
  for(Py_ssize_t i=0; i<network_length; i++){
    //Get item from network list, get first item from tuple, convert to string, capitalize.
    net_names[i]=toupper( PyString_AsString( PyTuple_GetItem( PySequence_GetItem( py_network, i ), 0))[0]);
    //Get item from network list, get second item from tuple, convert to float
    net_hor_dists[i]=PyFloat_AsDouble(PyTuple_GetItem(PySequence_GetItem(py_network, i),1));
  }


  printf("Building Network...\n");
  LALDetector network[(size_t)network_length];
  str2network(network, (size_t)network_length, net_names);
  printf("Done.\n");

  npy_intp distbins[1]={distances_length};
  PyObject *out_array=PyArray_SimpleNew(1,distbins,NPY_DOUBLE);
  printf("Simulating %li runs.\n", samples);
  //double efficiency[(size_t)distances_length];
  simulate((double*)PyArray_DATA(out_array), network, (size_t)network_length, net_hor_dists, jet_semiangle_deg, m1_min, m1_max, m2_min, m2_max, rho_thresh, distances, (size_t)distances_length, samples);

  return out_array;

}

const char docstring[] = "Simulates a detector network using\n"
  "single threshold detection of waveforms\n"
  "assuming X^2 SNR.\n"
  "========================\n"
  "       Parameters\n"
  "========================\n"
  "Network: Sequence of (Name, Horizon_distance[1.4,1.4 masses]) tuples.\n"
  "jet_semiangle_deg: Maximum jet angle iota (0,90)\n"
  "m1_min: Min mass for m1.\n"
  "m1_max: Max mass for m1.\n"
  "m2_min: Min mass for m2.\n"
  "m2_max: Max mass for m2.\n"
  "distances: Array of distance values to evaluate efficiency at.\n"
  "rho_thresh: SNR cutoff threshold.\n"
  "N: Sample size.\n";

static struct PyMethodDef methods[] = {
  {"simulate", pylalsimulate, METH_VARARGS, docstring},
  {NULL, NULL, 0, NULL}/* Sentinel */
};

void initcbc_network_efficiency(void) {
  (void) Py_InitModule3("pylal.cbc_network_efficiency", methods, docstring);
  import_array();
}
