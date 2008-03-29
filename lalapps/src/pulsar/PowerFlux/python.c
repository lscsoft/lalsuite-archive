#include <Python.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fenv.h>
#include <unistd.h>
#include <time.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "global.h"
#include "hookup.h"
#include "rastermagic.h"
#include "cmdline.h"
#include "fine_grid.h"
#include "lines.h"
#include "grid.h"
#include "polarization.h"
#include "statistics.h"
#include "dataset.h"
#include "candidates.h"
#include "util.h"
#include "jobs.h"

FILE *LOG=NULL, *FILE_LOG=NULL;

char *output_dir=NULL;
int first_bin;
int nsegments;
int useful_bins, nbins,side_cut;

char *earth_ephemeris=NULL, *sun_ephemeris=NULL;
double resolution; /* this is actual resolution, not the resolution argument passed on command line */
int fake_injection=0;

int no_am_response;
int do_CutOff=1;

SKY_GRID *fine_grid=NULL;
SKY_SUPERGRID *super_grid=NULL, *proto_super_grid=NULL;
SKY_GRID *patch_grid=NULL;

char *subinstance_name="python";

double band_axis_norm;
double band_axis[3];

double spindown;
extern INT64 spindown_start;

struct gengetopt_args_info args_info;

char *cmdline=NULL;

extern DATASET *datasets;
extern int d_free;

void *do_alloc(long a, long b)
{
void *r;
int i=0;
if(a<1)a=1;
if(b<1)b=1;
r=calloc(a,b);
while(r==NULL){
	fprintf(stderr,"Could not allocate %ld chunks of %ld bytes each (%ld bytes total), current memory usage %ld\n",a,b,a*b, MEMUSAGE);
	if(i>10)exit(-1);
	condor_safe_sleep(10);
	r=calloc(a,b);
	i++;
	}
return r;
}


static PyObject * powerflux_init(PyObject *self, PyObject *args)
{
char *cmd;
int i;

if (!PyArg_ParseTuple(args, "s", &cmd))
        return NULL;
fprintf(stderr, "Initializing PowerFlux\n");

if(cmdline!=NULL)free(cmdline);
cmdline=strdup(cmd);

fprintf(stderr, "%s\n", cmdline);

if(cmdline_parser_string(cmdline, &args_info, "powerflux"))return NULL;
for(i=0;i<args_info.config_given;i++)
	if(cmdline_parser_configfile(args_info.config_arg[i], &args_info, 0, 0, 0))return NULL;

gsl_rng_env_setup();
gsl_set_error_handler_off();

init_threads(args_info.num_threads_arg);
init_jobs();
init_hookup();
init_statistics();
tabulate_hann_filter7();
init_power_cache();
do_CutOff=args_info.do_cutoff_arg;

if(args_info.earth_ephemeris_given){
	earth_ephemeris=args_info.earth_ephemeris_arg;
	} else {
	earth_ephemeris=do_alloc(strlen(args_info.ephemeris_path_arg)+20,1);
	sprintf(earth_ephemeris,"%s/earth00-04.dat",args_info.ephemeris_path_arg);
	}
	
if(args_info.sun_ephemeris_given){
	sun_ephemeris=args_info.sun_ephemeris_arg;
	} else {
	sun_ephemeris=do_alloc(strlen(args_info.ephemeris_path_arg)+20,1);
	sprintf(sun_ephemeris,"%s/sun00-04.dat",args_info.ephemeris_path_arg);
	}
	
init_ephemeris();

resolution=(4500.0)/(args_info.first_bin_arg+args_info.nbins_arg/2);

resolution*=args_info.skymap_resolution_ratio_arg;
if(args_info.skymap_resolution_given){
	resolution=args_info.skymap_resolution_arg;
	}

fprintf(stderr, "resolution : %f\n", resolution);
if(!strcasecmp("sin_theta", args_info.sky_grid_arg)){
	patch_grid=make_sin_theta_grid(resolution*args_info.fine_factor_arg);
	proto_super_grid=make_sin_theta_supergrid(patch_grid, args_info.fine_factor_arg);
	} else
if(!strcasecmp("plain_rectangular", args_info.sky_grid_arg)){
	patch_grid=make_rect_grid(ceil(2.0*M_PI/(resolution*args_info.fine_factor_arg)), ceil(M_PI/(resolution*args_info.fine_factor_arg)));
	proto_super_grid=make_rect_supergrid(patch_grid, args_info.fine_factor_arg, args_info.fine_factor_arg);
	} else
if(!strcasecmp("arcsin", args_info.sky_grid_arg)){
	patch_grid=make_arcsin_grid(ceil(2.0*M_PI/(resolution*args_info.fine_factor_arg)), ceil(M_PI/(resolution*args_info.fine_factor_arg)));
	proto_super_grid=make_rect_supergrid(patch_grid, args_info.fine_factor_arg, args_info.fine_factor_arg);
	} else {
	fprintf(stderr,"*** Unknown sky grid type: \"%s\"\n", args_info.sky_grid_arg);
	exit(-1);
	}
fine_grid=proto_super_grid->super_grid;

fprintf(stderr,"fine grid: max_n_ra=%d max_n_dec=%d\n", 
	fine_grid->max_n_ra, fine_grid->max_n_dec);

side_cut=args_info.side_cut_arg;
if(!args_info.side_cut_given){
	double max_spindown;
	
	max_spindown=fabs(args_info.spindown_start_arg+(args_info.spindown_count_arg-1)*args_info.spindown_step_arg);
	if(fabs(args_info.spindown_start_arg)>max_spindown)max_spindown=fabs(args_info.spindown_start_arg);
	/* determine side cut from resolution, 6.0 factor is empirical */
	/* also add in spindown contribution - for now just plan for 4 months of data */
	side_cut=260+ceil(M_PI/resolution)/6.0+ceil((1800.0)*max_spindown*args_info.expected_timebase_arg*3600*24*31);
	/* round it up to a multiple of 450 */
	side_cut=450*ceil(side_cut/450.0);
	}
fprintf(stderr,"side_cut=%d\n", side_cut);
first_bin=args_info.first_bin_arg-side_cut;
useful_bins=args_info.nbins_arg;
nbins=args_info.nbins_arg+2*side_cut;

if(args_info.fake_freq_given) {
	fake_injection=1;

   	fprintf(stderr,"fake signal injection: yes\n");
	fprintf(stderr,"fake psi : %f\n", args_info.fake_psi_arg);
	fprintf(stderr,"fake phi : %f\n", args_info.fake_phi_arg);
	fprintf(stderr,"fake iota : %f\n", args_info.fake_iota_arg);
	fprintf(stderr,"fake ra : %f\n", args_info.fake_ra_arg);
	fprintf(stderr,"fake dec: %f\n", args_info.fake_dec_arg);
	fprintf(stderr,"fake spindown: %g\n", args_info.fake_spindown_arg);
	fprintf(stderr,"fake strain: %g\n", args_info.fake_strain_arg);
	fprintf(stderr,"fake frequency: %f\n", args_info.fake_freq_arg);
	fprintf(stderr,"fake reference time: %f\n", args_info.fake_ref_time_arg);
	}

no_am_response=args_info.no_am_response_arg;

init_polarizations0();

/* INPUT stage */

if(args_info.dataset_given) {
	fprintf(stderr, "Loading data from dataset %s\n", args_info.dataset_arg);
	load_dataset_from_file(args_info.dataset_arg);
	nsegments=total_segments();
	}

if(!args_info.spindown_start_time_given) {
	spindown_start=min_gps();
	} else {
	spindown_start=args_info.spindown_start_time_arg;
	}

post_init_datasets();

return Py_BuildValue("i", 0);
}


static double get_double(PyObject *obj, char *field)
{
PyObject *f;
double ans;
f=PyObject_GetAttrString(obj, field);
if(f==NULL) {
	fprintf(stderr, "Could not retrieve field %s from object %p\n", field, obj);
	return(NAN);
	}
ans=PyFloat_AsDouble(f);
Py_DECREF(f);
return(ans);
}

int set_double(PyObject *obj,char *field, double value)
{
return(PyObject_SetAttrString(obj, field, PyFloat_FromDouble(value)));
}

#define RETRIEVE_D(a)	{ cand.a=get_double(ocand, #a); }
#define STORE_D(a)	{ set_double(ocand, #a, cand.a); }

static PyObject * powerflux_compute_scores(PyObject *self, PyObject *args)
{
CANDIDATE cand;
PyObject *ocand;

memset(&cand, 0, sizeof(cand));

ocand=PyTuple_GetItem(args, 0);
if (ocand==NULL)return NULL;

RETRIEVE_D(frequency);
RETRIEVE_D(ra);
RETRIEVE_D(dec);
RETRIEVE_D(spindown);
RETRIEVE_D(psi);
RETRIEVE_D(iota);

compute_scores(&cand, 0);

STORE_D(coherence_score);
STORE_D(chi_sq);
STORE_D(power_cor);
STORE_D(snr);
STORE_D(strain);
STORE_D(strain_err);
STORE_D(total_weight);
STORE_D(f_max);
STORE_D(ifo_freq);
STORE_D(ifo_freq_sd);

return Py_BuildValue("i", 0);
}

static PyObject * powerflux_compute_matched_snr(PyObject *self, PyObject *args)
{
CANDIDATE cand;
PyObject *ocand;

memset(&cand, 0, sizeof(cand));

ocand=PyTuple_GetItem(args, 0);
if (ocand==NULL)return NULL;

RETRIEVE_D(frequency);
RETRIEVE_D(ra);
RETRIEVE_D(dec);
RETRIEVE_D(spindown);
RETRIEVE_D(psi);
RETRIEVE_D(iota);

single_pass_compute_matched_snr(&cand);

STORE_D(coherence_score);
STORE_D(chi_sq);
STORE_D(power_cor);
STORE_D(snr);
STORE_D(strain);
STORE_D(strain_err);
STORE_D(total_weight);
STORE_D(f_max);
STORE_D(ifo_freq);
STORE_D(ifo_freq_sd);

return Py_BuildValue("i", 0);
}

static PyObject * powerflux_compute_simple_snr(PyObject *self, PyObject *args)
{
CANDIDATE cand;
PyObject *ocand;

memset(&cand, 0, sizeof(cand));

ocand=PyTuple_GetItem(args, 0);
if (ocand==NULL)return NULL;

RETRIEVE_D(frequency);
RETRIEVE_D(ra);
RETRIEVE_D(dec);
RETRIEVE_D(spindown);
RETRIEVE_D(psi);
RETRIEVE_D(iota);

single_pass_compute_simple_snr(&cand);

STORE_D(coherence_score);
STORE_D(chi_sq);
STORE_D(power_cor);
STORE_D(snr);
STORE_D(strain);
STORE_D(strain_err);
STORE_D(total_weight);
STORE_D(f_max);
STORE_D(ifo_freq);
STORE_D(ifo_freq_sd);

return Py_BuildValue("i", 0);
}

static PyObject * powerflux_get_datasets(PyObject *self, PyObject *args)
{
PyObject *dataset_list;
PyObject *new_weighted_mean_list;
PyObject *TMedians_list;
PyObject *FMedians_list;
PyObject *gps_list;
PyObject *det_vel_list;
PyObject *veto_list;
long a, min_gps, max_gps;
int i,j;
DATASET *d;

dataset_list=PyList_New(d_free);

for(j=0;j<d_free;j++) {
	d=&(datasets[j]);

	new_weighted_mean_list=PyList_New(d->nbins);
	for(i=0;i<d->nbins;i++)PyList_SetItem(new_weighted_mean_list, i, PyFloat_FromDouble(d->new_weighted_mean[i]));

	FMedians_list=PyList_New(d->nbins);
	for(i=0;i<d->nbins;i++)PyList_SetItem(FMedians_list, i, PyFloat_FromDouble(d->FMedians[i]));

	TMedians_list=PyList_New(d->free);

	for(i=0;i<d->free;i++)PyList_SetItem(TMedians_list, i, PyFloat_FromDouble(d->TMedians[i]));

	gps_list=PyList_New(d->free);
	if(d->free<1) {
		max_gps=-1;
		min_gps=-1;
		} else {
		max_gps=d->gps[0];
		min_gps=d->gps[0];
		}
	for(i=0;i<d->free;i++) {
		a=d->gps[i];
		if(a>max_gps)max_gps=a;
		if(a<min_gps)min_gps=a;
		PyList_SetItem(gps_list, i, PyInt_FromLong(a));
		}

	det_vel_list=PyList_New(d->free);
	for(i=0;i<d->free;i++)PyList_SetItem(det_vel_list, i, Py_BuildValue("(fff)", d->detector_velocity[3*i], d->detector_velocity[3*i+1], d->detector_velocity[3*i+2]));

	veto_list=PyList_New(d->free);
	for(i=0;i<d->free;i++)PyList_SetItem(veto_list, i, PyInt_FromLong(d->sft_veto[i]));

	PyList_SetItem(dataset_list,j,Py_BuildValue("{s:s,s:s,s:i,s:i,s:i,s:f,s:i,s:i,s:i,s:i,s:f,s:f,s:N,s:N,s:N,s:N,s:N,s:N}",
		"name",	d->name,
		"detector", d->detector,
		"sft_count", d->free,
		"nbins", d->nbins,
		"first_bin", d->first_bin,
		"coherence_time", d->coherence_time,
		"gps_start", d->gps_start,
		"gps_stop", d->gps_stop,
		"min_gps", min_gps,
		"max_gps", max_gps,
		"TMedians", d->TMedian,
		"weight", d->weight,
		"new_weighted_means", new_weighted_mean_list,
		"FMedians", FMedians_list,
		"TMedians", TMedians_list,
		"gps", gps_list,
		"detector_velocity", det_vel_list,
		"sft_veto", veto_list));
	}
return dataset_list;
}

static PyObject * powerflux_set_veto(PyObject *self, PyObject *args)
{
long gps_start, gps_stop;
int dataset, flag;
int i, j;
int count;
DATASET *d;

if (!PyArg_ParseTuple(args, "llii", &gps_start, &gps_stop, &dataset, &flag))
        return NULL;

count=0;
for(j=0;j<d_free;j++) {
	if( (dataset>=0) && (dataset!=j))continue;
	d=&(datasets[j]);

	for(i=0;i<d->free;i++)
		if( (d->gps[i]>=gps_start) && (d->gps[i]<=gps_stop)) {
			d->sft_veto[i]=flag;
			count++;
			}
	}

return Py_BuildValue("i", count);
}


static PyObject * powerflux_get_gps(PyObject *self, PyObject *args)
{
PyObject *gps_list;
int nsegments;
int i,j,k;
DATASET *d;

nsegments=0;
for(j=0;j<d_free;j++) {
	d=&(datasets[j]);
	nsegments+=d->free;
	}

gps_list=PyList_New(nsegments);

i=0;
for(j=0;j<d_free;j++) {
	d=&(datasets[j]);
	/* process single SFTs */
	for(k=0;k<d->free;k++) {
		PyList_SetItem(gps_list, i, PyInt_FromLong(d->gps[k]));
		i++;
		}
	}
return gps_list;
}

static PyObject * powerflux_get_power_sum(PyObject *self, PyObject *args)
{
PyObject *power_list;
int nsegments;
int i,j,k;
DATASET *d;
CANDIDATE cand;
PyObject *ocand;
POLARIZATION *pl;
float a, f_plus, f_cross, f_plus_sq, f_cross_sq, doppler;
double f, x, y, power, weight, response;
float mismatch;
int signal_bin;
float filter[7];
float e[26];
float a_plus, a_cross, a_plus_sq, a_cross_sq, c_proj, p_proj;
float *bin_re, *bin_im;
double mult;

memset(&cand, 0, sizeof(cand));

ocand=PyTuple_GetItem(args, 0);
if (ocand==NULL)return NULL;

RETRIEVE_D(frequency);
RETRIEVE_D(ra);
RETRIEVE_D(dec);
RETRIEVE_D(spindown);
RETRIEVE_D(psi);
RETRIEVE_D(iota);

precompute_am_constants(e, cand.ra, cand.dec);
p_proj=cos(2*cand.psi);
c_proj=sin(2*cand.psi);

a=cos(cand.iota);
a_plus=(1.0+a*a)/2.0;
a_cross=a;

/* precompute squares */
a_plus_sq=a_plus*a_plus;
a_cross_sq=a_cross*a_cross;

/* correction factor to convert power into (approx) strain units */
mult=args_info.strain_norm_factor_arg/(1800.0*16384.0);
mult=mult*mult;

nsegments=0;
for(j=0;j<d_free;j++) {
	d=&(datasets[j]);
	nsegments+=d->free;
	}

power_list=PyList_New(nsegments);

i=0;
for(j=0;j<d_free;j++) {
	d=&(datasets[j]);
	pl=&(d->polarizations[0]);
	/* process single SFTs */
	for(k=0;k<d->free;k++) {
		if(d->sft_veto[k]) {
			PyList_SetItem(power_list, i, Py_BuildValue("(fff)", 0.0, 0.0, 0.0));
			i++;
			continue;
			}


		f_plus=F_plus_coeff(k, e, pl->AM_coeffs);
		f_cross=F_plus_coeff(k, e, pl->conjugate->AM_coeffs);

// 		f_plus_sq=f_plus*f_plus;
// 		f_cross_sq=f_cross*f_cross;
// 
// 		a=f_plus*a_plus+f_cross*a_cross;
// 		response[k]=0.25*((f_plus_sq+f_cross_sq)*(a_plus_sq+a_cross_sq)+(f_plus_sq-f_cross_sq)*(a_plus_sq-a_cross_sq)*p_proj+2.0*f_plus*f_cross*(a_plus_sq-a_cross_sq))*c_proj;

		a=f_plus*p_proj+f_cross*c_proj;
		f_cross=f_cross*p_proj-f_plus*c_proj;
		f_plus=a;

		f_plus_sq=f_plus*f_plus;
		f_cross_sq=f_cross*f_cross;

		//response[k]=0.25*((f_plus_sq+f_cross_sq)*(a_plus_sq+a_cross_sq)+(f_plus_sq-f_cross_sq)*(a_plus_sq-a_cross_sq));
		response=f_plus_sq*a_plus_sq+f_cross_sq*a_cross_sq;

		weight=d->expTMedians[k]*d->weight*response;

		doppler=e[0]*d->detector_velocity[3*k+0]+
			e[1]*d->detector_velocity[3*k+1]+
			e[2]*d->detector_velocity[3*k+2];

		f=cand.frequency+cand.frequency*doppler+cand.spindown*(d->gps[k]-spindown_start+d->coherence_time*0.5);

		signal_bin=rintf(1800.0*f-first_bin);
		mismatch=1800.0*f-first_bin-signal_bin;

		if( (signal_bin<0) || (signal_bin>=d->nbins)) {
			PyList_SetItem(power_list, i, Py_BuildValue("(fff)", 0.0, 0.0, 0.0));
			i++;
			continue;
			}

		bin_re=&(d->re[k*nbins+signal_bin]);
		bin_im=&(d->im[k*nbins+signal_bin]);

		fill_hann_filter7(filter, mismatch);

		x=bin_re[-3]*filter[0]+bin_re[-2]*filter[1]+bin_re[-1]*filter[2]+bin_re[0]*filter[3]+bin_re[1]*filter[4]+bin_re[2]*filter[5]+bin_re[3]*filter[6];
		y=bin_im[-3]*filter[0]+bin_im[-2]*filter[1]+bin_im[-1]*filter[2]+bin_im[0]*filter[3]+bin_im[1]*filter[4]+bin_im[2]*filter[5]+bin_im[3]*filter[6];

		power=x*x+y*y;


		PyList_SetItem(power_list, i, Py_BuildValue("(fff)", weight, power*mult, f));
		i++;
		}
	}
return power_list;
}

static PyObject * powerflux_get_power_hist(PyObject *self, PyObject *args)
{
PyObject *power_list;
int i,j,k;
DATASET *d;
CANDIDATE cand;
PyObject *ocand;
POLARIZATION *pl;
float a, f_plus, f_cross, f_plus_sq, f_cross_sq, doppler;
double f, x, y, power, weight, response;
float mismatch;
int signal_bin;
float filter[7];
float e[26];
float a_plus, a_cross, a_plus_sq, a_cross_sq, c_proj, p_proj;
float *bin_re, *bin_im;
double mult;
double *accum_weight;
double *accum_power;
int *count;

memset(&cand, 0, sizeof(cand));

ocand=PyTuple_GetItem(args, 0);
if (ocand==NULL)return NULL;

RETRIEVE_D(frequency);
RETRIEVE_D(ra);
RETRIEVE_D(dec);
RETRIEVE_D(spindown);
RETRIEVE_D(psi);
RETRIEVE_D(iota);

precompute_am_constants(e, cand.ra, cand.dec);
p_proj=cos(2*cand.psi);
c_proj=sin(2*cand.psi);

a=cos(cand.iota);
a_plus=(1.0+a*a)/2.0;
a_cross=a;

/* precompute squares */
a_plus_sq=a_plus*a_plus;
a_cross_sq=a_cross*a_cross;

/* correction factor to convert power into (approx) strain units */
mult=args_info.strain_norm_factor_arg/(1800.0*16384.0);
mult=mult*mult;

accum_weight=do_alloc(nbins, sizeof(*accum_weight));
accum_power=do_alloc(nbins, sizeof(*accum_power));
count=do_alloc(nbins, sizeof(*count));

for(i=0;i<nbins;i++) {
	accum_weight[i]=0.0;
	accum_power[i]=0.0;
	count[i]=0;
	}

i=0;
for(j=0;j<d_free;j++) {
	d=&(datasets[j]);
	pl=&(d->polarizations[0]);
	/* process single SFTs */
	for(k=0;k<d->free;k++) {
		if(d->sft_veto[k]) {
			continue;
			}


		f_plus=F_plus_coeff(k, e, pl->AM_coeffs);
		f_cross=F_plus_coeff(k, e, pl->conjugate->AM_coeffs);

// 		f_plus_sq=f_plus*f_plus;
// 		f_cross_sq=f_cross*f_cross;
// 
// 		a=f_plus*a_plus+f_cross*a_cross;
// 		response[k]=0.25*((f_plus_sq+f_cross_sq)*(a_plus_sq+a_cross_sq)+(f_plus_sq-f_cross_sq)*(a_plus_sq-a_cross_sq)*p_proj+2.0*f_plus*f_cross*(a_plus_sq-a_cross_sq))*c_proj;

		a=f_plus*p_proj+f_cross*c_proj;
		f_cross=f_cross*p_proj-f_plus*c_proj;
		f_plus=a;

		f_plus_sq=f_plus*f_plus;
		f_cross_sq=f_cross*f_cross;

		//response[k]=0.25*((f_plus_sq+f_cross_sq)*(a_plus_sq+a_cross_sq)+(f_plus_sq-f_cross_sq)*(a_plus_sq-a_cross_sq));
		response=f_plus_sq*a_plus_sq+f_cross_sq*a_cross_sq;

		weight=d->expTMedians[k]*d->weight*response;

		doppler=e[0]*d->detector_velocity[3*k+0]+
			e[1]*d->detector_velocity[3*k+1]+
			e[2]*d->detector_velocity[3*k+2];

		f=cand.frequency+cand.frequency*doppler+cand.spindown*(d->gps[k]-spindown_start+d->coherence_time*0.5);

		signal_bin=rintf(1800.0*f-first_bin);
		mismatch=1800.0*f-first_bin-signal_bin;

		if( (signal_bin<0) || (signal_bin>=d->nbins)) {
			continue;
			}

		bin_re=&(d->re[k*nbins+signal_bin]);
		bin_im=&(d->im[k*nbins+signal_bin]);

		fill_hann_filter7(filter, mismatch);

		x=bin_re[-3]*filter[0]+bin_re[-2]*filter[1]+bin_re[-1]*filter[2]+bin_re[0]*filter[3]+bin_re[1]*filter[4]+bin_re[2]*filter[5]+bin_re[3]*filter[6];
		y=bin_im[-3]*filter[0]+bin_im[-2]*filter[1]+bin_im[-1]*filter[2]+bin_im[0]*filter[3]+bin_im[1]*filter[4]+bin_im[2]*filter[5]+bin_im[3]*filter[6];

		power=x*x+y*y;

		accum_weight[signal_bin]+=weight*response;
		accum_power[signal_bin]+=weight*power;
		count[signal_bin]++;
		}
	}

power_list=PyList_New(nbins);

for(i=0;i<nbins;i++) {
	PyList_SetItem(power_list, i, Py_BuildValue("(ffi)", accum_weight[i], accum_power[i], count[i]));
	}

free(accum_weight);
free(accum_power);
free(count);

return power_list;
}

static PyObject * powerflux_get_power_sum_stats(PyObject *self, PyObject *args)
{
int nsegments;
int i,j,k, m;
DATASET *d;
CANDIDATE cand;
PyObject *ocand;
POLARIZATION *pl;
float a, f_plus, f_cross, f_plus_sq, f_cross_sq, doppler;
double f, x, y, power, weight, response;
float mismatch;
int signal_bin;
float filter[7];
float e[26];
float a_plus, a_cross, a_plus_sq, a_cross_sq, c_proj, p_proj;
float *bin_re, *bin_im;
double mult;
float mean;
double sums[5]={0,0,0, 0, 0};

memset(&cand, 0, sizeof(cand));

ocand=PyTuple_GetItem(args, 0);
if (ocand==NULL)return NULL;

RETRIEVE_D(frequency);
RETRIEVE_D(ra);
RETRIEVE_D(dec);
RETRIEVE_D(spindown);
RETRIEVE_D(psi);
RETRIEVE_D(iota);

precompute_am_constants(e, cand.ra, cand.dec);
p_proj=cos(2*cand.psi);
c_proj=sin(2*cand.psi);

a=cos(cand.iota);
a_plus=(1.0+a*a)/2.0;
a_cross=a;

/* precompute squares */
a_plus_sq=a_plus*a_plus;
a_cross_sq=a_cross*a_cross;

/* correction factor to convert power into (approx) strain units */
mult=args_info.strain_norm_factor_arg/(1800.0*16384.0);
mult=mult*mult;

nsegments=0;
for(j=0;j<d_free;j++) {
	d=&(datasets[j]);
	nsegments+=d->free;
	}


i=0;
for(j=0;j<d_free;j++) {
	d=&(datasets[j]);
	pl=&(d->polarizations[0]);
	/* process single SFTs */
	for(k=0;k<d->free;k++) {
		if(d->sft_veto[k]) {
			i++;
			continue;
			}


		f_plus=F_plus_coeff(k, e, pl->AM_coeffs);
		f_cross=F_plus_coeff(k, e, pl->conjugate->AM_coeffs);

// 		f_plus_sq=f_plus*f_plus;
// 		f_cross_sq=f_cross*f_cross;
// 
// 		a=f_plus*a_plus+f_cross*a_cross;
// 		response[k]=0.25*((f_plus_sq+f_cross_sq)*(a_plus_sq+a_cross_sq)+(f_plus_sq-f_cross_sq)*(a_plus_sq-a_cross_sq)*p_proj+2.0*f_plus*f_cross*(a_plus_sq-a_cross_sq))*c_proj;

		a=f_plus*p_proj+f_cross*c_proj;
		f_cross=f_cross*p_proj-f_plus*c_proj;
		f_plus=a;

		f_plus_sq=f_plus*f_plus;
		f_cross_sq=f_cross*f_cross;

		//response[k]=0.25*((f_plus_sq+f_cross_sq)*(a_plus_sq+a_cross_sq)+(f_plus_sq-f_cross_sq)*(a_plus_sq-a_cross_sq));
		response=f_plus_sq*a_plus_sq+f_cross_sq*a_cross_sq;

		weight=d->expTMedians[k]*d->weight*response;

		doppler=e[0]*d->detector_velocity[3*k+0]+
			e[1]*d->detector_velocity[3*k+1]+
			e[2]*d->detector_velocity[3*k+2];

		f=cand.frequency+cand.frequency*doppler+cand.spindown*(d->gps[k]-spindown_start+d->coherence_time*0.5);

		signal_bin=rintf(1800.0*f-first_bin);
		mismatch=1800.0*f-first_bin-signal_bin;

		if( (signal_bin<0) || (signal_bin>=d->nbins)) {
			i++;
			continue;
			}

		bin_re=&(d->re[k*nbins+signal_bin]);
		bin_im=&(d->im[k*nbins+signal_bin]);

		fill_hann_filter7(filter, mismatch);

		x=bin_re[-3]*filter[0]+bin_re[-2]*filter[1]+bin_re[-1]*filter[2]+bin_re[0]*filter[3]+bin_re[1]*filter[4]+bin_re[2]*filter[5]+bin_re[3]*filter[6];
		y=bin_im[-3]*filter[0]+bin_im[-2]*filter[1]+bin_im[-1]*filter[2]+bin_im[0]*filter[3]+bin_im[1]*filter[4]+bin_im[2]*filter[5]+bin_im[3]*filter[6];

		power=x*x+y*y;

		mean=0.0;
		for(m=10;m<30;m++) {
			x=bin_re[m-3]*filter[0]+bin_re[m-2]*filter[1]+bin_re[m-1]*filter[2]+bin_re[m+0]*filter[3]+bin_re[m+1]*filter[4]+bin_re[m+2]*filter[5]+bin_re[m+3]*filter[6];
			y=bin_im[m-3]*filter[0]+bin_im[m-2]*filter[1]+bin_im[m-1]*filter[2]+bin_im[m+0]*filter[3]+bin_im[m+1]*filter[4]+bin_im[m+2]*filter[5]+bin_im[m+3]*filter[6];

			mean+=x*x+y*y;
			
			x=bin_re[-m-3]*filter[0]+bin_re[-m-2]*filter[1]+bin_re[-m-1]*filter[2]+bin_re[-m+0]*filter[3]+bin_re[-m+1]*filter[4]+bin_re[-m+2]*filter[5]+bin_re[-m+3]*filter[6];
			y=bin_im[-m-3]*filter[0]+bin_im[-m-2]*filter[1]+bin_im[-m-1]*filter[2]+bin_im[-m+0]*filter[3]+bin_im[-m+1]*filter[4]+bin_im[-m+2]*filter[5]+bin_im[-m+3]*filter[6];

			mean+=x*x+y*y;
			}
		mean*=0.025;

		/* weight accumulation */
		sums[0]+=weight*response;
		sums[1]+=mult*(power-mean)*weight;


		/* statistical sums */
		sums[2]+=sums[0]*sums[0]*weight*response;
		sums[3]+=sums[0]*sums[1]*weight*response;
		sums[4]+=sums[1]*sums[1]*weight*response;

		i++;
		}
	}

sums[2]/=sums[0];
sums[3]/=sums[0];
sums[4]/=sums[0];

mult=- sums[3]/sums[2];

//fprintf(stderr, "%g %g %g %g %g mult=%g\n", sums[0], sums[1], sums[2], sums[3], sums[4], mult);

cand.total_weight=sums[0];
cand.power_cor=mult;
x=sqrt(sums[4]+2*sums[3]*mult+sums[2]*mult*mult);
cand.strain_err=x;
cand.snr=sums[1]/x;
if(sums[1]>0)cand.strain=sqrt(sums[1]/sums[0]);

STORE_D(coherence_score);
STORE_D(chi_sq);
STORE_D(power_cor);
STORE_D(snr);
STORE_D(strain);
STORE_D(strain_err);
STORE_D(total_weight);
STORE_D(f_max);
STORE_D(ifo_freq);
STORE_D(ifo_freq_sd);

return Py_BuildValue("i", 0);
}

static PyMethodDef PowerFluxMethods[] = {
    	{"init",  powerflux_init, METH_VARARGS, "Initialize powerflux"},
    	{"compute_scores",  powerflux_compute_scores, METH_VARARGS, "Compute candidate scores"},
    	{"compute_matched_snr",  powerflux_compute_matched_snr, METH_VARARGS, "Compute matched snr"},
    	{"compute_single_snr",  powerflux_compute_simple_snr, METH_VARARGS, "Compute simple snr"},
    	{"get_datasets",  powerflux_get_datasets, METH_VARARGS, "Get information about datasets loaded"},
    	{"set_veto",  powerflux_set_veto, METH_VARARGS, "Set veto flag for range of GPS values"},
    	{"get_gps",  powerflux_get_gps, METH_VARARGS, "Get list of start times of SFTs"},
    	{"get_power_sum",  powerflux_get_power_sum, METH_VARARGS, "Get power sum"},
    	{"get_power_hist",  powerflux_get_power_hist, METH_VARARGS, "Get power distribution across frequency bins at the detector"},
    	{"get_power_sum_stats",  powerflux_get_power_sum_stats, METH_VARARGS, "Get power sum with stats"},
    	{NULL, NULL, 0, NULL}        /* Sentinel */
	};

void initpowerflux(void)
{
FILE_LOG=fopen("/dev/null", "w");
LOG=fopen("/dev/null", "w");
if(FILE_LOG==NULL) {
	FILE_LOG=stdout;
	}
if(LOG==NULL) {
	LOG=stdout;
	}
(void) Py_InitModule("powerflux", PowerFluxMethods);
}

