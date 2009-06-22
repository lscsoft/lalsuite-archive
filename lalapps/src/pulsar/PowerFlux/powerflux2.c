/*
*  Copyright (C) 2007 Vladimir Dergachev
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
#include "lines.h"
#include "grid.h"
#include "polarization.h"
#include "statistics.h"
#include "dataset.h"
#include "candidates.h"
#include "power_sum_stats.h"
#include "outer_loop.h"
#include "util.h"
#include "jobs.h"
#include "skymarks.h"

extern int ntotal_polarizations, nlinear_polarizations;
extern POLARIZATION *polarizations;

FILE *LOG=NULL, *FILE_LOG=NULL, *DATA_LOG=NULL;
time_t start_time, end_time, stage_time;

extern INT64 spindown_start;

extern DATASET *datasets;
extern int d_free;

int nsegments;
int useful_bins, nbins,side_cut;

/* TODO: this variable is not used in 2.0 */
int stored_fine_bins=0;

double average_det_velocity[3];
float det_vel_ra, det_vel_dec;
double orbital_axis[3];
double band_axis_norm;
double band_axis[3];
float band_axis_ra, band_axis_dec;
double large_S=-1;

int first_bin;

char *sky_marks=NULL;
int sky_marks_free=0;
int sky_marks_size=0;
SKYMARK *compiled_skymarks=NULL;

int subinstance;
char *subinstance_name;

int do_CutOff=1;

struct gengetopt_args_info args_info;

double spindown;

char *earth_ephemeris=NULL, *sun_ephemeris=NULL;
double resolution; /* this is actual resolution, not the resolution argument passed on command line */
int fake_injection=0;

int no_am_response;

SKY_GRID *fine_grid=NULL;
SKY_SUPERGRID *super_grid=NULL, *proto_super_grid=NULL;
SKY_GRID *patch_grid=NULL;

char *output_dir;


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

static float compute_median(float *firstbin, int step, int count)
{
float *tmp;
int i;
tmp=aligned_alloca(count*sizeof(float));
for(i=0;i<count;i++)tmp[i]=firstbin[i*step];
sort_floats(tmp, count);
if(!(count & 1))return (tmp[count>>1]+tmp[(count>>1)-1])/2.0;
return tmp[count>>1];
}


static double exponential_distribution(double x, double lambda)
{
if(x<=0.0)return 0.0;
return (1-exp(-lambda*x));
}

void nonparametric(float *firstbin, int step, int count, float *median, float *ks_test)
{
float *tmp;
double a,b,lambda;
int i;
tmp=aligned_alloca(count*sizeof(float));
for(i=0;i<count;i++)tmp[i]=firstbin[i*step];
sort_floats(tmp, count);
if(count & 1)*median=(tmp[count>>1]+tmp[(count>>1)+1])/2.0;
	else *median=tmp[count>>1];
if(*median==0)return;
lambda=M_LN2/(*median);
a=-2.0;
for(i=0;i<count;i++){
	b=(((i+1)*1.0)/count)-exponential_distribution(tmp[i], lambda);
	if(a<b)a=b;
	b=exponential_distribution(tmp[i], lambda)-((i*1.0)/count);
	if(a<b)a=b;
	}
*ks_test=a;
}

void wrap_up(void)
{
time(&end_time);
fprintf(stderr, "exit memory: %g MB\n", (MEMUSAGE*10.0/(1024.0*1024.0))/10.0);
fprintf(LOG,"seconds elapsed: %ld\n",end_time-start_time);
fprintf(stderr,"seconds elapsed: %ld\n",end_time-start_time);
fclose(LOG);
fclose(FILE_LOG);
fclose(DATA_LOG);
}

int main(int argc, char *argv[])
{
RGBPic *p;
PLOT *plot;
char s[20000];
int i, j;
struct rlimit rl;

/* INIT stage */

time(&start_time);

fedisableexcept(FE_ALL_EXCEPT);

fprintf(stderr, "Initial memory: %g MB\n", (MEMUSAGE*10.0/(1024.0*1024.0))/10.0);

if(getrlimit(RLIMIT_AS, &rl)<0) {
	perror("Could not obtain virtual memory limit:");
	} else {
	fprintf(stderr, "Virtual memory limits soft=%ld hard=%ld\n", rl.rlim_cur, rl.rlim_max);
	}

if(getrlimit(RLIMIT_DATA, &rl)<0) {
	perror("Could not obtain data segment limit:");
	} else {
	fprintf(stderr, "Data segment limits soft=%ld hard=%ld\n", rl.rlim_cur, rl.rlim_max);
	}

if(cmdline_parser(argc, argv, &args_info))exit(-1);

for(i=0;i<args_info.config_given;i++)
	if(cmdline_parser_configfile(args_info.config_arg[i], &args_info, 0, 0, 0))exit(-1);

if(!args_info.dataset_given){
	fprintf(stderr,"** You must specify dataset description file (--dataset)\n");
	exit(-1);
	}
if(!args_info.ephemeris_path_given &&
	!(args_info.earth_ephemeris_given && args_info.sun_ephemeris_given)){
	fprintf(stderr,"** You must specify path to ephemeris files (--detreponse-path or --XXX-ephemeris)\n");
	exit(-1);
	}
if(!args_info.first_bin_given){
	fprintf(stderr,"** You must specify first bin to analyze (--first-bin)\n");
	exit(-1);
	}

if(!args_info.dataset_given && !args_info.detector_given){
	fprintf(stderr,"** You must specify detector (--detector)\n");
	exit(-1);
	}

gsl_rng_env_setup();
gsl_set_error_handler_off();

/* create output directories if not present */
if(args_info.output_given){
	mkdir(args_info.output_arg, 0777);
	output_dir=do_alloc(strlen(args_info.output_arg)+30, sizeof(*output_dir));
	sprintf(output_dir, "%s/%d-%f/", args_info.output_arg, args_info.first_bin_arg,args_info.first_bin_arg/1800.0);
	} else {
	output_dir=do_alloc(30, sizeof(*output_dir));
	sprintf(output_dir, "%d-%f/", args_info.first_bin_arg,args_info.first_bin_arg/1800.0);
	}
mkdir(output_dir, 0777);

if(args_info.dump_points_arg){
	snprintf(s, 20000, "%s/points", output_dir);
	mkdir(s, 0777);
	}


snprintf(s,20000,"%s/powerflux.log", output_dir);
LOG=fopen(s,"w");

while(LOG==NULL) {
	i=errno;
	fprintf(stderr, "Could not open log file \"%s\": %s\n", s, strerror(errno));
	condor_safe_sleep(60);
	LOG=fopen(s, "w");
	}

snprintf(s,20000,"%s/file.log", output_dir);
FILE_LOG=fopen(s,"w");

snprintf(s,20000,"%s/data.log", output_dir);
DATA_LOG=fopen(s,"w");
setbuffer(DATA_LOG, do_alloc(1024*1024*32, 1), 1024*1024*32);

if(args_info.label_given){
	fprintf(LOG, "label: \"%s\"\n", args_info.label_arg);
	fprintf(stderr, "label: \"%s\"\n", args_info.label_arg);
	}
	
if(gethostname(s, 19999)>=0){
	fprintf(stderr, "Running on %s\n", s);
	fprintf(LOG, "node: %s\n", s);
	} else {
	fprintf(stderr, "Could not obtain hostname\n");
	fprintf(LOG, "node: unknown\n");
	}

for(i=0;i<args_info.config_given;i++) {
	fprintf(LOG, "using config file: %s\n", args_info.config_arg[i]);
	}

init_threads(args_info.num_threads_arg);
init_jobs();
init_hookup();
init_statistics();
tabulate_hann_filter7();
init_power_sum_stats();
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

/* PREP1 stage */

fprintf(stderr,	"Initializing sky grids\n");

	/* Compute resolution from frequency.
	
	   The resolution needs to be fine enough so that 
	   the doppler shifts stay within 1/1800 Hz within 
	   single sky bin on a fine grid.
	  
	   This is usually much finer than the resolution needed
	   to resolve amplitude modulation.
	   
	   The variation of the Doppler shift depending on the sky
	   position is determined by the angle between average velocity
	   vector and the sky position.
	   
	   Thus VARIATION=CONST * AVERAGE_VELOCITY_VECTOR * FREQUENCY * RESOLUTION
	   
	   Since we need the condition to hold during the entire run, we use
	   
	   RESOLUTION = VARIATION/(CONST * MAX_AVERAGE_VELOCITY_VECTOR * FREQUENCY)
	   
	   We can now take both the VARIATION and FREQUENCY to be measured in 
	   frequency bins. 
	   
	   The value VARIATION/(CONST * MAX_AVERAGE_VELOCITY_VECTOR) can be
	   determined empirically from a few bands of PowerFlux.
	 */
	 
resolution=(4500.0)/(args_info.first_bin_arg+args_info.nbins_arg/2);
/* AM response is computed on a patch grid - this can be too coarse at low frequencies */
//if(resolution*args_info.fine_factor_arg> 0.025) resolution=0.025/args_info.fine_factor_arg;

fprintf(LOG,"resolution (auto) : %f\n", resolution);


resolution*=args_info.skymap_resolution_ratio_arg;
if(args_info.skymap_resolution_given){
	resolution=args_info.skymap_resolution_arg;
	}

fprintf(LOG,"resolution : %f\n", resolution);
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

fprintf(LOG, "sky map orientation: %s\n", args_info.skymap_orientation_arg);

if(!strcasecmp("ecliptic", args_info.skymap_orientation_arg)){
	rotate_grid_xy(patch_grid, -M_PI/2.0);
	rotate_grid_xy(fine_grid, -M_PI/2.0);

	rotate_grid_xz(patch_grid, -M_PI*23.44/180.0);
	rotate_grid_xz(fine_grid, -M_PI*23.44/180.0);

	rotate_grid_xy(patch_grid, M_PI/2.0);
	rotate_grid_xy(fine_grid, M_PI/2.0);
	} else
if(!strcasecmp("band_axis", args_info.skymap_orientation_arg)){
	rotate_grid_xy(patch_grid, -band_axis_ra);
	rotate_grid_xy(fine_grid, -band_axis_ra);

	rotate_grid_xz(patch_grid, -band_axis_dec+M_PI/2.0);
	rotate_grid_xz(fine_grid, -band_axis_dec+M_PI/2.0);

	rotate_grid_xy(patch_grid, band_axis_ra);
	rotate_grid_xy(fine_grid, band_axis_ra);
	}

fprintf(stderr, "Full grid memory: %g MB\n", (MEMUSAGE*10.0/(1024.0*1024.0))/10.0);
fprintf(LOG, "Full grid memory: %g MB\n", (MEMUSAGE*10.0/(1024.0*1024.0))/10.0);

free_values(fine_grid);

// fprintf(stderr, "No values grid memory: %g MB\n", (MEMUSAGE*10.0/(1024.0*1024.0))/10.0);
// fprintf(LOG, "No values memory: %g MB\n", (MEMUSAGE*10.0/(1024.0*1024.0))/10.0);

	/* Determine side_cut - amount of extra bins to load to accomodate 
	   Doppler shifts and spindown.

	   Since sky positions directly ahead and directly benhind average velocity
	   vector are the extrema of Doppler shifts the side_cut can be chosen 
	   large enough to be independent of the duration of the run, but without
	   loading an excessively large amount of data.
	   
	   Since the Doppler shifts do not change by more than 1/1800 Hz when moving
	   the distance of resolution radians on the sky the total variation is less than
	   M_PI/resolution. We use this as an upper bounds which turns out to be reasonably 
	   precise.
	   
	   Earth speed is approx v=30km/sec c=300000 km/sec
	   v/c=1e-4
	   sidecut = 1e-4 * freq_in_bins
	   
	   Compare to formula used:
	   
	   side_cut = M_PI/(6.0*resolution) = M_PI*freq_in_bins/(6*4500) = 1.1e-4 * freq_in_bins
	*/

side_cut=args_info.side_cut_arg;
if(!args_info.side_cut_given){
	double max_spindown;
	
	max_spindown=fabs(args_info.spindown_start_arg+(args_info.spindown_count_arg-1)*args_info.spindown_step_arg);
	if(fabs(args_info.spindown_start_arg)>max_spindown)max_spindown=fabs(args_info.spindown_start_arg);
	/* determine side cut from resolution, 6.0 factor is empirical */
	/* also add in spindown contribution - for now just plan for 4 months of data */
	side_cut=260+ceil(M_PI/resolution)/6.0+ceil((1800.0)*max_spindown*args_info.expected_timebase_arg*3600*24*31);
	/* round it up to a multiple of 450 */
/*	side_cut=450*ceil(side_cut/450.0);*/
	}
fprintf(stderr,"side_cut=%d\n", side_cut);
first_bin=args_info.first_bin_arg-side_cut;
useful_bins=args_info.nbins_arg;
nbins=args_info.nbins_arg+2*side_cut;

if(fine_grid->max_n_dec<800) {
	p=make_RGBPic(fine_grid->max_n_ra*(800/fine_grid->max_n_dec)+140, fine_grid->max_n_dec*(800/fine_grid->max_n_dec));
	} else
	p=make_RGBPic(fine_grid->max_n_ra+140, fine_grid->max_n_dec);

plot=make_plot(p->width, p->height);

#if 0  /* Debugging.. - do not remove unless structures and algorithm change */
{
	float *tmp;
	tmp=do_alloc(fine_grid->npoints, sizeof(*tmp));

	for(i=0;i<fine_grid->npoints;i++)tmp[i]=super_grid->list_map[i];
	plot_grid_f(p, fine_grid, tmp,1);
	RGBPic_dump_png("list_map.png", p);
	dump_ints("list_map.dat", super_grid->list_map, super_grid->super_grid->npoints, 1);

	for(i=0;i<fine_grid->npoints;i++)tmp[i]=super_grid->reverse_map[i];
	plot_grid_f(p, fine_grid, tmp,1);
	RGBPic_dump_png("reverse_map.png", p);
	dump_ints("reverse_map.dat", super_grid->reverse_map, super_grid->super_grid->npoints, 1);

	for(i=0;i<patch_grid->npoints;i++)tmp[i]=super_grid->first_map[i];
	plot_grid_f(p, patch_grid, tmp,1);
	RGBPic_dump_png("first_map.png", p);
	dump_ints("first_map.dat", super_grid->first_map, patch_grid->npoints, 1);
}


dump_floats("e0.dat",patch_grid->e[0],patch_grid->npoints,1);
dump_floats("e1.dat",patch_grid->e[1],patch_grid->npoints,1);
dump_floats("e2.dat",patch_grid->e[2],patch_grid->npoints,1);
dump_floats("e3.dat",patch_grid->e[3],patch_grid->npoints,1);
dump_floats("e4.dat",patch_grid->e[4],patch_grid->npoints,1);
dump_floats("e5.dat",patch_grid->e[5],patch_grid->npoints,1);
#endif

no_am_response=args_info.no_am_response_arg;


fprintf(LOG,"powerflux : %s\n",VERSION);
if(no_am_response){
	fprintf(LOG,"no_am_response : true\n");
	fprintf(stderr,"NO_AM_RESPONSE flag passed\n");
	}
fprintf(LOG,"firstbin  : %d\n",first_bin);
fprintf(LOG,"band start: %g Hz\n",first_bin/1800.0);
fprintf(LOG,"nbins     : %d\n",nbins);
fprintf(LOG,"side_cut  : %d\n",side_cut);
fprintf(LOG,"useful bins : %d\n",useful_bins);
fprintf(LOG,"useful band start: %g Hz\n",(first_bin+side_cut)/1800.0);
fprintf(LOG,"averaging mode: %s\n", args_info.averaging_mode_arg);

fprintf(LOG,"patch_type: %s\n", patch_grid->name);
fprintf(LOG,"patch_grid: %dx%d\n", patch_grid->max_n_ra, patch_grid->max_n_dec);
fprintf(LOG,"patch_grid npoints : %d\n", patch_grid->npoints);
fprintf(LOG,"fine_factor: %d\n", args_info.fine_factor_arg);
fprintf(LOG,"fine_type : %s\n",fine_grid->name);
fprintf(LOG,"fine_grid npoints  : %d\n", fine_grid->npoints);
fprintf(LOG,"fine_grid : %dx%d\n", fine_grid->max_n_ra, fine_grid->max_n_dec);
if(args_info.dataset_given)fprintf(LOG,"dataset: %s\n",args_info.dataset_arg);
fprintf(LOG,"first spindown: %g\n", args_info.spindown_start_arg);
fprintf(LOG,"spindown step : %g\n", args_info.spindown_step_arg);
fprintf(LOG,"spindown count: %d\n", args_info.spindown_count_arg);
fprintf(LOG,"niota: %d\n", args_info.niota_arg);
fprintf(LOG,"npsi: %d\n", args_info.npsi_arg);
fprintf(LOG,"nfshift: %d\n", args_info.nfshift_arg);
fprintf(LOG,"summing step: %d\n", args_info.summing_step_arg);
fprintf(LOG,"max_first_shift: %d\n", args_info.max_first_shift_arg);
fprintf(LOG,"Doppler multiplier: %g\n", args_info.doppler_multiplier_arg);
fprintf(LOG,"orientation: %g\n", args_info.orientation_arg);
fprintf(LOG,"make cutoff: %s\n",do_CutOff ? "yes (unused)" : "no" );
fprintf(LOG, "weight cutoff fraction: %g\n", args_info.weight_cutoff_fraction_arg);
fprintf(LOG, "per dataset weight cutoff fraction: %g\n", args_info.per_dataset_weight_cutoff_fraction_arg);
fprintf(LOG, "noise level: %s\n", args_info.tmedian_noise_level_arg ? "TMedian" : "in_place_sd");
fprintf(LOG, "skymarks: %s\n", args_info.fine_grid_skymarks_arg ? "spindown_independent" : "spindown_dependent");

fprintf(LOG, "subtract background: %s\n", args_info.subtract_background_arg ? "yes" : "no");
fflush(LOG);

/* we do not use precomputed polarization arrays in powerflux2 
 init_polarizations0(); */
ntotal_polarizations=0;

/* do we need to inject fake signal ? */
if(args_info.fake_freq_given) {
	fake_injection=1;

   	fprintf(LOG,"fake signal injection: yes\n");
	fprintf(LOG,"fake psi : %f\n", args_info.fake_psi_arg);
	fprintf(LOG,"fake phi : %f\n", args_info.fake_phi_arg);
	fprintf(LOG,"fake iota : %f\n", args_info.fake_iota_arg);
	fprintf(LOG,"fake ra : %f\n", args_info.fake_ra_arg);
	fprintf(LOG,"fake dec: %f\n", args_info.fake_dec_arg);
	fprintf(LOG,"fake spindown: %g\n", args_info.fake_spindown_arg);
	fprintf(LOG,"fake strain: %g\n", args_info.fake_strain_arg);
	fprintf(LOG,"fake frequency: %f\n", args_info.fake_freq_arg);
	fprintf(LOG,"fake reference time: %f\n", args_info.fake_ref_time_arg);
	
   	} else {
   	fprintf(LOG,"fake signal injection: none\n");
	}

init_datasets();
test_datasets();

/* INPUT stage */

if(args_info.dataset_given) {
	fprintf(stderr, "Loading data from dataset %s\n", args_info.dataset_arg);
	load_dataset_from_file(args_info.dataset_arg);
	nsegments=total_segments();
	if(nsegments==vetoed_segments()) {
		fprintf(LOG, "All SFTs vetoed, aborting!\n");
		fprintf(stderr, "All SFTs vetoed, aborting!\n");
		exit(-1);
		}
	}

/* This diagnostics should be moved into dataset.c when tested */
for(i=0;i<d_free;i++) {
	fprintf(LOG, "FMedians: \"%s\" \"%s\"", args_info.label_arg, datasets[i].name);
	for(j=0;j<datasets[i].nbins;j++)fprintf(LOG, " %g", datasets[i].TMedians[j]);
	fprintf(LOG, "\n");
	}

for(i=0;i<d_free;i++) {
	fprintf(LOG, "new_weighted_mean: \"%s\" \"%s\"", args_info.label_arg, datasets[i].name);
	for(j=0;j<datasets[i].nbins;j++)fprintf(LOG, " %g", datasets[i].new_weighted_mean[j]);
	fprintf(LOG, "\n");
	}

if(args_info.sky_marks_file_given) {
	FILE *f;
	int a;
	f=fopen(args_info.sky_marks_file_arg, "r");
	if(f==NULL) {
		perror("Could not open skymarks file:");
		exit(-1);
		}
	sky_marks_size=10000;
	sky_marks=do_alloc(sky_marks_size, sizeof(char));
	sky_marks_free=0;

	while(1) {
		a=fread(&(sky_marks[sky_marks_free]), 1, 10000, f);
		if(a>0)sky_marks_free+=a;
		if(a<10000)break;

		if(sky_marks_free+10000>=sky_marks_size) {
			char *ptr;
			sky_marks_size*=2;
			ptr=do_alloc(sky_marks_size, sizeof(char));
			memcpy(ptr, sky_marks, sky_marks_free);
			free(sky_marks);
			sky_marks=ptr;
			}
		}
	fclose(f);
	fprintf(LOG, "sky marks file: \"%s\"\n", args_info.sky_marks_file_arg);
	}
		

time(&stage_time);
fprintf(LOG, "input complete: %d\n", (int)(stage_time-start_time));
fprintf(stderr, "input complete: %d\n", (int)(stage_time-start_time));

if(args_info.dump_data_given) {
	dump_datasets(args_info.dump_data_arg);
	}


if(nsegments==0){
	fprintf(stderr,"ERROR: no input data found !\n");
	return -1;
	}
post_init_datasets();

if(args_info.dump_sftv2_given) {
	sftv2_dump_datasets(args_info.dump_sftv2_arg);
	}

fprintf(LOG,"nsegments : %d\n", nsegments);

fprintf(LOG,"first gps : %lld\n", min_gps());
fprintf(LOG,"last gps  : %lld\n", max_gps());

if(!args_info.spindown_start_time_given){
	spindown_start=min_gps();
	} else {
	spindown_start=args_info.spindown_start_time_arg;
	}
fprintf(LOG, "spindown start time: %lld\n", spindown_start);

fflush(LOG);


/* DIAG2 stage */
fprintf(stderr,"Computing detector speed\n");
datasets_average_detector_speed(average_det_velocity);
fprintf(LOG,"average detector velocity: %g %g %g\n", 
	average_det_velocity[0],
	average_det_velocity[1],
	average_det_velocity[2]);
det_vel_dec=atan2f(average_det_velocity[2], 
	sqrt(average_det_velocity[0]*average_det_velocity[0]+average_det_velocity[1]*average_det_velocity[1]));
det_vel_ra=atan2f(average_det_velocity[1], average_det_velocity[0]);
if(det_vel_ra<0)det_vel_ra+=2.0*M_PI;
fprintf(LOG,"average detector velocity RA (degrees) : %f\n", det_vel_ra*180.0/M_PI);
fprintf(LOG,"average detector velocity DEC (degrees): %f\n", det_vel_dec*180.0/M_PI);
fprintf(stderr,"average detector velocity RA (degrees) : %f\n", det_vel_ra*180.0/M_PI);
fprintf(stderr,"average detector velocity DEC (degrees): %f\n", det_vel_dec*180.0/M_PI);

orbital_axis[0]=0.0;
orbital_axis[1]=-sin(M_PI*23.44/180.0);
orbital_axis[2]=cos(M_PI*23.44/180.0);

/* crossproduct gives the vector perpedicular to both the average doppler shift and
  orbital axis */
band_axis[0]=orbital_axis[1]*average_det_velocity[2]-orbital_axis[2]*average_det_velocity[1];
band_axis[1]=orbital_axis[2]*average_det_velocity[0]-orbital_axis[0]*average_det_velocity[2];
band_axis[2]=orbital_axis[0]*average_det_velocity[1]-orbital_axis[1]*average_det_velocity[0];

/* Normalize */
band_axis_norm=sqrt(band_axis[0]*band_axis[0]+
	band_axis[1]*band_axis[1]+
	band_axis[2]*band_axis[2]);
/* replace 0.0 with something more reasonable later */
if(band_axis_norm<=0.0){
	band_axis[0]=0.0;
	band_axis[1]=0.0;
	band_axis[2]=1.0;
	} else {
	band_axis[0]/=band_axis_norm;
	band_axis[1]/=band_axis_norm;
	band_axis[2]/=band_axis_norm;
	}
/* 
  Normalize band_axis_norm so it matches the definition of \vec{u}
*/
band_axis_norm*=2.0*M_PI/(365.0*24.0*3600.0);

fprintf(LOG, "auto band axis norm: %g\n", band_axis_norm);
fprintf(LOG, "maximum S contribution from Doppler shifts: %g\n", band_axis_norm*(first_bin+nbins*0.5)/1800.0);

if(args_info.band_axis_norm_given) {
	band_axis_norm=args_info.band_axis_norm_arg;
	}

fprintf(LOG, "actual band axis norm: %g\n", band_axis_norm);

large_S=6.0/(1800.0*(max_gps()-min_gps()+1800.0));

fprintf(LOG, "auto large S: %g\n", large_S);

if(args_info.large_S_given){
	large_S=args_info.large_S_arg;
	}
fprintf(LOG, "large S: %g\n", large_S);
	

fprintf(LOG,"auto band axis: %g %g %g\n", 
	band_axis[0],
	band_axis[1],
	band_axis[2]);
fprintf(stderr,"auto band axis: %g %g %g\n", 
	band_axis[0],
	band_axis[1],
	band_axis[2]);

fprintf(LOG, "band_axis: %s\n", args_info.band_axis_arg);
fprintf(stderr, "band_axis: %s\n", args_info.band_axis_arg);
if(!strcasecmp(args_info.band_axis_arg, "equatorial")){
	band_axis[0]=0.0;
	band_axis[1]=0.0;
	band_axis[2]=1.0;	
	} else
if(!strncasecmp(args_info.band_axis_arg, "explicit", 8)){
	int q;
	q=sscanf(args_info.band_axis_arg+8, "(%lf,%lf,%lf)", 
		&(band_axis[0]),
		&(band_axis[1]),
		&(band_axis[2]));
	if(q!=3){
		fprintf(stderr,"Warning ! not all explicit band axis values were assigned. Format error ?\n");
		fprintf(LOG,"Warning ! not all explicit band axis values were assigned. Format error ?\n");
		}
	}
	
fprintf(LOG,"actual band axis: %g %g %g\n", 
	band_axis[0],
	band_axis[1],
	band_axis[2]);
fprintf(stderr,"actual band axis: %g %g %g\n", 
	band_axis[0],
	band_axis[1],
	band_axis[2]);

band_axis_dec=atan2f(band_axis[2], 
	sqrt(band_axis[0]*band_axis[0]+band_axis[1]*band_axis[1]));
band_axis_ra=atan2f(band_axis[1], band_axis[0]);

if(band_axis_ra<0)band_axis_ra+=2.0*M_PI;
fprintf(stderr,"band axis RA (degrees) : %f\n", band_axis_ra*180.0/M_PI);
fprintf(stderr,"band axis DEC (degrees): %f\n", band_axis_dec*180.0/M_PI);
fprintf(LOG,"band axis RA (degrees) : %f\n", band_axis_ra*180.0/M_PI);
fprintf(LOG,"band axis DEC (degrees): %f\n", band_axis_dec*180.0/M_PI);

/* now that we have new grid positions plot them */

plot_grid_f(p, patch_grid, patch_grid->latitude,1);
RGBPic_dump_png("patch_latitude.png", p);
dump_floats("patch_latitude.dat", patch_grid->latitude, patch_grid->npoints, 1);

plot_grid_f(p, patch_grid, patch_grid->longitude,1);
RGBPic_dump_png("patch_longitude.png", p);
dump_floats("patch_longitude.dat", patch_grid->longitude, patch_grid->npoints, 1);

plot_grid_f(p, fine_grid, fine_grid->latitude,1);
RGBPic_dump_png("fine_latitude.png", p);
dump_floats("fine_latitude.dat", fine_grid->latitude, fine_grid->npoints, 1);

plot_grid_f(p, fine_grid, fine_grid->longitude,1);
RGBPic_dump_png("fine_longitude.png", p);
dump_floats("fine_longitude.dat", fine_grid->longitude, fine_grid->npoints, 1);


/* COMP3 stage */

if(args_info.no_decomposition_arg){
	fprintf(stderr,"Exiting as requested (--no-decomposition=1\n");
	fprintf(LOG,"Exiting as requested (--no-decomposition=1\n");
	wrap_up();
	exit(0);
	}

/* DIAG4 stage */


if(args_info.no_demodulation_arg){
	fprintf(stderr,"Exiting as requested (--no-demodulation=1\n");
	fprintf(LOG,"Exiting as requested (--no-demodulation=1\n");
	wrap_up();
	exit(0);	
	}

/* PREP5 stage */

output_datasets_info();

/* TODO: this is here so we can fallback on old skymark code (for comparison). It would be nice to remove this altogether */
spindown=args_info.spindown_start_arg+0.5*(args_info.spindown_count_arg-1)*args_info.spindown_step_arg;

if(sky_marks!=NULL) {
	process_marks(fine_grid, sky_marks, sky_marks_free);
	compiled_skymarks=compile_marks(fine_grid, sky_marks, sky_marks_free);

	process_marks(patch_grid, "disk \"sky\" \"\" 0 0 100", 21);
	propagate_far_points_from_super_grid(patch_grid, proto_super_grid);
	//process_marks(patch_grid, sky_marks, sky_marks_free);
	} else {
	/* assign bands */
	if(!strcasecmp("S", args_info.skyband_method_arg)) {
		S_assign_bands(patch_grid, args_info.nskybands_arg, large_S, spindown, (first_bin+nbins*0.5)/1800.0);
		S_assign_bands(fine_grid, args_info.nskybands_arg, large_S, spindown, (first_bin+nbins*0.5)/1800.0);
		} else 
	if(!strcasecmp("angle", args_info.skyband_method_arg)) {
		angle_assign_bands(patch_grid, args_info.nskybands_arg);
		angle_assign_bands(fine_grid, args_info.nskybands_arg);
		} else {
		fprintf(stderr, "*** ERROR: unknown band assigment method \"%s\"\n",
			args_info.skyband_method_arg);
		fprintf(LOG, "*** ERROR: unknown band assigment method \"%s\"\n",
			args_info.skyband_method_arg);
		exit(-1);
		}
	}
/* mask points if requested */

if(args_info.focus_ra_given && 
	args_info.focus_dec_given && 
	args_info.focus_radius_given) {
	fprintf(LOG, "focus ra    : %f\n", args_info.focus_ra_arg);
	fprintf(LOG, "focus dec   : %f\n", args_info.focus_dec_arg);
	fprintf(LOG, "focus radius: %f\n", args_info.focus_radius_arg);
	mask_far_points(fine_grid, args_info.focus_ra_arg, args_info.focus_dec_arg, args_info.focus_radius_arg);
	propagate_far_points_from_super_grid(patch_grid, proto_super_grid);
	}

if(args_info.only_large_cos_given) {
	fprintf(LOG, "only large cos level: %f\n", args_info.only_large_cos_arg);
	mask_small_cos(fine_grid, band_axis[0], band_axis[1], band_axis[3], args_info.only_large_cos_arg);
	propagate_far_points_from_super_grid(patch_grid, proto_super_grid);
	}

super_grid=reduced_supergrid(proto_super_grid);
fine_grid=super_grid->super_grid;

print_grid_statistics(LOG, "", fine_grid);
precompute_values(fine_grid);
precompute_values(patch_grid);
verify_dataset_whole_sky_AM_response();

plot_grid_f(p, patch_grid, patch_grid->band_f, 1);
snprintf(s, 20000, "bands.png");
RGBPic_dump_png(s, p);
snprintf(s, 20000, "bands.dat");
dump_ints(s, patch_grid->band, patch_grid->npoints, 1);
fflush(LOG);


power_cache_selftest();
power_sum_stats_selftest();

/* MAIN LOOP stage */
time(&stage_time);
fprintf(LOG, "outer_loop_start: %d\n", (int)(stage_time-start_time));
fprintf(stderr, "outer_loop_start: %d\n", (int)(stage_time-start_time));

outer_loop();

fflush(LOG);

/*	fine_grid_free_arrays();*/
fine_grid=proto_super_grid->super_grid;
free_supergrid(super_grid);
super_grid=NULL;

fflush(LOG);

wrap_up();
return 0;
}
