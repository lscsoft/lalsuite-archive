#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <time.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "global.h"
#include "hookup.h"
#include "rastermagic.h"
#include "cmdline.h"
#include "fine_grid.h"
#include "lines.h"
#include "grid.h"
#include "polarization.h"

extern int npolarizations;
extern POLARIZATION *polarizations;

FILE *LOG=NULL, *FILE_LOG=NULL;
time_t start_time, end_time;

float *power=NULL; /* nsegments, nbins */
INT64 *gps=NULL;
long nsegments, nbins,side_cut;
long useful_bins;

float *det_velocity=NULL;
double average_det_velocity[3];
float det_vel_ra, det_vel_dec;
double orbital_axis[3];
double band_axis[3];
float band_axis_ra, band_axis_dec;

long first_bin;


float *TMedians=NULL,*FMedians=NULL, *expTMedians=NULL, *hours=NULL,*frequencies=NULL,*ks_test=NULL,*median=NULL,
      *max_residuals, *min_residuals;
float TMedian,expTMedian;
float *tm=NULL;
double *mean=NULL,*weighted_mean=NULL,*weight=NULL,*sigma=NULL,*freq_d=NULL,*hours_d=NULL;


int do_CutOff=1;

struct gengetopt_args_info args_info;

double spindown;
double orientation;
char *earth_ephemeris=NULL, *sun_ephemeris=NULL;
double resolution; /* this is actual resolution, not the resolution argument passed on command line */
int fake_injection=0;

int no_am_response;

SKY_GRID *fine_grid=NULL;
SKY_SUPERGRID *super_grid=NULL;
SKY_GRID *patch_grid=NULL;

SKY_GRID_TYPE *AM_coeffs_plus=NULL,*AM_coeffs_cross=NULL;
long AM_coeffs_size=0;

char *output_dir;


void *do_alloc(long a, long b)
{
void *r;
int i=0;
r=calloc(a,b);
while(r==NULL){
	fprintf(stderr,"Could not allocate %ld chunks of %ld bytes each (%ld bytes total)\n",a,b,a*b);
	if(i>10)exit(-1);
	sleep(10);
	r=calloc(a,b);
	i++;
	}
return r;
}

int float_cmp(float *a, float *b)
{
if(*a<*b)return -1;
if(*a>*b)return 1;
return 0;
}

static float compute_median(float *firstbin, long step, long count)
{
float *tmp;
long i;
tmp=alloca(count*sizeof(float));
for(i=0;i<count;i++)tmp[i]=firstbin[i*step];
qsort(tmp,count,sizeof(float),float_cmp);
if(!(count & 1))return (tmp[count>>1]+tmp[(count>>1)-1])/2.0;
return tmp[count>>1];
}

void compute_noise_curves(void)
{
float *tmp;
float *p,*t;
float a;
long i,j;
float b, b_initial;
tmp=do_alloc(nsegments*nbins,sizeof(float));
TRACE("Computing log10")
for(i=0;i<nsegments;i++){
	t=&(tmp[i*nbins]);
	p=&(power[i*nbins]);
	for(j=0;j<nbins;j++){
		t[j]=log10(p[j]);
		}
	}
TRACE("Computed log10")
TRACE("Factoring log power")
b=0;
for(i=0;i<nsegments;i++){
	a=compute_median(tmp+i*nbins,1,nbins);
	TMedians[i]=a;
	b+=a*a;
	t=&(tmp[i*nbins]);
	for(j=0;j<nbins;j++){
		t[j]-=a;
		}
	}
for(j=0;j<nbins;j++){
	a=compute_median(tmp+j,nbins,nsegments);
	FMedians[j]=a;
	b+=a*a;
	t=&(tmp[j]);
	for(i=0;i<nsegments;i++){
		t[i*nbins]-=a;
		}
	}
b_initial=b;
while(b>0){
	b=0;
	for(i=0;i<nsegments;i++){
		a=compute_median(tmp+i*nbins,1,nbins);
		TMedians[i]+=a;
		b+=a*a;
		t=&(tmp[i*nbins]);
		for(j=0;j<nbins;j++){
			t[j]-=a;
			}
		}
	for(j=0;j<nbins;j++){
		a=compute_median(tmp+j,nbins,nsegments);
		FMedians[j]+=a;
		b+=a*a;
		t=&(tmp[j]);
		for(i=0;i<nsegments;i++){
			t[i*nbins]-=a;
			}
		}
	fprintf(stderr,"%g\n",b);
	/* if one of r nbins or nsegments are even break because of order of
	   magnitude, do not try to solve exactly */
	if(!((nsegments &1)&&(nbins&1)) && (b<(b_initial*(1E-16))))break;
	}
TRACE("Factored log power")

for(i=0;i<nsegments;i++){
	max_residuals[i]=tmp[i*nbins];
	min_residuals[i]=tmp[i*nbins];
	for(j=1;j<nbins;j++){
		if(max_residuals[i]<tmp[i*nbins+j]){
			max_residuals[i]=tmp[i*nbins+j];
			}
		if(min_residuals[i]>tmp[i*nbins+j]){
			min_residuals[i]=tmp[i*nbins+j];
			}
		}
	}

free(tmp);
TMedian=compute_median(TMedians,1,nsegments);
fprintf(LOG,"Median noise level (TMedian): %g\n",TMedian);
fflush(LOG);
}

/* Note: this messes with tm variable,
  tm[i]=exp(log(10)*TMedians[i])/Modulation[i] */
float FindCutOff(float *tm)
{
long i;
double sum,sum_squared,mean,sigma;
double best_cutoff,smallest_sigma;
long best_i;
qsort(tm,nsegments,sizeof(float),float_cmp);
sum=0;
sum_squared=0;
best_i=0;
smallest_sigma=10*tm[0];
best_cutoff=tm[0];
for(i=0;i<nsegments;i++){
	sum+=tm[i];
	sum_squared+=tm[i]*tm[i];
	mean=sum/(i+1);
	sigma=10*sqrt((sum_squared))/(i+1);
	
	if(sigma<smallest_sigma){
		smallest_sigma=sigma;
		best_i=i;
		best_cutoff=tm[i];
		}
	}
//fprintf(stderr,"Cutoff: i=%ld sigma=%g cutoff=%g (%g,%g,..)\n",best_i, smallest_sigma, best_cutoff,tm[0],tm[1]);
return best_cutoff;
}

static double exponential_distribution(double x, double lambda)
{
if(x<=0.0)return 0.0;
return (1-exp(-lambda*x));
}

void nonparametric(float *firstbin, long step, long count, float *median, float *ks_test)
{
float *tmp;
double a,b,lambda;
long i;
tmp=alloca(count*sizeof(float));
for(i=0;i<count;i++)tmp[i]=firstbin[i*step];
qsort(tmp,count,sizeof(float),float_cmp);
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
fprintf(LOG,"seconds elapsed: %ld\n",end_time-start_time);
fprintf(stderr,"seconds elapsed: %ld\n",end_time-start_time);
fclose(LOG);
fclose(FILE_LOG);
}

void inject_fake_signal(void)
{
long i, bin;
float plus, cross, a, b, e[3], fake_power, factor;
gsl_rng *rng=NULL;
double average_freq, weight, mixed, plus_sq, cross_sq;

rng=gsl_rng_alloc(gsl_rng_default);
orientation=args_info.fake_orientation_arg;

fake_power=args_info.fake_strain_arg*args_info.fake_strain_arg*1800.0*16384.0*1800.0*16384.0;

average_freq=0.0;
weight=0.0;
mixed=0.0;
for(i=0;i<nsegments;i++){
	get_AM_response(gps[i]+900, args_info.fake_dec_arg, args_info.fake_ra_arg,
		&plus, &cross);
	/* effective power, ignoring phase, and taking into account windowing and
	   fourier coefficient (we are keeping only half of them)*/
//	a=plus*plus*fake_power*0.7*0.5;
	mixed+=plus*cross;
	plus_sq+=plus*plus;
	cross_sq+=cross*cross;
	factor=plus;
	a=factor*factor*fake_power*0.7*0.5;
	/* effective frequency */
	e[0]=cos(args_info.fake_dec_arg)*cos(args_info.fake_ra_arg);
	e[1]=cos(args_info.fake_dec_arg)*sin(args_info.fake_ra_arg);
	e[2]=sin(args_info.fake_dec_arg);
	b=args_info.fake_spindown_arg*(gps[i]-gps[0])+
	  args_info.fake_freq_arg*(1.0+(det_velocity[3*i]*e[0]+
		det_velocity[3*i+1]*e[1]+
		det_velocity[3*i+2]*e[2]));
	average_freq+=b*a;
	weight+=a;
	bin=rint(b*1800.0);
	if(bin<first_bin)continue;
	if(bin>=first_bin+nbins)continue;
	power[i*nbins+bin-first_bin]+=a+2.0*sqrt(a*power[i*nbins+bin-first_bin])*cos(gsl_rng_uniform(rng)*M_PI);
	}

average_freq/=weight;
fprintf(LOG,"average effective fake frequency: %f\n", average_freq);
fprintf(LOG,"mixed: %g\n",mixed/nsegments);
fprintf(LOG,"plus_sq: %g\n",plus_sq/nsegments);
fprintf(LOG,"cross_sq: %g\n",cross_sq/nsegments);

/* restore old value */
orientation=args_info.orientation_arg;
gsl_rng_free(rng);
}

int main(int argc, char *argv[])
{
RGBPic *p;
PLOT *plot;
long i,j,m,count;
float CutOff,b;
double a,w,a1,a2;
char s[20000];

time(&start_time);

if(cmdline_parser(argc, argv, &args_info))exit(-1);
if(args_info.config_given)
	if(cmdline_parser_configfile(args_info.config_arg, &args_info, 0))exit(-1);

if(!args_info.input_given){
	fprintf(stderr,"** You must specify patch to input files (--input)\n");
	exit(-1);
	}
if(!args_info.ephemeris_path_given &&
	!(args_info.earth_ephemeris_given && args_info.sun_ephemeris_given)){
	fprintf(stderr,"** You must specify patch to ephemeris files (--detreponse-path or --XXX-ephemeris)\n");
	exit(-1);
	}
if(!args_info.first_bin_given){
	fprintf(stderr,"** You must specify first bin to analyze (--first-bin)\n");
	exit(-1);
	}
if(!args_info.detector_given){
	fprintf(stderr,"** You must specify detector (--detector)\n");
	exit(-1);
	}

gsl_rng_env_setup();

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

snprintf(s,20000,"%s/powerflux.log", output_dir);
LOG=fopen(s,"w");

snprintf(s,20000,"%s/file.log", output_dir);
FILE_LOG=fopen(s,"w");

if(gethostname(s, 19999)>=0){
	fprintf(stderr, "Running on %s\n", s);
	fprintf(LOG, "node: %s\n", s);
	} else {
	fprintf(stderr, "Could not obtain hostname\n");
	fprintf(LOG, "node: unknown\n");
	}

init_hookup();

do_CutOff=args_info.do_cutoff_arg;
spindown=args_info.spindown_arg;
orientation=args_info.orientation_arg;

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

fprintf(stderr,	"Initializing sky grids\n");
resolution=(4500.0*args_info.resolution_ratio_arg)/(args_info.first_bin_arg+args_info.nbins_arg/2);
if(!strcasecmp("sin_theta", args_info.sky_grid_arg)){
	patch_grid=make_sin_theta_grid(resolution*args_info.fine_factor_arg);
	super_grid=make_sin_theta_supergrid(patch_grid, args_info.fine_factor_arg);
	} else
if(!strcasecmp("plain_rectangular", args_info.sky_grid_arg)){
	patch_grid=make_rect_grid(ceil(2.0*M_PI/(resolution*args_info.fine_factor_arg)), ceil(M_PI/(resolution*args_info.fine_factor_arg)));
	super_grid=make_rect_supergrid(patch_grid, args_info.fine_factor_arg, args_info.fine_factor_arg);
	} else
if(!strcasecmp("arcsin", args_info.sky_grid_arg)){
	patch_grid=make_arcsin_grid(ceil(2.0*M_PI/(resolution*args_info.fine_factor_arg)), ceil(M_PI/(resolution*args_info.fine_factor_arg)));
	super_grid=make_rect_supergrid(patch_grid, args_info.fine_factor_arg, args_info.fine_factor_arg);
	} else {
	fprintf(stderr,"*** Unknown sky grid type: \"%s\"\n", args_info.sky_grid_arg);
	exit(-1);
	}
fine_grid=super_grid->super_grid;

fprintf(stderr,"fine grid: max_n_ra=%d max_n_dec=%d\n", 
	fine_grid->max_n_ra, fine_grid->max_n_dec);


side_cut=args_info.side_cut_arg;
if(!args_info.side_cut_given){
	/* determine side cut from resolution, 6.0 factor is empirical */
	/* also add in spindown contribution - for now just plan for 4 months of data */
	side_cut=10+ceil(M_PI/resolution)/6.0+ceil(fabs((1800.0)*args_info.spindown_arg*10368000));
	}
fprintf(stderr,"side_cut=%d\n", side_cut);
first_bin=args_info.first_bin_arg-side_cut;
useful_bins=args_info.nbins_arg;
nbins=args_info.nbins_arg+2*side_cut;

if(fine_grid->max_n_dec<800){
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
fprintf(LOG,"firstbin  : %ld\n",first_bin);
fprintf(LOG,"band start: %g Hz\n",first_bin/1800.0);
fprintf(LOG,"nbins     : %ld\n",nbins);
fprintf(LOG,"side_cut  : %ld\n",side_cut);
fprintf(LOG,"useful bins : %ld\n",useful_bins);
fprintf(LOG,"useful band start: %g Hz\n",(first_bin+side_cut)/1800.0);
fprintf(LOG,"fine_factor: %ld\n", args_info.fine_factor_arg);
fprintf(LOG,"resolution : %f\n", resolution);

fprintf(LOG,"patch_type: %s\n", patch_grid->name);
fprintf(LOG,"patch_grid: %ldx%ld\n", patch_grid->max_n_ra, patch_grid->max_n_dec);
fprintf(LOG,"patch_grid npoints : %ld\n", patch_grid->npoints);
fprintf(LOG,"fine_type : %s\n",fine_grid->name);
fprintf(LOG,"fine_grid npoints  : %ld\n", fine_grid->npoints);
fprintf(LOG,"fine_grid : %ldx%ld\n", fine_grid->max_n_ra, fine_grid->max_n_dec);
fprintf(LOG,"input_data: %s\n",args_info.input_arg);
fprintf(LOG,"spindown  : %g\n", spindown);
fprintf(LOG,"orientation: %g\n", orientation);
fprintf(LOG,"make cutoff: %s\n",do_CutOff ? "yes" : "no" );
fflush(LOG);

read_directory(args_info.input_arg,1,4000, first_bin, nbins,
	&nsegments, &power, &gps);
if(nsegments==0){
	fprintf(stderr,"ERROR: no input data found !\n");
	return -1;
	}
fprintf(LOG,"nsegments : %ld\n",nsegments);
fprintf(LOG,"first gps : %lld\n",gps[0]);
fprintf(LOG,"last gps  : %lld\n",gps[nsegments-1]);
fflush(LOG);

fprintf(stderr,"Computing detector speed\n");
det_velocity=do_alloc(3*nsegments, sizeof(*det_velocity));
average_det_velocity[0]=0.0;
average_det_velocity[1]=0.0;
average_det_velocity[2]=0.0;
for(i=0;i<nsegments;i++){
	/* middle of the 30min interval */
	get_detector_vel(gps[i]+900,&(det_velocity[3*i]));
	average_det_velocity[0]+=det_velocity[3*i];
	average_det_velocity[1]+=det_velocity[3*i+1];
	average_det_velocity[2]+=det_velocity[3*i+2];
	}
average_det_velocity[0]/=nsegments;
average_det_velocity[1]/=nsegments;
average_det_velocity[2]/=nsegments;
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
a=sqrt(band_axis[0]*band_axis[0]+
	band_axis[1]*band_axis[1]+
	band_axis[2]*band_axis[2]);
/* replace 0.0 with something more reasonable later */
if(a<=0.0){
	band_axis[0]=0.0;
	band_axis[1]=0.0;
	band_axis[2]=1.0;
	} else {
	band_axis[0]/=a;
	band_axis[1]/=a;
	band_axis[2]/=a;
	}

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

/* assign bands */
assign_bands(patch_grid, args_info.nbands_arg);
assign_bands(fine_grid, args_info.nbands_arg);

if(args_info.focus_ra_given && 
   args_info.focus_dec_given && 
   args_info.focus_radius_given){
   	fprintf(LOG, "focus ra    : %f\n", args_info.focus_ra_arg);
   	fprintf(LOG, "focus dec   : %f\n", args_info.focus_dec_arg);
   	fprintf(LOG, "focus radius: %f\n", args_info.focus_radius_arg);
   	mask_far_points(fine_grid, args_info.focus_ra_arg, args_info.focus_dec_arg, args_info.focus_radius_arg);
	propagate_far_points_from_super_grid(patch_grid, super_grid);
   	}

if(args_info.only_large_cos_given){
	fprintf(LOG, "only large cos level: %f\n", args_info.only_large_cos_arg);
   	mask_small_cos(fine_grid, band_axis[0], band_axis[1], band_axis[3], args_info.only_large_cos_arg);
	propagate_far_points_from_super_grid(patch_grid, super_grid);
	}

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

plot_grid_f(p, fine_grid, fine_grid->band_f,1);
RGBPic_dump_png("bands.png", p);
dump_ints("bands.dat", fine_grid->band, fine_grid->npoints, 1);

/* do we need to inject fake signal ? */
if(args_info.fake_ra_given ||
   args_info.fake_dec_given ||
   args_info.fake_strain_given ||
   args_info.fake_orientation_given ||
   args_info.fake_spindown_given){

	fake_injection=1;
	if(!args_info.fake_freq_given){
		args_info.fake_freq_arg=(first_bin+nbins/2)/1800.0;
		}


   	fprintf(LOG,"fake signal injection: yes\n");
	fprintf(LOG,"fake ra : %f\n", args_info.fake_ra_arg);
	fprintf(LOG,"fake dec: %f\n", args_info.fake_dec_arg);
	fprintf(LOG,"fake orientation: %f\n", args_info.fake_orientation_arg);
	fprintf(LOG,"fake spindown: %g\n", args_info.fake_spindown_arg);
	fprintf(LOG,"fake strain: %g\n", args_info.fake_strain_arg);
	fprintf(LOG,"fake frequency: %f\n", args_info.fake_freq_arg);
	
	inject_fake_signal();
   	} else {
   	fprintf(LOG,"fake signal injection: none\n");
	}
	
/* compute time from the start of the run */
hours=do_alloc(nsegments, sizeof(*hours));
hours_d=do_alloc(nsegments, sizeof(*hours_d));
for(i=0;i<nsegments;i++){
	hours[i]=(1.0*(gps[i]-gps[0]))/3600.0;
	hours_d[i]=(1.0*(gps[i]-gps[0]))/3600.0;
	}
/* compute frequency array */
frequencies=do_alloc(nbins, sizeof(*frequencies));
freq_d=do_alloc(nbins, sizeof(*freq_d));
for(i=0;i<nbins;i++){
	freq_d[i]=(1.0*(first_bin+i))/1800.0;
	frequencies[i]=(1.0*(first_bin+i))/1800.0;
	}

mean=do_alloc(nbins, sizeof(*mean));
weighted_mean=do_alloc(nbins, sizeof(*weighted_mean));
weight=do_alloc(nbins, sizeof(*weight));
sigma=do_alloc(nbins, sizeof(*sigma));
median=do_alloc(nbins, sizeof(*median));
ks_test=do_alloc(nbins, sizeof(*ks_test));

/* first pass - raw statistics */
fprintf(stderr,"Computing background statistics\n");

for(i=0;i<nbins;i++){
	mean[i]=0;
	sigma[i]=0;
	}
for(i=0;i<nsegments;i++){
	for(j=0;j<nbins;j++){
		a=power[i*nbins+j];
		mean[j]+=a;
		sigma[j]+=a*a;
		}
	}
for(i=0;i<nbins;i++){
	mean[i]/=nsegments;
	sigma[i]=sqrt((sigma[i]-nsegments*mean[i]*mean[i])/(nsegments-1));
	nonparametric(power+i,nbins,nsegments,&(median[i]),&(ks_test[i]));
	}
	
adjust_plot_limits_d(plot, freq_d, mean, nbins, 1, 1, 1);
draw_grid(p, plot, 0, 0);
draw_points_d(p, plot, COLOR(255,0,0), freq_d, mean, nbins, 1, 1);
RGBPic_dump_png("mean.png", p);
dump_doubles("mean.dat", mean, nbins, 1);

adjust_plot_limits_d(plot, freq_d, sigma, nbins, 1, 1, 1);
draw_grid(p, plot, 0, 0);
draw_points_d(p, plot, COLOR(255,0,0), freq_d, sigma, nbins, 1, 1);
RGBPic_dump_png("sigma.png", p);

adjust_plot_limits_f(plot, frequencies, median, nbins, 1, 1, 1);
draw_grid(p, plot, 0, 0);
draw_points_f(p, plot, COLOR(255,0,0), frequencies, median, nbins, 1, 1);
RGBPic_dump_png("median.png", p);
dump_floats("median.dat", median, nbins, 1);

adjust_plot_limits_f(plot, frequencies, ks_test, nbins, 1, 1, 1);
draw_grid(p, plot, 0, 0);
draw_points_f(p, plot, COLOR(255,0,0), frequencies, ks_test, nbins, 1, 1);
RGBPic_dump_png("ks_test.png", p);
dump_floats("ks_test.dat", ks_test, nbins, 1);

if(args_info.no_decomposition_arg){
	fprintf(stderr,"Exiting as requested (--no-decomposition=1\n");
	fprintf(LOG,"Exiting as requested (--no-decomposition=1\n");
	wrap_up();
	exit(0);
	}

/* allocate TMedians, FMedians */
TMedians=do_alloc(nsegments, sizeof(*TMedians));
expTMedians=do_alloc(nsegments, sizeof(*expTMedians));
tm=do_alloc(nsegments, sizeof(*tm));
FMedians=do_alloc(nbins, sizeof(*FMedians));
max_residuals=do_alloc(nsegments, sizeof(*max_residuals));
min_residuals=do_alloc(nsegments, sizeof(*max_residuals));

/* decompose noise into FMedians, TMedians and residuals */
compute_noise_curves();

for(i=0;i<nsegments;i++){
	expTMedians[i]=exp(-M_LN10*2.0*(TMedians[i]-TMedian));
	}
expTMedian=exp(-M_LN10*2.0*TMedian);

adjust_plot_limits_f(plot, hours, TMedians, nsegments, 1, 1, 1);
draw_grid(p, plot, 0, 0);
draw_points_f(p, plot, COLOR(255,0,0), hours, TMedians, nsegments, 1, 1);
RGBPic_dump_png("TMedians.png", p);
dump_floats("TMedians.dat", TMedians, nsegments, 1);

adjust_plot_limits_f(plot, hours, max_residuals, nsegments, 1, 1, 1);
draw_grid(p, plot, 0, 0);
draw_points_f(p, plot, COLOR(255,0,0), hours, max_residuals, nsegments, 1, 1);
RGBPic_dump_png("max_residuals.png", p);
dump_floats("max_residuals.dat", max_residuals, nsegments, 1);

adjust_plot_limits_f(plot, hours, min_residuals, nsegments, 1, 1, 1);
draw_grid(p, plot, 0, 0);
draw_points_f(p, plot, COLOR(255,0,0), hours, min_residuals, nsegments, 1, 1);
RGBPic_dump_png("min_residuals.png", p);
dump_floats("min_residuals.dat", min_residuals, nsegments, 1);

adjust_plot_limits_f(plot, frequencies, FMedians,nbins, 1, 1, 1);
draw_grid(p, plot, 0, 0);
draw_points_f(p, plot, COLOR(255,0,0), frequencies, FMedians, nbins, 1, 1);
RGBPic_dump_png("FMedians.png", p);
dump_floats("FMedians.dat", FMedians, nbins, 1);

/* compute CutOff for non-demodulated analysis */
for(i=0;i<nsegments;i++)tm[i]=exp(M_LN10*TMedians[i]);
CutOff=log10(FindCutOff(tm));
fprintf(stderr,"CutOff=%g\n",CutOff);

for(i=0;i<nbins;i++){
	weighted_mean[i]=0;
	weight[i]=0;
	}
for(i=0;i<nsegments;i++){
	w=expTMedians[i];
	for(j=0;j<nbins;j++){
		a=power[i*nbins+j];
		weight[j]+=w;
		weighted_mean[j]+=a*w;
		}
	}
for(i=0;i<nbins;i++){
	weighted_mean[i]/=weight[i];
	}

adjust_plot_limits_d(plot, freq_d, weighted_mean, nbins, 1, 1, 1);
draw_grid(p, plot, 0, 0);
draw_points_d(p, plot, COLOR(255,0,0), freq_d, weighted_mean, nbins, 1, 1);
RGBPic_dump_png("weighted_mean.png", p);
dump_doubles("weighted_mean.dat", weighted_mean, nbins, 1);

/* second pass - apply CutOff */

for(i=0;i<nbins;i++){
	mean[i]=0;
	sigma[i]=0;
	weighted_mean[i]=0;
	weight[i]=0;
	}
count=0;
for(i=0;i<nsegments;i++){
	if(TMedians[i]>=CutOff)continue;
	w=expTMedians[i];
	count++;
	for(j=0;j<nbins;j++){
		a=power[i*nbins+j];
		mean[j]+=a;
		sigma[j]+=a*a;
		weight[j]+=w;
		weighted_mean[j]+=a*w;
		}
	}
for(i=0;i<nbins;i++){
	mean[i]/=count;
	sigma[i]=sqrt((sigma[i]-count*mean[i]*mean[i])/(count-1));
	weighted_mean[i]/=weight[i];
	}
adjust_plot_limits_d(plot, freq_d, mean, nbins, 1, 1, 1);
draw_grid(p, plot, 0, 0);
draw_points_d(p, plot, COLOR(255,0,0), freq_d, mean, nbins, 1, 1);
RGBPic_dump_png("new_mean.png", p);
dump_doubles("new_mean.dat", mean, nbins, 1);

adjust_plot_limits_d(plot, freq_d, weighted_mean, nbins, 1, 1, 1);
draw_grid(p, plot, 0, 0);
draw_points_d(p, plot, COLOR(255,0,0), freq_d, weighted_mean, nbins, 1, 1);
RGBPic_dump_png("new_weighted_mean.png", p);
dump_doubles("new_weighted_mean.dat", weighted_mean, nbins, 1);

adjust_plot_limits_d(plot, freq_d, sigma, nbins, 1, 1, 1);
draw_grid(p, plot, 0, 0);
draw_points_d(p, plot, COLOR(255,0,0), freq_d, sigma, nbins, 1, 1);
RGBPic_dump_png("new_sigma.png", p);

fprintf(stderr,"Checking background lines\n");
detect_background_lines(weighted_mean);

if(args_info.no_demodulation_arg){
	fprintf(stderr,"Exiting as requested (--no-demodulation=1\n");
	fprintf(LOG,"Exiting as requested (--no-demodulation=1\n");
	wrap_up();
	exit(0);	
	}

get_whole_sky_AM_response(gps, nsegments, &AM_coeffs_plus, &AM_coeffs_cross, &AM_coeffs_size);

init_polarizations();

#if 0 /* verify patch grid */
{
float *patch_cross=NULL, *patch_plus=NULL;

fprintf(stderr,"Computing patch grid\n");
fprintf(LOG,"patch grid size: %f KB\n",
	2*patch_grid->npoints*nsegments*sizeof(*patch_cross)/(1024.0));
fflush(LOG);
patch_cross=do_alloc(patch_grid->npoints*nsegments, sizeof(*patch_cross));
patch_plus=do_alloc(patch_grid->npoints*nsegments, sizeof(*patch_plus));

for(i=0;i<nsegments;i++){
	for(j=0;j<patch_grid->npoints;j++){
		get_AM_response(gps[i]+900, 
			patch_grid->latitude[j], 
			patch_grid->longitude[j], 
			&(patch_plus[j+i*patch_grid->npoints]), &(patch_cross[j+i*patch_grid->npoints]));
		fprintf(stderr,"%g ", patch_plus[j+i*patch_grid->npoints]-AM_response(i, patch_grid, j, AM_coeffs_plus));
		fprintf(stderr,"%g ", patch_cross[j+i*patch_grid->npoints]-AM_response(i, patch_grid, j, AM_coeffs_cross));

		patch_plus[j+i*patch_grid->npoints]-=AM_response(i, patch_grid, j, AM_coeffs_plus);
		patch_cross[j+i*patch_grid->npoints]-=AM_response(i, patch_grid, j, AM_coeffs_cross);

		}
	}
dump_floats("patch_plus0.dat",patch_plus,patch_grid->npoints,1);
dump_floats("patch_cross0.dat",patch_cross,patch_grid->npoints,1);

/* just testing.. */
adjust_plot_limits_f(plot, hours, patch_plus, nsegments, 1, patch_grid->npoints, 1);
draw_grid(p, plot, 0, 0);
draw_points_f(p, plot, COLOR(255,0,0), hours, patch_plus, nsegments, 1, patch_grid->npoints);
RGBPic_dump_png("patch_plus0.png", p);

/* just testing.. */
plot_grid_f(p, patch_grid, patch_plus,1);
RGBPic_dump_png("patch_plus_gps0.png", p);

plot_grid_f(p, patch_grid, patch_cross,1);
RGBPic_dump_png("patch_cross_gps0.png", p);
}
#endif		

fprintf(stderr,"Computing cutoff values\n");
/* compute CutOff values for each patch */

for(i=0;i<patch_grid->npoints;i++){
	
	for(m=0;m<npolarizations;m++){
		if(patch_grid->band[i]<0){
			polarizations[m].patch_CutOff[i]=0.0;
			continue;
			}
		for(j=0;j<nsegments;j++)tm[j]=1.0/(sqrt(expTMedians[j])*AM_response(j, patch_grid, i, polarizations[m].AM_coeffs));
		polarizations[m].patch_CutOff[i]=FindCutOff(tm);
		}
	}

for(m=0;m<npolarizations;m++){
	plot_grid_f(p, patch_grid, polarizations[m].patch_CutOff,1);
	snprintf(s,20000,"patch_CutOff_%s.png",polarizations[m].name);
	RGBPic_dump_png(s, p);
	}
	
#if 0
/* free patch modulations data - it can be quite large, esp for S3 */
free(patch_plus);
patch_plus=NULL;
free(patch_cross);
patch_cross=NULL;
#endif

fine_grid_stage();

wrap_up();
return 0;
}
