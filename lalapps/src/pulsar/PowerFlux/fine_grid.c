#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <string.h>

#include "global.h"
#include "rastermagic.h"
#include "cmdline.h"
#include "hookup.h"
#include "grid.h"
#include "polarization.h"
#include "statistics.h"

extern POLARIZATION *polarizations;
extern int npolarizations;

extern float *power;
extern float *TMedians,*expTMedians;
extern float TMedian;

extern struct gengetopt_args_info args_info;

extern long nsegments, nbins, first_bin, side_cut, useful_bins;

extern INT64 *gps;

extern int lines_list[];

extern SKY_GRID *fine_grid, *patch_grid;
extern SKY_SUPERGRID *super_grid;
extern float *det_velocity;

extern int do_CutOff;
extern int fake_injection;

extern FILE *LOG;

extern double spindown;
extern double resolution;

extern float *frequencies;

extern INT64 *gps;

SUM_TYPE normalizing_weight;

long stored_fine_bins=0;

long max_shift=0, min_shift=0;


SUM_TYPE  *low_ul, *low_ul_freq;  /* lower (over polarizations) accumulation limits */
SUM_TYPE *skymap_low_ul, *skymap_low_ul_freq; /* skymaps */
SUM_TYPE *spectral_plot_low_ul; /* spectral plots */

/* the value 1.22 is obtained by S.R script. It is valid for single-bin processing
   the value 0.88 is obtained by the same script. It is valid for 3 averaged neighboring bins */

float quantile2std=1.22;

/* We divide by 0.7 to compensate for spreading due to Hanning window */
/* We also divide by 0.85 to compensate for non bin centered signals */
/* for 3 averaged neighbouring bins that value should be 3.0/0.85 */
float upper_limit_comp=1.0/(0.7*0.85);

/* single bin version */

void (*process_patch)(POLARIZATION *pol, long pi, long k, float CutOff);

static void process_patch1(POLARIZATION *pol, long pi, long k, float CutOff)
{
long i,kk,b,b0,b1,n,offset;
long bin_shift;
SUM_TYPE mod;
SUM_TYPE a,w,w2;
SUM_TYPE *sum,*sq_sum;
float *p;
float doppler;

CutOff=2*CutOff; /* weighted sum can benefit from more SFTs */


for(i=0,kk=super_grid->first_map[pi];kk>=0;kk=super_grid->list_map[kk],i++)
		{

		if(fine_grid->band[kk]<0)continue;
		
		mod=1.0/(AM_response(k, fine_grid, kk, pol->AM_coeffs)); 

		if(do_CutOff && (mod>CutOff))continue;
		
		/* this assumes that both relfreq and spindown are small enough that the bin number
		   down not change much within the band - usually true */
		doppler=fine_grid->e[0][kk]*det_velocity[3*k+0]+
			fine_grid->e[1][kk]*det_velocity[3*k+1]+
			fine_grid->e[2][kk]*det_velocity[3*k+2];

		bin_shift=-rint((first_bin+nbins*0.5)*doppler+1800.0*spindown*(gps[k]-gps[0]));

		if(bin_shift>max_shift)max_shift=bin_shift;
		if(bin_shift<min_shift)min_shift=bin_shift;
		b0=side_cut-bin_shift;
		b1=(nbins-side_cut)-bin_shift;
		if((b0<0)||(b1>nbins)){
			fprintf(stderr,"Working frequency range obscured by bin_shift shift: bin_shift=%ld kk=%ld i=%ld pi=%ld\n",
				bin_shift, kk, i, pi);
			exit(-1);
			}
		
		/* prime array pointers */
		#ifdef WEIGHTED_SUM
		w2=expTMedians[k]/mod;
		w=w2/mod;
		pol->skymap.total_weight[kk]+=w;
		#else
		pol->skymap.total_count[kk]++;
		#endif
		
		sum=&(pol->fine_grid_sum[b0-side_cut+bin_shift+useful_bins*i]);
		#ifdef COMPUTE_SIGMA
		sq_sum=&(pol->fine_grid_sq_sum[b0-side_cut+bin_shift+useful_bins*i]);
		#endif
		p=&(power[k*nbins+b0]);
		
		/* cycle over bins */
		for(b=b0;b<b1;b++){			    
			#ifdef WEIGHTED_SUM
			a=(*p)*w2;
			(*sum)+=a;
			#else
			a=(*p)*mod;
			(*sum)+=a;
			#ifdef COMPUTE_SIGMA
			(*sq_sum)+=a*a;
			sq_sum++;
			#endif			
			#endif
		
			p++;
			sum++;
			}
		/* subtract lines */
		for(n=0;(lines_list[n]>=0)&&(n<5);n++){
			b=lines_list[n];
			if(b<b0)continue;
			if(b>=b1)continue;
			offset=b-side_cut+bin_shift+useful_bins*i;
			/* prime array pointers */
			sum=&(pol->fine_grid_sum[offset]);
			p=&(power[k*nbins+b]);

			#ifdef WEIGHTED_SUM
			pol->fine_grid_weight[offset]+=w;
			#else
			pol->fine_grid_count[offset]++;
			#endif
			
			#ifdef WEIGHTED_SUM
			a=(*p)*w2;
			(*sum)-=a;
			#else
			a=(*p)*mod;
			(*sum)-=a;
			#ifdef COMPUTE_SIGMA
			pol->fine_grid_sq_sum[offset]-=a*a;
			#endif			
			#endif
			}
		}

}

/* three bin version */

static void process_patch3(POLARIZATION *pol, long pi, long k, float CutOff)
{
long i,kk,b,b0,b1,n,offset;
long bin_shift;
SUM_TYPE mod;
SUM_TYPE a,w,w2;
SUM_TYPE *sum,*sq_sum;
float *p;
float doppler;

CutOff=2*CutOff/3.0; /* weighted sum can benefit from more SFTs */


for(i=0,kk=super_grid->first_map[pi];kk>=0;kk=super_grid->list_map[kk],i++)
		{
		if(fine_grid->band[kk]<0)continue;

		mod=1.0/(AM_response(k, fine_grid, kk, pol->AM_coeffs)); 

		if(do_CutOff && (mod>CutOff))continue;
		
		/* this assumes that both relfreq and spindown are small enough that the bin number
		   down not change much within the band - usually true */
		doppler=fine_grid->e[0][kk]*det_velocity[3*k+0]+
			fine_grid->e[1][kk]*det_velocity[3*k+1]+
			fine_grid->e[2][kk]*det_velocity[3*k+2];

		bin_shift=-rint((first_bin+nbins*0.5)*doppler+1800.0*spindown*(gps[k]-gps[0]));

		if(bin_shift>max_shift)max_shift=bin_shift;
		if(bin_shift<min_shift)min_shift=bin_shift;
		b0=side_cut-bin_shift;
		b1=(nbins-side_cut)-bin_shift;
		if((b0<1)||(b1>(nbins-1))){
			fprintf(stderr,"Working frequency range obscured by bin_shift shift: bin_shift=%ld kk=%ld i=%ld pi=%ld\n",
				bin_shift, kk, i, pi);
			exit(-1);
			}
		
		/* prime array pointers */
		#ifdef WEIGHTED_SUM
		w2=expTMedians[k]/mod;
		w=w2/(3.0*mod);
		pol->skymap.total_weight[kk]+=w;
		#else
		pol->skymap.total_count[kk]++;
		#endif
		
		sum=&(pol->fine_grid_sum[b0-side_cut+bin_shift+useful_bins*i]);
		#ifdef COMPUTE_SIGMA
		sq_sum=&(pol->fine_grid_sq_sum[b0-side_cut+bin_shift+useful_bins*i]);
		#endif
		p=&(power[k*nbins+b0]);
		
		/* cycle over bins */
		for(b=b0;b<b1;b++){			    
			#ifdef WEIGHTED_SUM
			a=(p[-1]+p[0]+p[1])*w2;
			(*sum)+=a;
			#else
			a=(*p)*mod;
			(*sum)+=a;
			#ifdef COMPUTE_SIGMA
			(*sq_sum)+=a*a;
			sq_sum++;
			#endif			
			#endif
		
			p++;
			sum++;
			}
		/* subtract lines */
		for(n=0;(lines_list[n]>=0)&&(n<5);n++){
			b=lines_list[n];
			if(b<b0)continue;
			if(b>=b1)continue;
			offset=b-side_cut+bin_shift+useful_bins*i;
			/* prime array pointers */
			sum=&(pol->fine_grid_sum[offset]);
			p=&(power[k*nbins+b]);

			#ifdef WEIGHTED_SUM
			pol->fine_grid_weight[offset]+=w;
			#else
			pol->fine_grid_count[offset]++;
			#endif
			
			#ifdef WEIGHTED_SUM
			a=(p[-1]+p[0]+p[1])*w2;
			(*sum)-=a;
			#else
			a=(p[-1]+p[0]+p[1])*mod;
			(*sum)-=a;
			#ifdef COMPUTE_SIGMA
			pol->fine_grid_sq_sum[offset]-=a*a;
			#endif			
			#endif
			}
		}

}

void dump_pic(char *file, double *z)
{
RGBPic *p;

if(!clear_name_png(file))return;

p=make_RGBPic(fine_grid->max_n_ra+140, fine_grid->max_n_dec);

plot_grid_d(p, fine_grid, z, 1);
RGBPic_dump_png(file, p);

free_RGBPic(p);
}

int float_cmp(float *a, float *b);

static float upper_limit90(float dx)
{
if(dx>=0.0)return dx+1.64;
return dx+sqrt(dx*dx+1.64);
}

static inline float upper_limit95(float dx)
{
if(dx>=-0.1)return dx+1.96;
dx+=0.1;
return dx+sqrt(dx*dx+1.86);
}

static inline float lower_limit95(float dx)
{
if(dx<1.7)return 0;
if(dx>=2.1)return (0.455+0.766*(dx-2.1));
if(dx<1.8)return 0.16;
if(dx<1.9)return 0.26;
if(dx<2.0)return 0.35;
return 0.455;
}

void make_limits(POLARIZATION *pol, long pi)
{
SUM_TYPE M,Q80,Q20,S,dx;
SUM_TYPE a,b,c;
SUM_TYPE *tmp=NULL;
long i,j,k,offset,band;
NORMAL_STATS nstats;

/* allocate on stack, for speed */
tmp=alloca(useful_bins*sizeof(*tmp));

/* sort to compute robust estimates */
nstats.flag=STAT_FLAG_INPLACE_SORT_DATA
	| STAT_FLAG_ESTIMATE_MEAN
	| STAT_FLAG_ESTIMATE_SIGMA;

if(args_info.ks_test_arg){
	nstats.flag|=STAT_FLAG_ESTIMATE_KS_LEVEL
		| STAT_FLAG_COMPUTE_KS_TEST;
	}

for(i=0,offset=super_grid->first_map[pi];offset>=0;offset=super_grid->list_map[offset],i++){

        band=fine_grid->band[offset];
	if(band<0)continue;

	/* figure out offset to put results into */

	memcpy(tmp,pol->fine_grid_sum+i*useful_bins,useful_bins*sizeof(*tmp));
	/* compute correlations - for diagnostics only */
	a=0.0;
	b=0.0;
	for(j=0;j<useful_bins-2;j++){
		pol->skymap.cor1[offset]+=(tmp[j]*normalizing_weight)*(tmp[j+1]*normalizing_weight);
		pol->skymap.cor2[offset]+=(tmp[j]*normalizing_weight)*(tmp[j+2]*normalizing_weight);
		b+=(tmp[j]*normalizing_weight);
		a+=(tmp[j]*normalizing_weight)*(tmp[j]*normalizing_weight);
		}
	c=(b-tmp[0]-tmp[1]+tmp[useful_bins-2]+tmp[useful_bins-1]);
	pol->skymap.cor2[offset]=(pol->skymap.cor2[offset]-b*c/(useful_bins-2))/
		(sqrt((a-b*b/(useful_bins-2))*
		(a-tmp[0]*tmp[0]-tmp[1]*tmp[1]
			+tmp[useful_bins-2]*tmp[useful_bins-2]+tmp[useful_bins-1]*tmp[useful_bins-1]
			-c*c/(useful_bins-2))));
			
	b+=tmp[useful_bins-2];
	a+=tmp[useful_bins-2]*tmp[useful_bins-2];
	pol->skymap.cor1[offset]+=tmp[useful_bins-2]*tmp[useful_bins-1];
	c=b-tmp[0]+tmp[useful_bins-1];
	pol->skymap.cor1[offset]=(pol->skymap.cor1[offset]-b*c/(useful_bins-1))/
		(sqrt(
		(a-b*b/(useful_bins-1))*
		(a-tmp[0]*tmp[0]+tmp[useful_bins-1]*tmp[useful_bins-1]-c*c/(useful_bins-1))
		));
	
	/* Output point data if requested */
	if(args_info.dump_points_arg){
		char s[20000];
		snprintf(s, 20000, "points/%s_%d.png", pol->name, offset);
		if(clear_name_png(s)){
			RGBPic *p;
			PLOT *plot;
			float *freq_f;

			freq_f=&(frequencies[side_cut]);
		
			if(fine_grid->max_n_dec<800)
				p=make_RGBPic(fine_grid->max_n_ra*(800/fine_grid->max_n_dec)+140, fine_grid->max_n_dec*(800/fine_grid->max_n_dec));
				else 
				p=make_RGBPic(fine_grid->max_n_ra+140, fine_grid->max_n_dec);	

			plot=make_plot(p->width, p->height);


			adjust_plot_limits_f(plot, freq_f, tmp, useful_bins, 1, 1, 1);
			draw_grid(p, plot, 0, 0);
			draw_points_f(p, plot, COLOR(255,0,0), freq_f, tmp, useful_bins, 1, 1);
			RGBPic_dump_png(s, p);

			free_plot(plot);
			free_RGBPic(p);
			}
			
		snprintf(s, 20000, "points/%s_%d.dat", pol->name, offset);
		dump_floats(s, tmp, useful_bins, 1);
		}
	
	compute_normal_stats(tmp, useful_bins, &nstats);

	pol->skymap.ks_test[offset]=nstats.ks_test;
	pol->skymap.ks_count[offset]=nstats.ks_count;
	
	M=nstats.mean;
	S=nstats.sigma;
	
	pol->skymap.M_map[offset]=M;
	pol->skymap.S_map[offset]=S;
	pol->skymap.max_upper_limit[offset]=0;
	
	for(k=0;k<useful_bins;k++){
		dx=(pol->fine_grid_sum[i*useful_bins+k]-M)/S;		
		a=upper_limit95(dx)*S;
		if(a>pol->skymap.max_upper_limit[offset]){
			pol->skymap.max_upper_limit[offset]=a;
			pol->skymap.freq_map[offset]=(first_bin+side_cut+k)/1800.0;
			}
		if(a>pol->spectral_plot.max_upper_limit[k+band*useful_bins]){
			pol->spectral_plot.max_upper_limit[k+band*useful_bins]=a;
			pol->spectral_plot.ul_ra[k+band*useful_bins]=fine_grid->longitude[offset];
			pol->spectral_plot.ul_dec[k+band*useful_bins]=fine_grid->latitude[offset];
			}

		if(dx>pol->skymap.max_dx[offset]){
			pol->skymap.max_dx[offset]=dx;
			}
		if(dx>pol->spectral_plot.max_dx[k+band*useful_bins]){
			pol->spectral_plot.max_dx[k+band*useful_bins]=dx;
			pol->spectral_plot.dx_ra[k+band*useful_bins]=fine_grid->longitude[offset];
			pol->spectral_plot.dx_dec[k+band*useful_bins]=fine_grid->latitude[offset];
			}
		if(a<low_ul[i*useful_bins+k]){
			low_ul[i*useful_bins+k]=a;
			low_ul_freq[i*useful_bins+k]=(first_bin+side_cut+k)/1800.0;
			}

			
		a=lower_limit95(dx)*S;
		if(a>pol->skymap.max_lower_limit[offset]){
			pol->skymap.max_lower_limit[offset]=a;
			}
			
		if(lines_list[0]<0)a=0.0;
			else {
			#ifdef WEIGHTED_SUM
			a=pol->fine_grid_weight[i*useful_bins+k]/pol->skymap.total_weight[offset];
			#else
			a=pol->fine_grid_count[i*useful_bins+k]/pol->skymap.total_count[offset];
			#endif
			}
		if(a>pol->spectral_plot.max_mask_ratio[k+band*useful_bins]){
			pol->spectral_plot.max_mask_ratio[k+band*useful_bins]=a;
			}
		}
	}

}

void make_unified_limits(long pi)
{
long i, offset;
int band, k;
SUM_TYPE a,b,c;
for(i=0,offset=super_grid->first_map[pi];offset>=0;offset=super_grid->list_map[offset],i++){

        band=fine_grid->band[offset];
	if(band<0)continue;

	for(k=0;k<useful_bins;k++){
		a=low_ul[i*useful_bins+k];
		if(a>spectral_plot_low_ul[k+band*useful_bins]){
			spectral_plot_low_ul[k+band*useful_bins]=a;
			}
		if(a>skymap_low_ul[offset]){
			skymap_low_ul[offset]=a;
			skymap_low_ul_freq[offset]=low_ul_freq[i*useful_bins+k];
			}
		}
	}
}

void output_limits(POLARIZATION *pol)
{
char s[20000];
RGBPic *p;
PLOT *plot;
long i, max_dx_i, masked, k;
SUM_TYPE max_dx;
float *max_band, *masked_max_band;
long *max_band_arg, *masked_max_band_arg;
float *freq_f;
float max_ratio;
HISTOGRAM *hist;

freq_f=&(frequencies[side_cut]);

max_band=do_alloc(args_info.nbands_arg, sizeof(*max_band));
masked_max_band=do_alloc(args_info.nbands_arg, sizeof(*max_band));
max_band_arg=do_alloc(args_info.nbands_arg, sizeof(*max_band_arg));
masked_max_band_arg=do_alloc(args_info.nbands_arg, sizeof(*max_band_arg));
hist=new_histogram(args_info.hist_bins_arg, args_info.nbands_arg);

if(fine_grid->max_n_dec<800){
	p=make_RGBPic(fine_grid->max_n_ra*(800/fine_grid->max_n_dec)+140, fine_grid->max_n_dec*(800/fine_grid->max_n_dec));
	} else 
	p=make_RGBPic(fine_grid->max_n_ra+140, fine_grid->max_n_dec);	

plot=make_plot(p->width, p->height);

snprintf(s,19999,"%s_weight.png",pol->name);
if(clear_name_png(s)){
	plot_grid_f(p, fine_grid, pol->skymap.total_weight, 1);
	RGBPic_dump_png(s, p);
	}

snprintf(s,19999,"%s_cor1.png",pol->name);
if(clear_name_png(s)){
	plot_grid_f(p, fine_grid, pol->skymap.cor1, 1);
	RGBPic_dump_png(s, p);
	}

snprintf(s,19999,"%s_cor2.png",pol->name);
if(clear_name_png(s)){
	plot_grid_f(p, fine_grid, pol->skymap.cor2, 1);
	RGBPic_dump_png(s, p);
	}

if(args_info.ks_test_arg){
	snprintf(s,19999,"%s_ks_test.png",pol->name);
	if(clear_name_png(s)){
		plot_grid_f(p, fine_grid, pol->skymap.ks_test, 1);
		RGBPic_dump_png(s, p);
		}
	snprintf(s,19999,"%s_ks_test.dat",pol->name);
	dump_floats(s,pol->skymap.ks_test,fine_grid->npoints,1);
	compute_histogram_f(hist, pol->skymap.ks_test, fine_grid->band, fine_grid->npoints);
	snprintf(s,19999,"hist_%s_ks_test",pol->name);
	print_histogram(LOG, hist, s);
	
	snprintf(s,19999,"%s_ks_count.png",pol->name);
	if(clear_name_png(s)){
		plot_grid_f(p, fine_grid, pol->skymap.ks_count, 1);
		RGBPic_dump_png(s, p);
		}
	snprintf(s,19999,"%s_ks_count.dat",pol->name);
	dump_floats(s,pol->skymap.ks_count,fine_grid->npoints,1);
	}
	
snprintf(s,19999,"%s_max_upper_limit.png",pol->name);
if(clear_name_png(s)){
	plot_grid_f(p, fine_grid, pol->skymap.max_upper_limit, 1);
	RGBPic_dump_png(s, p);
	}

snprintf(s,19999,"%s_max_lower_limit.png",pol->name);
if(clear_name_png(s)){
	plot_grid_f(p, fine_grid, pol->skymap.max_lower_limit, 1);
	RGBPic_dump_png(s, p);
	}

snprintf(s,19999,"%s_arg_freq.png",pol->name);
if(clear_name_png(s)){
	plot_grid_f(p, fine_grid, pol->skymap.freq_map, 1);
	RGBPic_dump_png(s, p);
	}

snprintf(s,19999,"%s_arg_freq.dat",pol->name);
dump_floats(s,pol->skymap.freq_map,fine_grid->npoints,1);

snprintf(s,19999,"%s_max_dx.dat",pol->name);
dump_floats(s,pol->skymap.max_dx,fine_grid->npoints,1);

snprintf(s,19999,"%s_S_map.dat",pol->name);
dump_floats(s,pol->skymap.S_map,fine_grid->npoints,1);
	

for(i=0;i<fine_grid->npoints;i++){
	if(fine_grid->band[i]<0){
		pol->skymap.max_upper_limit[i]=-1.0;
		pol->skymap.max_lower_limit[i]=-1.0;
		continue;
		}
	pol->skymap.max_upper_limit[i]=sqrt(2.0*pol->skymap.max_upper_limit[i]*upper_limit_comp)/(1800.0*16384.0);
	/* lower limit is unchanged */
	pol->skymap.max_lower_limit[i]=sqrt(2.0*pol->skymap.max_lower_limit[i])/(1800.0*16384.0);
	}

/* output interesting points around fake injection */
if(fake_injection){
	float largest;
	double ds, best_ds;
	long best_i=-1, largest_i=-1;
	fprintf(LOG,"Interesting points: index longitude latitude pol max_dx upper_strain lower_strain freq\n");
	for(i=0;i<fine_grid->npoints;i++){
		/* e[2][i] is just cos of latitude */
		/* Approximate spherical distance */
		#if 1
		ds=sqr_f(fine_grid->latitude[i]-args_info.fake_dec_arg)+
			sqr_f((fine_grid->longitude[i]-args_info.fake_ra_arg)*fine_grid->e[3][i]);
		if(ds<9*resolution*resolution){
		#elif 0
		/* Simplified spherical distance - should be better, different from exact
		   by a monotonic function, but harder to find interesting points.  */
		ds=1.0-(fine_grid->e[2][i]*sin(args_info.fake_dec_arg)+
			fine_grid->e[3][i]*cos(args_info.fake_dec_arg)*
			cos(fine_grid->longitude[i]-args_info.fake_ra_arg));
		if(args_info.fake_dec_arg*fine_grid->latitude[i]<0)ds=1;
		if(ds<9*resolution){
		#else
		/* Exact spherical distance */
		ds=acos(fine_grid->e[2][i]*sin(args_info.fake_dec_arg)+
			fine_grid->e[3][i]*cos(args_info.fake_dec_arg)*
			cos(fine_grid->longitude[i]-args_info.fake_ra_arg));		
		if(args_info.fake_dec_arg*fine_grid->latitude[i]<0)ds=1;
		if(ds<3*resolution){
		#endif
		
		   	fprintf(LOG, "%d %f %f %s %f %g %g %f\n",
				i,
				fine_grid->longitude[i], fine_grid->latitude[i], 
				pol->name, pol->skymap.max_dx[i], 
				pol->skymap.max_upper_limit[i], pol->skymap.max_lower_limit[i],
				pol->skymap.freq_map[i]);

			if(largest_i<0){
				largest=pol->skymap.max_upper_limit[i];
				largest_i=i;
				} else 
			if(largest<pol->skymap.max_upper_limit[i]){
				largest=pol->skymap.max_upper_limit[i];
				largest_i=i;
				}
		   	}

		if(best_i<0){
			best_ds=ds;
			best_i=i;
			} else
		if(ds<best_ds){
			best_ds=ds;
			best_i=i;
			}
		
		}
	if(best_i>=0)
	fprintf(LOG, "closest: %f %f %s %f %g %g %f\n",fine_grid->longitude[best_i], fine_grid->latitude[best_i], 
		pol->name, pol->skymap.max_dx[best_i], 
		pol->skymap.max_upper_limit[best_i], pol->skymap.max_lower_limit[best_i], pol->skymap.freq_map[best_i]);
	if(largest_i>=0)
	fprintf(LOG, "largest: %f %f %s %f %g %g %f\n",fine_grid->longitude[largest_i], fine_grid->latitude[largest_i], 
		pol->name, pol->skymap.max_dx[largest_i], 
		pol->skymap.max_upper_limit[largest_i], pol->skymap.max_lower_limit[largest_i], pol->skymap.freq_map[largest_i]);
	}



snprintf(s,19999,"%s_max_strain.dat",pol->name);
dump_floats(s,pol->skymap.max_upper_limit,fine_grid->npoints,1);

max_dx=0.0;
max_dx_i=0;
masked=0;
for(i=0;i<args_info.nbands_arg;i++){
	max_band[i]=-1.0;
	masked_max_band[i]=-1.0;
	max_band_arg[i]=-1;
	masked_max_band_arg[i]=-1;
	}
	
for(i=0;i<fine_grid->npoints;i++){
	k=fine_grid->band[i];
	if(k<0)continue;

	if(pol->skymap.max_upper_limit[i]>max_band[k]){
		max_band[k]=pol->skymap.max_upper_limit[i];
		max_band_arg[k]=i;
		}

	if(pol->skymap.max_sub_weight[i]>=pol->skymap.total_weight[i]*(1-args_info.small_weight_ratio_arg)){
		pol->skymap.max_upper_limit[i]=0.0;
		pol->skymap.max_lower_limit[i]=0.0;
		pol->skymap.max_dx[i]=0.0;
		masked++;
		}
	if(pol->skymap.max_dx[i]>max_dx){
		max_dx=pol->skymap.max_dx[i];
		max_dx_i=i;
		}

	if(pol->skymap.max_upper_limit[i]>masked_max_band[k]){
		masked_max_band[k]=pol->skymap.max_upper_limit[i];
		masked_max_band_arg[k]=i;
		}
	}
fprintf(LOG, "masked: %s %ld\n", pol->name, masked);
fprintf(LOG, "largest signal: longitude latitude pol max_dx upper_strain lower_strain freq\n");	
fprintf(LOG, "max_dx: %f %f %s %f %g %g %f\n",fine_grid->longitude[max_dx_i], fine_grid->latitude[max_dx_i], 
				pol->name, pol->skymap.max_dx[max_dx_i], 
				pol->skymap.max_upper_limit[max_dx_i], 
				pol->skymap.max_lower_limit[max_dx_i], pol->skymap.freq_map[max_dx_i]);

fprintf(LOG, "max/masked band format: band_num longitude latitude pol max_dx upper_strain freq\n");
for(i=0;i<args_info.nbands_arg;i++){

	fprintf(LOG, "max_band: %ld %f %f %s %f %g %f\n", i, fine_grid->longitude[max_band_arg[i]], fine_grid->latitude[max_band_arg[i]], 
				pol->name, pol->skymap.max_dx[max_band_arg[i]], 
				max_band[i], 
				pol->skymap.freq_map[max_band_arg[i]]);

	fprintf(LOG, "masked_max_band: %ld %f %f %s %f %g %f\n", i, fine_grid->longitude[masked_max_band_arg[i]], fine_grid->latitude[masked_max_band_arg[i]], 
				pol->name, pol->skymap.max_dx[masked_max_band_arg[i]], 
				masked_max_band[i], 
				pol->skymap.freq_map[masked_max_band_arg[i]]);

	snprintf(s,19999,"%s_max_upper_limit_band_%d.png",pol->name, i);
	if(clear_name_png(s)){
		adjust_plot_limits_f(plot, freq_f, &(pol->spectral_plot.max_upper_limit[i*useful_bins]), useful_bins, 1, 1, 1);
		draw_grid(p, plot, 0, 0);
		draw_points_f(p, plot, COLOR(255,0,0), freq_f, &(pol->spectral_plot.max_upper_limit[i*useful_bins]), useful_bins, 1, 1);
		RGBPic_dump_png(s, p);
		}
	snprintf(s,19999,"%s_max_upper_limit_band_%d.dat",pol->name, i);
	dump_floats(s, &(pol->spectral_plot.max_upper_limit[i*useful_bins]), useful_bins, 1);
	
	snprintf(s,19999,"%s_max_dx_band_%d.png",pol->name, i);
	if(clear_name_png(s)){
		adjust_plot_limits_f(plot, freq_f, &(pol->spectral_plot.max_dx[i*useful_bins]), useful_bins, 1, 1, 1);
		draw_grid(p, plot, 0, 0);
		draw_points_f(p, plot, COLOR(255,0,0), freq_f, &(pol->spectral_plot.max_dx[i*useful_bins]), useful_bins, 1, 1);
		RGBPic_dump_png(s, p);
		}
	snprintf(s,19999,"%s_max_dx_band_%d.dat",pol->name, i);
	dump_floats(s, &(pol->spectral_plot.max_dx[i*useful_bins]), useful_bins, 1);
	
	if(lines_list[0]>=0){
		snprintf(s,19999,"%s_max_mask_ratio_band_%d.png",pol->name, i);
		if(clear_name_png(s)){
			adjust_plot_limits_f(plot, freq_f, &(pol->spectral_plot.max_mask_ratio[i*useful_bins]), useful_bins, 1, 1, 1);
			draw_grid(p, plot, 0, 0);
			draw_points_f(p, plot, COLOR(255,0,0), freq_f, &(pol->spectral_plot.max_mask_ratio[i*useful_bins]), useful_bins, 1, 1);
			RGBPic_dump_png(s, p);
			}
		snprintf(s,19999,"%s_max_mask_ratio_band_%d.dat",pol->name, i);
		dump_floats(s, &(pol->spectral_plot.max_mask_ratio[i*useful_bins]), useful_bins, 1);
		}
	max_ratio=pol->spectral_plot.max_mask_ratio[i*useful_bins];
	for(k=1;k<useful_bins;k++)
		if(max_ratio<pol->spectral_plot.max_mask_ratio[i*useful_bins+k]){
			max_ratio=pol->spectral_plot.max_mask_ratio[i*useful_bins+k];
			}
	fprintf(LOG,"max_ratio: %d %s %f\n", i, pol->name, max_ratio);
        /* old 
	fprintf(LOG, "max_band: %ld %s %g\n", i, pol->name, max_band[i]);
	fprintf(LOG, "masked_max_band: %ld %s %g\n", i, pol->name, masked_max_band[i]);
	*/
	}

for(i=0;i<args_info.nbands_arg*useful_bins;i++){
	pol->spectral_plot.max_upper_limit[i]=sqrt(2.0*pol->spectral_plot.max_upper_limit[i]*upper_limit_comp)/(1800.0*16384.0);
	}

for(i=0;i<args_info.nbands_arg;i++){
	snprintf(s,19999,"%s_max_upper_strain_band_%d.png",pol->name, i);
	if(clear_name_png(s)){
		adjust_plot_limits_f(plot, freq_f, &(pol->spectral_plot.max_upper_limit[i*useful_bins]), useful_bins, 1, 1, 1);
		draw_grid(p, plot, 0, 0);
		draw_points_f(p, plot, COLOR(255,0,0), freq_f, &(pol->spectral_plot.max_upper_limit[i*useful_bins]), useful_bins, 1, 1);
		RGBPic_dump_png(s, p);
		}
	snprintf(s,19999,"%s_max_upper_strain_band_%d.dat",pol->name, i);
	dump_floats(s, &(pol->spectral_plot.max_upper_limit[i*useful_bins]), useful_bins, 1);
	}
	
snprintf(s,19999,"%s_max_upper_strain.png",pol->name);
if(clear_name_png(s)){
	plot_grid_f(p, fine_grid, pol->skymap.max_upper_limit, 1);
	RGBPic_dump_png(s, p);
	}
compute_histogram_f(hist, pol->skymap.max_upper_limit, fine_grid->band, fine_grid->npoints);
snprintf(s,19999,"hist_%s_max_upper_strain",pol->name);
print_histogram(LOG, hist, s);

snprintf(s,19999,"%s_max_lower_strain.png",pol->name);
if(clear_name_png(s)){
	plot_grid_f(p, fine_grid, pol->skymap.max_lower_limit, 1);
	RGBPic_dump_png(s, p);
	}

snprintf(s,19999,"%s_max_dx.png",pol->name);
if(clear_name_png(s)){
	plot_grid_f(p, fine_grid, pol->skymap.max_dx, 1);
	RGBPic_dump_png(s, p);
	}

snprintf(s,19999,"%s_M_map.png",pol->name);
if(clear_name_png(s)){
	plot_grid_f(p, fine_grid, pol->skymap.M_map, 1);
	RGBPic_dump_png(s, p);
	}

snprintf(s,19999,"%s_S_map.png",pol->name);
if(clear_name_png(s)){
	plot_grid_f(p, fine_grid, pol->skymap.S_map, 1);
	RGBPic_dump_png(s, p);
	}

fflush(LOG);
free_histogram(hist);
free_plot(plot);
free_RGBPic(p);
free(max_band);
free(masked_max_band);
free(max_band_arg);
free(masked_max_band_arg);
}

void output_unified_limits(void)
{
RGBPic *p;
PLOT *plot;
long i,k;
char s[20000];
SUM_TYPE *skymap_high_ul, *skymap_high_ul_freq;
SUM_TYPE *spectral_plot_high_ul;
float *freq_f;
SUM_TYPE max_high_ul, max_low_ul;
int max_high_ul_i, max_low_ul_i;
HISTOGRAM *hist;

freq_f=&(frequencies[side_cut]);

skymap_high_ul=do_alloc(fine_grid->npoints, sizeof(*skymap_high_ul));
skymap_high_ul_freq=do_alloc(fine_grid->npoints, sizeof(*skymap_high_ul_freq));
spectral_plot_high_ul=do_alloc(useful_bins*args_info.nbands_arg, sizeof(*spectral_plot_high_ul));
hist=new_histogram(args_info.hist_bins_arg, args_info.nbands_arg);

if(fine_grid->max_n_dec<800){
	p=make_RGBPic(fine_grid->max_n_ra*(800/fine_grid->max_n_dec)+140, fine_grid->max_n_dec*(800/fine_grid->max_n_dec));
	} else 
	p=make_RGBPic(fine_grid->max_n_ra+140, fine_grid->max_n_dec);	

plot=make_plot(p->width, p->height);

max_high_ul_i=-1;
max_low_ul_i=-1;
for(i=0;i<fine_grid->npoints;i++){
	if(fine_grid->band[i]<0){
		skymap_low_ul[i]=-1.0;
		skymap_high_ul[i]=-1.0;
		skymap_high_ul_freq[i]=-1.0;
		continue;
		}
	skymap_low_ul[i]=sqrt(2.0*skymap_low_ul[i]*upper_limit_comp)/(1800.0*16384.0);
	
	skymap_high_ul[i]=polarizations[0].skymap.max_upper_limit[i];
	skymap_high_ul_freq[i]=polarizations[0].skymap.freq_map[i];
	for(k=1;k<npolarizations;k++){	
		if(skymap_high_ul[i]<polarizations[k].skymap.max_upper_limit[i]){
			skymap_high_ul[i]=polarizations[k].skymap.max_upper_limit[i];
			skymap_high_ul_freq[i]=polarizations[k].skymap.freq_map[i];
			}
		}
	if(max_high_ul_i<0){
		max_high_ul_i=i;
		max_high_ul=skymap_high_ul[i];
		} else {
		if(max_high_ul<skymap_high_ul[i]){
			max_high_ul_i=i;
			max_high_ul=skymap_high_ul[i];
			}
		}
	if(max_low_ul_i<0){
		max_low_ul_i=i;
		max_low_ul=skymap_low_ul[i];
		} else {
		if(max_low_ul<skymap_low_ul[i]){
			max_low_ul_i=i;
			max_low_ul=skymap_low_ul[i];
			}
		}
	}
if(max_high_ul_i>=0){
	fprintf(LOG, "max_high_ul legend: RA DEC high_ul freq\n");
	fprintf(LOG, "max_high_ul: %f %f %g %f\n", 
		fine_grid->longitude[max_high_ul_i],
		fine_grid->latitude[max_high_ul_i],
		max_high_ul,
		skymap_high_ul_freq[max_high_ul_i]
		);
	}
if(max_low_ul_i>=0){
	fprintf(LOG, "max_low_ul legend: RA DEC low_ul freq\n");
	fprintf(LOG, "max_low_ul: %f %f %g %f\n", 
		fine_grid->longitude[max_low_ul_i],
		fine_grid->latitude[max_low_ul_i],
		max_low_ul,
		skymap_low_ul_freq[max_low_ul_i]
		);
	}

if(clear_name_png("low_ul.png")){
	plot_grid_f(p, fine_grid, skymap_low_ul, 1);
	RGBPic_dump_png("low_ul.png", p);
	}
dump_floats("low_ul.dat", skymap_low_ul, fine_grid->npoints, 1);
compute_histogram_f(hist, skymap_low_ul, fine_grid->band, fine_grid->npoints);
print_histogram(LOG, hist, "hist_low_ul");

if(clear_name_png("high_ul.png")){
	plot_grid_f(p, fine_grid, skymap_high_ul, 1);
	RGBPic_dump_png("high_ul.png", p);
	}
dump_floats("high_ul.dat", skymap_high_ul, fine_grid->npoints, 1);
compute_histogram_f(hist, skymap_high_ul, fine_grid->band, fine_grid->npoints);
print_histogram(LOG, hist, "hist_high_ul");


for(i=0;i<useful_bins*args_info.nbands_arg;i++){
	spectral_plot_low_ul[i]=sqrt(2.0*spectral_plot_low_ul[i]*upper_limit_comp)/(1800.0*16384.0);
	
	spectral_plot_high_ul[i]=polarizations[0].spectral_plot.max_upper_limit[i];
	for(k=1;k<npolarizations;k++){	
		if(spectral_plot_high_ul[i]<polarizations[k].spectral_plot.max_upper_limit[i])spectral_plot_high_ul[i]=polarizations[k].spectral_plot.max_upper_limit[i];
		}
	}

fprintf(LOG,"band upper limits: band UL freq\n");

for(i=0;i<args_info.nbands_arg;i++){

	max_high_ul_i=0;
	max_high_ul=spectral_plot_high_ul[i*useful_bins];
	for(k=1;k<useful_bins;k++){
		if(max_high_ul<spectral_plot_high_ul[i*useful_bins+k]){
			max_high_ul_i=k;
			max_high_ul=spectral_plot_high_ul[i*useful_bins+k];
			}
		}
	fprintf(LOG, "max_high_ul_band: %d %g %f\n", 
		i, max_high_ul, freq_f[max_high_ul_i]);
	
	max_low_ul_i=0;
	max_low_ul=spectral_plot_low_ul[i*useful_bins];
	for(k=1;k<useful_bins;k++){
		if(max_low_ul<spectral_plot_low_ul[i*useful_bins+k]){
			max_low_ul_i=k;
			max_low_ul=spectral_plot_low_ul[i*useful_bins+k];
			}
		}
	fprintf(LOG, "max_low_ul_band: %d %g %f\n", 
		i, max_low_ul, freq_f[max_low_ul_i]);
	
	snprintf(s,19999,"low_band_%d_ul.png", i);
	if(clear_name_png(s)){
		adjust_plot_limits_f(plot, freq_f, &(spectral_plot_low_ul[i*useful_bins]), useful_bins, 1, 1, 1);
		draw_grid(p, plot, 0, 0);
		draw_points_f(p, plot, COLOR(255,0,0), freq_f, &(spectral_plot_low_ul[i*useful_bins]), useful_bins, 1, 1);
		RGBPic_dump_png(s, p);
		}
	snprintf(s,19999,"low_band_%d_ul.dat", i);
	dump_floats(s, &(spectral_plot_low_ul[i*useful_bins]), useful_bins, 1);

	snprintf(s,19999,"high_band_%d_ul.png", i);
	if(clear_name_png(s)){
		adjust_plot_limits_f(plot, freq_f, &(spectral_plot_high_ul[i*useful_bins]), useful_bins, 1, 1, 1);
		draw_grid(p, plot, 0, 0);
		draw_points_f(p, plot, COLOR(255,0,0), freq_f, &(spectral_plot_high_ul[i*useful_bins]), useful_bins, 1, 1);
		RGBPic_dump_png(s, p);
		}
	snprintf(s,19999,"high_band_%d_ul.dat", i);
	dump_floats(s, &(spectral_plot_high_ul[i*useful_bins]), useful_bins, 1);
	}

free_histogram(hist);
free_plot(plot);
free_RGBPic(p);
free(skymap_high_ul);
}

void compute_mean(long pi)
{
SUM_TYPE a,c;
long i,k,m;
long offset;
for(k=0,offset=super_grid->first_map[pi];offset>=0;offset=super_grid->list_map[offset],k++){
	if(fine_grid->band[offset]<0)continue;
	
	for(m=0;m<npolarizations;m++){
		polarizations[m].skymap.max_sub_weight[offset]=0.0;
	
		for(i=0;i<useful_bins;i++){

			#ifdef WEIGHTED_SUM		
			if(polarizations[m].fine_grid_weight[i+k*useful_bins]>polarizations[m].skymap.max_sub_weight[offset]){
				polarizations[m].skymap.max_sub_weight[offset]=polarizations[m].fine_grid_weight[i+k*useful_bins];
				}
	
			c=(polarizations[m].skymap.total_weight[offset]-polarizations[m].fine_grid_weight[i+k*useful_bins]);
			#else
			c=(polarizations[m].skymap.total_count[offset]-polarizations[m].fine_grid_count[i+k*useful_bins]);
			#endif
			if(c>0){
				a=polarizations[m].fine_grid_sum[i+k*useful_bins];
				#ifdef COMPUTE_SIGMA
				b=polarizations[m].fine_grid_sq_sum[i+k*useful_bins];
				#endif
		
				#ifdef WEIGHTED_SUM		
				polarizations[m].fine_grid_sum[i+k*useful_bins]=a/c;
				#else
				polarizations[m].fine_grid_sum[i+k*useful_bins]=a/c;
				#endif
				#ifdef COMPUTE_SIGMA
				polarizations[m].fine_grid_sq_sum[i+k*useful_bins]=sqrt((b*polarizations[m].fine_grid_count[i+k*useful_bins]-a)/c);
				#endif
				}
			}
		}
	}
}

void compute_mean_no_lines(long pi)
{
SUM_TYPE a,c;
long i,k,offset,m;
for(k=0,offset=super_grid->first_map[pi];offset>=0;offset=super_grid->list_map[offset],k++){
	if(fine_grid->band[offset]<0)continue;

	for(m=0;m<npolarizations;m++){
		polarizations[m].skymap.max_sub_weight[offset]=0.0;
		
		for(i=0;i<useful_bins;i++){
		
			#ifdef WEIGHTED_SUM		
			c=(polarizations[m].skymap.total_weight[offset]);
			#else
			c=(polarizations[m].skymap.total_count[offset]);
			#endif
			if(c>0){
				a=polarizations[m].fine_grid_sum[i+k*useful_bins];
				#ifdef COMPUTE_SIGMA
				b=polarizations[m].fine_grid_sq_sum[i+k*useful_bins];
				#endif
		
				#ifdef WEIGHTED_SUM		
				polarizations[m].fine_grid_sum[i+k*useful_bins]=a/c;
				#else
				polarizations[m].fine_grid_sum[i+k*useful_bins]=a/c;
				#endif
				#ifdef COMPUTE_SIGMA
				polarizations[m].fine_grid_sq_sum[i+k*useful_bins]=sqrt((b*polarizations[m].fine_grid_count[i+k*useful_bins]-a)/c);
				#endif
				}
			}
		}
	}
}


void fine_grid_stage(void)
{
long pi,i,k,m;
double a,b;

normalizing_weight=exp(-M_LN10*TMedian);

stored_fine_bins=super_grid->max_npatch;

allocate_polarization_arrays();
	
low_ul=do_alloc(stored_fine_bins*useful_bins, sizeof(*low_ul));
low_ul_freq=do_alloc(stored_fine_bins*useful_bins, sizeof(*low_ul_freq));
skymap_low_ul=do_alloc(fine_grid->npoints, sizeof(*skymap_low_ul));
skymap_low_ul_freq=do_alloc(fine_grid->npoints, sizeof(*skymap_low_ul_freq));
spectral_plot_low_ul=do_alloc(useful_bins*args_info.nbands_arg, sizeof(*spectral_plot_low_ul));

for(i=0;i<fine_grid->npoints;i++){
	skymap_low_ul[i]=-1.0;
	skymap_low_ul_freq[i]=-1.0;	
	}
for(i=0;i<useful_bins*args_info.nbands_arg;i++){
	spectral_plot_low_ul[i]=-1.0;
	}
	
/* see comments above variables */
if(args_info.three_bins_arg){
	quantile2std=0.88;
	upper_limit_comp=3.0/(0.85);
	process_patch=process_patch3;
	fprintf(LOG,"mode: 3 bins\n");
	} else 
	{
	quantile2std=1.22;
	upper_limit_comp=1.0/(0.7*0.85);
	process_patch=process_patch1;
	fprintf(LOG,"mode: 1 bin\n");
	}

fprintf(stderr,"Main loop: %ld patches to process.\n", patch_grid->npoints);
for(pi=0;pi<patch_grid->npoints;pi++){
	if(patch_grid->band[pi]<0)continue;
	
	clear_accumulation_arrays();	

	/* process single patch */
	for(k=0;k<nsegments;k++){
		a=expTMedians[k];
		
		for(m=0;m<npolarizations;m++){
			b=polarizations[m].patch_CutOff[pi];
			/* process polarization */
			if(!do_CutOff || (b*a*AM_response(k, patch_grid, pi, polarizations[m].AM_coeffs)<4))
				process_patch(&(polarizations[m]), pi,k,b*sqrt(a));
			}
		}
		
	/* compute means */
	if(lines_list[0]<0)compute_mean_no_lines(pi);
		else compute_mean(pi);
		
	for(i=0;i<stored_fine_bins*useful_bins;i++){
		low_ul[i]=1.0e23; /* Sufficiently large number, 
		                     even for SFTs done with 
				     make_fake_data */
		}
	/* compute upper limits */
	for(i=0;i<npolarizations;i++){
		make_limits(&(polarizations[i]), pi);
		}
	make_unified_limits(pi);

	if(! (pi % 100))fprintf(stderr,"%ld ",pi);
	}
fprintf(stderr,"%ld ",pi);

fprintf(LOG,"Maximum bin shift: %ld\n", max_shift);
fprintf(LOG,"Minimum bin shift: %ld\n", min_shift);
fflush(LOG);

fprintf(stderr,"\nMaximum bin shift is %ld\n", max_shift);
fprintf(stderr,"Minimum bin shift is %ld\n", min_shift);

for(i=0;i<npolarizations;i++){
	output_limits(&(polarizations[i]));
	}
	
output_unified_limits();
}

