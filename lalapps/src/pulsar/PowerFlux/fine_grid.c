#include <stdio.h>
#include <stdlib.h>
/* We need this define to get NAN values */
#define __USE_ISOC99
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
extern int nlinear_polarizations, ntotal_polarizations;

extern float *power;
extern float *TMedians,*expTMedians;
extern float TMedian;

extern struct gengetopt_args_info args_info;

extern int nsegments, nbins, first_bin, side_cut, useful_bins;

INT64 *gps=NULL;
INT64 spindown_start;

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

extern int subinstance;
extern char *subinstance_name;

SUM_TYPE normalizing_weight;

int stored_fine_bins=0;

int max_shift=0, min_shift=0;


SUM_TYPE  *circ_ul, *circ_ul_freq;  /* lower (over polarizations) accumulation limits */
SUM_TYPE *skymap_circ_ul, *skymap_circ_ul_freq; /* skymaps */
SUM_TYPE *spectral_plot_circ_ul; /* spectral plots */

/* the value 1.22 is obtained by S.R script. It is valid for single-bin processing
   the value 0.88 is obtained by the same script. It is valid for 3 averaged neighboring bins */

float quantile2std=1.22;

float upper_limit_comp;
float lower_limit_comp;

/* single bin version */

void (*process_patch)(POLARIZATION *pol, int pi, int k, float CutOff);

static void process_patch1(POLARIZATION *pol, int pi, int k, float CutOff)
{
int i,kk,b,b0,b1,n,offset;
int bin_shift;
SUM_TYPE mod;
SUM_TYPE a,w,w2;
SUM_TYPE *sum,*sq_sum;
float *p;
float doppler;
float f_plus, f_cross, beta1, beta2;

CutOff=2*CutOff; /* weighted sum can benefit from more SFTs */


for(i=0,kk=super_grid->first_map[pi];kk>=0;kk=super_grid->list_map[kk],i++)
		{

		if(fine_grid->band[kk]<0)continue;
		
		/* Get amplitude response */
		f_plus=F_plus(k, fine_grid, kk, pol->AM_coeffs);
		f_cross=F_plus(k, fine_grid, kk, pol->conjugate->AM_coeffs);


		mod=1.0/(pol->plus_factor*f_plus*f_plus+pol->cross_factor*f_cross*f_cross); 

		if(do_CutOff && (mod>CutOff))continue;

		/* this assumes that both relfreq and spindown are small enough that the bin number
		   down not change much within the band - usually true */
		doppler=fine_grid->e[0][kk]*det_velocity[3*k+0]+
			fine_grid->e[1][kk]*det_velocity[3*k+1]+
			fine_grid->e[2][kk]*det_velocity[3*k+2];

		bin_shift=-rint((first_bin+nbins*0.5)*doppler+1800.0*spindown*(gps[k]-spindown_start));

		if(bin_shift>max_shift)max_shift=bin_shift;
		if(bin_shift<min_shift)min_shift=bin_shift;
		b0=side_cut-bin_shift;
		b1=(nbins-side_cut)-bin_shift;
		if((b0<0)||(b1>nbins)){
			fprintf(stderr,"Working frequency range obscured by bin_shift shift: bin_shift=%d kk=%d i=%d pi=%d\n",
				bin_shift, kk, i, pi);
			exit(-1);
			}
		
		if(args_info.compute_betas_arg){
		        beta1=f_cross*f_plus*mod;	
			beta2=(-pol->cross_factor*f_plus*f_plus+pol->plus_factor*f_cross*f_cross)*mod;
			}

		/* prime array pointers */
		#ifdef WEIGHTED_SUM
		w2=expTMedians[k]/mod;
		w=w2/mod;
		pol->skymap.total_weight[kk]+=w;
	        
		if(args_info.compute_betas_arg){
			pol->skymap.beta1[kk]+=w*beta1;
			pol->skymap.beta2[kk]+=w*beta2;
			}
		#else
		pol->skymap.total_count[kk]++;

		if(args_info.compute_betas_arg){
			pol->skymap.beta1[kk]+=beta1;
			pol->skymap.beta2[kk]+=beta2;
			}
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

static void process_patch3(POLARIZATION *pol, int pi, int k, float CutOff)
{
int i,kk,b,b0,b1,n,offset;
int bin_shift;
SUM_TYPE mod;
SUM_TYPE a,w,w2;
SUM_TYPE *sum,*sq_sum;
float *p;
float doppler;
float f_plus, f_cross, beta1, beta2;

CutOff=2*CutOff/3.0; /* weighted sum can benefit from more SFTs */


for(i=0,kk=super_grid->first_map[pi];kk>=0;kk=super_grid->list_map[kk],i++)
		{
		if(fine_grid->band[kk]<0)continue;

		/* Get amplitude response */
		f_plus=F_plus(k, fine_grid, kk, pol->AM_coeffs);
		f_cross=F_plus(k, fine_grid, kk, pol->conjugate->AM_coeffs);


		mod=1.0/(pol->plus_factor*f_plus*f_plus+pol->cross_factor*f_cross*f_cross); 

		if(do_CutOff && (mod>CutOff))continue;

		
		/* this assumes that both relfreq and spindown are small enough that the bin number
		   down not change much within the band - usually true */
		doppler=fine_grid->e[0][kk]*det_velocity[3*k+0]+
			fine_grid->e[1][kk]*det_velocity[3*k+1]+
			fine_grid->e[2][kk]*det_velocity[3*k+2];

		bin_shift=-rint((first_bin+nbins*0.5)*doppler+1800.0*spindown*(gps[k]-spindown_start));

		if(bin_shift>max_shift)max_shift=bin_shift;
		if(bin_shift<min_shift)min_shift=bin_shift;
		b0=side_cut-bin_shift;
		b1=(nbins-side_cut)-bin_shift;
		if((b0<1)||(b1>(nbins-1))){
			fprintf(stderr,"Working frequency range obscured by bin_shift shift: bin_shift=%d kk=%d i=%d pi=%d\n",
				bin_shift, kk, i, pi);
			exit(-1);
			}
		
		if(args_info.compute_betas_arg){
		        beta1=f_cross*f_plus*mod;	
			beta2=(-pol->cross_factor*f_plus*f_plus+pol->plus_factor*f_cross*f_cross)*mod;
			}

		/* prime array pointers */
		#ifdef WEIGHTED_SUM
		w2=expTMedians[k]/mod;
		w=w2/(3.0*mod);
		pol->skymap.total_weight[kk]+=w;

		if(args_info.compute_betas_arg){
			pol->skymap.beta1[kk]+=w*beta1;
			pol->skymap.beta2[kk]+=w*beta2;
			}
		#else
		pol->skymap.total_count[kk]++;

		if(args_info.compute_betas_arg){
			pol->skymap.beta1[kk]+=beta1;
			pol->skymap.beta2[kk]+=beta2;
			}
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

void make_limits(POLARIZATION *pol, int pi)
{
SUM_TYPE M,Q80,Q20,S,dx;
SUM_TYPE a,b,c;
SUM_TYPE *tmp=NULL;
int i,j,k,offset,band;
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
		snprintf(s, 20000, "points/%s%s_%d.png", subinstance_name, pol->name, offset);
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
			
		snprintf(s, 20000, "points/%s%s_%d.dat", subinstance_name, pol->name, offset);
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
			
		/* circ_ul describes limit on circularly polarized signals */
		/* this formula should only use linear polarizations */
		if(args_info.compute_betas_arg && (pol->cross_factor==0)){
			a*=1.0/(1.0+pol->skymap.beta2[offset]);
			if(a<circ_ul[i*useful_bins+k]){
				circ_ul[i*useful_bins+k]=a;
				circ_ul_freq[i*useful_bins+k]=(first_bin+side_cut+k)/1800.0;
				}
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

void make_unified_limits(int pi)
{
int i, offset;
int band, k;
SUM_TYPE a,b,c;
for(i=0,offset=super_grid->first_map[pi];offset>=0;offset=super_grid->list_map[offset],i++){

        band=fine_grid->band[offset];
	if(band<0)continue;

	for(k=0;k<useful_bins;k++){
		a=circ_ul[i*useful_bins+k];
		if(a>spectral_plot_circ_ul[k+band*useful_bins]){
			spectral_plot_circ_ul[k+band*useful_bins]=a;
			}
		if(a>skymap_circ_ul[offset]){
			skymap_circ_ul[offset]=a;
			skymap_circ_ul_freq[offset]=circ_ul_freq[i*useful_bins+k];
			}
		}
	}
}

void output_limits(POLARIZATION *pol)
{
char s[20000];
RGBPic *p;
PLOT *plot;
int i, max_dx_i, largest_i, masked, k;
SUM_TYPE max_dx, largest;
float *max_band, *masked_max_band;
int *max_band_arg, *masked_max_band_arg;
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

#define OUTPUT_SKYMAP(format, field)	{\
	snprintf(s,19999, "%s" format ".png", subinstance_name, pol->name); \
	if(clear_name_png(s)){ \
		plot_grid_f(p, fine_grid, pol->skymap.field, 1); \
		RGBPic_dump_png(s, p); \
		} \
	snprintf(s,19999, "%s" format ".dat", subinstance_name, pol->name); \
	dump_floats(s, pol->skymap.field, fine_grid->npoints, 1); \
	}

OUTPUT_SKYMAP("%s_weight", total_weight);

if(args_info.compute_betas_arg){
	OUTPUT_SKYMAP("%s_beta1", beta1);
	OUTPUT_SKYMAP("%s_beta2", beta2);
	}

OUTPUT_SKYMAP("%s_cor1", cor1);
OUTPUT_SKYMAP("%s_cor2", cor2);

if(args_info.ks_test_arg){
	OUTPUT_SKYMAP("%s_ks_test", ks_test);

	compute_histogram_f(hist, pol->skymap.ks_test, fine_grid->band, fine_grid->npoints);
	snprintf(s,19999,"%shist_%s_ks_test", subinstance_name, pol->name);
	print_histogram(LOG, hist, s);
	
	OUTPUT_SKYMAP("%s_ks_count", ks_count);
	}
	
OUTPUT_SKYMAP("%s_max_upper_limit", max_upper_limit);
OUTPUT_SKYMAP("%s_max_lower_limit", max_lower_limit);
OUTPUT_SKYMAP("%s_arg_freq", freq_map);

snprintf(s,19999,"%s%s_max_dx.dat", subinstance_name, pol->name);
dump_floats(s,pol->skymap.max_dx,fine_grid->npoints,1);

snprintf(s,19999,"%s%s_S_map.dat", subinstance_name, pol->name);
dump_floats(s,pol->skymap.S_map,fine_grid->npoints,1);
	

for(i=0;i<fine_grid->npoints;i++){
	if(fine_grid->band[i]<0){
		pol->skymap.max_upper_limit[i]=-1.0;
		pol->skymap.max_lower_limit[i]=-1.0;
		continue;
		}
	pol->skymap.max_upper_limit[i]=sqrt(pol->skymap.max_upper_limit[i])*upper_limit_comp;
	pol->skymap.max_lower_limit[i]=sqrt(pol->skymap.max_lower_limit[i])*lower_limit_comp;
	}

/* output interesting points around fake injection */
if(fake_injection){
	double ds, best_ds;
	int best_i=-1;
	fprintf(LOG,"Interesting points: index longitude latitude pol max_dx upper_strain lower_strain freq beta1 beta2\n");
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
		
		   	fprintf(LOG, "%d %f %f %s %f %g %g %f %f %f\n",
				i,
				fine_grid->longitude[i], fine_grid->latitude[i], 
				pol->name, pol->skymap.max_dx[i], 
				pol->skymap.max_upper_limit[i], pol->skymap.max_lower_limit[i],
				pol->skymap.freq_map[i], 
				args_info.compute_betas_arg?pol->skymap.beta1[i]:NAN,
				args_info.compute_betas_arg?pol->skymap.beta2[i]:NAN);

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
	fprintf(LOG, "i_closest: %d %f %f %s %f %g %g %f %f %f\n", best_i, fine_grid->longitude[best_i], fine_grid->latitude[best_i], 
		pol->name, pol->skymap.max_dx[best_i], 
		pol->skymap.max_upper_limit[best_i], pol->skymap.max_lower_limit[best_i], pol->skymap.freq_map[best_i],
		args_info.compute_betas_arg?pol->skymap.beta1[best_i]:NAN, args_info.compute_betas_arg?pol->skymap.beta2[best_i]:NAN);
	if(largest_i>=0)
	fprintf(LOG, "i_largest: %d %f %f %s %f %g %g %f %f %f\n", largest_i, fine_grid->longitude[largest_i], fine_grid->latitude[largest_i], 
		pol->name, pol->skymap.max_dx[largest_i], 
		pol->skymap.max_upper_limit[largest_i], pol->skymap.max_lower_limit[largest_i], pol->skymap.freq_map[largest_i],
		args_info.compute_betas_arg?pol->skymap.beta1[largest_i]:NAN, 
		args_info.compute_betas_arg?pol->skymap.beta2[largest_i]:NAN);
	}



snprintf(s,19999,"%s%s_max_strain.dat", subinstance_name, pol->name);
dump_floats(s,pol->skymap.max_upper_limit,fine_grid->npoints,1);

max_dx=0.0;
max_dx_i=0;
largest=0.0;
largest_i=0;
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

	if(pol->skymap.max_upper_limit[i]>largest){
		largest=pol->skymap.max_upper_limit[i];
		largest_i=i;
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
fprintf(LOG, "masked: %s %d\n", pol->name, masked);
fprintf(LOG, "strongest signal: longitude latitude pol max_dx upper_strain lower_strain freq beta1 beta2\n");	
fprintf(LOG, "max_dx: %f %f %s %f %g %g %f %f %f\n",fine_grid->longitude[max_dx_i], fine_grid->latitude[max_dx_i], 
				pol->name, pol->skymap.max_dx[max_dx_i], 
				pol->skymap.max_upper_limit[max_dx_i], 
				pol->skymap.max_lower_limit[max_dx_i], 
				pol->skymap.freq_map[max_dx_i],
				args_info.compute_betas_arg?pol->skymap.beta1[max_dx_i]:NAN, 
				args_info.compute_betas_arg?pol->skymap.beta2[max_dx_i]:NAN);

fprintf(LOG, "largest signal: longitude latitude pol max_dx upper_strain lower_strain freq beta1 beta2\n");	
fprintf(LOG, "largest: %f %f %s %f %g %g %f %f %f\n",fine_grid->longitude[largest_i], fine_grid->latitude[largest_i], 
				pol->name, pol->skymap.max_dx[largest_i], 
				pol->skymap.max_upper_limit[largest_i], 
				pol->skymap.max_lower_limit[largest_i], 
				pol->skymap.freq_map[largest_i],
				args_info.compute_betas_arg?pol->skymap.beta1[largest_i]:NAN, 
				args_info.compute_betas_arg?pol->skymap.beta2[largest_i]:NAN);

fprintf(LOG, "max/masked band format: band_num longitude latitude pol max_dx upper_strain freq beta1 beta2\n");
for(i=0;i<args_info.nbands_arg;i++){
	if(max_band_arg[i]<0){
		fprintf(LOG, "max_band: %d NAN NAN %s NAN NAN NAN NAN NAN\n", i, pol->name); 
		fprintf(LOG, "masked_max_band: %d NAN NAN %s NAN NAN NAN NAN NAN\n", i, pol->name);
		fprintf(LOG,"max_ratio: %d %s NAN\n", i, pol->name);
		continue;
		}

	fprintf(LOG, "max_band: %d %f %f %s %f %g %f %f %f\n", i, fine_grid->longitude[max_band_arg[i]], fine_grid->latitude[max_band_arg[i]], 
				pol->name, pol->skymap.max_dx[max_band_arg[i]], 
				max_band[i], 
				pol->skymap.freq_map[max_band_arg[i]],
				args_info.compute_betas_arg?pol->skymap.beta1[max_band_arg[i]]:NAN, 
				args_info.compute_betas_arg?pol->skymap.beta2[max_band_arg[i]]:NAN);

	fprintf(LOG, "masked_max_band: %d %f %f %s %f %g %f %f %f\n", i, fine_grid->longitude[masked_max_band_arg[i]], fine_grid->latitude[masked_max_band_arg[i]], 
				pol->name, pol->skymap.max_dx[masked_max_band_arg[i]], 
				masked_max_band[i], 
				pol->skymap.freq_map[masked_max_band_arg[i]],
				args_info.compute_betas_arg?pol->skymap.beta1[masked_max_band_arg[i]]:NAN, 
				args_info.compute_betas_arg?pol->skymap.beta2[masked_max_band_arg[i]]:NAN);

	snprintf(s,19999,"%s%s_max_upper_limit_band_%d.png", subinstance_name, pol->name, i);
	if(clear_name_png(s)){
		adjust_plot_limits_f(plot, freq_f, &(pol->spectral_plot.max_upper_limit[i*useful_bins]), useful_bins, 1, 1, 1);
		draw_grid(p, plot, 0, 0);
		draw_points_f(p, plot, COLOR(255,0,0), freq_f, &(pol->spectral_plot.max_upper_limit[i*useful_bins]), useful_bins, 1, 1);
		RGBPic_dump_png(s, p);
		}
	snprintf(s,19999,"%s%s_max_upper_limit_band_%d.dat", subinstance_name, pol->name, i);
	dump_floats(s, &(pol->spectral_plot.max_upper_limit[i*useful_bins]), useful_bins, 1);
	
	snprintf(s,19999,"%s%s_max_dx_band_%d.png", subinstance_name, pol->name, i);
	if(clear_name_png(s)){
		adjust_plot_limits_f(plot, freq_f, &(pol->spectral_plot.max_dx[i*useful_bins]), useful_bins, 1, 1, 1);
		draw_grid(p, plot, 0, 0);
		draw_points_f(p, plot, COLOR(255,0,0), freq_f, &(pol->spectral_plot.max_dx[i*useful_bins]), useful_bins, 1, 1);
		RGBPic_dump_png(s, p);
		}
	snprintf(s,19999,"%s%s_max_dx_band_%d.dat", subinstance_name, pol->name, i);
	dump_floats(s, &(pol->spectral_plot.max_dx[i*useful_bins]), useful_bins, 1);
	
	if(lines_list[0]>=0){
		snprintf(s,19999,"%s%s_max_mask_ratio_band_%d.png", subinstance_name, pol->name, i);
		if(clear_name_png(s)){
			adjust_plot_limits_f(plot, freq_f, &(pol->spectral_plot.max_mask_ratio[i*useful_bins]), useful_bins, 1, 1, 1);
			draw_grid(p, plot, 0, 0);
			draw_points_f(p, plot, COLOR(255,0,0), freq_f, &(pol->spectral_plot.max_mask_ratio[i*useful_bins]), useful_bins, 1, 1);
			RGBPic_dump_png(s, p);
			}
		snprintf(s,19999,"%s%s_max_mask_ratio_band_%d.dat", subinstance_name, pol->name, i);
		dump_floats(s, &(pol->spectral_plot.max_mask_ratio[i*useful_bins]), useful_bins, 1);
		}
	max_ratio=pol->spectral_plot.max_mask_ratio[i*useful_bins];
	for(k=1;k<useful_bins;k++)
		if(max_ratio<pol->spectral_plot.max_mask_ratio[i*useful_bins+k]){
			max_ratio=pol->spectral_plot.max_mask_ratio[i*useful_bins+k];
			}
	fprintf(LOG,"max_ratio: %d %s %f\n", i, pol->name, max_ratio);
        /* old 
	fprintf(LOG, "max_band: %d %s %g\n", i, pol->name, max_band[i]);
	fprintf(LOG, "masked_max_band: %d %s %g\n", i, pol->name, masked_max_band[i]);
	*/
	}

for(i=0;i<args_info.nbands_arg*useful_bins;i++){
	pol->spectral_plot.max_upper_limit[i]=sqrt(pol->spectral_plot.max_upper_limit[i])*upper_limit_comp;
	}

for(i=0;i<args_info.nbands_arg;i++){
	snprintf(s,19999,"%s%s_max_upper_strain_band_%d.png", subinstance_name, pol->name, i);
	if(clear_name_png(s)){
		adjust_plot_limits_f(plot, freq_f, &(pol->spectral_plot.max_upper_limit[i*useful_bins]), useful_bins, 1, 1, 1);
		draw_grid(p, plot, 0, 0);
		draw_points_f(p, plot, COLOR(255,0,0), freq_f, &(pol->spectral_plot.max_upper_limit[i*useful_bins]), useful_bins, 1, 1);
		RGBPic_dump_png(s, p);
		}
	snprintf(s,19999,"%s%s_max_upper_strain_band_%d.dat", subinstance_name, pol->name, i);
	dump_floats(s, &(pol->spectral_plot.max_upper_limit[i*useful_bins]), useful_bins, 1);
	}
	
snprintf(s,19999,"%s%s_max_upper_strain.png", subinstance_name, pol->name);
if(clear_name_png(s)){
	plot_grid_f(p, fine_grid, pol->skymap.max_upper_limit, 1);
	RGBPic_dump_png(s, p);
	}
compute_histogram_f(hist, pol->skymap.max_upper_limit, fine_grid->band, fine_grid->npoints);
snprintf(s,19999,"%shist_%s_max_upper_strain", subinstance_name, pol->name);
print_histogram(LOG, hist, s);

snprintf(s,19999,"%s%s_max_lower_strain.png", subinstance_name, pol->name);
if(clear_name_png(s)){
	plot_grid_f(p, fine_grid, pol->skymap.max_lower_limit, 1);
	RGBPic_dump_png(s, p);
	}

snprintf(s,19999,"%s%s_max_dx.png", subinstance_name, pol->name);
if(clear_name_png(s)){
	plot_grid_f(p, fine_grid, pol->skymap.max_dx, 1);
	RGBPic_dump_png(s, p);
	}

snprintf(s,19999,"%s%s_M_map.png", subinstance_name, pol->name);
if(clear_name_png(s)){
	plot_grid_f(p, fine_grid, pol->skymap.M_map, 1);
	RGBPic_dump_png(s, p);
	}

snprintf(s,19999,"%s%s_S_map.png", subinstance_name, pol->name);
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
int i,k;
char s[20000];
SUM_TYPE *skymap_high_ul, *skymap_high_ul_freq;
SUM_TYPE *spectral_plot_high_ul;
float *freq_f;
SUM_TYPE max_high_ul, max_circ_ul;
int max_high_ul_i, max_circ_ul_i;
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
max_circ_ul_i=-1;
for(i=0;i<fine_grid->npoints;i++){
	if(fine_grid->band[i]<0){
		skymap_circ_ul[i]=-1.0;
		skymap_high_ul[i]=-1.0;
		skymap_high_ul_freq[i]=-1.0;
		continue;
		}
	skymap_circ_ul[i]=sqrt(skymap_circ_ul[i])*upper_limit_comp;
	
	skymap_high_ul[i]=polarizations[0].skymap.max_upper_limit[i];
	skymap_high_ul_freq[i]=polarizations[0].skymap.freq_map[i];
	for(k=1;k<ntotal_polarizations;k++){	
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
	if(max_circ_ul_i<0){
		max_circ_ul_i=i;
		max_circ_ul=skymap_circ_ul[i];
		} else {
		if(max_circ_ul<skymap_circ_ul[i]){
			max_circ_ul_i=i;
			max_circ_ul=skymap_circ_ul[i];
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
if(max_circ_ul_i>=0){
	fprintf(LOG, "max_circ_ul legend: RA DEC circ_ul freq\n");
	fprintf(LOG, "max_circ_ul: %f %f %g %f\n", 
		fine_grid->longitude[max_circ_ul_i],
		fine_grid->latitude[max_circ_ul_i],
		max_circ_ul,
		skymap_circ_ul_freq[max_circ_ul_i]
		);
	}

if(args_info.compute_betas_arg){
	if(clear_name_png("circ_ul.png")){
		plot_grid_f(p, fine_grid, skymap_circ_ul, 1);
		RGBPic_dump_png("circ_ul.png", p);
		}
	dump_floats("circ_ul.dat", skymap_circ_ul, fine_grid->npoints, 1);
	compute_histogram_f(hist, skymap_circ_ul, fine_grid->band, fine_grid->npoints);
	print_histogram(LOG, hist, "hist_circ_ul");
	}

if(clear_name_png("high_ul.png")){
	plot_grid_f(p, fine_grid, skymap_high_ul, 1);
	RGBPic_dump_png("high_ul.png", p);
	}
dump_floats("high_ul.dat", skymap_high_ul, fine_grid->npoints, 1);
compute_histogram_f(hist, skymap_high_ul, fine_grid->band, fine_grid->npoints);
print_histogram(LOG, hist, "hist_high_ul");


for(i=0;i<useful_bins*args_info.nbands_arg;i++){
	spectral_plot_circ_ul[i]=sqrt(spectral_plot_circ_ul[i])*upper_limit_comp;
	
	spectral_plot_high_ul[i]=polarizations[0].spectral_plot.max_upper_limit[i];
	for(k=1;k<ntotal_polarizations;k++){	
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
	
	max_circ_ul_i=0;
	max_circ_ul=spectral_plot_circ_ul[i*useful_bins];
	for(k=1;k<useful_bins;k++){
		if(max_circ_ul<spectral_plot_circ_ul[i*useful_bins+k]){
			max_circ_ul_i=k;
			max_circ_ul=spectral_plot_circ_ul[i*useful_bins+k];
			}
		}
	fprintf(LOG, "max_circ_ul_band: %d %g %f\n", 
		i, max_circ_ul, freq_f[max_circ_ul_i]);
	
	snprintf(s,19999,"%slow_band_%d_ul.png", subinstance_name, i);
	if(clear_name_png(s)){
		adjust_plot_limits_f(plot, freq_f, &(spectral_plot_circ_ul[i*useful_bins]), useful_bins, 1, 1, 1);
		draw_grid(p, plot, 0, 0);
		draw_points_f(p, plot, COLOR(255,0,0), freq_f, &(spectral_plot_circ_ul[i*useful_bins]), useful_bins, 1, 1);
		RGBPic_dump_png(s, p);
		}
	snprintf(s,19999,"%slow_band_%d_ul.dat", subinstance_name, i);
	dump_floats(s, &(spectral_plot_circ_ul[i*useful_bins]), useful_bins, 1);

	snprintf(s,19999,"%shigh_band_%d_ul.png", subinstance_name, i);
	if(clear_name_png(s)){
		adjust_plot_limits_f(plot, freq_f, &(spectral_plot_high_ul[i*useful_bins]), useful_bins, 1, 1, 1);
		draw_grid(p, plot, 0, 0);
		draw_points_f(p, plot, COLOR(255,0,0), freq_f, &(spectral_plot_high_ul[i*useful_bins]), useful_bins, 1, 1);
		RGBPic_dump_png(s, p);
		}
	snprintf(s,19999,"%shigh_band_%d_ul.dat", subinstance_name, i);
	dump_floats(s, &(spectral_plot_high_ul[i*useful_bins]), useful_bins, 1);
	}

free_histogram(hist);
free_plot(plot);
free_RGBPic(p);
free(skymap_high_ul);
}

void compute_mean(int pi)
{
SUM_TYPE a,c;
int i,k,m;
int offset;
for(k=0,offset=super_grid->first_map[pi];offset>=0;offset=super_grid->list_map[offset],k++){
	if(fine_grid->band[offset]<0)continue;
	
	for(m=0;m<ntotal_polarizations;m++){
		polarizations[m].skymap.max_sub_weight[offset]=0.0;

		#ifdef WEIGHTED_SUM		
		c=(polarizations[m].skymap.total_weight[offset]);
		#else
		c=(polarizations[m].skymap.total_count[offset]);
		#endif

		if(c>0){	
			if(args_info.compute_betas_arg){
			    polarizations[m].skymap.beta1[offset]/=c;
			    polarizations[m].skymap.beta2[offset]/=c;
			    }

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
		
				polarizations[m].fine_grid_sum[i+k*useful_bins]=a/c;

				#ifdef COMPUTE_SIGMA
				polarizations[m].fine_grid_sq_sum[i+k*useful_bins]=sqrt((b*polarizations[m].fine_grid_count[i+k*useful_bins]-a)/c);
				#endif
				}
			}
		   }
		}
	}
}

void compute_mean_no_lines(int pi)
{
SUM_TYPE a,c;
int i,k,offset,m;
for(k=0,offset=super_grid->first_map[pi];offset>=0;offset=super_grid->list_map[offset],k++){
	if(fine_grid->band[offset]<0)continue;

	for(m=0;m<ntotal_polarizations;m++){
		polarizations[m].skymap.max_sub_weight[offset]=0.0;
		
		#ifdef WEIGHTED_SUM		
		c=(polarizations[m].skymap.total_weight[offset]);
		#else
		c=(polarizations[m].skymap.total_count[offset]);
		#endif


		if(c>0){
			if(args_info.compute_betas_arg){
				polarizations[m].skymap.beta1[offset]/=c;
				polarizations[m].skymap.beta2[offset]/=c;
				}

			for(i=0;i<useful_bins;i++){
				a=polarizations[m].fine_grid_sum[i+k*useful_bins];
				#ifdef COMPUTE_SIGMA
				b=polarizations[m].fine_grid_sq_sum[i+k*useful_bins];
				#endif
		
				polarizations[m].fine_grid_sum[i+k*useful_bins]=a/c;

				#ifdef COMPUTE_SIGMA
				polarizations[m].fine_grid_sq_sum[i+k*useful_bins]=sqrt((b*polarizations[m].fine_grid_count[i+k*useful_bins]-a)/c);
				#endif
				}
			}
		}
	}
}

void init_fine_grid_stage(void)
{

normalizing_weight=exp(-M_LN10*TMedian);

stored_fine_bins=super_grid->max_npatch;

allocate_polarization_arrays();
	
circ_ul=do_alloc(stored_fine_bins*useful_bins, sizeof(*circ_ul));
circ_ul_freq=do_alloc(stored_fine_bins*useful_bins, sizeof(*circ_ul_freq));
skymap_circ_ul=do_alloc(fine_grid->npoints, sizeof(*skymap_circ_ul));
skymap_circ_ul_freq=do_alloc(fine_grid->npoints, sizeof(*skymap_circ_ul_freq));
spectral_plot_circ_ul=do_alloc(useful_bins*args_info.nbands_arg, sizeof(*spectral_plot_circ_ul));

if(args_info.three_bins_arg){
	quantile2std=0.88;
	process_patch=process_patch3;
	fprintf(LOG,"mode: 3 bins\n");
	} else 
	{
	quantile2std=1.22;
	process_patch=process_patch1;
	fprintf(LOG,"mode: 1 bin\n");
	}

if(!strcasecmp("Hann", args_info.upper_limit_comp_arg)){
	if(args_info.three_bins_arg){
		/* 3 bins should contain the entire signal, regardless
		   of positioning */
		upper_limit_comp=sqrt(3.0);
		} else 
		{
		/* 0.85 is a ratio between amplitude of 
		   half-bin centered signal and bin centered signal
		   *amplitude*

		   */
		upper_limit_comp=1.0/0.85;
		}
	} else {
	upper_limit_comp=atof(args_info.upper_limit_comp_arg);
	}
fprintf(LOG, "upper limit compensation factor: %8f\n", upper_limit_comp);

	/* Extra factor to convert to amplitude from RMS power */
upper_limit_comp*=sqrt(2.0);
	/* Extra factor to convert to strain from raw SFT units */
upper_limit_comp/=(1800.0*16384.0);
	/* Extra factor to account for the fact that only half of SFT
	   coefficients is stored */
upper_limit_comp*=sqrt(2.0);

if(!strcasecmp("Hann", args_info.lower_limit_comp_arg)){
	if(args_info.three_bins_arg){
		lower_limit_comp=sqrt(3.0);
		} else 
		{
		lower_limit_comp=1.0;
		}
	} else {
	lower_limit_comp=atof(args_info.lower_limit_comp_arg);
	}
fprintf(LOG, "lower limit compensation factor: %8f\n", lower_limit_comp);

	/* Extra factor to  convert to amplitude from RMS power */
lower_limit_comp*=sqrt(2.0);
	/* Extra factor to convert to strain from raw SFT units */
lower_limit_comp/=(1800.0*16384.0);
	/* Extra factor to account for the fact that only half of SFT
	   coefficients is stored */
lower_limit_comp*=sqrt(2.0);

}

void fine_grid_stage(void)
{
int pi,i,k,m;
double a,b;

clear_polarization_arrays();

min_shift=0;
max_shift=0;

for(i=0;i<fine_grid->npoints;i++){
	skymap_circ_ul[i]=-1.0;
	skymap_circ_ul_freq[i]=-1.0;	
	}
for(i=0;i<useful_bins*args_info.nbands_arg;i++){
	spectral_plot_circ_ul[i]=-1.0;
	}
	
fprintf(stderr,"Main loop: %d patches to process.\n", patch_grid->npoints);
for(pi=0;pi<patch_grid->npoints;pi++){
	if(patch_grid->band[pi]<0)continue;
	
	clear_accumulation_arrays();	

	/* process single patch */
	for(k=0;k<nsegments;k++){
		a=expTMedians[k];
		
		for(m=0;m<ntotal_polarizations;m++){
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
		circ_ul[i]=1.0e23; /* Sufficiently large number, 
		                     even for SFTs done with 
				     make_fake_data */
		}
	/* compute upper limits */
	for(i=0;i<ntotal_polarizations;i++){
		make_limits(&(polarizations[i]), pi);
		}
	make_unified_limits(pi);

	if(! (pi % 100))fprintf(stderr,"%d ",pi);
	}
fprintf(stderr,"%d ",pi);

fprintf(LOG,"Maximum bin shift: %d\n", max_shift);
fprintf(LOG,"Minimum bin shift: %d\n", min_shift);
fflush(LOG);

fprintf(stderr,"\nMaximum bin shift is %d\n", max_shift);
fprintf(stderr,"Minimum bin shift is %d\n", min_shift);

for(i=0;i<ntotal_polarizations;i++){
	output_limits(&(polarizations[i]));
	}
	
output_unified_limits();
}

