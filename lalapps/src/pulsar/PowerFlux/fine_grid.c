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

extern INT64 *gps;

SUM_TYPE normalizing_weight;

long stored_fine_bins=0;

long max_shift=0, min_shift=0;



/* the value 1.22 is obtained by S.R script. It is valid for single-bin processing
   the value 0.88 is obtained by the same script. It is valid for 3 averaged neighboring bins */

float quantile2std=1.22;

/* We divide by 0.7 to compensate for spreading due to Hanning window */
/* We also divide by 0.85 to compensate for non bin centered signals */
/* for 3 averaged neighbouring bins that value should be 3.0/0.85 */
float upper_limit_comp=1/(0.7*0.85);

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
		pol->total_weight[kk]+=w;
		#else
		pol->total_count[kk]++;
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
		pol->total_weight[kk]+=w;
		#else
		pol->total_count[kk]++;
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

static float upper_limit95(float dx)
{
if(dx>=-0.1)return dx+1.96;
dx+=0.1;
return dx+sqrt(dx*dx+1.86);
}

static float lower_limit95(float dx)
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
long i,j,k,offset;

/* allocate on stack, for speed */
tmp=alloca(useful_bins*sizeof(*tmp));

for(i=0,offset=super_grid->first_map[pi];offset>=0;offset=super_grid->list_map[offset],i++){
	/* figure out offset to put results into */
		
	memcpy(tmp,pol->fine_grid_sum+i*useful_bins,useful_bins*sizeof(*tmp));
	/* compute correlations - for diagnostics only */
	a=0.0;
	b=0.0;
	for(j=0;j<useful_bins-2;j++){
		pol->cor1[offset]+=(tmp[j]*normalizing_weight)*(tmp[j+1]*normalizing_weight);
		pol->cor2[offset]+=(tmp[j]*normalizing_weight)*(tmp[j+2]*normalizing_weight);
		b+=(tmp[j]*normalizing_weight);
		a+=(tmp[j]*normalizing_weight)*(tmp[j]*normalizing_weight);
		}
	c=(b-tmp[0]-tmp[1]+tmp[useful_bins-2]+tmp[useful_bins-1]);
	pol->cor2[offset]=(pol->cor2[offset]-b*c/(useful_bins-2))/
		(sqrt((a-b*b/(useful_bins-2))*
		(a-tmp[0]*tmp[0]-tmp[1]*tmp[1]
			+tmp[useful_bins-2]*tmp[useful_bins-2]+tmp[useful_bins-1]*tmp[useful_bins-1]
			-c*c/(useful_bins-2))));
			
	b+=tmp[useful_bins-2];
	a+=tmp[useful_bins-2]*tmp[useful_bins-2];
	pol->cor1[offset]+=tmp[useful_bins-2]*tmp[useful_bins-1];
	c=b-tmp[0]+tmp[useful_bins-1];
	pol->cor1[offset]=(pol->cor1[offset]-b*c/(useful_bins-1))/
		(sqrt(
		(a-b*b/(useful_bins-1))*
		(a-tmp[0]*tmp[0]+tmp[useful_bins-1]*tmp[useful_bins-1]-c*c/(useful_bins-1))
		));
	
	/* sort to compute robust estimates */
	qsort(tmp,useful_bins,sizeof(*tmp),float_cmp);
	/* median */
	M=tmp[useful_bins/2];
	/* 0.8 quantile */
	Q80=tmp[(useful_bins*4)/5];
	/* 0.2 quantile */
	Q20=tmp[useful_bins/5];
	S=(Q80-Q20)/quantile2std;
	pol->M_map[offset]=M;
	pol->S_map[offset]=S;
	pol->max_upper_limit[offset]=0;
	for(k=0;k<useful_bins;k++){
		dx=(pol->fine_grid_sum[i*useful_bins+k]-M)/S;		
		a=upper_limit95(dx)*S;
		if(a>pol->max_upper_limit[offset]){
			pol->max_upper_limit[offset]=a;
			pol->freq_map[offset]=(first_bin+side_cut+k)/1800.0;
			}
		if(dx>pol->max_dx[offset]){
			pol->max_dx[offset]=dx;
			}
			
		a=lower_limit95(dx)*S;
		if(a>pol->max_lower_limit[offset]){
			pol->max_lower_limit[offset]=a;
			}
		}
	}

}

void output_limits(POLARIZATION *pol)
{
char s[20000];
RGBPic *p;
long i, max_dx_i, masked, k;
SUM_TYPE max_dx;
float *max_band, *masked_max_band;
long *max_band_arg, *masked_max_band_arg;

max_band=do_alloc(args_info.dec_bands_arg, sizeof(*max_band));
masked_max_band=do_alloc(args_info.dec_bands_arg, sizeof(*max_band));
max_band_arg=do_alloc(args_info.dec_bands_arg, sizeof(*max_band_arg));
masked_max_band_arg=do_alloc(args_info.dec_bands_arg, sizeof(*max_band_arg));

if(fine_grid->max_n_dec<800){
	p=make_RGBPic(fine_grid->max_n_ra*(800/fine_grid->max_n_dec)+140, fine_grid->max_n_dec*(800/fine_grid->max_n_dec));
	} else 
	p=make_RGBPic(fine_grid->max_n_ra+140, fine_grid->max_n_dec);	

snprintf(s,19999,"%s_weight.png",pol->name);
plot_grid_f(p, fine_grid, pol->total_weight, 1);
RGBPic_dump_png(s, p);

snprintf(s,19999,"%s_cor1.png",pol->name);
plot_grid_f(p, fine_grid, pol->cor1, 1);
RGBPic_dump_png(s, p);

snprintf(s,19999,"%s_cor2.png",pol->name);
plot_grid_f(p, fine_grid, pol->cor2, 1);
RGBPic_dump_png(s, p);

snprintf(s,19999,"%s_max_upper_limit.png",pol->name);
plot_grid_f(p, fine_grid, pol->max_upper_limit, 1);
RGBPic_dump_png(s, p);

snprintf(s,19999,"%s_max_lower_limit.png",pol->name);
plot_grid_f(p, fine_grid, pol->max_lower_limit, 1);
RGBPic_dump_png(s, p);

snprintf(s,19999,"%s_arg_freq.png",pol->name);
plot_grid_f(p, fine_grid, pol->freq_map, 1);
RGBPic_dump_png(s, p);

snprintf(s,19999,"%s_arg_freq.dat",pol->name);
dump_floats(s,pol->freq_map,fine_grid->npoints,1);

snprintf(s,19999,"%s_max_dx.dat",pol->name);
dump_floats(s,pol->max_dx,fine_grid->npoints,1);

snprintf(s,19999,"%s_S_map.dat",pol->name);
dump_floats(s,pol->S_map,fine_grid->npoints,1);
	

for(i=0;i<fine_grid->npoints;i++){
	pol->max_upper_limit[i]=sqrt(2.0*pol->max_upper_limit[i]*upper_limit_comp)/(1800.0*16384.0);
	/* lower limit is unchanged */
	pol->max_lower_limit[i]=sqrt(2.0*pol->max_lower_limit[i])/(1800.0*16384.0);
	}

/* output interesting points around fake injection */
if(fake_injection){
	float largest;
	double ds, best_ds;
	long best_i=-1, largest_i=-1;
	fprintf(LOG,"Interesting points: longitude latitude pol max_dx upper_strain lower_strain freq\n");
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
		
		   	fprintf(LOG, "%f %f %s %f %g %g %f\n",fine_grid->longitude[i], fine_grid->latitude[i], 
				pol->name, pol->max_dx[i], 
				pol->max_upper_limit[i], pol->max_lower_limit[i],
				pol->freq_map[i]);

			if(largest_i<0){
				largest=pol->max_upper_limit[i];
				largest_i=i;
				} else 
			if(largest<pol->max_upper_limit[i]){
				largest=pol->max_upper_limit[i];
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
		pol->name, pol->max_dx[best_i], 
		pol->max_upper_limit[best_i], pol->max_lower_limit[best_i], pol->freq_map[best_i]);
	if(largest_i>=0)
	fprintf(LOG, "largest: %f %f %s %f %g %g %f\n",fine_grid->longitude[largest_i], fine_grid->latitude[largest_i], 
		pol->name, pol->max_dx[largest_i], 
		pol->max_upper_limit[largest_i], pol->max_lower_limit[largest_i], pol->freq_map[largest_i]);
	}



snprintf(s,19999,"%s_max_strain.dat",pol->name);
dump_floats(s,pol->max_upper_limit,fine_grid->npoints,1);

max_dx=0.0;
max_dx_i=0;
masked=0;
for(i=0;i<args_info.dec_bands_arg;i++){
	max_band[i]=-1.0;
	masked_max_band[i]=-1.0;
	max_band_arg[i]=-1;
	masked_max_band_arg[i]=-1;
	}
	
for(i=0;i<fine_grid->npoints;i++){
	k=floor((0.5+fine_grid->latitude[i]/M_PI)*args_info.dec_bands_arg);
	if(k==args_info.dec_bands_arg)k=args_info.dec_bands_arg-1;
	if(k<0)k=0;

	if(pol->max_upper_limit[i]>max_band[k]){
		max_band[k]=pol->max_upper_limit[i];
		max_band_arg[k]=i;
		}

	if(pol->max_sub_weight[i]>=pol->total_weight[i]*(1-args_info.small_weight_ratio_arg)){
		pol->max_upper_limit[i]=0.0;
		pol->max_lower_limit[i]=0.0;
		pol->max_dx[i]=0.0;
		masked++;
		}
	if(pol->max_dx[i]>max_dx){
		max_dx=pol->max_dx[i];
		max_dx_i=i;
		}

	if(pol->max_upper_limit[i]>masked_max_band[k]){
		masked_max_band[k]=pol->max_upper_limit[i];
		masked_max_band_arg[k]=i;
		}
	}
fprintf(LOG, "masked: %s %ld\n", pol->name, masked);
fprintf(LOG, "largest signal: longitude latitude pol max_dx upper_strain lower_strain freq\n");	
fprintf(LOG, "max_dx: %f %f %s %f %g %g %f\n",fine_grid->longitude[max_dx_i], fine_grid->latitude[max_dx_i], 
				pol->name, pol->max_dx[max_dx_i], 
				pol->max_upper_limit[max_dx_i], 
				pol->max_lower_limit[max_dx_i], pol->freq_map[max_dx_i]);

fprintf(LOG, "max/masked band format: band_num longitude latitude pol max_dx upper_strain freq\n");
for(i=0;i<args_info.dec_bands_arg;i++){

	fprintf(LOG, "max_band: %ld %f %f %s %f %g %f\n", i, fine_grid->longitude[max_band_arg[i]], fine_grid->latitude[max_band_arg[i]], 
				pol->name, pol->max_dx[max_band_arg[i]], 
				max_band[i], 
				pol->freq_map[max_band_arg[i]]);

	fprintf(LOG, "masked_max_band: %ld %f %f %s %f %g %f\n", i, fine_grid->longitude[masked_max_band_arg[i]], fine_grid->latitude[masked_max_band_arg[i]], 
				pol->name, pol->max_dx[masked_max_band_arg[i]], 
				masked_max_band[i], 
				pol->freq_map[masked_max_band_arg[i]]);

        /* old 
	fprintf(LOG, "max_band: %ld %s %g\n", i, pol->name, max_band[i]);
	fprintf(LOG, "masked_max_band: %ld %s %g\n", i, pol->name, masked_max_band[i]);
	*/
	}
	
snprintf(s,19999,"%s_max_upper_strain.png",pol->name);
plot_grid_f(p, fine_grid, pol->max_upper_limit, 1);
RGBPic_dump_png(s, p);

snprintf(s,19999,"%s_max_lower_strain.png",pol->name);
plot_grid_f(p, fine_grid, pol->max_lower_limit, 1);
RGBPic_dump_png(s, p);

snprintf(s,19999,"%s_max_dx.png",pol->name);
plot_grid_f(p, fine_grid, pol->max_dx, 1);
RGBPic_dump_png(s, p);

snprintf(s,19999,"%s_M_map.png",pol->name);
plot_grid_f(p, fine_grid, pol->M_map, 1);
RGBPic_dump_png(s, p);

snprintf(s,19999,"%s_S_map.png",pol->name);
plot_grid_f(p, fine_grid, pol->S_map, 1);
RGBPic_dump_png(s, p);

free_RGBPic(p);
free(max_band);
free(masked_max_band);
free(max_band_arg);
free(masked_max_band_arg);
}

void compute_mean(long pi)
{
SUM_TYPE a,b,c;
long i,k,m;
long offset;
for(k=0,offset=super_grid->first_map[pi];offset>=0;offset=super_grid->list_map[offset],k++){
	for(m=0;m<npolarizations;m++){
		polarizations[m].max_sub_weight[offset]=0.0;
	
		for(i=0;i<useful_bins;i++){

			#ifdef WEIGHTED_SUM		
			if(polarizations[m].fine_grid_weight[i+k*useful_bins]>polarizations[m].max_sub_weight[offset]){
				polarizations[m].max_sub_weight[offset]=polarizations[m].fine_grid_weight[i+k*useful_bins];
				}
	
			c=(polarizations[m].total_weight[offset]-polarizations[m].fine_grid_weight[i+k*useful_bins]);
			#else
			c=(polarizations[m].total_count[offset]-polarizations[m].fine_grid_count[i+k*useful_bins]);
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
SUM_TYPE a,b,c;
long i,k,offset,m;
for(k=0,offset=super_grid->first_map[pi];offset>=0;offset=super_grid->list_map[offset],k++){
	for(m=0;m<npolarizations;m++){
		polarizations[m].max_sub_weight[offset]=0.0;
		
		for(i=0;i<useful_bins;i++){
		
			#ifdef WEIGHTED_SUM		
			c=(polarizations[m].total_weight[offset]);
			#else
			c=(polarizations[m].total_count[offset]);
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
long pi,i,j,k,kk,m;
double a,b;
long total;

normalizing_weight=exp(-M_LN10*TMedian);

stored_fine_bins=super_grid->max_npatch;

allocate_polarization_arrays();
	
/* see comments above variables */
if(args_info.three_bins_arg){
	quantile2std=0.88;
	upper_limit_comp=3.0/(0.85);
	process_patch=process_patch3;
	fprintf(LOG,"mode: 3 bins\n");
	} else 
	{
	quantile2std=1.22;
	upper_limit_comp=1/(0.7*0.85);
	process_patch=process_patch1;
	fprintf(LOG,"mode: 1 bin\n");
	}

fprintf(stderr,"Main loop\n");
for(pi=0;pi<patch_grid->npoints;pi++){
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
	/* compute upper limits */
	for(i=0;i<npolarizations;i++){
		make_limits(&(polarizations[i]), pi);
		}

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
}

