#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "global.h"
#include "polarization.h"
#include "cmdline.h"


extern struct gengetopt_args_info args_info;
extern FILE *LOG;
extern SKY_GRID_TYPE *AM_coeffs_plus,*AM_coeffs_cross;
extern long AM_coeffs_size;

extern SKY_GRID *fine_grid, *patch_grid;
extern int lines_list[];

extern long stored_fine_bins;
extern long useful_bins;

int ntotal_polarizations=-1, nlinear_polarizations=-1;
POLARIZATION *polarizations=NULL;


void init_polarizations(void)
{
int i,k, nderived=1;
float a;

nlinear_polarizations=args_info.nlinear_polarizations_arg;;
ntotal_polarizations=nlinear_polarizations+nderived;

fprintf(LOG,"nlinear_polarizations: %d\n", nlinear_polarizations);
fprintf(LOG,"ntotal_polarizations : %d\n", ntotal_polarizations);
if(nlinear_polarizations<2){
	fprintf(stderr,"*** ERROR: number of polarizations is less than 2, aborting\n");
	exit(-1);
	}
if(nlinear_polarizations & 1){
	fprintf(stderr,"*** ERROR: number of polarizations is not even, aborting\n");
	exit(-1);
	}
polarizations=do_alloc(ntotal_polarizations, sizeof(*polarizations));
memset(polarizations, 0, (ntotal_polarizations)*sizeof(*polarizations));

fprintf(stderr,"Initializing polarizations:\n");

polarizations[0].orientation=0;
polarizations[0].name="plus";
polarizations[0].conjugate=&(polarizations[1]);
polarizations[0].plus_proj=1.0;
polarizations[0].cross_proj=0.0;
polarizations[0].AM_coeffs=AM_coeffs_plus;
polarizations[0].plus_factor=1.0;
polarizations[0].cross_factor=0.0;
fprintf(stderr,"\t%s 0.0\n",polarizations[0].name);

polarizations[1].orientation=M_PI/4.0;
polarizations[1].name="cross";
polarizations[1].conjugate=&(polarizations[0]);
polarizations[1].plus_proj=0.0;
polarizations[1].cross_proj=1.0;
polarizations[1].AM_coeffs=AM_coeffs_cross;
polarizations[1].plus_factor=1.0;
polarizations[1].cross_factor=0.0;
fprintf(stderr,"\t%s %g\n",polarizations[1].name, M_PI/4.0);

for(i=2;i<nlinear_polarizations;i++){
	polarizations[i].name=do_alloc(16,sizeof(char));
	polarizations[i].conjugate=&(polarizations[i ^ 1]);
	a=(i>>1)*M_PI/(2.0*nlinear_polarizations);
	if(i & 1){
		a+=M_PI/4.0;
		snprintf(polarizations[i].name,16,"pi_%d_%d",((i+nlinear_polarizations)>>1), 2*nlinear_polarizations);
		} else {
		snprintf(polarizations[i].name,16,"pi_%d_%d", (i>>1), 2*nlinear_polarizations);
		}
	polarizations[i].orientation=a;
	polarizations[i].plus_proj=cos(2*a);
	polarizations[i].cross_proj=sin(2*a);
	polarizations[i].plus_factor=1.0;
	polarizations[i].cross_factor=0.0;
	fprintf(stderr,"\t%s a=%g\n",polarizations[i].name, a);
	
	polarizations[i].AM_coeffs=do_alloc(AM_coeffs_size,sizeof(*(polarizations[i].AM_coeffs)));
	for(k=0;k<AM_coeffs_size;k++){
		polarizations[i].AM_coeffs[k]=AM_coeffs_plus[k]*polarizations[i].plus_proj+
			AM_coeffs_cross[k]*polarizations[i].cross_proj;
		}
	}

polarizations[nlinear_polarizations+0].orientation=0;
polarizations[nlinear_polarizations+0].name="circular";
polarizations[nlinear_polarizations+0].conjugate=&(polarizations[1]);
polarizations[nlinear_polarizations+0].plus_proj=1.0;
polarizations[nlinear_polarizations+0].cross_proj=0.0;
polarizations[nlinear_polarizations+0].AM_coeffs=AM_coeffs_plus;
polarizations[nlinear_polarizations+0].plus_factor=1.0;
polarizations[nlinear_polarizations+0].cross_factor=1.0;
fprintf(stderr,"\t%s %f %f\n",polarizations[nlinear_polarizations+0].name,
	polarizations[nlinear_polarizations+0].plus_factor, 
	polarizations[nlinear_polarizations+0].cross_factor);


for(i=0;i<ntotal_polarizations;i++){
	polarizations[i].patch_CutOff=do_alloc(patch_grid->npoints,sizeof(*polarizations[i].patch_CutOff));
	}
fflush(LOG);
}

void allocate_polarization_arrays(void)
{
long total,i,k;
total=ntotal_polarizations*sizeof(*polarizations[0].fine_grid_sum);

#ifdef COMPUTE_SIGMA
total+=ntotal_polarizations*sizeof(*polarizations[0].fine_grid_sq_sum);
#endif

if(lines_list[0]>=0){
	#ifdef WEIGHTED_SUM
	total+=ntotal_polarizations*sizeof(*polarizations[0].fine_grid_sum);
	#else
	total+=ntotal_polarizations*sizeof(*polarizations[0].fine_grid_count);
	#endif
	}

if(args_info.compute_betas_arg){
	fprintf(LOG,"Compute betas: true\n");
	} else {
	fprintf(LOG,"Compute betas: false\n");
	}


fprintf(stderr, "Allocating accumulation arrays, (%.1f KB total)\n", 
	(stored_fine_bins*useful_bins*total)/(1024.0));
fprintf(LOG, "Accumulation set size: %f KB\n", 
	(stored_fine_bins*useful_bins*total)/(1024.0));
	
fprintf(stderr, "Skymap arrays size: %.1f MB\n", ntotal_polarizations*(11.0+2.0*args_info.compute_betas_arg)*fine_grid->npoints*sizeof(SUM_TYPE)/(1024.0*1024.0));
fprintf(LOG, "Skymap arrays size: %f MB\n", ntotal_polarizations*(11.0+2.0*args_info.compute_betas_arg)*fine_grid->npoints*sizeof(SUM_TYPE)/(1024.0*1024.0));

fprintf(stderr, "Spectral plot arrays size: %.1f KB\n", ntotal_polarizations*7.0*useful_bins*args_info.nbands_arg*sizeof(SUM_TYPE)/1024.0);
fprintf(LOG, "Spectral plot arrays size: %f KB\n", ntotal_polarizations*7.0*useful_bins*args_info.nbands_arg*sizeof(SUM_TYPE)/1024.0);
		
for(i=0;i<ntotal_polarizations;i++){
	/* Accumulation arrays */
	polarizations[i].fine_grid_sum=do_alloc(stored_fine_bins*useful_bins,sizeof(*polarizations[i].fine_grid_sum));
	#ifdef COMPUTE_SIGMA
	polarizations[i].fine_grid_sq_sum=do_alloc(stored_fine_bins*useful_bins,sizeof(*polarizations[i].fine_grid_sq_sum));
	#endif


	#ifdef WEIGHTED_SUM
	if(lines_list[0]>=0){
		polarizations[i].fine_grid_weight=do_alloc(stored_fine_bins*useful_bins,sizeof(*polarizations[i].fine_grid_weight));
		}
	
	polarizations[i].skymap.total_weight=do_alloc(fine_grid->npoints,sizeof(*polarizations[i].fine_grid_weight));
	#else
	if(lines_list[0]>=0){
		polarizations[i].fine_grid_count=do_alloc(stored_fine_bins*useful_bins,sizeof(*polarizations[i].fine_grid_count));
		}
	
	polarizations[i].skymap.total_count=do_alloc(fine_grid->npoints,sizeof(*polarizations[i].fine_grid_count));
	#endif

	/* Output arrrays - skymaps*/
	polarizations[i].skymap.max_sub_weight=do_alloc(fine_grid->npoints, sizeof(SUM_TYPE));
	polarizations[i].skymap.max_dx=do_alloc(fine_grid->npoints, sizeof(SUM_TYPE));
	polarizations[i].skymap.M_map=do_alloc(fine_grid->npoints, sizeof(SUM_TYPE));
	polarizations[i].skymap.S_map=do_alloc(fine_grid->npoints, sizeof(SUM_TYPE));
	polarizations[i].skymap.max_upper_limit=do_alloc(fine_grid->npoints, sizeof(SUM_TYPE));
	polarizations[i].skymap.max_lower_limit=do_alloc(fine_grid->npoints, sizeof(SUM_TYPE));
	polarizations[i].skymap.freq_map=do_alloc(fine_grid->npoints, sizeof(SUM_TYPE));
	polarizations[i].skymap.cor1=do_alloc(fine_grid->npoints, sizeof(SUM_TYPE));
	polarizations[i].skymap.cor2=do_alloc(fine_grid->npoints, sizeof(SUM_TYPE));
	polarizations[i].skymap.ks_test=do_alloc(fine_grid->npoints, sizeof(SUM_TYPE));
	polarizations[i].skymap.ks_count=do_alloc(fine_grid->npoints, sizeof(SUM_TYPE));

	if(args_info.compute_betas_arg){
		polarizations[i].skymap.beta1=do_alloc(fine_grid->npoints, sizeof(SUM_TYPE));
		polarizations[i].skymap.beta2=do_alloc(fine_grid->npoints, sizeof(SUM_TYPE));
		} else {
		polarizations[i].skymap.beta1=NULL;
		polarizations[i].skymap.beta2=NULL;
		}

	for(k=0;k<fine_grid->npoints;k++){
		polarizations[i].skymap.max_dx[k]=-1.0;
		polarizations[i].skymap.M_map[k]=-1.0;
		polarizations[i].skymap.S_map[k]=-1.0;
		polarizations[i].skymap.max_upper_limit[k]=-1.0;
		polarizations[i].skymap.max_lower_limit[k]=-1.0;
		polarizations[i].skymap.freq_map[k]=-1.0;
		polarizations[i].skymap.cor1[k]=-1.0;
		polarizations[i].skymap.cor2[k]=-1.0;
		polarizations[i].skymap.ks_test[k]=-1.0;
		polarizations[i].skymap.ks_count[k]=-1.0;
		#ifdef WEIGHTED_SUM
		polarizations[i].skymap.total_weight[k]=0.0;
		#else
		polarizations[i].skymap.total_count[k]=0;
		#endif

		if(args_info.compute_betas_arg){
			polarizations[i].skymap.beta1[k]=0.0;
			polarizations[i].skymap.beta2[k]=0.0;
			}
		}
	/* Output arrays - spectral plots */
	polarizations[i].spectral_plot.max_upper_limit=do_alloc(useful_bins*args_info.nbands_arg, sizeof(SUM_TYPE));
	polarizations[i].spectral_plot.ul_dec=do_alloc(useful_bins*args_info.nbands_arg, sizeof(SUM_TYPE));
	polarizations[i].spectral_plot.ul_ra=do_alloc(useful_bins*args_info.nbands_arg, sizeof(SUM_TYPE));
	polarizations[i].spectral_plot.max_dx=do_alloc(useful_bins*args_info.nbands_arg, sizeof(SUM_TYPE));
	polarizations[i].spectral_plot.dx_dec=do_alloc(useful_bins*args_info.nbands_arg, sizeof(SUM_TYPE));
	polarizations[i].spectral_plot.dx_ra=do_alloc(useful_bins*args_info.nbands_arg, sizeof(SUM_TYPE));
	polarizations[i].spectral_plot.max_mask_ratio=do_alloc(useful_bins*args_info.nbands_arg, sizeof(SUM_TYPE));

	for(k=0;k<useful_bins*args_info.nbands_arg;k++){
		polarizations[i].spectral_plot.max_upper_limit[k]=-1.0;
		polarizations[i].spectral_plot.ul_dec[k]=-10.0;
		polarizations[i].spectral_plot.ul_ra[k]=-10.0;
		polarizations[i].spectral_plot.max_dx[k]=-1.0;
		polarizations[i].spectral_plot.dx_dec[k]=-10.0;
		polarizations[i].spectral_plot.dx_ra[k]=-10.0;
		polarizations[i].spectral_plot.max_mask_ratio[k]=-1.0;
		}
	}

}

void clear_accumulation_arrays(void)
{
long i,k;

for(k=0;k<ntotal_polarizations;k++){
	for(i=0;i<stored_fine_bins*useful_bins;i++){
		polarizations[k].fine_grid_sum[i]=0.0;
			#ifdef COMPUTE_SIGMA
		polarizations[k].fine_grid_sq_sum[i]=0.0;
		#endif
			if(lines_list[0]>=0){
			#ifdef WEIGHTED_SUM	
			polarizations[k].fine_grid_weight[i]=0.0;
			#else	
			polarizations[k].fine_grid_count[i]=0;
			#endif
			}
		}
	}
}
