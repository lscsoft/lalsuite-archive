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

int npolarizations=2;
POLARIZATION *polarizations=NULL;


void init_polarizations(void)
{
int i,k;
float a;
npolarizations=args_info.npolarizations_arg;
fprintf(LOG,"npolarizations: %d\n", npolarizations);
if(npolarizations<2){
	fprintf(stderr,"*** ERROR: number of polarizations is less than 2, aborting\n");
	}
polarizations=do_alloc(npolarizations, sizeof(*polarizations));
memset(polarizations, 0, npolarizations*sizeof(*polarizations));

fprintf(stderr,"Initializing polarizations:\n");

polarizations[0].orientation=0;
polarizations[0].name="plus";
polarizations[0].plus_proj=1.0;
polarizations[0].cross_proj=0.0;
polarizations[0].AM_coeffs=AM_coeffs_plus;
fprintf(stderr,"\t%s 0.0\n",polarizations[0].name);

polarizations[1].orientation=M_PI/4.0;
polarizations[1].name="cross";
polarizations[1].plus_proj=0.0;
polarizations[1].cross_proj=1.0;
polarizations[1].AM_coeffs=AM_coeffs_cross;
fprintf(stderr,"\t%s %g\n",polarizations[1].name, M_PI/4.0);

for(i=2;i<npolarizations;i++){
	a=(i-1)*M_PI/(4.0*(npolarizations-1));
	polarizations[i].orientation=a;
	polarizations[i].plus_proj=cos(2*a);
	polarizations[i].cross_proj=sin(2*a);
	polarizations[i].name=do_alloc(16,sizeof(char));
	snprintf(polarizations[i].name,16,"pol_%d_%d",i-1, 4*(npolarizations-1));
	fprintf(stderr,"\t%s a=%g\n",polarizations[i].name, a);
	
	polarizations[i].AM_coeffs=do_alloc(AM_coeffs_size,sizeof(*(polarizations[i].AM_coeffs)));
	for(k=0;k<AM_coeffs_size;k++){
		polarizations[i].AM_coeffs[k]=AM_coeffs_plus[k]*polarizations[i].plus_proj+
			AM_coeffs_cross[k]*polarizations[i].cross_proj;
		}
	}

for(i=0;i<npolarizations;i++){
	polarizations[i].patch_CutOff=do_alloc(patch_grid->npoints,sizeof(*polarizations[i].patch_CutOff));
	}
fflush(LOG);
}

void allocate_polarization_arrays(void)
{
long total,i,k;
total=npolarizations*sizeof(*polarizations[0].fine_grid_sum);

#ifdef COMPUTE_SIGMA
total+=npolarizations*sizeof(*polarizations[0].fine_grid_sq_sum);
#endif

if(lines_list[0]>=0){
	#ifdef WEIGHTED_SUM
	total+=npolarizations*sizeof(*polarizations[0].fine_grid_sum);
	#else
	total+=npolarizations*sizeof(*polarizations[0].fine_grid_count);
	#endif
	}

fprintf(stderr,"Allocating accumulation arrays, (%.1f KB total)\n", 
	(stored_fine_bins*useful_bins*total)/(1024.0));
fprintf(LOG," accumulation set size: %f KB\n", 
	(stored_fine_bins*useful_bins*total)/(1024.0));
		
for(i=0;i<npolarizations;i++){
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

for(k=0;k<npolarizations;k++){
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
