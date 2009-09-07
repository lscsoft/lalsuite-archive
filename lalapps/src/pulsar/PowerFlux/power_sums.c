#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/* We need this define to get NAN values */
#define __USE_ISOC99
#include <math.h>

#include "global.h"
#include "power_cache.h"
#include "power_sums.h"
#include "dataset.h"
#include "grid.h"
#include "skymarks.h"
#include "summing_context.h"
#include "cmdline.h"

extern struct gengetopt_args_info args_info;

extern SKY_GRID *fine_grid, *patch_grid;
extern SKY_SUPERGRID *super_grid;

extern int first_bin, nbins, useful_bins;
extern SKYMARK *compiled_skymarks;

extern DATASET *datasets;

INT64 spindown_start;

void generate_patch_templates(int pi, POWER_SUM **ps, int *count)
{
POWER_SUM *p;
int i, j, k, kk;
int skyband;
int fshift_count=args_info.nfshift_arg; /* number of frequency offsets */

float e[GRID_E_COUNT];
float patch_e[GRID_E_COUNT];

*count=0;

p=do_alloc(super_grid->max_npatch*args_info.spindown_count_arg*fshift_count, sizeof(*p));
*ps=p;

for(i=0;i<args_info.spindown_count_arg;i++) {
	for(kk=super_grid->first_map[pi];kk>=0;kk=super_grid->list_map[kk]) {
		for(k=0;k<GRID_E_COUNT;k++) {
			e[k]=fine_grid->e[k][kk];
			patch_e[k]=patch_grid->e[k][pi];
			}

		/* TODO - this effectively requires skybands to not depend on spindown it would be nice if that was not so */
		if(args_info.fine_grid_skymarks_arg)skyband=fine_grid->band[kk];
			else skyband=mark_sky_point(compiled_skymarks, kk, fine_grid->longitude[kk], fine_grid->latitude[kk], e, args_info.spindown_start_arg+i*args_info.spindown_step_arg);
		if(skyband<0)continue;

		for(j=0;j<fshift_count;j++) {
			p->freq_shift=args_info.frequency_offset_arg+j/(1800.0*fshift_count);
			p->spindown=args_info.spindown_start_arg+i*args_info.spindown_step_arg;
			p->ra=fine_grid->longitude[kk];
			p->dec=fine_grid->latitude[kk];
			p->min_gps=-1;
			p->max_gps=-1;
			p->patch_ra=patch_grid->longitude[pi];
			p->patch_dec=patch_grid->latitude[pi];

			memcpy(p->e, e, GRID_E_COUNT*sizeof(float));
			memcpy(p->patch_e, patch_e, GRID_E_COUNT*sizeof(float));

			/* TODO - this effectively requires skybands do not depend on spindown it would be nice if that was not so */			
			p->skyband=skyband;

			p->pps=allocate_partial_power_sum_F(useful_bins);
			zero_partial_power_sum_F(p->pps);

			(*count)++;
			p++;
			}
		}
	}
}

void clone_templates(POWER_SUM *ps, int count, POWER_SUM **ps_out)
{
int i, k;
*ps_out=do_alloc(count, sizeof(POWER_SUM));

for(i=0;i<count;i++) {
	(*ps_out)[i].freq_shift=ps[i].freq_shift;
	(*ps_out)[i].spindown=ps[i].spindown;
	(*ps_out)[i].ra=ps[i].ra;
	(*ps_out)[i].dec=ps[i].dec;
	(*ps_out)[i].min_gps=ps[i].min_gps;
	(*ps_out)[i].max_gps=ps[i].max_gps;
	(*ps_out)[i].patch_ra=ps[i].patch_ra;
	(*ps_out)[i].patch_dec=ps[i].patch_dec;

	for(k=0;k<GRID_E_COUNT;k++) {
		(*ps_out)[i].e[k]=ps[i].e[k];
		(*ps_out)[i].patch_e[k]=ps[i].patch_e[k];
		}

	(*ps_out)[i].skyband=ps[i].skyband;

	(*ps_out)[i].pps=allocate_partial_power_sum_F(useful_bins);
	zero_partial_power_sum_F((*ps_out)[i].pps);
	}
}

void free_templates(POWER_SUM *ps, int count)
{
int i;
for(i=0;i<count;i++) {
	free_partial_power_sum_F(ps[i].pps);
	ps[i].pps=NULL;
	}
free(ps);
}

void accumulate_power_sums(SUMMING_CONTEXT *ctx, POWER_SUM *ps, int count, double gps_start, double gps_stop, int veto_mask)
{
int segment_count;
SEGMENT_INFO *si, *si_local;
POWER_SUM *ps_local;
DATASET *d;
POLARIZATION *pl;
int gps_step=args_info.summing_step_arg;
int i, j, k;
float min_shift, max_shift, a;
double gps_idx;
float center_frequency=(first_bin+nbins*0.5);
int group_count=24;
/* for work with Doppler shifts sidereal day is best */
//#define SCALER (GROUP_COUNT/(23.0*3600.0+56.0*60.0+4.0))
float group_scaler=group_count; /* this determines sub-bin resolution in group formation */
SEGMENT_INFO **groups;
int *group_segment_count;
float avg_spindown=args_info.spindown_start_arg+0.5*args_info.spindown_step_arg*(args_info.spindown_count_arg-1);

float *patch_e=ps[0].patch_e; /* set of coefficients for this patch, used for amplitude response and bin shift estimation */

for(gps_idx=gps_start; gps_idx<gps_stop; gps_idx+=gps_step) {

	si=find_segments(gps_idx, gps_idx+gps_step, veto_mask, &segment_count);
	if(segment_count<1) {
		free(si);
		continue;
		}

	/* This assumes that we are patch bound !! *** WARNING ***
	TODO: make this assumption automatic in the data layout.
	 */
	si_local=si;
	min_shift=1000000; /* we should never load this many bins */
	max_shift=-1000000;
	for(j=0;j<segment_count;j++) {
// 		si_local->ra=ps[0].patch_ra;
// 		si_local->dec=ps[0].patch_dec;
// 		memcpy(si_local->e, ps[0].patch_e, GRID_E_COUNT*sizeof(SKY_GRID_TYPE));

		d=&(datasets[si_local->dataset]);
		pl=&(d->polarizations[0]);

// 		si_local->f_plus=F_plus_coeff(si_local->segment,  patch_e, pl->AM_coeffs);
// 		si_local->f_cross=F_plus_coeff(si_local->segment,  patch_e, pl->conjugate->AM_coeffs);

		si_local->f_plus=F_plus_coeff(si_local->segment,  patch_e, d->AM_coeffs_plus);
		si_local->f_cross=F_plus_coeff(si_local->segment,  patch_e, d->AM_coeffs_cross);


		a=center_frequency*(float)args_info.doppler_multiplier_arg*(patch_e[0]*si_local->detector_velocity[0]
						+patch_e[1]*si_local->detector_velocity[1]
						+patch_e[2]*si_local->detector_velocity[2])
			+si_local->coherence_time*avg_spindown*(float)(si_local->gps-spindown_start);
		if(a<min_shift)min_shift=a;
		if(a>max_shift)max_shift=a;
		si_local++;
		}

	//group_count=ceil(max_shift-min_shift);
	if(group_count>200) {
		fprintf(stderr, "Warning group count too large: %d\n", group_count);
		group_count=200;
		}

	group_segment_count=do_alloc(group_count, sizeof(*group_segment_count));
	groups=do_alloc(group_count, sizeof(*groups));

	for(k=0;k<group_count;k++) {
		group_segment_count[k]=0;
		groups[k]=do_alloc(segment_count, sizeof(SEGMENT_INFO));
		}

	/* group segments into bunches with similar shifts - mostly by sidereal time
           this way there is larger correllation of frequency shifts during summing and better use of power cache */
	si_local=si;
	for(j=0;j<segment_count;j++) {
		a=(center_frequency*(float)args_info.doppler_multiplier_arg*(patch_e[0]*si_local->detector_velocity[0]
						+patch_e[1]*si_local->detector_velocity[1]
						+patch_e[2]*si_local->detector_velocity[2])
			+si_local->coherence_time*avg_spindown*(float)(si_local->gps-spindown_start));
		//a*=0.25;
		k=floorf((a-floorf(a))*group_count);
		if(k<0)k=0;
		if(k>=group_count)k=group_count-1;

		memcpy(&(groups[k][group_segment_count[k]]), si_local, sizeof(SEGMENT_INFO));
		group_segment_count[k]++;
		
		si_local++;
		}

// 	for(k=0;k<GROUP_COUNT;k++) {
// 		fprintf(stderr, "group %d has %d segments\n", k, group_segment_count[k]);
// 		}

	/* loop over groups */

	for(k=0;k<group_count;k++) {
 		//fprintf(stderr, "group %d has %d segments\n", k, group_segment_count[k]);
		if(group_segment_count[k]<1)continue;
		ctx->reset_cache(ctx, group_segment_count[k], count);
	
		/* loop over templates */
		ps_local=ps;
		for(i=0;i<count;i++) {
			/* fill in segment info appropriate to this template */
			si_local=groups[k];
			for(j=0;j<group_segment_count[k];j++) {
	// 			si[j].ra=ps[i].patch_ra;
	// 			si[j].dec=ps[i].patch_dec;
	// 			memcpy(si[j].e, ps[i].patch_e, GRID_E_COUNT*sizeof(SKY_GRID_TYPE));
	
				si_local->bin_shift=si_local->coherence_time*(ps_local->freq_shift+ps_local->spindown*(float)(si_local->gps-spindown_start))+
					center_frequency*(float)args_info.doppler_multiplier_arg*(ps_local->e[0]*si_local->detector_velocity[0]
						+ps_local->e[1]*si_local->detector_velocity[1]
						+ps_local->e[2]*si_local->detector_velocity[2]);
				si_local++;
				}
	
			//accumulate_single_bin_power_sum_cached1(groups[k], group_segment_count[k], ps_local->pps);
			ctx->accumulate_power_sum_cached(ctx, groups[k], group_segment_count[k], ps_local->pps);
			ps_local++;
			}
		}
	for(k=0;k<group_count;k++) {
		free(groups[k]);
		}
	free(groups);
	free(group_segment_count);
	free(si);
	}

for(i=0;i<count;i++) {
	if(ps[i].min_gps<0 || ps[i].min_gps>gps_start)ps[i].min_gps=gps_start;
	if(ps[i].max_gps<0 || ps[i].max_gps<gps_stop)ps[i].max_gps=gps_stop;
	}
}
