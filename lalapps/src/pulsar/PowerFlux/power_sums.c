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
#include "cmdline.h"

extern struct gengetopt_args_info args_info;

extern SKY_GRID *fine_grid, *patch_grid;
extern SKY_SUPERGRID *super_grid;

extern int first_bin, nbins;
extern SKYMARK *compiled_skymarks;

INT64 spindown_start;

void generate_patch_templates(int pi, POWER_SUM **ps, int *count)
{
POWER_SUM *p;
int i, j, k, kk;
int skyband;
int fshift_count=args_info.nfshift_arg; /* number of frequency offsets */

*count=0;

p=do_alloc(super_grid->max_npatch*args_info.spindown_count_arg*fshift_count, sizeof(*p));
*ps=p;

for(i=0;i<args_info.spindown_count_arg;i++) {
	for(kk=super_grid->first_map[pi];kk>=0;kk=super_grid->list_map[kk]) {
		/* TODO - this effectively requires skybands to not depend on spindown it would be nice if that was not so */
		if(args_info.fine_grid_skymarks_arg)skyband=fine_grid->band[kk];
			else skyband=mark_sky_point(compiled_skymarks, kk, fine_grid->longitude[kk], fine_grid->latitude[kk], args_info.spindown_start_arg+i*args_info.spindown_step_arg);
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

			for(k=0;k<GRID_E_COUNT;k++) {
				p->e[k]=fine_grid->e[k][kk];
				p->patch_e[k]=patch_grid->e[k][pi];
				}

			/* TODO - this effectively requires skybands do not depend on spindown it would be nice if that was not so */			
			p->skyband=skyband;

			p->pps=allocate_partial_power_sum();
			zero_partial_power_sum(p->pps);

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
	(*ps_out)[i].patch_ra=ps[i].patch_ra;
	(*ps_out)[i].patch_dec=ps[i].patch_dec;

	for(k=0;k<GRID_E_COUNT;k++) {
		(*ps_out)[i].e[k]=ps[i].e[k];
		(*ps_out)[i].patch_e[k]=ps[i].patch_e[k];
		}

	(*ps_out)[i].skyband=ps[i].skyband;

	(*ps_out)[i].pps=allocate_partial_power_sum();
	zero_partial_power_sum((*ps_out)[i].pps);
	}
}

void free_templates(POWER_SUM *ps, int count)
{
int i;
for(i=0;i<count;i++) {
	free_partial_power_sum(ps[i].pps);
	ps[i].pps=NULL;
	}
free(ps);
}

void accumulate_power_sums(POWER_SUM *ps, int count, double gps_start, double gps_stop, int veto_mask)
{
int segment_count;
SEGMENT_INFO *si;
int gps_step=24*3600*5; /* set gps step to 5 days */
int i, j;
double gps_idx;
double center_frequency=(first_bin+nbins*0.5)/1800.0;

for(gps_idx=gps_start; gps_idx<gps_stop; gps_idx+=gps_step) {
	si=find_segments(gps_idx, gps_idx+gps_step, veto_mask, &segment_count);
	if(segment_count<1) {
		free(si);
		continue;
		}

	/* loop over templates */
	for(i=0;i<count;i++) {
		/* fill in segment info appropriate to this template */
		for(j=0;j<segment_count;j++) {
			si[j].ra=ps[i].patch_ra;
			si[j].dec=ps[i].patch_dec;
			memcpy(si[j].e, ps[i].patch_e, GRID_E_COUNT*sizeof(SKY_GRID_TYPE));

			si[j].freq_shift=ps[i].freq_shift+ps[i].spindown*(si[j].gps-spindown_start)+
				center_frequency*args_info.doppler_multiplier_arg*(ps[i].e[0]*si[j].detector_velocity[0]
					+ps[i].e[1]*si[j].detector_velocity[1]
					+ps[i].e[2]*si[j].detector_velocity[2]);
			}

		accumulate_single_bin_power_sum_cached1(si, segment_count, ps[i].pps);
		}
	free(si);
	}

for(i=0;i<count;i++) {
	if(ps[i].min_gps<0 || ps[i].min_gps>gps_start)ps[i].min_gps=gps_start;
	if(ps[i].max_gps<0 || ps[i].max_gps<gps_stop)ps[i].max_gps=gps_stop;
	}
}
