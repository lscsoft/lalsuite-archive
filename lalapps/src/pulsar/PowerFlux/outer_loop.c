#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/* We need this define to get NAN values */
#define __USE_ISOC99
#include <math.h>
#include <time.h>

#include "global.h"
#include "power_cache.h"
#include "power_sums.h"
#include "power_sum_stats.h"
#include "dataset.h"
#include "grid.h"
#include "cmdline.h"
#include "rastermagic.h"
#include "hookup.h"
#include "outer_loop.h"

extern struct gengetopt_args_info args_info;

extern SKY_GRID *fine_grid, *patch_grid;

extern FILE * DATA_LOG, * LOG, *FILE_LOG;

extern int first_bin, side_cut;

extern DATASET *datasets;
extern int d_free;

int data_log_index=0;

int write_data_log_header=1;

typedef struct {
	unsigned int veto_mask;
	char *name;
	} VETO_INFO;

/* we should have less than 5 detectors. Also, veto has only 8 bits some of which are used already - rework the scheme when one needs to accomodate more detectors */
VETO_INFO veto_info[4];
int veto_free=0;

void assign_veto(void)
{
int i,k, m;
memset(veto_info, 0, 10*sizeof(*veto_info));
veto_free=0;
for(k=0;k<d_free;k++) {
	for(m=0;veto_info[m].veto_mask!=0;m++) {
		if(!strcmp(veto_info[m].name, datasets[k].detector))break;
		if(m>=4) {
			fprintf(stderr, "*** INTERNAL ERROR: PowerFlux cannot handle this many detectors\n");
			exit(-1);
			}
		}
	if(veto_info[m].veto_mask==0) {
		veto_free=m+1;
		veto_info[m].name=strdup(datasets[k].detector);
		veto_info[m].veto_mask=(1<<(4+m));
		}
	for(i=0;i<datasets[k].free;i++) {
		datasets[k].sft_veto[i]=(datasets[k].sft_veto[i] & ~ ((1<<4)-1)) | veto_info[m].veto_mask;
		}
	}
}

void log_extremes(EXTREME_INFO *ei, int pi, POWER_SUM **ps, int nchunks, int count)
{
PARTIAL_POWER_SUM_F *pps;
POWER_SUM_STATS pstats, pstats_accum;
int i, k;
int highest_ul_idx=0;
int highest_circ_ul_idx=0;
int highest_snr_idx=0;
int skyband;

pps=allocate_partial_power_sum_F();

for(i=0;i<count;i++) {
	zero_partial_power_sum_F(pps);
	for(k=0;k<nchunks;k++) {
		accumulate_partial_power_sum_F1(pps, (ps[k][i].pps));
		}
	power_sum_stats(pps, &(pstats));

	skyband=ps[0][i].skyband;

	#define FILL_EXTRA_PARAMS(target) {\
		target.ra=ps[0][i].ra; \
		target.dec=ps[0][i].dec; \
		target.spindown=ps[0][i].spindown; \
		target.frequency=ps[0][i].freq_shift+((target).bin+first_bin+side_cut)/1800.0; \
		}

	#define FILL_POINT_STATS(target, source)	{\
		memcpy(&(target), &(source), sizeof(target)); \
		FILL_EXTRA_PARAMS(target); \
		}

	if(ei->band_info[skyband].max_weight<0) {
		memcpy(&(ei->band_info[skyband]), &pstats, sizeof(pstats));
		FILL_EXTRA_PARAMS(ei->band_info[skyband].highest_ul);
		FILL_EXTRA_PARAMS(ei->band_info[skyband].highest_circ_ul);
		FILL_EXTRA_PARAMS(ei->band_info[skyband].highest_snr);
		FILL_EXTRA_PARAMS(ei->band_info[skyband].highest_ks);
		} else {

		if(pstats.highest_ul.ul>ei->band_info[skyband].highest_ul.ul) {
			FILL_POINT_STATS(ei->band_info[skyband].highest_ul, pstats.highest_ul);
			}
	
		if(pstats.highest_circ_ul.ul>ei->band_info[skyband].highest_circ_ul.ul) {
			FILL_POINT_STATS(ei->band_info[skyband].highest_circ_ul, pstats.highest_circ_ul);
			}
	
		if(pstats.highest_snr.snr>ei->band_info[skyband].highest_snr.snr) {
			FILL_POINT_STATS(ei->band_info[skyband].highest_snr, pstats.highest_snr);
			}
	
		if(pstats.highest_ks.ks_value>ei->band_info[skyband].highest_ks.ks_value) {
			FILL_POINT_STATS(ei->band_info[skyband].highest_ks, pstats.highest_ks);
			}
	
		if(pstats.max_weight>ei->band_info[skyband].max_weight) {
			ei->band_info[skyband].max_weight=pstats.max_weight;
			}
	
		if(pstats.min_weight<ei->band_info[skyband].min_weight) {
			ei->band_info[skyband].min_weight=pstats.min_weight;
			}
	
		if(pstats.max_weight_loss_fraction>ei->band_info[skyband].max_weight_loss_fraction) {
			ei->band_info[skyband].max_weight_loss_fraction=pstats.max_weight_loss_fraction;
			}

		}


	if(i==0) {
		memcpy(&pstats_accum, &pstats, sizeof(pstats));
		continue;
		}

	if(pstats.highest_ul.ul>pstats_accum.highest_ul.ul) {
		memcpy(&pstats_accum.highest_ul, &pstats.highest_ul, sizeof(pstats.highest_ul));
		highest_ul_idx=i;
		}

	if(pstats.highest_circ_ul.ul>pstats_accum.highest_circ_ul.ul) {
		memcpy(&pstats_accum.highest_circ_ul, &pstats.highest_circ_ul, sizeof(pstats.highest_circ_ul));
		highest_circ_ul_idx=i;
		}

	if(pstats.highest_snr.snr>pstats_accum.highest_snr.snr) {
		memcpy(&pstats_accum.highest_snr, &pstats.highest_snr, sizeof(pstats.highest_snr));
		highest_snr_idx=i;
		}

	if(pstats.highest_ks.ks_value>pstats_accum.highest_ks.ks_value) {
		memcpy(&pstats_accum.highest_ks, &pstats.highest_ks, sizeof(pstats.highest_ks));
		}

	if(pstats.max_weight>pstats_accum.max_weight) {
		pstats_accum.max_weight=pstats.max_weight;
		}

	if(pstats.min_weight<pstats_accum.min_weight) {
		pstats_accum.min_weight=pstats.min_weight;
		}

	if(pstats.max_weight_loss_fraction>pstats_accum.max_weight_loss_fraction) {
		pstats_accum.max_weight_loss_fraction=pstats.max_weight_loss_fraction;
		}
	}
free_partial_power_sum_F(pps);

if(write_data_log_header) {
	write_data_log_header=0;
	fprintf(DATA_LOG, "kind index set pi pps_count first_bin min_gps max_gps skyband frequency spindown ra dec iota psi snr ul ll M S ks_value ks_count frequency_bin max_weight weight_loss_fraction max_ks_value max_weight_loss_fraction\n");
	}

/* now that we know extreme points go and characterize them */
#define WRITE_POINT(psum, pstat, kind)	{\
	fprintf(DATA_LOG, "%s %d %s %d %d %d %lf %lf %d %lf %lg %lf %lf %lf %lf %lf %lg %lg %lg %lg %lf %d %d %lg %lf %lf %lf\n", \
		kind, \
		data_log_index, \
		ei->name, \
		pi, \
		count, \
		first_bin+side_cut, \
		psum.min_gps, \
		psum.max_gps, \
		psum.skyband, \
		(pstat.bin+first_bin+side_cut)/1800.0+psum.freq_shift, \
		psum.spindown, \
		psum.ra, \
		psum.dec, \
		pstat.iota, \
		pstat.psi, \
		pstat.snr, \
		pstat.ul, \
		pstat.ll, \
		pstat.M, \
		pstat.S, \
		pstat.ks_value, \
		pstat.ks_count, \
		pstat.bin, \
		pstat.max_weight, \
		pstat.weight_loss_fraction, \
		pstats_accum.highest_ks.ks_value, \
		pstats_accum.max_weight_loss_fraction \
		); data_log_index++; }

if(args_info.output_cache_arg) {
	WRITE_POINT(ps[0][highest_ul_idx], pstats_accum.highest_ul, "ul");
	WRITE_POINT(ps[0][highest_circ_ul_idx], pstats_accum.highest_circ_ul, "circ");
	}
if(pstats_accum.highest_snr.snr>args_info.min_candidate_snr_arg)WRITE_POINT(ps[0][highest_snr_idx], pstats_accum.highest_snr, "snr");

#define FILL_SKYMAP(skymap, value)	if(ei->skymap!=NULL)ei->skymap[pi]=value;

FILL_SKYMAP(ul_skymap, pstats_accum.highest_ul.ul);
FILL_SKYMAP(ul_freq_skymap, (pstats_accum.highest_ul.bin)/1800.0+ps[0][highest_ul_idx].freq_shift);

FILL_SKYMAP(circ_ul_skymap, pstats_accum.highest_circ_ul.ul);
FILL_SKYMAP(circ_ul_freq_skymap, (pstats_accum.highest_circ_ul.bin)/1800.0+ps[0][highest_circ_ul_idx].freq_shift);

FILL_SKYMAP(snr_skymap, pstats_accum.highest_snr.snr);
FILL_SKYMAP(snr_ul_skymap, pstats_accum.highest_snr.ul);
FILL_SKYMAP(snr_freq_skymap, (pstats_accum.highest_snr.bin)/1800.0+ps[0][highest_snr_idx].freq_shift);

FILL_SKYMAP(max_weight_skymap, pstats_accum.max_weight);
FILL_SKYMAP(min_weight_skymap, pstats_accum.min_weight);
FILL_SKYMAP(weight_loss_fraction_skymap, pstats_accum.max_weight_loss_fraction);

FILL_SKYMAP(ks_skymap, pstats_accum.highest_ks.ks_value);
}

static char s[20000];

EXTREME_INFO * allocate_extreme_info(char *name)
{
EXTREME_INFO *ei;
int i;

ei=do_alloc(1, sizeof(*ei));
memset(ei, 0, sizeof(*ei));

ei->name=strdup(name);

if(args_info.compute_skymaps_arg) {
	ei->ul_skymap=do_alloc(patch_grid->npoints, sizeof(float));
	ei->ul_freq_skymap=do_alloc(patch_grid->npoints, sizeof(float));
	ei->circ_ul_skymap=do_alloc(patch_grid->npoints, sizeof(float));
	ei->circ_ul_freq_skymap=do_alloc(patch_grid->npoints, sizeof(float));
	ei->snr_skymap=do_alloc(patch_grid->npoints, sizeof(float));
	ei->snr_ul_skymap=do_alloc(patch_grid->npoints, sizeof(float));
	ei->snr_freq_skymap=do_alloc(patch_grid->npoints, sizeof(float));
	ei->max_weight_skymap=do_alloc(patch_grid->npoints, sizeof(float));
	ei->min_weight_skymap=do_alloc(patch_grid->npoints, sizeof(float));
	ei->weight_loss_fraction_skymap=do_alloc(patch_grid->npoints, sizeof(float));
	ei->ks_skymap=do_alloc(patch_grid->npoints, sizeof(float));
	}

ei->band_info=do_alloc(fine_grid->nbands, sizeof(*ei->band_info));
memset(ei->band_info, 0, fine_grid->nbands*sizeof(*ei->band_info));

for(i=0;i<fine_grid->nbands;i++) {
	ei->band_info[i].max_weight=-1;
	}

return ei;
}

void free_extreme_info(EXTREME_INFO *ei)
{

free(ei->name);
free(ei->band_info);

#define FREE(x)  {if(ei->x!=NULL)free(ei->x); ei->x=NULL; }

FREE(ul_skymap);
FREE(ul_freq_skymap);
FREE(circ_ul_skymap);
FREE(circ_ul_freq_skymap);
FREE(snr_skymap);
FREE(snr_ul_skymap);
FREE(snr_freq_skymap);
FREE(max_weight_skymap);
FREE(min_weight_skymap);
FREE(weight_loss_fraction_skymap);
FREE(ks_skymap);

free(ei);
}

void output_extreme_info(RGBPic *p, EXTREME_INFO *ei)
{
int skyband;

fprintf(LOG, "kind set first_bin skyband frequency spindown ra dec iota psi snr ul ll M S ks_value ks_count frequency_bin max_weight weight_loss_fraction max_ks_value max_weight_loss_fraction\n");

/* now that we know extreme points go and characterize them */
#define WRITE_SKYBAND_POINT(pstat, kind)	\
	fprintf(LOG, "band_info: %s %d %s %s %d %lf %lg %lf %lf %lf %lf %lf %lg %lg %lg %lg %lf %d %d %lg %lf %lf %lf\n", \
		kind, \
		skyband, \
		fine_grid->band_name[skyband], \
		ei->name, \
		first_bin+side_cut, \
		pstat.frequency, \
		pstat.spindown, \
		pstat.ra, \
		pstat.dec, \
		pstat.iota, \
		pstat.psi, \
		pstat.snr, \
		pstat.ul, \
		pstat.ll, \
		pstat.M, \
		pstat.S, \
		pstat.ks_value, \
		pstat.ks_count, \
		pstat.bin, \
		pstat.max_weight, \
		pstat.weight_loss_fraction, \
		ei->band_info[skyband].highest_ks.ks_value, \
		ei->band_info[skyband].max_weight_loss_fraction \
		); 

for(skyband=0;skyband<fine_grid->nbands;skyband++) {
	WRITE_SKYBAND_POINT(ei->band_info[skyband].highest_ul, "ul");
	WRITE_SKYBAND_POINT(ei->band_info[skyband].highest_circ_ul, "circ");
	WRITE_SKYBAND_POINT(ei->band_info[skyband].highest_snr, "snr");

	/* old-style output for compatibility, but we have to attach ei->name to distinguish different sets */
	fprintf(LOG, "max_high_ul_band: %d %lg %lf %lg %s\n", skyband, 
		ei->band_info[skyband].highest_ul.ul,
		ei->band_info[skyband].highest_ul.frequency,
		ei->band_info[skyband].highest_ul.spindown,
		ei->name);

	fprintf(LOG, "max_circ_ul_band: %d %lg %lf %lg %s\n", skyband, 
		ei->band_info[skyband].highest_circ_ul.ul,
		ei->band_info[skyband].highest_circ_ul.frequency,
		ei->band_info[skyband].highest_circ_ul.spindown,
		ei->name);

	fprintf(LOG, "max_dx_band: %d %s %lf %lf %lg %lf %lf %lf %lf %s\n", skyband,
		fine_grid->band_name[skyband],
		ei->band_info[skyband].highest_snr.snr,
		ei->band_info[skyband].highest_snr.frequency,
		ei->band_info[skyband].highest_snr.spindown,
		ei->band_info[skyband].highest_snr.ra,
		ei->band_info[skyband].highest_snr.dec,
		ei->band_info[skyband].highest_snr.iota,
		ei->band_info[skyband].highest_snr.psi,
		ei->name
		);

	}

#define OUTPUT_SKYMAP(array, tag)	{ \
	if(ei->array!=NULL) { \
		fprintf(stderr, "\t%s_%s\n", ei->name, tag); \
		snprintf(s,19999, "%s_%s.png", ei->name, tag); \
		if(clear_name_png(s)) { \
			plot_grid_f(p, patch_grid, ei->array, 1); \
			RGBPic_dump_png(s, p); \
			} \
		snprintf(s,19999, "%s_%s.dat", ei->name, tag); \
		dump_floats(s, ei->array, patch_grid->npoints, 1); \
		} \
	}

OUTPUT_SKYMAP(ul_skymap, "upper_limit");
OUTPUT_SKYMAP(ul_freq_skymap, "upper_limit_frequency");
OUTPUT_SKYMAP(circ_ul_skymap, "circular_upper_limit");
OUTPUT_SKYMAP(circ_ul_freq_skymap, "circular_upper_limit_frequency");
OUTPUT_SKYMAP(snr_skymap, "snr");
OUTPUT_SKYMAP(snr_ul_skymap, "snr_upper_limit");
OUTPUT_SKYMAP(snr_freq_skymap, "snr_frequency");
OUTPUT_SKYMAP(max_weight_skymap, "max_weight");
OUTPUT_SKYMAP(min_weight_skymap, "min_weight");
OUTPUT_SKYMAP(weight_loss_fraction_skymap, "weight_loss_fraction");
OUTPUT_SKYMAP(ks_skymap, "ks_value");
}

void outer_loop(void)
{
int pi, i, k, m;
int nchunks;
int nei;
POWER_SUM **ps, **ps_tmp;
EXTREME_INFO **ei;
int ps_tmp_len;
int count;
double gps_start=min_gps();
double gps_stop=max_gps()+1;
time_t start_time, end_time;
RGBPic *p;
PLOT *plot;

assign_veto();

nchunks=args_info.nchunks_arg*veto_free;
ps=do_alloc(nchunks, sizeof(*ps));
ps_tmp=do_alloc(nchunks, sizeof(*ps));

ei=do_alloc(args_info.nchunks_arg*(args_info.nchunks_arg-1)*(veto_free+1), sizeof(*ei));

fprintf(LOG, "nchunks: %d\n", args_info.nchunks_arg);
fprintf(LOG, "veto_free: %d\n", veto_free);

nei=0;

for(i=0;i<args_info.nchunks_arg;i++)
	for(k=1;k<=args_info.nchunks_arg-i;k++)
		for(m=-1;m<veto_free;m++) {
			if(m<0) {
				if(veto_free<=1)continue; /* if there is only one detector no reason to compute "all" twice */
				snprintf(s, 19999, "%d_%d_all", i, i+k-1);
				} else {
				snprintf(s, 19999, "%d_%d_%s", i, i+k-1, veto_info[m].name);
				}
			ei[nei]=allocate_extreme_info(s);
			ei[nei]->first_chunk=i;
			ei[nei]->nchunks=k;
			ei[nei]->veto_num=m;
			nei++;
			}

fprintf(LOG, "nei: %d\n", nei);

time(&start_time);

fprintf(stderr, "%d patches to process\n", patch_grid->npoints);
for(pi=0;pi<patch_grid->npoints;pi++) {
	fprintf(stderr, "%d\n", pi);
	generate_patch_templates(pi, &(ps[0]), &count);

	if(count<1)continue;

	for(i=1;i<nchunks;i++) {
		clone_templates(ps[0], count, &(ps[i]));
		}
	for(i=0;i<args_info.nchunks_arg;i++) {
		for(k=0;k<veto_free;k++) {
			accumulate_power_sums(ps[i*veto_free+k], count, gps_start+i*(gps_stop-gps_start)/args_info.nchunks_arg, gps_start+(i+1)*(gps_stop-gps_start)/args_info.nchunks_arg, veto_info[k].veto_mask);
			}
		}

	/* find largest strain and largest SNR candidates for this patch */
	for(i=0;i<nei;i++) {
		ps_tmp_len=0;
		for(k=0;k<ei[i]->nchunks;k++) {
			if(ei[i]->veto_num<0) {
				for(m=0;m<veto_free;m++) {
					ps_tmp[ps_tmp_len]=ps[k*veto_free+m];
					ps_tmp_len++;
					}
				} else {
				ps_tmp[ps_tmp_len]=ps[k*veto_free+ei[i]->veto_num];
				ps_tmp_len++;
				}
			}
		log_extremes(ei[i], pi, ps_tmp, ps_tmp_len, count);
		}

	for(i=0;i<nchunks;i++) {
		free_templates(ps[i], count);
		ps[i]=NULL;
		}
	}
time(&end_time);
if(end_time<start_time)end_time=start_time;
fprintf(stderr, "Patch speed: %f\n", patch_grid->npoints/(1.0*(end_time-start_time+1.0)));
fprintf(LOG, "Patch speed: %f\n", patch_grid->npoints/(1.0*(end_time-start_time+1.0)));
free(ps);
print_cache_stats();

fflush(DATA_LOG);
fflush(LOG);

if(patch_grid->max_n_dec<800){
	p=make_RGBPic(patch_grid->max_n_ra*(800/patch_grid->max_n_dec)+140, patch_grid->max_n_dec*(800/patch_grid->max_n_dec));
	} else {
	p=make_RGBPic(patch_grid->max_n_ra+140, patch_grid->max_n_dec);
	}

plot=make_plot(p->width, p->height);

fprintf(stderr, "%d points written into the data log\n", data_log_index);
fprintf(LOG, "%d points written into the data log\n", data_log_index);

fprintf(stderr, "Writing skymaps\n");

for(i=0;i<nei;i++)
	output_extreme_info(p, ei[i]);
}
