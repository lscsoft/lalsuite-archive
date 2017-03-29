#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
/* We need this define to get NAN values */
//#define __USE_ISOC99
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
#include "summing_context.h"
#include "jobs.h"

#define DIVERT_TEMPLATE_MARKED  (1<<0)
#define DIVERT_TEMPLATE_MARKED_UL  (1<<2)

#define DIVERT_TEMPLATE_WRITTEN (1<<7)

extern struct gengetopt_args_info args_info;

extern SKY_GRID *fine_grid, *patch_grid;

extern FILE * DATA_LOG, * LOG, *FILE_LOG, *DIVERT_LOG, *INPUT_TEMPLATE_LOG;

extern int first_bin, side_cut, nsegments, useful_bins;

extern DATASET *datasets;
extern int d_free;

extern int input_templates;
MUTEX input_template_mutex;

int data_log_index=0;

int write_data_log_header=1;
int write_diverted_log_header=1;

int total_diverted_count=0;

int divert_buffer_free=0;
int divert_buffer_size=0;
double divert_buffer_ul_threshold=0.0;
DIVERTED_ENTRY *divert_buffer=NULL;
float *divert_sort_buffer=NULL;
MUTEX divert_buffer_mutex;

int nei;
EXTREME_INFO **ei=NULL;
int nchunks;

typedef struct {
	unsigned int veto_mask;
	char *name;
	} VETO_INFO;

/* we should have less than 5 detectors. Also, veto has only 8 bits some of which are used already - rework the scheme when one needs to accomodate more detectors */
VETO_INFO veto_info[4];
int veto_free=0;

void assign_detector_veto(void)
{
int i,k, m;
memset(veto_info, 0, 4*sizeof(*veto_info));

fprintf(LOG, "split_ifos: %s\n", args_info.split_ifos_arg ? "yes" : "no");
if(!args_info.split_ifos_arg) {
	/* do not split detectors, one veto entry only */
	veto_free=1;
	veto_info[0].veto_mask=(((1<<4)-1)<<4);
	return;
	}

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
		fprintf(stderr, "veto_info: %d \"%s\" %02x\n", m, veto_info[m].name, veto_info[m].veto_mask);
		fprintf(LOG, "veto_info: %d \"%s\" %02x\n", m, veto_info[m].name, veto_info[m].veto_mask);
		}
	for(i=0;i<datasets[k].free;i++) {
		datasets[k].sft_veto[i]=(datasets[k].sft_veto[i] & ~ (((1<<4)-1)<<4) ) | veto_info[m].veto_mask;
		}
	}
}

typedef struct {
	double weight;
	int dataset;
	int segment;
	} WEIGHT_INFO;

int weight_cmp(WEIGHT_INFO *w1, WEIGHT_INFO *w2)
{
if(w1->weight<w2->weight)return -1;
if(w1->weight>w2->weight)return 1;
return 0;
}

/* instead of having per-skypatch 

 we simply discard a SFTs which contribute less than a certain fraction of weight in each dataset */
void assign_cutoff_veto(void)
{
int i,k,m;
WEIGHT_INFO *w;
double total_weight, accum_weight, a;
DATASET *d;

w=do_alloc(nsegments, sizeof(*w));
m=0;
total_weight=0.0;
for(k=0;k<d_free;k++) {
	d=&(datasets[k]);
	for(i=0;i<d->free;i++) {
		if(d->sft_veto[i] && d->sft_veto[i]!=3)continue;

		w[m].dataset=k;
		w[m].segment=i;
		a=d->expTMedians[i]*d->weight;
		total_weight+=a;
		w[m].weight=a;
		m++;
		}
	}

if(m<1) {
	fprintf(stderr, "No clean SFTs found, skipping veto\n");
	fprintf(LOG, "No clean SFTs found, skipping veto\n");
	return;
	}

//fprintf(stderr, "%d %d\n", m, nsegments);
qsort(w, m, sizeof(*w), weight_cmp);
//fprintf(stderr, "%g %g %g ... %g %g %g\n", w[0].weight, w[1].weight, w[2].weight, w[m-3].weight, w[m-2].weight, w[m-1].weight);accum_weight=0;

accum_weight=0;
for(i=0;i<m;i++) {
	accum_weight+=w[i].weight;
	if(accum_weight>args_info.weight_cutoff_fraction_arg*total_weight)break;
	datasets[w[i].dataset].sft_veto[w[i].segment]=3;
	}
accum_weight-=w[i].weight;
fprintf(stderr, "Vetoed %d sfts (out of %d, %f ratio) with %g weight out of %g total weight (%f fraction)\n",
	i, nsegments, (i*1.0)/nsegments, accum_weight, total_weight, accum_weight/total_weight);
fprintf(LOG, "Vetoed %d sfts (out of %d, %f ratio) with %g weight out of %g total weight (%f fraction)\n",
	i, nsegments, (i*1.0)/nsegments,  accum_weight, total_weight, accum_weight/total_weight);

for(k=0;k<d_free;k++) {
	d=&(datasets[k]);
	m=0;
	for(i=0;i<d->free;i++) {
		if(d->sft_veto[i])continue;
		m++;
		}
	fprintf(stderr, "Dataset %s has %d clean sfts (%f ratio)\n", d->name, m, (1.0*m)/d->free);
	}
free(w);
}

void assign_per_dataset_cutoff_veto(void)
{
int i,k,m;
WEIGHT_INFO *w;
double total_weight, accum_weight, a;
DATASET *d;

w=do_alloc(nsegments, sizeof(*w));
for(k=0;k<d_free;k++) {
	d=&(datasets[k]);
	m=0;
	total_weight=0.0;
	for(i=0;i<d->free;i++) {
		if(d->sft_veto[i] && d->sft_veto[i]!=3)continue;

		w[m].dataset=k;
		w[m].segment=i;
		a=d->expTMedians[i]*d->weight;
		total_weight+=a;
		w[m].weight=a;
		m++;
		}

	if(m<1)continue;

	//fprintf(stderr, "%d %d\n", m, nsegments);
	qsort(w, m, sizeof(*w), weight_cmp);
	//fprintf(stderr, "%g %g %g ... %g %g %g\n", w[0].weight, w[1].weight, w[2].weight, w[m-3].weight, w[m-2].weight, w[m-1].weight);accum_weight=0;
	
	accum_weight=0;
	for(i=0;i<m;i++) {
		accum_weight+=w[i].weight;
		if(accum_weight>args_info.per_dataset_weight_cutoff_fraction_arg*total_weight)break;
		datasets[w[i].dataset].sft_veto[w[i].segment]=3;
		}
	accum_weight-=w[i].weight;
	fprintf(stderr, "Dataset %s: vetoed %d sfts (out of %d, %f ratio) with %g weight out of %g total weight (%f fraction)\n", d->name,
		i, d->free, (i*1.0)/d->free, accum_weight, total_weight, accum_weight/total_weight);
	fprintf(LOG, "Dataset %s: vetoed %d sfts (out of %d, %f ratio) with %g weight out of %g total weight (%f fraction)\n", d->name,
		i, d->free, (i*1.0)/d->free,  accum_weight, total_weight, accum_weight/total_weight);
	}

for(k=0;k<d_free;k++) {
	d=&(datasets[k]);
	m=0;
	for(i=0;i<d->free;i++) {
		if(d->sft_veto[i])continue;
		m++;
		}
	fprintf(stderr, "Dataset %s has %d clean sfts (%f ratio)\n", d->name, m, (1.0*m)/d->free);
	fprintf(LOG, "Dataset %s has %d clean sfts (%f ratio)\n", d->name, m, (1.0*m)/d->free);
	}
free(w);
}

void replay_diverted_entry(DIVERTED_ENTRY *de, EXTREME_INFO *ei)
{
int skyband;
int fshift_count=args_info.nfshift_arg;
int pi=de->pi;

skyband=de->ti.skyband;

thread_mutex_lock(&(ei->mutex));

if(de->pstats_ul.max_weight_loss_fraction>=1) {
	ei->band_masked_count[skyband]++;
	return;
	}
if(de->pstats_snr.max_weight_loss_fraction>=1) {
	ei->band_masked_count[skyband]++;
	return;
	}
ei->band_valid_count[skyband]++;

#define RFILL_EXTRA_PARAMS(target) {\
	target.ra=de->ti.ra; \
	target.dec=de->ti.dec; \
	target.spindown=de->ti.spindown; \
	target.fdotdot=de->ti.fdotdot; \
	target.freq_modulation_freq=de->ti.freq_modulation_freq; \
	target.freq_modulation_depth=de->ti.freq_modulation_depth; \
	target.freq_modulation_phase=de->ti.freq_modulation_phase; \
	target.frequency=((double)(de->ti.snr_subbin)/fshift_count+first_bin+side_cut)/args_info.sft_coherence_time_arg; \
	}

#define RFILL_POINT_STATS(target, source)	{\
	memcpy(&(target), &(source), sizeof(target)); \
	RFILL_EXTRA_PARAMS(target); \
	}

/* Note: approximation - we only update MAX and MIN for highest_ul and highest_snr, not all freq_shifts as is for non-diverted templates */
/* This introduces a possible non-repeatability in multi-threaded execution affecting secondary quantities like max_m4 */
	
#define RUPDATE_MAX(target, field) {\
	if(de->pstats_ul.field>target.field) { \
		target.field=de->pstats_ul.field;\
		} \
	if(de->pstats_snr.field>target.field) { \
		target.field=de->pstats_snr.field;\
		} \
	}

#define RUPDATE_MIN(target, field) {\
	if(de->pstats_ul.field<target.field) { \
		target.field=de->pstats_ul.field;\
		} \
	if(de->pstats_snr.field<target.field) { \
		target.field=de->pstats_snr.field;\
		} \
	}

if(ei->band_info[skyband].max_weight<0) {
	memcpy(&(ei->band_info[skyband]), &(de->pstats_ul), sizeof(de->pstats_ul));
	RFILL_EXTRA_PARAMS(ei->band_info[skyband].highest_ul);
	RFILL_EXTRA_PARAMS(ei->band_info[skyband].highest_circ_ul);
	RFILL_EXTRA_PARAMS(ei->band_info[skyband].highest_snr);
	RFILL_EXTRA_PARAMS(ei->band_info[skyband].highest_ks);
	} else {

	ei->band_info[skyband].ntemplates+=de->pstats_ul.ntemplates*fshift_count;

	if(de->pstats_ul.highest_ul.ul>ei->band_info[skyband].highest_ul.ul) {
		RFILL_POINT_STATS(ei->band_info[skyband].highest_ul, de->pstats_ul.highest_ul);
		}

		/* circ_ul and highest_ks fields are not updated by this code */
		
// 	if(pstats->highest_circ_ul.ul>ei->band_info[skyband].highest_circ_ul.ul) {
// 		RFILL_POINT_STATS(ei->band_info[skyband].highest_circ_ul, pstats->highest_circ_ul);
// 		}

	if(de->pstats_snr.highest_snr.snr>ei->band_info[skyband].highest_snr.snr) {
		RFILL_POINT_STATS(ei->band_info[skyband].highest_snr, de->pstats_snr.highest_snr);
		}

// 	if(pstats->highest_ks.ks_value>ei->band_info[skyband].highest_ks.ks_value) {
// 		RFILL_POINT_STATS(ei->band_info[skyband].highest_ks, pstats->highest_ks);
// 		}

	RUPDATE_MAX(ei->band_info[skyband], max_weight);
	RUPDATE_MIN(ei->band_info[skyband], min_weight);
	RUPDATE_MAX(ei->band_info[skyband], max_weight_loss_fraction);

	RUPDATE_MAX(ei->band_info[skyband], max_m1_neg);
	RUPDATE_MIN(ei->band_info[skyband], min_m1_neg);
	RUPDATE_MAX(ei->band_info[skyband], max_m3_neg);
	RUPDATE_MIN(ei->band_info[skyband], min_m3_neg);
	RUPDATE_MAX(ei->band_info[skyband], max_m4);
	RUPDATE_MIN(ei->band_info[skyband], min_m4);
	}

thread_mutex_unlock(&(ei->mutex));

/* During main run operation the skymaps are not enabled (NULL) to save space and computation */

#define FILL_SKYMAP(skymap, value)	if(ei->skymap!=NULL)ei->skymap[pi]=value;

if(ei->ul_skymap!=NULL && (de->pstats_ul.highest_ul.ul>ei->ul_skymap[pi])) {
	FILL_SKYMAP(ul_skymap, de->pstats_ul.highest_ul.ul);
	FILL_SKYMAP(ul_freq_skymap, (de->pstats_ul.highest_ul.bin)/args_info.sft_coherence_time_arg+(double)de->highest_ul_j/fshift_count);
	}

/* Not filled in properly
FILL_SKYMAP(circ_ul_skymap, pstats_accum.highest_circ_ul.ul);
FILL_SKYMAP(circ_ul_freq_skymap, (pstats_accum.highest_circ_ul.bin)/args_info.sft_coherence_time_arg+ps[0][highest_circ_ul_idx].freq_shift);
*/

if(ei->snr_skymap!=NULL && (de->pstats_snr.highest_snr.snr>ei->snr_skymap[pi])) {
	FILL_SKYMAP(snr_skymap, de->pstats_snr.highest_snr.snr);
	FILL_SKYMAP(snr_ul_skymap, de->pstats_snr.highest_snr.ul);
	FILL_SKYMAP(snr_freq_skymap, (de->pstats_snr.highest_snr.bin)/args_info.sft_coherence_time_arg+(double)de->highest_snr_j/fshift_count);
	}

// Not filled in properly
// FILL_SKYMAP(max_weight_skymap, pstats_accum.max_weight);
// FILL_SKYMAP(min_weight_skymap, pstats_accum.min_weight);
// FILL_SKYMAP(weight_loss_fraction_skymap, pstats_accum.max_weight_loss_fraction);
// 
// FILL_SKYMAP(ks_skymap, pstats_accum.highest_ks.ks_value);
}

void trim_divert_buffer(void)
{
int i,j;
if(divert_buffer_free<args_info.divert_ul_count_limit_arg)return;

for(i=0;i<divert_buffer_free;i++)divert_sort_buffer[i]=divert_buffer[i].ti.ul;
sort_floats(divert_sort_buffer, divert_buffer_free);
fprintf(stderr, "trim_divert_buffer 1: %g %g %g\n", divert_sort_buffer[0], divert_sort_buffer[1], divert_sort_buffer[2]);
divert_buffer_ul_threshold=divert_sort_buffer[divert_buffer_free-args_info.divert_ul_count_limit_arg];

j=0;
for(i=0;i<divert_buffer_free;i++) {
	if(divert_buffer[i].ti.ul>divert_buffer_ul_threshold) {
		if(j<i) {
			divert_buffer[j]=divert_buffer[i];
			}
		j++;
		continue;
		}
	
	replay_diverted_entry(&(divert_buffer[i]), ei[divert_buffer[i].ei_idx]);
	
	}
divert_buffer_free=j;
}

MUTEX data_logging_mutex;

void pstats_log_extremes(SUMMING_CONTEXT *ctx, POWER_SUM_STATS *tmp_pstat, POWER_SUM **ps, int count, EXTREME_INFO *ei, int pi)
{
int i,j, i_start, highest_ul_j, highest_snr_j;
POWER_SUM_STATS pstats, pstats_accum, *tp;
int highest_ul_idx=0;
int highest_circ_ul_idx=0;
int highest_snr_idx=0;
int skyband;
char *diverted;
double average_S;
TEMPLATE_INFO ti;
int fshift_count=args_info.nfshift_arg; /* number of frequency offsets */

diverted=ctx->diverted;

memset(&pstats_accum, 0, sizeof(pstats_accum));
pstats_accum.max_weight=-1;

average_S=tmp_pstat[0].highest_circ_ul.S;
for(i=1;i<count;i++)average_S+=tmp_pstat[i].highest_circ_ul.S;
average_S/=count;

/* process SNR thresholds first */
if(total_diverted_count<args_info.divert_count_limit_arg) {
	for(i=0;i<count;i++) {
		if(diverted[i])continue;
		
		tp=&(tmp_pstat[i]);
		
		if(tp->max_weight_loss_fraction>=1) {
			continue;
			}
		
		if(args_info.divert_snr_arg>0 && (tp->highest_snr.snr>args_info.divert_snr_arg)) {
			
			/* divert all frequencies in the template */
			i_start=fshift_count*(i/fshift_count);
			for(j=0;j<fshift_count;j++)
				diverted[i_start+j]|=DIVERT_TEMPLATE_MARKED;
			
			continue;
			}
			
		}
	}

/* Process upper limit toplist */
if((ei->last_chunk-ei->first_chunk+1==args_info.nchunks_arg) && ((ei->veto_num==-1) || (veto_free<=1)))
for(i=0;i<count;i++) {
	if(diverted[i])continue;
	
	tp=&(tmp_pstat[i]);
	
	if(tp->max_weight_loss_fraction>=1) {
		continue;
		}
	
	if(((args_info.divert_ul_arg>0 && (tp->highest_ul.ul>args_info.divert_ul_arg)) ||
	(args_info.divert_ul_rel_arg>0 && (tp->highest_ul.ul>args_info.divert_ul_rel_arg*average_S))) && 
		(tp->highest_ul.ul>divert_buffer_ul_threshold)
		) {
		
		
		/* divert all frequencies in the template */
		i_start=fshift_count*(i/fshift_count);
		for(j=0;j<fshift_count;j++)
			diverted[i_start+j]|=DIVERT_TEMPLATE_MARKED_UL;
		}
	}
	
	
thread_mutex_lock(&(ei->mutex));
for(i=0;i<count;i++) {
	skyband=ps[0][i].skyband;
	
	if(diverted[i] & DIVERT_TEMPLATE_MARKED) {
		ei->band_diverted_count[skyband]++;
		continue;
		}

	if((ei->last_chunk-ei->first_chunk+1==args_info.nchunks_arg) && ((ei->veto_num==-1) || (veto_free<=1)) && (diverted[i] & DIVERT_TEMPLATE_MARKED_UL)) {
		continue;
		}
		
	memcpy(&pstats, &(tmp_pstat[i]), sizeof(POWER_SUM_STATS));
	
	if(pstats.max_weight_loss_fraction>=1) {
		ei->band_masked_count[skyband]++;
		continue;
		}
	ei->band_valid_count[skyband]++;

	#define FILL_EXTRA_PARAMS(target) {\
		target.ra=ps[0][i].ra; \
		target.dec=ps[0][i].dec; \
		target.spindown=ps[0][i].spindown; \
		target.fdotdot=ps[0][i].fdotdot; \
		target.freq_modulation_freq=ps[0][i].freq_modulation_freq; \
		target.freq_modulation_depth=ps[0][i].freq_modulation_depth; \
		target.freq_modulation_phase=ps[0][i].freq_modulation_phase; \
		target.frequency=(double)ps[0][i].freq_shift+((target).bin+first_bin+side_cut)/args_info.sft_coherence_time_arg; \
		}

	#define FILL_POINT_STATS(target, source)	{\
		memcpy(&(target), &(source), sizeof(target)); \
		FILL_EXTRA_PARAMS(target); \
		}

	#define UPDATE_MAX(target, field) {\
		if(pstats.field>target.field) { \
			target.field=pstats.field;\
			} \
		}

	#define UPDATE_MIN(target, field) {\
		if(pstats.field<target.field) { \
			target.field=pstats.field;\
			} \
		}

	if(ei->band_info[skyband].max_weight<0) {
		memcpy(&(ei->band_info[skyband]), &pstats, sizeof(pstats));
		FILL_EXTRA_PARAMS(ei->band_info[skyband].highest_ul);
		FILL_EXTRA_PARAMS(ei->band_info[skyband].highest_circ_ul);
		FILL_EXTRA_PARAMS(ei->band_info[skyband].highest_snr);
		FILL_EXTRA_PARAMS(ei->band_info[skyband].highest_ks);
		} else {

		ei->band_info[skyband].ntemplates+=pstats.ntemplates;

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
	
		UPDATE_MAX(ei->band_info[skyband], max_weight);
		UPDATE_MIN(ei->band_info[skyband], min_weight);
		UPDATE_MAX(ei->band_info[skyband], max_weight_loss_fraction);

		UPDATE_MAX(ei->band_info[skyband], max_m1_neg);
		UPDATE_MIN(ei->band_info[skyband], min_m1_neg);
		UPDATE_MAX(ei->band_info[skyband], max_m3_neg);
		UPDATE_MIN(ei->band_info[skyband], min_m3_neg);
		UPDATE_MAX(ei->band_info[skyband], max_m4);
		UPDATE_MIN(ei->band_info[skyband], min_m4);

		}


	/* No need to fill extra parameters as results are printed in this function */
	if(pstats_accum.max_weight<0) {
		memcpy(&pstats_accum, &pstats, sizeof(pstats));
		continue;
		}

	pstats_accum.ntemplates+=pstats.ntemplates;

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

	UPDATE_MAX(pstats_accum, max_weight);
	UPDATE_MIN(pstats_accum, min_weight);
	UPDATE_MAX(pstats_accum, max_weight_loss_fraction);

	UPDATE_MAX(pstats_accum, max_m1_neg);
	UPDATE_MIN(pstats_accum, min_m1_neg);
	UPDATE_MAX(pstats_accum, max_m3_neg);
	UPDATE_MIN(pstats_accum, min_m3_neg);
	UPDATE_MAX(pstats_accum, max_m4);
	UPDATE_MIN(pstats_accum, min_m4);

	}
//free_partial_power_sum_F(pps);
thread_mutex_unlock(&(ei->mutex));

thread_mutex_lock(&data_logging_mutex);

if(write_diverted_log_header) {
	write_diverted_log_header=0;
	fprintf(LOG, "diverted_log_structure_size: %ld\n", sizeof(ti));
	}

for(i=0;i<count;i+=fshift_count) {
	if((diverted[i] & DIVERT_TEMPLATE_WRITTEN))continue;
	if(!(diverted[i] & DIVERT_TEMPLATE_MARKED))continue;
	
	tp=&(tmp_pstat[i]);
	
	ti.skyband=ps[0][i].skyband;
	
	ti.snr=tp->highest_snr.snr;
	ti.snr_subbin=tp->highest_snr.bin*fshift_count;
	highest_snr_j=0;
	
	for(j=1;j<fshift_count;j++)
		if(ti.snr<tp[j].highest_snr.snr) {
			ti.snr=tp[j].highest_snr.snr;
			ti.snr_subbin=tp[j].highest_snr.bin*fshift_count+j;
			highest_snr_j=j;
			}
			
	ti.ul=tp->highest_ul.ul;
	ti.circ_ul=tp->highest_circ_ul.ul;
	
	ti.freq_modulation_freq=ps[0][i].freq_modulation_freq;
	ti.freq_modulation_phase=ps[0][i].freq_modulation_phase;
	ti.freq_modulation_depth=ps[0][i].freq_modulation_depth;
	
	ti.spindown=ps[0][i].spindown;
	ti.fdotdot=ps[0][i].fdotdot;
	ti.ra=ps[0][i].ra;
	ti.dec=ps[0][i].dec;
	
	/* This just shows the start of the band to analyze 
	 * The actual frequency where the maximum is achieved can be different for SNR and UL.
	 */
	/* ti.frequency=(double)ps[0][i].freq_shift+ti.first_bin/args_info.sft_coherence_time_arg; */
	
	ti.first_chunk=ei->first_chunk;
	ti.last_chunk=ei->last_chunk;
	ti.veto_num=ei->veto_num;
	

	if(diverted[i]==DIVERT_TEMPLATE_MARKED_UL) {
		highest_ul_j=0;
		for(j=1;j<fshift_count;j++)
			if(ti.ul<tp[j].highest_ul.ul) {
				ti.ul=tp[j].highest_ul.ul;
				highest_ul_j=j;
				}
	
		
		
		thread_mutex_lock(&divert_buffer_mutex);
		
		if(divert_buffer_free>=divert_buffer_size)trim_divert_buffer();
		
		memcpy(&(divert_buffer[divert_buffer_free].ti), &ti, sizeof(ti));
		memcpy(&(divert_buffer[divert_buffer_free].pstats_ul), &(tp[highest_ul_j]), sizeof(*tp));
		memcpy(&(divert_buffer[divert_buffer_free].pstats_snr), &(tp[highest_snr_j]), sizeof(*tp));
		divert_buffer[divert_buffer_free].ei_idx=ei->idx;
		divert_buffer[divert_buffer_free].pi=pi;
		divert_buffer[divert_buffer_free].highest_snr_j=highest_snr_j;
		divert_buffer[divert_buffer_free].highest_ul_j=highest_ul_j;
		divert_buffer_free++;		
		
		thread_mutex_unlock(&divert_buffer_mutex);
		diverted[i]|=DIVERT_TEMPLATE_WRITTEN;
		total_diverted_count++;
		continue;
		}
	
	fwrite(&ti, sizeof(ti), 1, DIVERT_LOG);
	diverted[i]|=DIVERT_TEMPLATE_WRITTEN;
	total_diverted_count++;
	}

if(write_data_log_header) {
	write_data_log_header=0;
	/* we write this into the main log file so that data.log files can simply be concatenated together */
	fprintf(LOG, "data_log: kind label index set pi pps_count template_count first_bin min_gps max_gps skyband frequency spindown fdotdot freq_modulation_freq freq_modulation_depth freq_modulation_phase ra dec iota psi snr ul ll M S ks_value ks_count m1_neg m3_neg m4 frequency_bin max_weight weight_loss_fraction max_ks_value max_m1_neg min_m1_neg max_m3_neg min_m3_neg max_m4 min_m4 max_weight_loss_fraction\n");
	}

/* now that we know extreme points go and characterize them */
#define WRITE_POINT(psum, pstat, kind)	{\
	fprintf(DATA_LOG, "%s \"%s\" %d %s %d %d %d %d %lf %lf %d %lf %lg %lg %lg %lg %lg %lf %lf %lf %lf %lf %lg %lg %lg %lg %lf %d %lf %lf %lf %d %lg %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", \
		kind, \
		args_info.label_arg, \
		data_log_index, \
		ei->name, \
		pi, \
		count, \
		pstats_accum.ntemplates, \
		first_bin+side_cut, \
		psum.min_gps, \
		psum.max_gps, \
		psum.skyband, \
		(pstat.bin+first_bin+side_cut)/args_info.sft_coherence_time_arg+psum.freq_shift, \
		psum.spindown, \
		psum.fdotdot, \
		psum.freq_modulation_freq, \
		psum.freq_modulation_depth, \
		psum.freq_modulation_phase, \
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
		pstat.m1_neg, \
		pstat.m3_neg, \
		pstat.m4, \
		pstat.bin, \
		pstat.max_weight, \
		pstat.weight_loss_fraction, \
		pstats_accum.highest_ks.ks_value, \
		pstats_accum.max_m1_neg, \
		pstats_accum.min_m1_neg, \
		pstats_accum.max_m3_neg, \
		pstats_accum.min_m3_neg, \
		pstats_accum.max_m4, \
		pstats_accum.min_m4, \
		pstats_accum.max_weight_loss_fraction \
		); data_log_index++; }

if(args_info.output_cache_arg) {
	WRITE_POINT(ps[0][highest_ul_idx], pstats_accum.highest_ul, "ul");
	WRITE_POINT(ps[0][highest_circ_ul_idx], pstats_accum.highest_circ_ul, "circ");
	}
if((pstats_accum.highest_snr.snr>args_info.min_candidate_snr_arg) &
   (pstats_accum.highest_snr.bin>=args_info.tail_veto_arg) &
   (pstats_accum.highest_snr.bin<(useful_bins-args_info.tail_veto_arg)))WRITE_POINT(ps[0][highest_snr_idx], pstats_accum.highest_snr, "snr");

thread_mutex_unlock(&data_logging_mutex);

FILL_SKYMAP(ul_skymap, pstats_accum.highest_ul.ul);
FILL_SKYMAP(ul_freq_skymap, (pstats_accum.highest_ul.bin)/args_info.sft_coherence_time_arg+ps[0][highest_ul_idx].freq_shift);

FILL_SKYMAP(circ_ul_skymap, pstats_accum.highest_circ_ul.ul);
FILL_SKYMAP(circ_ul_freq_skymap, (pstats_accum.highest_circ_ul.bin)/args_info.sft_coherence_time_arg+ps[0][highest_circ_ul_idx].freq_shift);

FILL_SKYMAP(snr_skymap, pstats_accum.highest_snr.snr);
FILL_SKYMAP(snr_ul_skymap, pstats_accum.highest_snr.ul);
FILL_SKYMAP(snr_freq_skymap, (pstats_accum.highest_snr.bin)/args_info.sft_coherence_time_arg+ps[0][highest_snr_idx].freq_shift);

FILL_SKYMAP(max_weight_skymap, pstats_accum.max_weight);
FILL_SKYMAP(min_weight_skymap, pstats_accum.min_weight);
FILL_SKYMAP(weight_loss_fraction_skymap, pstats_accum.max_weight_loss_fraction);

FILL_SKYMAP(ks_skymap, pstats_accum.highest_ks.ks_value);
}

void log_extremes(SUMMING_CONTEXT *ctx, POWER_SUM_STATS *tmp_pstat, EXTREME_INFO *ei, int pi, POWER_SUM **ps, int nchunks, int count)
{
PARTIAL_POWER_SUM_F *pps=ctx->log_extremes_pps;
int i, k;

//fprintf(stderr, "count=%d\n", count);
//pps=allocate_partial_power_sum_F(useful_bins, 1);

for(i=0;i<count;i++) {
	zero_partial_power_sum_F(pps);
	for(k=0;k<nchunks;k++) {
#if MANUAL_SSE
		sse_accumulate_partial_power_sum_F(pps, (ps[k][i].pps));
#else
		accumulate_partial_power_sum_F(pps, (ps[k][i].pps));
#endif
		}
	power_sum_stats(pps, &(tmp_pstat[i]));
	
	if(args_info.dump_power_sums_arg) {
		thread_mutex_lock(&data_logging_mutex);
		
		fprintf(DATA_LOG, "power_sum %s %d %d %lf %lf %lf %lg ", ei->name, pi, first_bin+side_cut, ps[0][i].ra, ps[0][i].dec, ps[0][i].freq_shift, ps[0][i].spindown);
		dump_partial_power_sum_F(DATA_LOG, pps);
		fprintf(DATA_LOG, "\n");
		
		thread_mutex_unlock(&data_logging_mutex);
		}

	}
	
pstats_log_extremes(ctx, tmp_pstat, ps, count, ei, pi);
}

static char s[20000];

EXTREME_INFO * allocate_extreme_info(char *name)
{
EXTREME_INFO *ei;
int i;

ei=do_alloc(1, sizeof(*ei));
memset(ei, 0, sizeof(*ei));

ei->name=strdup(name);

thread_mutex_init(&(ei->mutex));

if(args_info.compute_skymaps_arg) {
	if(input_templates>=0) {
		fprintf(stderr, "compute-skymaps is incompatible with binary-template-file\n");
		exit(-1);
		}
	
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

ei->band_valid_count=do_alloc(fine_grid->nbands, sizeof(*ei->band_valid_count));
ei->band_masked_count=do_alloc(fine_grid->nbands, sizeof(*ei->band_masked_count));
ei->band_diverted_count=do_alloc(fine_grid->nbands, sizeof(*ei->band_diverted_count));
memset(ei->band_valid_count, 0, fine_grid->nbands*sizeof(*ei->band_valid_count));
memset(ei->band_masked_count, 0, fine_grid->nbands*sizeof(*ei->band_masked_count));
memset(ei->band_diverted_count, 0, fine_grid->nbands*sizeof(*ei->band_diverted_count));

for(i=0;i<fine_grid->nbands;i++) {
	ei->band_info[i].max_weight=-1;
	}

return ei;
}

void free_extreme_info(EXTREME_INFO *ei)
{

free(ei->name);
free(ei->band_info);
free(ei->band_valid_count);
free(ei->band_masked_count);

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

fprintf(LOG, "tag: kind label skyband skyband_name set first_bin frequency spindown fdotdot freq_modulation_freq freq_modulation_depth freq_modulation_phase ra dec iota psi snr ul ll M S ks_value ks_count m1_neg m3_neg m4 frequency_bin max_weight weight_loss_fraction max_ks_value max_m1_neg min_m1_neg max_m3_neg min_m3_neg max_m4 min_m4 max_weight_loss_fraction valid_count masked_count diverted_count template_count\n");

/* now that we know extreme points go and characterize them */
#define WRITE_SKYBAND_POINT(pstat, kind)	\
	fprintf(LOG, "band_info: %s \"%s\" %d %s %s %d %lf %lg %lg %lg %lg %lg %lf %lf %lf %lf %lf %lg %lg %lg %lg %lf %d %lf %lf %lf %d %lg %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %d %d\n", \
		kind, \
		args_info.label_arg, \
		skyband, \
		fine_grid->band_name[skyband], \
		ei->name, \
		first_bin+side_cut, \
		pstat.frequency, \
		pstat.spindown, \
		pstat.fdotdot, \
		pstat.freq_modulation_freq, \
		pstat.freq_modulation_depth, \
		pstat.freq_modulation_phase, \
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
		pstat.m1_neg, \
		pstat.m3_neg, \
		pstat.m4, \
		pstat.bin, \
		pstat.max_weight, \
		pstat.weight_loss_fraction, \
		ei->band_info[skyband].highest_ks.ks_value, \
		ei->band_info[skyband].max_m1_neg, \
		ei->band_info[skyband].min_m1_neg, \
		ei->band_info[skyband].max_m3_neg, \
		ei->band_info[skyband].min_m3_neg, \
		ei->band_info[skyband].max_m4, \
		ei->band_info[skyband].min_m4, \
		ei->band_info[skyband].max_weight_loss_fraction, \
		ei->band_valid_count[skyband], \
		ei->band_masked_count[skyband], \
		ei->band_diverted_count[skyband], \
		ei->band_info[skyband].ntemplates \
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

void create_segments(EXTREME_INFO ***out_ei, int *out_nei)
{
int i, k, m, nei;
EXTREME_INFO **ei;
ei=do_alloc(args_info.nchunks_arg*(args_info.nchunks_arg+1)*(veto_free+1)/2, sizeof(*ei));

fprintf(LOG, "nchunks: %d\n", args_info.nchunks_arg);
fprintf(LOG, "nchunks refinement: %d\n", args_info.nchunks_refinement_arg);
fprintf(LOG, "min nchunks: %d\n", args_info.min_nchunks_arg);
fprintf(LOG, "veto_free: %d\n", veto_free);

nei=0;

for(i=0;i<args_info.nchunks_arg;i+=args_info.nchunks_refinement_arg)
	for(k=0;k< args_info.nchunks_arg-i;k+=args_info.nchunks_refinement_arg) {
		if(k+1<args_info.min_nchunks_arg)continue;
		
		for(m=-1;m<veto_free;m++) {
			if(m<0) {
				if((veto_free<=1) && args_info.split_ifos_arg)continue; /* if there is only one detector no reason to compute "all" twice */
				snprintf(s, 19999, "%d_%d_all", i/args_info.nchunks_refinement_arg, (i+k)/args_info.nchunks_refinement_arg);
				} else {
				if(!args_info.split_ifos_arg)continue; /* combine data from all detectors */
				snprintf(s, 19999, "%d_%d_%s", i/args_info.nchunks_refinement_arg, (i+k)/args_info.nchunks_refinement_arg, veto_info[m].name);
				}
			ei[nei]=allocate_extreme_info(s);
			ei[nei]->first_chunk=i;
			ei[nei]->last_chunk=i+k;
			ei[nei]->veto_num=m;
			ei[nei]->idx=nei;
			nei++;
			}
		}

*out_nei=nei;
*out_ei=ei;
}

SUMMING_CONTEXT **summing_contexts=NULL;
struct {
	POWER_SUM **ps;
	POWER_SUM **ps_tmp;
	float *temp;
	} *cruncher_contexts=NULL;
int n_contexts=0;


double gps_start;
double gps_stop;

#if MANUAL_SSE
#define MODE(a)	(args_info.sse_arg ? (sse_ ## a) : (a) )
#else
#define MODE(a)	(a)
#endif

extern ALIGNMENT_COEFFS *alignment_grid;
extern int alignment_grid_free;

void log_extremes_viterbi(SUMMING_CONTEXT *ctx, int pi, POWER_SUM **ps, int count)
{
int i,j,k,m,r, lm;
POINT_STATS pst;
POWER_SUM_STATS *stats;
float *tmp;
float *tmp2a, *tmp2b, *tmp_min_weight, *tmp_max_weight;
float *p1, *p2, *p3, *p4, *p5;
float min_weight, max_weight, inv_weight, weight, *total_weight;
int fshift_count=args_info.nfshift_arg; /* number of frequency offsets */
int shift;
long tmp_stride=(useful_bins+(ALIGNMENT-1)) & (~(ALIGNMENT-1));
long tmp_size;

if(args_info.filter_lines_arg) {
	/* To fix this we would need to pass and process per-frequency weight arrays. This needs to be done carefully to maintain efficiency 
	 * The code will work as is if this is disabled, using an approximation to true weight. 
	 * But the checks that enough weight was accumulated to compute power will not work */
	fprintf(stderr, "*** ERROR: viterbi filtering is incompatible with filter-lines=1\n");
	exit(-1);
	}

	/* size of tmp array */
tmp_size=args_info.nchunks_arg*veto_free*fshift_count*tmp_stride*sizeof(float);
	/* size of tmp2 arrays */
tmp_size+=2*(tmp_stride*fshift_count)*sizeof(float);
	/* size of stats array */
tmp_size+=nei*count*sizeof(*stats)+ALIGNMENT;
	/* sizes of tmp_min_weight and tmp_max_weight arrays */
tmp_size+=2*args_info.nchunks_arg*veto_free*fshift_count*sizeof(float)+2*ALIGNMENT;
	/* sizes of total_weight array */
tmp_size+=fshift_count*sizeof(float)+2*ALIGNMENT;
	/* size of diverted array */
tmp_size+=count*sizeof(*ctx->diverted)+2*ALIGNMENT;

if(ctx->log_extremes_pstats_scratch_size<tmp_size) {
	free(ctx->log_extremes_pstats_scratch);

	ctx->log_extremes_pstats_scratch_size=tmp_size;

	ctx->log_extremes_pstats_scratch=do_alloc(1, tmp_size);
	
	p1=(float *)ctx->log_extremes_pstats_scratch;
	PRAGMA_IVDEP
	for(i=0;i<(ctx->log_extremes_pstats_scratch_size/sizeof(*p1));i++)p1[i]=NAN;
	
	fprintf(stderr, "Expanded log_extremes_pstats_scratch to %f MB nchunks=%d veto_free=%d count=%d nei=%d\n", ctx->log_extremes_pstats_scratch_size*1e-6, args_info.nchunks_arg, veto_free, count, nei);
	}

p1=(float *)ctx->log_extremes_pstats_scratch;
ctx->diverted=p1; p1=ALIGN_POINTER(p1+((count+3)>>2));
tmp=p1; p1=ALIGN_POINTER(p1+args_info.nchunks_arg*veto_free*fshift_count*tmp_stride);
tmp2a=p1; p1=ALIGN_POINTER(p1+tmp_stride*fshift_count);
tmp2b=p1; p1=ALIGN_POINTER(p1+tmp_stride*fshift_count);
stats=(POWER_SUM_STATS *)p1; p1=ALIGN_POINTER(p1+((nei*count*sizeof(*stats)+3)>>2));
tmp_min_weight=p1; p1=ALIGN_POINTER(p1+args_info.nchunks_arg*veto_free*fshift_count);
tmp_max_weight=p1; p1=ALIGN_POINTER(p1+args_info.nchunks_arg*veto_free*fshift_count);
total_weight=p1; p1=ALIGN_POINTER(p1+fshift_count);

/* Check that size was computed accurately */
if(((char *)p1)-ctx->log_extremes_pstats_scratch>ctx->log_extremes_pstats_scratch_size) {
	fprintf(stderr, "*** ERROR: log_extremes_pstats_scratch_size=%ld but need %ld memory\n", 
		ctx->log_extremes_pstats_scratch_size, ((char *)p1)-ctx->log_extremes_pstats_scratch);
	exit(-1);
	}
	
memset(ctx->diverted, 0, count*sizeof(*ctx->diverted));

for(i=0;i<nei;i++) {
	for(j=0;j<count;j++)
		prepare_power_sum_stats(&stats[i*count+j]);
	}

for(lm=0;lm<alignment_grid_free;lm++) {


	for(j=0;j<count;j+=fshift_count) {
		
	for(i=0;i<args_info.nchunks_arg;i++) {
		for(k=0;k<veto_free;k++) {
			for(shift=0;shift<fshift_count;shift++) {
 				MODE(compute_power)(ps[i*veto_free+k][j+shift].pps, &(alignment_grid[lm]), &(tmp[((i*veto_free+k)*fshift_count+shift)*tmp_stride]), &(tmp_min_weight[(i*veto_free+k)*fshift_count+shift]), &(tmp_max_weight[(i*veto_free+k)*fshift_count+shift]));
				
				}
			}
		}
		
	for(i=0;i<nei;i++) {
		max_weight=0;
		min_weight=0;
		memset(tmp2a, 0, sizeof(*tmp2a)*tmp_stride*fshift_count);
		memset(total_weight, 0, sizeof(*total_weight)*fshift_count);
		
		for(k=ei[i]->first_chunk;k<=ei[i]->last_chunk;k++) {
			/* Accumulate tmp2 from tmp using Viterbi-like algorithm */
			/* We need to do all sub-bin shifts in one go */
			if(k>ei[i]->first_chunk) {
				for(shift=0;shift<fshift_count;shift++) {
					p2=&(tmp2a[tmp_stride*shift]);
					
					if(shift>0)p3=&(tmp2a[tmp_stride*(shift-1)]);
						else p3=&(tmp2a[tmp_stride*(shift+fshift_count-1)-1]);
						
					if(shift<fshift_count-1)p4=&(tmp2a[tmp_stride*(shift+1)]);
						else p4=&(tmp2a[tmp_stride*(shift-fshift_count+1)+1]);
						
					p5=&(tmp2b[tmp_stride*shift]);
						
					if(shift==0)p5[0]=fmaxf(p2[0], p4[0]);
						else
						p5[0]=fmaxf(p2[0], fmaxf(p3[0], p4[0]));

					if(shift==fshift_count-1)p5[useful_bins-1]=fmaxf(p2[useful_bins-1], p3[useful_bins-1]);
						else
						p5[useful_bins-1]=fmaxf(p2[useful_bins-1], fmaxf(p3[useful_bins-1], p4[useful_bins-1]));
							
					PRAGMA_IVDEP
					for(r=1;r<useful_bins-1;r++) {
						p5[r]=fmaxf(p2[r], fmaxf(p3[r], p4[r]));
						}
				
					}
				p1=tmp2b;
				tmp2b=tmp2a;
				tmp2a=p1;
				}
				
			if(ei[i]->veto_num<0) {
				for(m=0;m<veto_free;m++) {
					/* Accumulate tmp2 from tmp incorporating all per-ifo pieces */
					/* We need to do all sub-bin shifts in one go */
					for(shift=0;shift<fshift_count;shift++) {
						/* It could be that the chunk is too small to contain any SFTs. 
						 * Skip and continue */
						if(tmp_max_weight[(k*veto_free+m)*fshift_count+shift]<=0)continue;
						
						p1=&(tmp[((k*veto_free+m)*fshift_count+shift)*tmp_stride]);
						p5=&(tmp2a[tmp_stride*shift]);
															
						/* This is at best an approximation when line filtering is on */
						weight=0.5*(tmp_min_weight[(k*veto_free+m)*fshift_count+shift]+tmp_max_weight[(k*veto_free+m)*fshift_count+shift]);
						total_weight[shift]+=weight;
						
						PRAGMA_IVDEP
						for(r=0;r<useful_bins;r++) {
							p5[r]+=p1[r]*weight;
							}
					
					
						min_weight+=tmp_min_weight[(k*veto_free+m)*fshift_count+shift];
						max_weight+=tmp_max_weight[(k*veto_free+m)*fshift_count+shift];
// 						if(min_weight>tmp_min_weight[(k*veto_free+m)*fshift_count+shift])min_weight=tmp_min_weight[(k*veto_free+m)*fshift_count+shift];
// 						if(max_weight<tmp_max_weight[(k*veto_free+m)*fshift_count+shift])max_weight=tmp_max_weight[(k*veto_free+m)*fshift_count+shift];
						}
					
					}
				} else {
				m=ei[i]->veto_num;
				
				/* Accumulate tmp2 from tmp incorporating all per-ifo pieces */
				/* We need to do all sub-bin shifts in one go */
				for(shift=0;shift<fshift_count;shift++) {
					/* It could be that the chunk is too small to contain any SFTs. 
						* Skip and continue */
					if(tmp_max_weight[(k*veto_free+m)*fshift_count+shift]<=0)continue;
					
					p1=&(tmp[((k*veto_free+m)*fshift_count+shift)*tmp_stride]);
					p5=&(tmp2a[tmp_stride*shift]);
														
					/* This is at best an approximation when line filtering is on */
					weight=0.5*(tmp_min_weight[(k*veto_free+m)*fshift_count+shift]+tmp_max_weight[(k*veto_free+m)*fshift_count+shift]);
					total_weight[shift]+=weight;
					
					PRAGMA_IVDEP
					for(r=0;r<useful_bins;r++) {
						p5[r]+=p1[r]*weight;
						}
				

					min_weight+=tmp_min_weight[(k*veto_free+m)*fshift_count+shift];
					max_weight+=tmp_max_weight[(k*veto_free+m)*fshift_count+shift];
// 					if(min_weight>tmp_min_weight[(k*veto_free+m)*fshift_count+shift])min_weight=tmp_min_weight[(k*veto_free+m)*fshift_count+shift];
// 					if(max_weight<tmp_max_weight[(k*veto_free+m)*fshift_count+shift])max_weight=tmp_max_weight[(k*veto_free+m)*fshift_count+shift];
					}
				}
			}
			
		for(shift=0;shift<fshift_count;shift++) {
			p5=&(tmp2a[tmp_stride*shift]);
			
			
			inv_weight=1.0f/total_weight[shift];
			PRAGMA_IVDEP
			for(r=0;r<useful_bins;r++) {
				p5[r]*=inv_weight;
				}
			MODE(compute_universal_statistics)(p5, min_weight, max_weight, &(alignment_grid[lm]), &pst);
			
			if(!isfinite(pst.snr)) {
				/* This is not fatal, but possibly needs checking */
				fprintf(stderr, "*** ERROR: non-finite max_dx pi=%d lm=%d i=%d j=%d shift=%d count=%d tmp_stride=%ld min_weight=%g max_weight=%g total_weight=%g veto_num=%d %d_%d\n", pi, lm, i, j, shift, count, tmp_stride, min_weight, max_weight, total_weight[shift], ei[i]->veto_num, ei[i]->first_chunk, ei[i]->last_chunk); 
// 				fprintf(stderr, "p5={");
// 				for(r=0;r<useful_bins;r++)fprintf(stderr, " %g", p5[r]);
// 				fprintf(stderr, "}\n");
				}
			
			update_power_sum_stats(&pst, &(alignment_grid[lm]), &(stats[i*count+j+shift]));
			}
		}
	}
			
}
/* find largest strain and largest SNR candidates for this patch and log info */
for(i=0;i<nei;i++) {
	pstats_log_extremes(ctx, &(stats[i*count]), ps, count, ei[i], pi);
	}
}


void outer_loop_cruncher(int thread_id, void *data)
{
int pi=(long)data;
SUMMING_CONTEXT *ctx=summing_contexts[thread_id+1];
int ps_tmp_len;
int i,k,m,count;
POWER_SUM **ps=cruncher_contexts[thread_id+1].ps;
POWER_SUM **ps_tmp=cruncher_contexts[thread_id+1].ps_tmp;
POWER_SUM_STATS *le_pstats;
TEMPLATE_INFO ti;

ctx->nchunks=nchunks;
ctx->power_sums_idx=0;

//fprintf(stderr, "%d ", pi);

if(input_templates>=0) {
	thread_mutex_lock(&input_template_mutex);
	
	fseek(INPUT_TEMPLATE_LOG, pi*sizeof(TEMPLATE_INFO), SEEK_SET);
	fread(&ti, sizeof(TEMPLATE_INFO), 1, INPUT_TEMPLATE_LOG);
	
	thread_mutex_unlock(&input_template_mutex);
	
	generate_followup_templates(ctx, &ti, &(ps[0]), &count);
	} else {
	generate_patch_templates(ctx, pi, &(ps[0]), &count);
	}

if(count<1) {
	//free(ps[0]);
	ps[0]=NULL;
	return;
	}

for(i=1;i<nchunks;i++) {
	ctx->power_sums_idx=i;
	clone_templates(ctx, ps[0], count, &(ps[i]));
	}
for(i=0;i<args_info.nchunks_arg;i++) {
	for(k=0;k<veto_free;k++) {
		ctx->accumulate_power_sums(ctx, ps[i*veto_free+k], count, gps_start+i*(gps_stop-gps_start)/args_info.nchunks_arg, gps_start+(i+1)*(gps_stop-gps_start)/args_info.nchunks_arg, veto_info[k].veto_mask);
		}
	}

if(args_info.viterbi_power_sums_arg) {
	log_extremes_viterbi(ctx, pi, ps, count);
	} else {
	/* find largest strain and largest SNR candidates for this patch */
	
	
	if(ctx->log_extremes_pstats_scratch_size<count*(sizeof(*le_pstats)+sizeof(*ctx->diverted))+2*ALIGNMENT) {
		free(ctx->log_extremes_pstats_scratch);

		ctx->log_extremes_pstats_scratch_size=count*(sizeof(*le_pstats)+sizeof(*ctx->diverted))+2*ALIGNMENT;

		ctx->log_extremes_pstats_scratch=do_alloc(count, sizeof(*le_pstats));
		}
	ctx->diverted=(char *)ctx->log_extremes_pstats_scratch;
	le_pstats=(POWER_SUM_STATS *)ALIGN_POINTER(ctx->diverted+count);
	
	memset(ctx->diverted, 0, count*sizeof(*ctx->diverted));
	
	for(i=0;i<nei;i++) {
		ps_tmp_len=0;
		for(k=ei[i]->first_chunk;k<=ei[i]->last_chunk;k++) {
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

		log_extremes(ctx, le_pstats, ei[i], pi, ps_tmp, ps_tmp_len, count);
		}
	}

for(i=0;i<nchunks;i++) {
	free_templates_ctx(ctx, ps[i], count);
	ps[i]=NULL;
	}
//fprintf(stderr, "pps_hits=%ld pps_misses=%ld pps_rollbacks=%ld\n", ctx->pps_hits, ctx->pps_misses, ctx->pps_rollbacks);
}

void outer_loop(void)
{
int pi, i, k;
time_t start_time, end_time;
RGBPic *p;
PLOT *plot;

thread_mutex_init(&data_logging_mutex);
thread_mutex_init(&input_template_mutex);
thread_mutex_init(&divert_buffer_mutex);

assign_per_dataset_cutoff_veto();
assign_cutoff_veto();
assign_detector_veto();

nchunks=args_info.nchunks_arg*veto_free;

create_segments(&ei, &nei);

n_contexts=get_max_threads();
summing_contexts=do_alloc(n_contexts, sizeof(*summing_contexts));
for(i=0;i<n_contexts;i++)
	summing_contexts[i]=create_summing_context();

fprintf(stderr, "veto_free=%d\n", veto_free);
cruncher_contexts=do_alloc(n_contexts, sizeof(*cruncher_contexts));
for(i=0;i<n_contexts;i++) {
	cruncher_contexts[i].ps=do_alloc(nchunks, sizeof(*cruncher_contexts[i].ps));
	cruncher_contexts[i].ps_tmp=do_alloc(nchunks, sizeof(*cruncher_contexts[i].ps_tmp));
	}

fprintf(LOG, "nei: %d\n", nei);

divert_buffer_size=args_info.divert_buffer_size_arg;
divert_buffer_free=0;
divert_buffer_ul_threshold=0.0;
fprintf(stderr, "Allocating %.1fMB for divert buffer\n", divert_buffer_size*sizeof(DIVERTED_ENTRY));
divert_buffer=do_alloc(divert_buffer_size, sizeof(*divert_buffer));
divert_sort_buffer=do_alloc(divert_buffer_size, sizeof(*divert_sort_buffer));

gps_start=min_gps();
gps_stop=max_gps()+1;

reset_jobs_done_ratio();

fprintf(stderr, "Outer loop iteration start memory: %g MB\n", (MEMUSAGE*10.0/(1024.0*1024.0))/10.0);
fprintf(LOG, "Outer loop iteration start memory: %g MB\n", (MEMUSAGE*10.0/(1024.0*1024.0))/10.0);

time(&start_time);

if(input_templates>=0) {
	fprintf(stderr, "%d patches to process\n", input_templates);
	for(pi=0;pi<input_templates;pi++) {
		submit_job(outer_loop_cruncher, (void *)((long)pi));
		}
	} else {
	fprintf(stderr, "%d patches to process\n", patch_grid->npoints);
	for(pi=0;pi<patch_grid->npoints;pi++) {
		submit_job(outer_loop_cruncher, (void *)((long)pi));
		}
	}
k=0;
while(do_single_job(-1)) {
	#if 0
	time(&end_time);
	if(end_time<start_time)end_time=start_time;
	fprintf(stderr, "%d (%f patches/sec)\n", k, k/(1.0*(end_time-start_time+1.0)));
	summing_contexts[0]->print_cache_stats(summing_contexts[0]);
	fprintf(stderr, "%d\n", pi);
	#endif

	k++;
	if(k > args_info.progress_update_interval_arg) {
		fprintf(stderr, "% 3.1f ", jobs_done_ratio()*100);
		k=0;
		}
	}
wait_for_all_done();
fprintf(stderr, "\n");

time(&end_time);
if(end_time<start_time)end_time=start_time;
fprintf(stderr, "Patch speed: %f\n", patch_grid->npoints/(1.0*(end_time-start_time+1.0)));
fprintf(LOG, "Patch speed: %f\n", patch_grid->npoints/(1.0*(end_time-start_time+1.0)));

for(i=0;i<n_contexts;i++) {
	summing_contexts[i]->print_cache_stats(summing_contexts[i]);
	free_summing_context(summing_contexts[i]);
	summing_contexts[i]=NULL;

	free(cruncher_contexts[i].ps);
	free(cruncher_contexts[i].ps_tmp);
	}
free(summing_contexts);
summing_contexts=NULL;
free(cruncher_contexts);
cruncher_contexts=NULL;


fflush(DATA_LOG);
fflush(LOG);

if(patch_grid->max_n_dec<800){
	p=make_RGBPic(patch_grid->max_n_ra*(800/patch_grid->max_n_dec)+140, patch_grid->max_n_dec*(800/patch_grid->max_n_dec));
	} else {
	p=make_RGBPic(patch_grid->max_n_ra+140, patch_grid->max_n_dec);
	}


fprintf(stderr, "%d points written into the data log\n", data_log_index);
fprintf(LOG, "%d points written into the data log\n", data_log_index);

fprintf(stderr, "total_diverted_count: %d\n", total_diverted_count);
fprintf(LOG, "total_diverted_count: %d\n", total_diverted_count);

fprintf(LOG, "divert_buffer_ul_threshold: %g\n", divert_buffer_ul_threshold);
fprintf(LOG, "divert_buffer_free: %g\n", divert_buffer_free);

for(i=0;i<divert_buffer_free;i++) {
	fwrite(&(divert_buffer[i].ti), sizeof(divert_buffer[i].ti), 1, DIVERT_LOG);
	}

if(total_diverted_count>=args_info.divert_count_limit_arg && n_contexts>1) {
	/* This is because we don't know which thread hits the limit first */
	fprintf(stderr, "*** WARNING: divert count limit exceeded, the results are non-repeatable when multiple threads are used\n");
	fprintf(LOG, "*** WARNING: divert count limit exceeded, the results are non-repeatable when multiple threads are used\n");
	}

fprintf(stderr, "Writing skymaps\n");

plot=make_plot(p->width, p->height);

for(i=0;i<nei;i++)
	output_extreme_info(p, ei[i]);

fprintf(stderr, "Outer loop extreme info done memory: %g MB\n", (MEMUSAGE*10.0/(1024.0*1024.0))/10.0);
fprintf(LOG, "Outer loop extreme info done memory: %g MB\n", (MEMUSAGE*10.0/(1024.0*1024.0))/10.0);

for(i=0;i<nei;i++)
	free_extreme_info(ei[i]);

free_plot(plot);
free_RGBPic(p);
}
