#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/* We need this define to get NAN values */
#define __USE_ISOC99
#include <math.h>

#include "global.h"
#include "power_cache.h"
#include "power_sums.h"
#include "power_sum_stats.h"
#include "dataset.h"
#include "statistics.h"
#include "grid.h"
#include "cmdline.h"

extern struct gengetopt_args_info args_info;
extern int nbins, first_bin, side_cut, useful_bins;

double upper_limit_comp=NAN, strain_comp=NAN; 

ALIGNMENT_COEFFS *alignment_grid=NULL;
int alignment_grid_free=0;
int alignment_grid_size=0;

extern FILE *LOG;

/* Include file with Feldman-Cousins upper limits directly
   so as to benefit from function inlining */

#include "fc.c"

void compute_alignment_coeffs1(ALIGNMENT_COEFFS *ac)
{
double a, a_plus_sq, a_cross_sq, cpsi, spsi;
a=cos(ac->iota);
a=a*a;

a_plus_sq=(1+a)*0.5;
a_plus_sq=a_plus_sq*a_plus_sq;
a_cross_sq=a;

cpsi=cos(2*ac->psi);
spsi=sin(2*ac->psi);

ac->pp=(cpsi*cpsi*a_plus_sq+spsi*spsi*a_cross_sq);
ac->pc=2*cpsi*spsi*(a_plus_sq-a_cross_sq);
ac->cc=(spsi*spsi*a_plus_sq+cpsi*cpsi*a_cross_sq);

ac->pppp=ac->pp*ac->pp;
ac->pppc=2*ac->pp*ac->pc;
ac->ppcc=2*ac->pp*ac->cc+ac->pc*ac->pc;
ac->pccc=2*ac->pc*ac->cc;
ac->cccc=ac->cc*ac->cc;
}

/* this is identical to the one above, but follows the formula in polarization.pdf */

void compute_alignment_coeffs(ALIGNMENT_COEFFS *ac)
{
double a, a_plus_sq, a_cross_sq, cpsi, spsi, asum, adiff;
a=cos(ac->iota);
a=a*a;

a_plus_sq=(1+a)*0.5;
a_plus_sq=a_plus_sq*a_plus_sq;
a_cross_sq=a;

cpsi=cos(4*ac->psi);
spsi=sin(4*ac->psi);

asum=0.25*(a_plus_sq+a_cross_sq);
adiff=0.25*(a_plus_sq-a_cross_sq);

ac->pp=(asum+adiff*cpsi);
ac->pc=2*adiff*spsi;
ac->cc=(asum-adiff*cpsi);

ac->pppp=ac->pp*ac->pp;
ac->pppc=2*ac->pp*ac->pc;
ac->ppcc=2*ac->pp*ac->cc+ac->pc*ac->pc;
ac->pccc=2*ac->pc*ac->cc;
ac->cccc=ac->cc*ac->cc;
}

void generate_alignment_grid(void)
{
int i, j, k;
int npsi=args_info.npsi_arg, niota=args_info.niota_arg;

alignment_grid_size=npsi*niota+1;
alignment_grid_free=npsi*niota+1;
alignment_grid=do_alloc(alignment_grid_size, sizeof(*alignment_grid));

/* First polarization is circular one */
alignment_grid[0].psi=0.0;
alignment_grid[0].iota=0.0;

k=1;
for(j=0;j<niota;j++)
	for(i=0;i<npsi;i++) {
		/* note: make better grid by not overcovering circular polarization neighbourhood */
		alignment_grid[k].psi=(0.5*M_PI*(i+0.5*(j&1)))/npsi;
		alignment_grid[k].iota=acos((1.0*j)/niota); 
		k++;
		}

for(k=0;k<alignment_grid_free;k++) {
	fprintf(LOG, "alignment entry %d: %f %f\n", k, alignment_grid[k].iota, alignment_grid[k].psi);
	compute_alignment_coeffs(&(alignment_grid[k]));
	}
fprintf(stderr, "Alignment grid size: %d\n", alignment_grid_free);
fprintf(LOG, "Alignment grid size: %d\n", alignment_grid_free);
}

/* an implementation of merge sort - this modifies input array */
void merge_sort_floats2(float *data, int count, int step)
{
int i, k, touched;
float a,b;
/* handcode cases of small number of items */
k=step;
if(k>=count)return;
k+=step;
if(k>=count) {
	a=data[0];
	b=data[step];
	if(a>b) {
		data[0]=b;
		data[step]=a;
		}
	return;
	}
merge_sort_floats2(data, count, step<<1);
merge_sort_floats2(data+step, count-step, step<<1);
touched=1;
/*while(touched) {
	touched=0;
	a=data[0];
	for(i=step;i<count;i+=step) {
		b=data[i];
		if(a>b) {
			data[i-step]=b;
			data[i]=a;
			touched=1;
			} else {
			a=b;
			}
		}
	}*/
while(touched) {
	touched=0;
	for(i=0;i<count-step;i+=step<<1) {
		a=data[i];
		b=data[i+step];
		if(a>b) {
			data[i]=b;
			data[i+step]=a;
			touched=1;
			}
		}
	for(i=step;i<count-step;i+=step<<1) {
		a=data[i];
		b=data[i+step];
		if(a>b) {
			data[i]=b;
			data[i+step]=a;
			touched=1;
			}
		}
	}
}

/* an implementation of merge sort - this modifies input array */
void merge_sort_floats(float *data, int count)
{
int i, j, k, m;
float a,b;
int step;
float *tmp;
/* handcode cases of small number of items */
if(count<=1)return;

a=data[0];
b=data[1];
if(a>b) {
	data[0]=b;
	data[1]=a;
	}
step=2;

tmp=alloca(count*sizeof(*tmp));
while(step<count) {
	k=step*2;
	if(k>count) {
		k=count;
		}
	merge_sort_floats(&(data[step]), k-step);

	a=data[0];
	b=data[step];
	for(i=0,j=step, m=0;(i<step) && (j<k);m++) {
		if(a<b) {
			tmp[m]=a;
			i++;
			a=data[i];
			} else {
			tmp[m]=b;
			j++;
			b=data[j];
			}
		}
	for(;i<step;i++,m++) tmp[m]=data[i];
	for(;j<k;j++,m++) tmp[m]=data[j];
	memcpy(data, tmp, k*sizeof(*tmp));
	step=k;
	}

}

/* an implementation of quick sort - this modifies input array */
void quick_sort_floats1(float *data, int count)
{
int i, j;
float a,b;
/* handcode cases of small number of items */
if(count<=1)return;
a=data[0];
if(count==2) {
	b=data[1];
	if(a>b) {
		data[0]=b;
		data[1]=a;
		}
	return;
	}

/* we do an average so that qsort performs well in the case of sloped array */
/*	for(i=0;i<count;i++)fprintf(stderr, "%g ", data[i]);
	fprintf(stderr, "\n");*/
i=0;
j=count-1;
//a=0.5*(a+data[j]);
while(i<j) {
	if(data[i]<=a) {
		i++;
		continue;
		}
	if(data[j]>=a) {
		j--;
		continue;
		}
	b=data[i];
	data[i]=data[j];
	data[j]=b;
	i++;
	j--;
	}
if(i==j) {
	if(data[j]>a)j--;
	}
if(data[j]<a)j++;
if(j==count) {
	b=data[0];
	data[0]=data[count-1];
	data[count-1]=b;
	quick_sort_floats1(data, count-1);
	return;
	}
if(j==0) {
	quick_sort_floats1(data+1, count-1);
	return;
	}
/*if(!j || j==count) {
	fprintf(stderr, "*** INTERNAL ERROR: recursion on quick sort count=%d j=%d a=%f\n", count, j, a);
	for(i=0;i<count;i++)fprintf(stderr, "%g ", data[i]);
	fprintf(stderr, "\n");
	exit(-1);
	}*/
quick_sort_floats1(data, j);
quick_sort_floats1(&(data[j]), count-j);
}

/* an optimized implementation of quick sort */
static inline void manual_sort_floats(float *data, int count)
{
float a,b,c;
if(count<=1)return;
a=data[0];
if(count==2) {
	b=data[1];
	if(a>b) {
		data[0]=b;
		data[1]=a;
		}
	return;
	}
if(count==3) {
	b=data[1];
	if(a>b) {
		b=a;
		a=data[1];
		}
	c=data[2];
	if(c<b) {
		if(c<a) {
			c=b;
			b=a;
			a=data[2];
			} else {
			c=b;
			b=data[2];
			}
		}
	data[0]=a;
	data[1]=b;
	data[2]=c;
	return;
	}
if(count==4) {
	}
return;
}

static inline int partition_floats(float *data, int count)
{
float a,b,c;
float *first=data;
float *last=&(data[count-1]);

a=*first;
first++;
//a=0.5*(a+data[j]);
while(first<last) {
	b=*first;
	if(b<=a) {
		first++;
		continue;
		}
	c=*last;
	if(c>=a) {
		last--;
		continue;
		}
	*first=c;
	*last=b;
	first++;
	last--;
	}
if(first==last) {
	if(*last>a)last--;
	}
if(*last<a)last++;
return (last-data);
}

/* an implementation of quick sort - this modifies input array */
void quick_sort_floats(float *data, int count)
{
int j;
float b;
int stack_len=0;

while(1) {
	if(count<4) {
		manual_sort_floats(data, count);
		return;
		} else {
		j=partition_floats(data, count);
		}
	if(j<0) return;
	if(j==count) {
		b=data[0];
		data[0]=data[count-1];
		data[count-1]=b;
		count--;
		continue;
		}
	if(j==0) {
		data++;
		count--;
		continue;
		}
/*if(!j || j==count) {
	fprintf(stderr, "*** INTERNAL ERROR: recursion on quick sort count=%d j=%d a=%f\n", count, j, a);
	for(i=0;i<count;i++)fprintf(stderr, "%g ", data[i]);
	fprintf(stderr, "\n");
	exit(-1);
	}*/
	if(j<4)manual_sort_floats(data, j);
		else quick_sort_floats(data, j);
	stack_len++;
	data+=j;
	count-=j;
	continue;
	}
}

void bucket_sort_floats(float *data, int count)
{
float *tmp[4];
int tmp_count[4];
float a,b,c, mult;
int i, k, m;
if(count<=16) {
	merge_sort_floats2(data, count, 1);
	return;
	}

// if(count<=2) {
// 	if(count<=1)return;
// 	a=data[0];
// 	b=data[1];
// 	if(a>b) {
// 		data[0]=b;
// 		data[1]=a;
// 		}
// 	return;
// 	}

for(k=0;k<4;k++) {
	tmp[k]=alloca(count*sizeof(**tmp));
	tmp_count[k]=0;
	}

a=data[0];
b=a;
for(i=1;i<count;i++) {
	c=data[i];
	if(c>b)b=c;
		else
	if(c<a)a=c;
	}
if(b<=a)b=a+1;
mult=4.0/(b-a);

for(i=0;i<count;i++) {
	c=data[i];
	k=((c-a)*mult);
	if(k>3)k=3;
	if(k<0)k=0;
	tmp[k][tmp_count[k]]=c;
	tmp_count[k]++;
	}
i=0;
for(k=0;k<4;k++) {
	if(tmp_count[k]==count) {
		merge_sort_floats2(data, count, 1);
		return;
		}
	bucket_sort_floats(tmp[k], tmp_count[k]);
	for(m=0;m<tmp_count[k];m++) {
		data[i]=tmp[k][m];
		i++;
		}
	}
}

int is_sorted(float *data, int count)
{
int i;
float a;
a=data[0];
for(i=1;i<count;i++) {
	if(a>data[i])return 0;
	}
return 1;
}

void point_power_sum_stats(PARTIAL_POWER_SUM_F *pps, ALIGNMENT_COEFFS *ag, POINT_STATS *pst)
{
int i;
float M, S, a, inv_S, inv_weight;
float *tmp=NULL;
NORMAL_STATS nstats;
float max_dx;
int max_dx_bin;
float weight, min_weight, max_weight;

/* allocate on stack, for speed */
tmp=aligned_alloca(useful_bins*sizeof(*tmp));

memset(&nstats, 0, sizeof(nstats));

/* sort to compute robust estimates */
nstats.flag= STAT_FLAG_ESTIMATE_MEAN
	| STAT_FLAG_ESTIMATE_SIGMA;

if(args_info.ks_test_arg){
	nstats.flag|=STAT_FLAG_ESTIMATE_KS_LEVEL
		| STAT_FLAG_COMPUTE_KS_TEST;
	}


if(pps->weight_arrays_non_zero) {
	max_weight=0;
	min_weight=1e50;

	if(!pps->collapsed_weight_arrays) {
		for(i=0;i<useful_bins;i++) {
			pps->weight_pppp[i]+=pps->c_weight_pppp;
			pps->weight_pppc[i]+=pps->c_weight_pppc;
			pps->weight_ppcc[i]+=pps->c_weight_ppcc;
			pps->weight_pccc[i]+=pps->c_weight_pccc;
			pps->weight_cccc[i]+=pps->c_weight_cccc;
			}
		pps->c_weight_pppp=0;
		pps->c_weight_pppc=0;
		pps->c_weight_ppcc=0;
		pps->c_weight_pccc=0;
		pps->c_weight_cccc=0;
		pps->collapsed_weight_arrays=1;
		}

	for(i=0;i<useful_bins;i++) {
		weight=(pps->weight_pppp[i]*ag->pppp+
			pps->weight_pppc[i]*ag->pppc+
			pps->weight_ppcc[i]*ag->ppcc+
			pps->weight_pccc[i]*ag->pccc+
			pps->weight_cccc[i]*ag->cccc);
	
		if(weight>max_weight)max_weight=weight;
		if(weight<min_weight)min_weight=weight;

		tmp[i]=(pps->power_pp[i]*ag->pp+pps->power_pc[i]*ag->pc+pps->power_cc[i]*ag->cc)/weight;
			
		}
	} else {
	weight=(pps->c_weight_pppp*ag->pppp+
		pps->c_weight_pppc*ag->pppc+
		pps->c_weight_ppcc*ag->ppcc+
		pps->c_weight_pccc*ag->pccc+
		pps->c_weight_cccc*ag->cccc);
	max_weight=weight;
	min_weight=weight;

	inv_weight=1.0/weight;

	for(i=0;i<useful_bins;i++) {
		tmp[i]=((float)pps->power_pp[i]*ag->pp+(float)pps->power_pc[i]*ag->pc+(float)pps->power_cc[i]*ag->cc)*inv_weight;
		}
	}

/* 0 weight can happen due to extreme line veto at low frequencies and small spindowns */
if(min_weight<=0) {
	pst->bin=0;
	pst->iota=ag->iota;
	pst->psi=ag->psi;
	
	/* convert to upper limit units */
	pst->S=-1;
	pst->M=-1;
	
	pst->ul=-1;
	pst->ll=-1;
	pst->centroid=-1;
	pst->snr=-1;
	
	pst->max_weight=-1;
	pst->weight_loss_fraction=1;
	pst->ks_value=-1;
	pst->ks_count=-1;
	return;
	}
	
/* find highest bin */
max_dx=tmp[0];
max_dx_bin=0;

for(i=1;i<useful_bins;i++){
	a=tmp[i];
	if(a>max_dx) {
		max_dx=a;
		max_dx_bin=i;
		}
	}

/* sort to compute statistics */
/*merge_sort_floats(tmp, useful_bins, 1);*/
/*bucket_sort_floats(tmp, useful_bins);*/
/* merge_sort_floats(tmp, useful_bins); */
quick_sort_floats(tmp, useful_bins);

if(!is_sorted(tmp, useful_bins)) {
	fprintf(stderr, "Internal error: incorrectly sorted array\n");
	for(i=0;i<useful_bins;i++)fprintf(stderr, " %f", tmp[i]);
	fprintf(stderr, "\n");
	exit(-1);
	}

compute_normal_stats(tmp, useful_bins, &nstats);

M=nstats.mean;
S=nstats.sigma;

inv_S=1.0/S;

/* convert to SNR from the highest power */
max_dx=(max_dx-M)*inv_S;

if(max_dx<=0) {
	/* In theory we could have max_dx=0 because the distribution is flat, but we really should not have this */
	fprintf(stderr, "***ERROR - max_dx<=0  max_dx=%g max_dx_bin=%d M=%g S=%g inv_S=%g\n",
			max_dx,
			max_dx_bin,
			M,
			S,
			inv_S);
	}

pst->bin=max_dx_bin;
pst->iota=ag->iota;
pst->psi=ag->psi;

/* convert to upper limit units */
pst->S=sqrt(S)*strain_comp;
pst->M=sqrt(M)*strain_comp;

pst->ul=sqrt(upper_limit95(max_dx)*S)*strain_comp*upper_limit_comp;
pst->ll=sqrt(lower_limit95(max_dx)*S)*strain_comp;
pst->centroid=sqrt(max_dx*S)*strain_comp*upper_limit_comp;
pst->snr=max_dx;

pst->max_weight=max_weight;
pst->weight_loss_fraction=(max_weight-min_weight)/max_weight;
pst->ks_value=nstats.ks_test;
pst->ks_count=nstats.ks_count;
}

void power_sum_stats(PARTIAL_POWER_SUM_F *pps, POWER_SUM_STATS *stats)
{
int k;
POINT_STATS pst;

stats->highest_ul.ul=-1;
stats->highest_circ_ul.ul=-1;
stats->highest_snr.snr=-1;
stats->highest_ks.ks_value=-1;
stats->highest_M.M=-1;
stats->highest_S.S=-1;
stats->max_weight_loss_fraction=-1;
stats->max_weight=-1;
stats->min_weight=1e50;

memset(&pst, 0, sizeof(pst));

for(k=0;k<alignment_grid_free;k++) {
	point_power_sum_stats(pps, &(alignment_grid[k]), &(pst));

	if(pst.snr>stats->highest_snr.snr) {
		memcpy(&(stats->highest_snr), &(pst), sizeof(pst));
		}

	if(pst.ul>stats->highest_ul.ul) {
		memcpy(&(stats->highest_ul), &(pst), sizeof(pst));
		}

	if(pst.ks_value>stats->highest_ks.ks_value) {
		memcpy(&(stats->highest_ks), &(pst), sizeof(pst));
		}

	if(pst.M>stats->highest_M.M) {
		memcpy(&(stats->highest_M), &(pst), sizeof(pst));
		}

	if(pst.S>stats->highest_S.S) {
		memcpy(&(stats->highest_S), &(pst), sizeof(pst));
		}

	/* Let us consider anything with iota < 1e-5 as circular. 
		In practice this should only be one point */
	if(alignment_grid[k].iota<1e-5) {
		memcpy(&(stats->highest_circ_ul), &(pst), sizeof(pst));
		}

	if(pst.max_weight>stats->max_weight)stats->max_weight=pst.max_weight;
	if(pst.max_weight<stats->min_weight)stats->min_weight=pst.max_weight;
	if(pst.weight_loss_fraction>stats->max_weight_loss_fraction)stats->max_weight_loss_fraction=pst.weight_loss_fraction;
	}
}

void init_power_sum_stats(void)
{
init_fc_ul();
init_fc_ll();
verify_limits();

/* Account for power loss due to Hann windowing */
upper_limit_comp=1.0/0.85; 


// // /*	/* Extra factor to convert to amplitude from RMS power */
// // strain_comp=sqrt(2.0);*/
	/* New AM response correctly computes expected power from h0 */
strain_comp=1.0;
	/* Extra factor to convert to strain from raw SFT units */
strain_comp/=(1800.0*16384.0);
	/* Extra factor to account for the fact that only half of SFT
	   coefficients is stored */
strain_comp*=sqrt(2.0);
	/* Revert strain normalization */
strain_comp*=args_info.strain_norm_factor_arg;

generate_alignment_grid();
}