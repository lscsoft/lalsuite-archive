#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <gsl/gsl_sf.h>


#include "global.h"
#include "statistics.h"

#define LEFT_NORMAL_BOUND   -4.0
#define RIGHT_NORMAL_BOUND   4.0
#define NORMAL_SAMPLE_COUNT  200
#define NORMAL_RANGE   (RIGHT_NORMAL_BOUND-LEFT_NORMAL_BOUND)
#define NORMAL_STEP    (NORMAL_RANGE/NORMAL_SAMPLE_COUNT)

STAT_TYPE normal_distribution_table[NORMAL_SAMPLE_COUNT];

static inline STAT_TYPE exact_normal_distribution(STAT_TYPE x)
{
return (1.0-gsl_sf_erf_Q(x));
}


STAT_TYPE sampled_normal_distribution(STAT_TYPE x)
{
int i;
if(x<LEFT_NORMAL_BOUND)return 0.0;
if(x>RIGHT_NORMAL_BOUND)return 1.0;
i=floor((x-LEFT_NORMAL_BOUND)/NORMAL_STEP-0.5);
if(i>=NORMAL_SAMPLE_COUNT)i=NORMAL_SAMPLE_COUNT-1;
return normal_distribution_table[i];
}

STAT_TYPE interpolated_normal_distribution(STAT_TYPE x)
{
int i;
STAT_TYPE dx;
if(x<LEFT_NORMAL_BOUND)return 0.0;
if(x>RIGHT_NORMAL_BOUND)return 1.0;
if(x<(LEFT_NORMAL_BOUND+NORMAL_STEP*0.5)){
	return normal_distribution_table[0]*(x-LEFT_NORMAL_BOUND)/(NORMAL_STEP*0.5);
	}
if(x>(RIGHT_NORMAL_BOUND-NORMAL_STEP*0.5)){
	return normal_distribution_table[NORMAL_SAMPLE_COUNT-1]*
		(RIGHT_NORMAL_BOUND-x)/(NORMAL_STEP*0.5)+
		(x+NORMAL_STEP*0.5-RIGHT_NORMAL_BOUND)/(NORMAL_STEP*0.5);
	}
i=floor((x-LEFT_NORMAL_BOUND-NORMAL_STEP*0.5)/NORMAL_STEP);
if(i>=NORMAL_SAMPLE_COUNT-1)i=NORMAL_SAMPLE_COUNT-2;
dx=(x-LEFT_NORMAL_BOUND)/NORMAL_STEP-(i+0.5);
return normal_distribution_table[i]*(1.0-dx)+normal_distribution_table[i+1]*dx;
}

#define normal_distribution interpolated_normal_distribution

extern FILE *LOG;

void init_statistics(void)
{
int i;
STAT_TYPE x,y,err;
for(i=0;i<NORMAL_SAMPLE_COUNT;i++){
	normal_distribution_table[i]=
		exact_normal_distribution((NORMAL_STEP/2.0)
				        +LEFT_NORMAL_BOUND
					+i*NORMAL_STEP);
	}
err=0.0;
/* do it the brute force way - much harder to screw up */
for(x=-10.0;x<10.0;x+=0.001){
	y=fabs(normal_distribution(x)-exact_normal_distribution(x));
	if(y>err)err=y;
	}
fprintf(stderr, "Normal distribution approximation error is %f\n", err);
fprintf(LOG, "Normal distribution approximation error: %f\n", err);
}

void compute_normal_sorted_stats(float *data, long count, NORMAL_STATS *stats)
{
long i, ks_count_plus, ks_count_minus, dir;
STAT_TYPE a,b,quantile2std, ks_level;
float mean, sigma, step;
double d;

if(count==0)return;

if(stats->flag & STAT_FLAG_ESTIMATE_MEAN){
	if(count & 1){
		stats->mean=data[(count>>1)];
		} else {
		stats->mean=(data[(count>>1)]+data[(count>>1)-1])/2.0;
		}
	mean=0.0;
	for(i=50;i<count-50;i++)mean+=data[i];
	mean/=(count-100);
	stats->mean=mean;
	}

if(stats->flag & STAT_FLAG_ESTIMATE_SIGMA){
	/* add a method to autocompute this value later */
	quantile2std=1.22;
	/* 1.62 appears to work better.. */
	quantile2std=1.62;
	/* 1.69 even better.. */
	quantile2std=1.69;

	stats->sigma=(data[(count*4)/5]-data[count/5])/quantile2std; 
	
	sigma=0.0;
	mean=stats->mean;
	for(i=50;i<count-50;i++){
		sigma+=(data[i]-mean)*(i-50.0)/400.0;
		}
	stats->sigma=sigma*5.0/400.0;
	}
	
if(stats->flag & STAT_FLAG_ESTIMATE_KS_LEVEL){
	stats->ks_level=0.022;
	}

if(0 && stats->flag & STAT_FLAG_AUTOFIT){
	mean=stats->mean;
	sigma=stats->sigma;
	dir=0;
	step=0.1;
	while(1) {
		a=-2.0;
		ks_count_plus=0;
		ks_count_minus=0;
		ks_level=stats->ks_level;
		for(i=0;i<count;i++){
			b=(((i+1)*1.0)/count)-normal_distribution((data[i]-mean)/sigma);
			//if(a<b)a=b;
			if(b>ks_level)ks_count_plus++;
			b=normal_distribution((data[i]-mean)/sigma)-((i*1.0)/count);
			//if(a<b)a=b;
			if(b>ks_level)ks_count_minus++;
			}
		if(ks_count_plus>(ks_count_minus+1)){
			if(dir!=1){
				step=step*0.5;
				dir=1;
                                if(step<1e-6)break;
				}
			mean-=sigma*step;
			continue;
			}
			else
		if((ks_count_plus+1)<ks_count_minus){
			if(dir!=-1){
				step=step*0.5;
				dir=-1;
				if(step<1e-6)break;
				}
			mean+=sigma*step;
			continue;
			}
		break;
		} 
//	fprintf(stderr, "{%d %d} ", ks_count_plus, ks_count_minus);
	stats->mean=mean;
	}

if(stats->flag & STAT_FLAG_COMPUTE_KS_TEST){
	a=-2.0;
	ks_count_plus=0;
	ks_count_minus=0;
	ks_level=stats->ks_level;
	mean=stats->mean;
	sigma=stats->sigma;
	for(i=0;i<count;i++){
		b=(((i+1)*1.0)/count)-normal_distribution((data[i]-mean)/sigma);
		if(a<b)a=b;
		if(b>ks_level)ks_count_plus++;
		b=normal_distribution((data[i]-mean)/sigma)-((i*1.0)/count);
		if(a<b)a=b;
		if(b>ks_level)ks_count_minus++;
		}
	stats->ks_count=ks_count_plus+ks_count_minus;
	stats->ks_test=a;
	}
}

static int float_cmp(float *a, float *b)
{
if(*a<*b)return -1;
if(*a>*b)return 1;
return 0;
}

void sort_floats(float *data, long count)
{
qsort(data, count, sizeof(float), float_cmp);
}


void compute_normal_stats(float *data, long count, NORMAL_STATS *stats)
{
float *tmp;

if(stats->flag & STAT_FLAG_SORT_DATA){

	if(count>=1000){
		tmp=do_alloc(count, sizeof(*tmp));
		} else {
		tmp=alloca(count*sizeof(*tmp));
		}
	memcpy(tmp, data, count*sizeof(*tmp));

	sort_floats(tmp, count);

	compute_normal_sorted_stats(tmp, count, stats);

	if(count>=1000){
		free(tmp);
		}
		
	} else 
if(stats->flag & STAT_FLAG_INPLACE_SORT_DATA){

	sort_floats(data, count);
	compute_normal_sorted_stats(data, count, stats);
	
	} else {
	
	compute_normal_sorted_stats(data, count, stats);
	
	}
}

HISTOGRAM * new_histogram(int nbins, int nbands)
{
HISTOGRAM *h;

h=do_alloc(1, sizeof(*h));
h->nbands=nbands;
h->nbins=nbins;
h->max=do_alloc(nbands, sizeof(*h->max));
h->min=do_alloc(nbands, sizeof(*h->min));
h->hist=do_alloc(nbands*nbins, sizeof(*h->hist));

return h;
}

void free_histogram(HISTOGRAM *h)
{
free(h->hist);
free(h->max);
free(h->min);
free(h);
}

void compute_histogram_f(HISTOGRAM *h, float *data, int *bands, long count)
{
int i,j,k;
float f;

for(i=0;i<h->nbands;i++){
	h->max[i]=-2e6;
	h->min[i]=-1e6;
	}
for(i=0;i<h->nbins*h->nbands;i++)h->hist[i]=0;
/* 1st pass */
/* computer histogram limits */
k=0; /* in case bands==NULL */
for(i=0;i<count;i++){
	if(bands!=NULL){
		k=bands[i];
		if((k<0) || (k>=h->nbands))continue;
		}
	f=data[i];
	if(h->min[k]>h->max[k]){
		h->min[k]=f;
		h->max[k]=f;
		continue;
		}
	if(f<h->min[k])h->min[k]=f;
	if(f>h->max[k])h->max[k]=f;
	}
/* make sure that max>min */
for(i=0;i<h->nbands;i++){
	if(h->min[i]==h->max[i])h->max[i]*=2.0;
	}
/* 2nd pass */
/* compute histogram */
for(i=0;i<count;i++){
	if(bands!=NULL){
		k=bands[i];
		if((k<0) || (k>=h->nbands))continue;
		}
	j=floor((h->nbins*(data[i]-h->min[k]))/(h->max[k]-h->min[k]));
	if(j>=h->nbins)j=h->nbins-1;
	h->hist[k*h->nbins+j]++;
	}
}

void print_histogram(FILE *f, HISTOGRAM *h, char *prefix)
{
int i,k;
fprintf(f,"histogram: band min max counts..\n");
for(k=0;k<h->nbands;k++){
	fprintf(f, "%s: %d %g %g", prefix, k, h->min[k], h->max[k]);
	for(i=0;i<h->nbins;i++)
		fprintf(f, " %ld", h->hist[k*h->nbins+i]);
	fprintf(f,"\n");
	}
}
