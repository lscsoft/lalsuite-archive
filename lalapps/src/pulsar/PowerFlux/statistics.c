#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gsl/gsl_sf.h>


#include "global.h"
#include "statistics.h"

static inline STAT_TYPE normal_distribution(STAT_TYPE x)
{
return (1.0-gsl_sf_erf_Q(x));
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
