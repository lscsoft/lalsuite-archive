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

void compute_normal_sorted_ks_test(float *data, long count, NORMAL_STATS *stats)
{
long i, ks_count;
STAT_TYPE a,b,quantile2std;

if(count==0)return;

if(count & 1){
	stats->mean=data[(count>>1)];
	} else {
	stats->mean=(data[(count>>1)]+data[(count>>1)-1])/2.0;
	}

/* add a method to autocompute this value later */
quantile2std=1.22;
stats->sigma=(data[(count*4)/5]-data[count/5])/quantile2std; 

a=-2.0;
ks_count=0;
for(i=0;i<count;i++){
	b=(((i+1)*1.0)/count)-normal_distribution((data[i]-stats->mean)/stats->sigma);
	if(a<b)a=b;
	if(b>0.07)ks_count++;
	b=normal_distribution((data[i]-stats->mean)/stats->sigma)-((i*1.0)/count);
	if(a<b)a=b;
	if(b>0.07)ks_count++;
	}
stats->ks_count=ks_count;
stats->ks_test=a;
}

static int float_cmp(float *a, float *b)
{
if(*a<*b)return -1;
if(*a>*b)return 1;
return 0;
}

void sort_floats(float *data, long count)
{
qsort(data,count,sizeof(float),float_cmp);
}


void compute_normal_ks_test(float *data, long count, NORMAL_STATS *stats)
{
float *tmp;

if(count>=1000){
	tmp=do_alloc(count, sizeof(*tmp));
	} else {
	tmp=alloca(count*sizeof(*tmp));
	}
memcpy(tmp, data, count*sizeof(*tmp));

sort_floats(tmp, count);

compute_normal_sorted_ks_test(tmp, count, stats);

if(count>=1000){
	free(tmp);
	}
}
