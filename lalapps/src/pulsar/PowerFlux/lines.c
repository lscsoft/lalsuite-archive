#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "global.h"
#include "lines.h"
#include "cmdline.h"

extern long nbins, first_bin;
extern FILE *LOG;

unsigned char *lines;
int lines_list[5]={-1,-1,-1,-1,-1};

extern struct gengetopt_args_info args_info;

int double_cmp(double *a, double *b)
{
if(*a<*b)return -1;
if(*a>*b)return 1;
return 0;
}

LINES_REPORT *make_lines_report(int x0, int x1, int nlines)
{
LINES_REPORT *lr;
int i;

lr=do_alloc(1, sizeof(*lr));
lr->x0=x0;
lr->x1=x1;
lr->lines=do_alloc(x1-x0+1, sizeof(*(lr->lines)));
/* zero out the bitmap.. in current implementation do_alloc
does this (as it calls calloc), but we might want to change
do_alloc implementation later */
memset(lr->lines, 0, (x1-x0+1)*sizeof(*(lr->lines)));

lr->nlines=nlines;
lr->lines_list=do_alloc(nlines, sizeof(*lr->lines_list));
for(i=0;i<nlines;i++)lr->lines_list[i]=-1; /* no line */

lr->nlittle=30;
return lr;
}

void free_lines_report(LINES_REPORT *lr)
{
free(lr->lines);
free(lr->lines_list);
}

#if 0 
/* old simpler code.. comment it out */
void detect_lines(double *mean)
{
double *tmp;
double median,q5,q30;
int i,j;
lines=do_alloc(nbins, sizeof(*lines));
tmp=do_alloc(nbins, sizeof(*tmp));
memcpy(tmp, mean, nbins*sizeof(*tmp));
qsort(tmp, nbins, sizeof(double), double_cmp);

median=tmp[nbins>>1];
q5=tmp[nbins-5];
q30=tmp[nbins-30];

for(i=0;i<nbins;i++){
	if(mean[i]<q5)lines[i]=0;
		else {
		if(mean[i]>(2*q30-median)){
			lines[i]=1;
			fprintf(stderr,"Line found in %d bin: \t%g Hz\n",
				i, (first_bin+i)/1800.0);
			for(j=0;j<5;j++){
				if(lines_list[j]<0)lines_list[j]=i;
				}
			}
		}
	}

free(tmp);
}
#endif

void detect_lines_d(double *z, LINES_REPORT *lr)
{
double *tmp;
double median, qlines, qmost, unit_strength;
int i, j, k, nb, vh_count=0;

nb=lr->x1-lr->x0+1;

if(nb<=lr->nlines){
	fprintf(stderr,"ERROR: No data to analyze.. x0=%d x1=%d nlines=%d\n",
		lr->x0, lr->x1, lr->nlines);
	return;
	}

if(nb<=lr->nlittle){
	fprintf(stderr,"ERROR: No data to analyze.. x0=%d x1=%d nlittle=%d\n",
		lr->x0, lr->x1, lr->nlittle);
	return;
	}

tmp=do_alloc(nb, sizeof(*tmp));
memcpy(tmp, z+lr->x0, nb*sizeof(*tmp));
qsort(tmp, nb, sizeof(*tmp), double_cmp);

median=tmp[nb>>1];
qlines=tmp[nb-lr->nlines];
qmost=tmp[nb-lr->nlittle];

unit_strength=qmost-median;
if(unit_strength<=0.0)unit_strength=1.0;

/* first pass - mark suspicious bins */
for(i=lr->x0;i<=lr->x1;i++){
	lr->lines[i]=0;
	if(z[i]>qmost)lr->lines[i]|=LINE_HIGH;
	if(z[i]>qlines)lr->lines[i]|=LINE_CANDIDATE;
	if(z[i]>(2*qmost-median)){
		lr->lines[i]|=LINE_VERY_HIGH;
		vh_count++;
		}
	}

/* second pass - mark clustered lines */
for(i=lr->x0;i<=lr->x1;i++){
	if(lr->lines[i] & LINE_HIGH){
		for(j=i;(j<=lr->x1)&&(lr->lines[j]&LINE_HIGH);j++);
		if((j-i)>=5){
			for(k=i;k<j;k++)lr->lines[j]|=LINE_CLUSTERED;
			} else
		if((i>lr->x0) && (i<lr->x1) && ((z[i]-median)>2*(z[i-1]-median)) && ((z[i]-median)>2*(z[i+1]-median)))
			lr->lines[i]|=LINE_ISOLATED;
		}
	}
	
/* third pass - decide what will be considered to be a line */
for(i=lr->x0+1;i<=lr->x1-1;i++){
	 if(lr->lines[i]){
	 	/* Note: strength is in units of unit_strength, 
		   which is something like 2-3 sigma, 
		   depending on what was passed in */
	 	fprintf(stderr, "line detected: bin=%d z=%g strength=%f flag=0x%08x\n", i, z[i], (z[i]-median)/unit_strength,lr->lines[i]);
	 	fprintf(LOG, "line detected: bin=%d z=%g strength=%f flag=0x%08x\n", i, z[i], (z[i]-median)/unit_strength, lr->lines[i]);
		}		
	 
	/* be very conservative mark only points which are LINE_VERY_HIGH
	   and have both neighbours  below 0.4 level of the center line (i.e. side lobes (due to Hann window)
	    of strictly bin-centered line */
	 if(((lr->lines[i] & (LINE_CANDIDATE | LINE_VERY_HIGH | LINE_CLUSTERED))==(LINE_CANDIDATE | LINE_VERY_HIGH))){
		lr->lines[i]|=LINE_YES;
		fprintf(stderr,"line detected: i=%d z[i]=%g\n", i, z[i]);
	
		for(j=0;j<lr->nlines;j++){
			/* this can happen if we are adding new lines */
			if(lr->lines_list[j]==i)break;
			/* empty slot */
			if(lr->lines_list[j]<0){
				lr->lines_list[j]=i;
				break;
				}
			}
		}
	  if(lr->lines[i] & LINE_CANDIDATE){
	  	fprintf(stderr,"i=%d %g %g %g\n", i,
			z[i-1]-median, 
			z[i]-median, 
			z[i+1]-median);
	  	}
	}
lr->median=median;
lr->qmost=qmost;
lr->qlines=qlines;
fprintf(stderr,"median=%g qmost=%g qlines=%g\n",
		 median, qmost, qlines);
fprintf(stderr,"bins marked \"very high\": %d\n", vh_count);
free(tmp);
}

void detect_lines_f(float *z, LINES_REPORT *lr)
{
float *tmp;
float median, qlines, qmost, unit_strength;
int i, j, k, nb;
int vh_count=0;

nb=lr->x1-lr->x0+1;

if(nb<=lr->nlines){
	fprintf(stderr,"ERROR: No data to analyze.. x0=%d x1=%d nlines=%d\n",
		lr->x0, lr->x1, lr->nlines);
	return;
	}

if(nb<=lr->nlittle){
	fprintf(stderr,"ERROR: No data to analyze.. x0=%d x1=%d nlittle=%d\n",
		lr->x0, lr->x1, lr->nlittle);
	return;
	}

tmp=do_alloc(nb, sizeof(*tmp));
memcpy(tmp, z+lr->x0, nb*sizeof(*tmp));
qsort(tmp, nb, sizeof(*tmp), double_cmp);

median=tmp[nb>>1];
qlines=tmp[nb-lr->nlines];
qmost=tmp[nb-lr->nlittle];

unit_strength=qmost-median;
if(unit_strength<=0.0)unit_strength=1.0;

/* first pass - mark suspision bins */
for(i=lr->x0;i<=lr->x1;i++){
	lr->lines[i]=0;
	if(z[i]>qmost)lr->lines[i]|=LINE_HIGH;
	if(z[i]>qlines)lr->lines[i]|=LINE_CANDIDATE;
	if(z[i]>(2*qmost-median)){
		lr->lines[i]|=LINE_VERY_HIGH;
		vh_count++;
		}
	}

/* second pass - mark clustered lines */
for(i=lr->x0;i<=lr->x1;i++){
	if(lr->lines[i] & LINE_HIGH){
		for(j=i;(j<=lr->x1)&&(lr->lines[j]&LINE_HIGH);j++);
		if((j-i)>=5){
			for(k=i;k<j;k++)lr->lines[j]|=LINE_CLUSTERED;
			} else
		if((i>lr->x0) && (i<lr->x1) && ((z[i]-median)>2*(z[i-1]-median)) && ((z[i]-median)>2*(z[i+1]-median)))
			lr->lines[i]|=LINE_ISOLATED;

		}
	}
	
/* third pass - decide what will be considered to be a line */
for(i=lr->x0+1;i<=lr->x1-1;i++){
	 if(lr->lines[i]){
	 	/* Note: strength is in units of unit_strength, 
		   which is something like 2-3 sigma, 
		   depending on what was passed in */
	 	fprintf(stderr, "line detected: bin=%d z=%g strength=%f flag=0x%08x\n", i, z[i], (z[i]-median)/unit_strength,lr->lines[i]);
	 	fprintf(LOG, "line detected: bin=%d z=%g strength=%f flag=0x%08x\n", i, z[i], (z[i]-median)/unit_strength, lr->lines[i]);
		}		

	/* be very conservative: mark only points which are LINE_VERY_HIGH
	   and have both neighbours  below 0.4 level of the center line (i.e. side lobes (due to Hann window)
	    of strictly bin-centered line */
	 if(((lr->lines[i] & (LINE_CANDIDATE | LINE_VERY_HIGH | LINE_CLUSTERED))==(LINE_CANDIDATE | LINE_VERY_HIGH))){
		lr->lines[i]|=LINE_YES;
		fprintf(stderr,"line detected: i=%d z[i]=%g\n", i, z[i]);
	
		for(j=0;j<lr->nlines;j++){
			/* this can happen if we are adding new lines */
			if(lr->lines_list[j]==i)break;
			/* empty slot */
			if(lr->lines_list[j]<0){
				lr->lines_list[j]=i;
				break;
				}
			}
		}
	}

lr->median=median;
lr->qmost=qmost;
lr->qlines=qlines;
fprintf(stderr,"median=%g qmost=%g qlines=%g\n", median, qmost, qlines);
fprintf(stderr,"bins marked \"very high\": %d\n", vh_count);

free(tmp);
}

void print_lines_report(FILE *f,LINES_REPORT *lr,char *tag)
{
int i;
fprintf(f,"%s median: %g\n",tag,lr->median);
fprintf(f,"%s qmost : %g\n",tag,lr->qmost);
fprintf(f,"%s qlines: %g\n",tag,lr->qlines);
fprintf(f,"%s cutoff: %g\n",tag,2*lr->qmost-lr->median);
for(i=0;i<lr->nlines;i++)
	if(lr->lines_list[i]>=0)
		fprintf(f,"%s line bin: %d\n",tag,lr->lines_list[i]);
for(i=0;i<lr->nlines;i++)
	if(lr->lines_list[i]>=0)
		fprintf(f,"%s line freq: %g Hz\n",tag,(first_bin+lr->lines_list[i])/1800.0);
}

void detect_background_lines(double *mean)
{
LINES_REPORT lr;
int i;
/* there is a reason why we use global arrays here: 
   It is important to avoid another inner cycle when accumulating
   mean. By having a fixed array with a fixed number of 
   background lines and a hardcoded if we have only one
   if statement which is almost never taken - a good situation for
   todays cpus. A for loop would be more expensive and I am not certain
   the compiler would be smart enough to optimize a fixed loop away
   
   Note: alternative approach would be to accumulate line detection value
   using a for loop with hardcoded number of lines. The hardcoding can be
   done as a #define so that the compiler sees a small fixed loop and unwraps it.

   The latter requires testing and will certainly not happen without optimization
   flags turned on.
    */
lr.nlines=5;
lr.lines_list=lines_list;
lines=do_alloc(nbins, sizeof(*lines));
lr.lines=lines;
lr.x0=0;
lr.x1=nbins-1;
lr.nlittle=30;
for(i=0;i<lr.nlines;i++)lines_list[i]=-1; /* no line */

detect_lines_d(mean, &lr);
print_lines_report(LOG,&lr,"background");

if(!args_info.filter_lines_arg){
	fprintf(stderr,"Background lines removal will not be performed\n");
	fprintf(LOG, "Background lines removal will not be performed\n");
	for(i=0;i<lr.nlines;i++)lines_list[i]=-1; /* no line */
	}
fflush(LOG);
}
