#include <stdio.h>
#include <stdlib.h>

#define __USE_ISOC99	1
#include <math.h>

#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/RealFFT.h>
#include <lal/Window.h>
#include <lal/BandPassTimeSeries.h>

#include "make_sft_global.h"
#include "op_method.h"

/* free refers to the first unoccupied entry */
COMPLEX8 *R_0=NULL, *H_0=NULL, *C_0=NULL;
long R_size=-1, R_free=0, C_size=-1, C_free=0, H_size=-1, H_free=0;
double R_rate, C_rate;
long cal_sft_length=0;

/* from make_sft.c */
extern long samples_per_second;
extern INT64 gps_start;
extern INT64 gps_end;
extern long total_samples;
extern int dont_fail_on_missing_calibration;
extern int dont_fail_on_zero_ds;

void load_response_file(char *name, COMPLEX8 **data, long *_size, long *_free, double *rate)
{
FILE *f;
unsigned char s[2000];
REAL4 magnitude, phase, freq;
int count;
COMPLEX8 *p;
if(*_size<0){
	*_size=cal_sft_length/2+1; /* a good value to start with - this is what current files use at the moment */
	*_free=0;
	*data=do_alloc(*_size, sizeof(**data));
	}
fault_file(name);
f=fopen(name, "r");
if(f==NULL){
	fprintf(stderr, "** Could not open \"%s\" for reading:", name);
	perror("");
	exit(-1);
	}
*_free=0;
count=0;
while(!feof(f)){
	fgets(s, 2000, f);
	if((s[0]=='#')||(s[0]=='%')||(s[0]==0))continue; /* skip comments and empty lines */
	/* we assume that entries in response files are from 0 up to whatever
	   last bin is there */
	sscanf(s, "%g %g %g", &freq, &magnitude, &phase);
	count++;
	if(*_free>=*_size){
		*_size*=2;
		p=do_alloc(*_size, sizeof(**data));
		if(_free>0)memcpy(p, (*data), *_free*sizeof(**data));
		free(*data);
		(*data)=p;
		p=NULL;
		}
	if(isnan(magnitude)||isnan(phase)){
		(*data)[*_free].re=0.0;
		(*data)[*_free].im=0.0;
		(*_free)++;
		fprintf(stderr,"** NaN in (*data) bin %ld\n", *_free);
		continue;
		}
	(*data)[*_free].re=magnitude*cos(phase);
	(*data)[*_free].im=magnitude*sin(phase);
	(*_free)++;
	}
*rate=(count-1)/freq;
fclose(f);
}

int resample_array(COMPLEX8 **data, long new_size, long entries, double A_factor, double B_factor)
{
COMPLEX8 *p;
long i,k;
REAL4 x,delta;
p=do_alloc(new_size, sizeof(*p));
for(i=0;i<new_size;i++){
	p[i].re=0.0;
	p[i].im=0.0;
	}
for(i=0;i<new_size;i++){
	x=(i*A_factor)/B_factor;
	k=x;
	delta=x-k;
	if(k<0){
		/* this could occur if response files do not start with 0 - but they do */
		fprintf(stderr,"** Unexpected situation when rescaling array: k<0\n");
		fprintf(stderr,"**    detail: k=%ld x=%g A=%g B=%g i=%ld\n", k, x, A_factor, B_factor, i);
		continue;
		}
	if(k>=entries-1){
		/* we run out of entries */
		break;
		}
	p[i].re=(*data)[k].re*(1.0-delta)+(*data)[k+1].re*delta;
	p[i].im=(*data)[k].im*(1.0-delta)+(*data)[k+1].im*delta;
	}
free(*data);
*data=p;
/* return number of valid entries assigned.. */
return i;
}

void load_R(char *name)
{
cal_sft_length=samples_per_second*60;
if(R_free>0){
	fprintf(stderr,"** ERROR: Attempt to override R file with file \"%s\"\n", name);
	exit(-1);
	}
fprintf(stderr,"Attempting to load R_file \"%s\"\n", name);
load_response_file(name, &R_0, &R_size, &R_free, &R_rate);
fprintf(stderr,"Loaded R_file: %d samples, rate=%f\n", R_free, R_rate);
}

void load_C(char *name)
{
cal_sft_length=samples_per_second*60;
if(C_free>0){
	fprintf(stderr,"** ERROR: Attempt to override C file with file \"%s\"\n", name);
	exit(-1);
	}
fprintf(stderr,"Attempting to load C_file \"%s\"\n", name);
load_response_file(name, &C_0, &C_size, &C_free, &C_rate);
fprintf(stderr,"Loaded C_file: %d samples, rate=%f\n", C_free, C_rate);
}

void post_init_response_files(void)
{
long i;
REAL4 ar,ai,br,bi;

if((R_free<=0)){
	if(!dont_fail_on_missing_calibration){
		fprintf(stderr, "** R constants were not loaded, exit as requested.\n");
		exit(-1);
		}
	fprintf(stderr, "** R constants were not loaded, no calibration will be performed.\n");
	return;
	}
if(C_free<=0){
	if(!dont_fail_on_missing_calibration){
		fprintf(stderr, "** C constants were not loaded, exit as requested.\n");
		exit(-1);
		}
	fprintf(stderr, "** C constants were not loaded, no calibration will be performed.\n");
	return;
	}

cal_sft_length=samples_per_second*60;
/* if we got extra entries truncate data */
if(R_free>C_free){
	fprintf(stderr, "** extra entries in R file or missing entries in C file\n");
	R_free=C_free;
	} else
if(R_free<C_free){
	fprintf(stderr, "** extra entries in C file or missing entries in R file\n");
	C_free=R_free;
	}
if(R_rate!=C_rate){
	fprintf(stderr, "** C_rate!=R_rate, diff=%g\n", R_rate-C_rate);
	if(fabs(R_rate-C_rate)>1e-6){
		fprintf(stderr, "*** Large mismatch in C_rate versus R_rate, aborting\n");
		exit(-1);
		}
	R_rate=C_rate;
	}

H_size=cal_sft_length/2+1;
if(H_size<R_free)H_size=R_free;
H_0=do_alloc(H_size, sizeof(*H_0));
H_free=R_free;
for(i=0;i<R_free;i++){
	ar=R_0[i].re;
	ai=R_0[i].im;
	br=C_0[i].re;
	bi=C_0[i].im;
	H_0[i].re=ar*br-ai*bi-1.0;
	H_0[i].im=ar*bi+ai*br;
	#if 0
	printf("R=[%g,%g] C=[%g,%g] H=[%g,%g]\n", ar, ai, br, bi, H_0[i].re, H_0[i].im);
	#endif
	}
C_free=resample_array(&C_0, cal_sft_length/2+1, C_free, C_rate, 60);
H_free=resample_array(&H_0, cal_sft_length/2+1, H_free, C_rate, 60);
C_size=cal_sft_length/2+1;
H_size=C_size;
}

typedef struct {
	/* the first sample this constant applies to */
	long samples_start;
	REAL4 alpha;
	REAL4 alphabeta;
	} ALPHABETA;

/* ab_ptr is used for scanning the array */	
ALPHABETA *ab_data=NULL;
long ab_free=0, ab_size=-1, ab_ptr=0;


void add_alpha_beta(INT64 gps, REAL4 alpha, REAL4 alphabeta)
{
ALPHABETA *p;

if(total_samples<0){
	fprintf(stderr, "** ALPHA_BETA_* commands should follow SAMPLES_PER_SECOND and DURATION commands\n");
	exit(-1);
	}
/* do not store alpha/beta values outside of processing segment */
if(gps>=(gps_end+60))return;
if(gps<(gps_start-60))return;

cal_sft_length=samples_per_second*60;
if(ab_free>=ab_size){
	ab_size=2*ab_size+10;
	p=do_alloc(ab_size, sizeof(*p));
	if(ab_free>0)memcpy(p, ab_data, ab_free*sizeof(*p));
	if(ab_data!=NULL)free(ab_data);
	ab_data=p;
	}
p=&(ab_data[ab_free]);
p->samples_start=(gps-gps_start)*samples_per_second;

fprintf(stderr,"Adding gps=%lld alpha=%g alphabeta=%g -> samples_start=%ld\n", gps, alpha, alphabeta, p->samples_start);

p->alpha=alpha;
p->alphabeta=alphabeta;
ab_free++;

if((alpha<0.01)||(alphabeta<0.01)){
	if(!dont_fail_on_missing_calibration){
		fprintf(stderr, "** Bogus alpha/beta constants found for sample %ld (GPS time %lld), exit as requested.\n", p->samples_start, gps);
		exit(-1);
		}
	}
}

void load_alpha_beta(char *name)
{
FILE *f;
unsigned char s[2000];
REAL4 alpha, alphabeta;
INT64 gps;

fprintf(stderr,"Attempting to load alpha/beta file \"%s\"\n", name);
fault_file(name);
f=fopen(name, "r");
if(f==NULL){
	fprintf(stderr, "** Could not open \"%s\" for reading:", name);
	perror("");
	exit(-1);
	}
while(!feof(f)){
	fgets(s, 2000, f);
	if((s[0]=='#')||(s[0]=='%')||(s[0]==0))continue; /* ignore comments and empty lines */
	sscanf(s, "%lld %g %g", &gps, &alphabeta, &alpha);
	add_alpha_beta(gps, alpha, alphabeta);
	}
fclose(f);
}

int ab_cmp(ALPHABETA *A, ALPHABETA *B)
{
if(A->samples_start>B->samples_start)return 1;
if(A->samples_start<B->samples_start)return -1;
return 0;
}

void post_init_alpha_beta(void)
{
if(ab_free<=0){
	if(!dont_fail_on_missing_calibration){
		fprintf(stderr, "** No relevant alpha/beta were loaded, exit as requested.\n");
		exit(-1);
		}
	fprintf(stderr, "** No relevant alpha/beta constants were loaded, assuming alpha=beta=1.0.\n");
	return;
	}
/* sort alpha/beta values .. They should be sorted already but one more time does not hurt */
qsort(ab_data, ab_free, sizeof(*ab_data), ab_cmp);
#if 1
	/* check that the order of the data is increasing with samples_start */
	fprintf(stderr,"ab_data[].samples_start= %ld %ld %ld\n", ab_data[0].samples_start, ab_data[1].samples_start, ab_data[2].samples_start);
#endif
/* set pointer to the first entry */
ab_ptr=0;
}

/* find pair of ALPHABETA structures which samples_start numbers bracket
   first_sample.

   For now, this is simple linear search, good for case of sequential reading.
   
   Also, 1800/60=30 - we should not have more than 31 ALPHABETA entries.   
   
  */
void position_alpha_beta_ptr(long first_sample)
{
if(ab_ptr>=(ab_free-1))ab_ptr=ab_free-2;
if(ab_ptr<0)ab_ptr=0;
while((ab_ptr<(ab_free-1))&& (ab_data[ab_ptr+1].samples_start<first_sample)){
	ab_ptr++;
	}
while((ab_ptr>0) && (ab_data[ab_ptr].samples_start>first_sample)){
	ab_ptr--;
	}
}

/* Note that, when processing interpolated data, the alpha and beta constants can very well be missing.
   However, at the moment, we are only dealing with contiguous segments */

void get_alpha_beta(long first_sample, REAL4 *alpha, REAL4 *alphabeta)
{
long x;
if(ab_free<=0){
	/* no constants have been loaded */
	*alpha=1.0;
	*alphabeta=1.0;
	return;
	}
if(ab_free==1){
	*alpha=ab_data->alpha;
	*alphabeta=ab_data->alphabeta;
	return;
	}
position_alpha_beta_ptr(first_sample);
/* Are we too far to the left ? */
if(first_sample<ab_data[ab_ptr].samples_start){
	if(!dont_fail_on_missing_calibration){
		fprintf(stderr, "** No alpha/beta constants found for sample %ld, exit as requested.\n", first_sample);
		exit(-1);
		}
	fprintf(stderr,"No alpha/beta constants found for sample %ld, using left tail\n", first_sample);
	*alpha=ab_data[ab_ptr].alpha;
	*alphabeta=ab_data[ab_ptr].alphabeta;
	return;
	}
/* Are we too far to the right ? */
x=ab_data[ab_ptr+1].samples_start-first_sample;
if(x<0){
	if(!dont_fail_on_missing_calibration){
		fprintf(stderr, "** No alpha/beta constants found for sample %ld, exit as requested.\n", first_sample);
		exit(-1);
		}
	fprintf(stderr,"No alpha/beta constants found for sample %ld, using right tail\n", first_sample);
	*alpha=ab_data[ab_ptr+1].alpha;
	*alphabeta=ab_data[ab_ptr+1].alphabeta;
	return;
	}
if(x>cal_sft_length)x=cal_sft_length;
/* interpolate alpha/beta using length of calibrated sft */
*alpha=(ab_data[ab_ptr].alpha*x+ab_data[ab_ptr+1].alpha*(cal_sft_length-x))/cal_sft_length;
*alphabeta=(ab_data[ab_ptr].alphabeta*x+ab_data[ab_ptr+1].alphabeta*(cal_sft_length-x))/cal_sft_length;
/* check for validity */
if((*alpha<0.01)||(*alphabeta<0.01)){
	if(!dont_fail_on_missing_calibration){
		fprintf(stderr, "** Bogus alpha/beta constants found for sample %ld, exit as requested.\n", first_sample);
		exit(-1);
		}
	}
}

void calibrate_fft(COMPLEX8Vector *fft, long first_sample)
{
REAL4 alpha, alphabeta,nr,ni,dr,di,ds,xi,xr,yi,yr;
int i;
if((R_free<=0)||(C_free<=0)){
	/* no calibration constants have been loaded.. nothing to do */
	return;
	}
get_alpha_beta(first_sample, &alpha, &alphabeta);
for(i=0;i<H_free;i++){
	/* compute denominator */
	dr=C_0[i].re*alpha;
	di=C_0[i].im*alpha;
	ds=dr*dr+di*di;
	/* the following is possibly bogus.. better checking for 0 ? 
	   note that we are supposed to be windowing the result anyway, 
	   so, presumably, we don't care about NaNs in odd places.
		*/
	if(ds<=0.0){
		if(!dont_fail_on_zero_ds){
			fprintf(stderr,"ds=0.0 for i=%d\n", i);
			exit(-1);
			}
		fft->data[i].re=0.0;
		fft->data[i].im=0.0;
		continue;
		}

	/* compute numerator */
	nr=H_0[i].re*alphabeta+1;
	ni=H_0[i].im*alphabeta;
	
	/* divide, i.e. multiply by C conjugate and divide by |C|^2 */
	xr=(nr*dr+ni*di)/ds;
	xi=(ni*dr-nr*di)/ds;
	
	/* copy old data. Counting on smart compiler to put this into registers */
	yr=fft->data[i].re;
	yi=fft->data[i].im;
	
	/* compute result */
	fft->data[i].re=yr*xr-yi*xi;
	fft->data[i].im=yr*xi+yi*xr;
	}
}

void compute_Phi(PHI_DATA3 *phi_data)
{
int i,k;
double a,b;
double *p;
if(ab_free<=0){
	fprintf(stderr, "** ERROR: no alpha/beta constants were loaded, cannot compute Phi\n");
	exit(-1);
	}
for(i=0;i<5;i++)phi_data->a[i]=0.0;
phi_data->a[0]=ab_free; /* N */
for(i=0;i<ab_free;i++){
	a=((2*ab_data[i].samples_start-total_samples)*1.0)/total_samples;
	b=a*a;
	phi_data->a[1]+=a;
	phi_data->a[2]+=b;
	phi_data->a[3]+=b*a;
	phi_data->a[4]+=b*b;
	}
p=phi_data->a;
phi_data->det=(p[0]*p[2]*p[4]-p[0]*p[3]*p[3]-p[1]*p[1]*p[4]+2.0*p[1]*p[2]*p[3]-p[2]*p[2]*p[2]);
if(phi_data->det==0.0){
	fprintf(stderr,"** ERROR: Determinant of Phi is zero\n");
	exit(-1);
	}
phi_data->inverse[0][0]=p[2]*p[4]-p[3]*p[3];
phi_data->inverse[0][1]=p[2]*p[3]-p[1]*p[4];
phi_data->inverse[0][2]=p[1]*p[3]-p[2]*p[2];
phi_data->inverse[1][0]=phi_data->inverse[0][1];
phi_data->inverse[1][1]=p[0]*p[4]-p[2]*p[2];
phi_data->inverse[1][2]=p[1]*p[2]-p[0]*p[3];
phi_data->inverse[2][0]=phi_data->inverse[0][2];
phi_data->inverse[2][1]=phi_data->inverse[1][2];
phi_data->inverse[2][2]=p[0]*p[2]-p[1]*p[1];
for(i=0;i<3;i++)
	for(k=0;k<3;k++)phi_data->inverse[i][k]/=phi_data->det;
fprintf(stderr,"Computed Phi: ab_free=%ld det=%g a[]=(%g,%g,%g,%g,%g)\n", 
	ab_free, phi_data->det,
	phi_data->a[0],
	phi_data->a[1],
	phi_data->a[2],
	phi_data->a[3],
	phi_data->a[4]);
}

void compute_phi_r(PHI_RESPONSE3 *phi_response, long n)
{
int i;
double alpha, alphabeta,a,b;
double dr,di,ds,nr,ni,xr,xi;
for(i=0;i<3;i++)phi_response->phi_r_re[i]=0.0;
for(i=0;i<3;i++)phi_response->phi_r_im[i]=0.0;
n=(n+15)/30; /* we only have data every 60 seconds, not every 1800 */
if((n)>=H_free)return; /* this actually hardcodes the length of the segment as 1800 sec */

for(i=0;i<ab_free;i++){
	/* load alpha/beta data */
	a=((2*ab_data[i].samples_start-total_samples)*1.0)/total_samples;
	b=a*a;
	alpha=ab_data[i].alpha;
	alphabeta=ab_data[i].alphabeta;

	/* compute denominator */
	dr=C_0[n].re*alpha;
	di=C_0[n].im*alpha;
	ds=dr*dr+di*di;
	/* the following is possibly bogus.. better checking for 0 ? 
	   note that we are supposed to be windowing the result anyway, 
	   so, presumably, we don't care about NaNs in odd places.
		*/
	if(ds<=0.0){
		for(i=0;i<3;i++)phi_response->phi_r_re[i]=0.0;
		for(i=0;i<3;i++)phi_response->phi_r_im[i]=0.0;
		fprintf(stderr,"ds=0.0 for n=%ld\n", n);
		return;		
		}

	/* compute numerator */
	nr=H_0[n].re*alphabeta+1;
	ni=H_0[n].im*alphabeta;
	
	/* divide, i.e. multiply by C conjugate and divide by |C|^2 */
	xr=(nr*dr+ni*di)/ds;
	xi=(ni*dr-nr*di)/ds;

	phi_response->phi_r_re[0]+=xr;	
	phi_response->phi_r_re[1]+=xr*a;	
	phi_response->phi_r_re[2]+=xr*b;	

	phi_response->phi_r_im[0]+=xi;	
	phi_response->phi_r_im[1]+=xi*a;	
	phi_response->phi_r_im[2]+=xi*b;	
	}
}

void compute_estimation_errors(PHI_DATA3 *phi, PHI_RESPONSE3 *phi_r, long nbin, double *max1, double *max3)
{
int i,j,k,n;
double a,b,c,d;
double dr,di,ds,nr,ni,xr,xi,alpha,alphabeta,est3_r,est3_i,est1_r,est1_i,err,mag;
double m[3];
*max3=0;
*max1=0;
n=(nbin+15)/30; /* we only have data every 60 seconds, not every 1800 */
if((n)>=H_free)return; /* this actually hardcodes the length of the segment as 1800 sec */
for(i=0;i<ab_free;i++){
	/* load alpha/beta data */
	a=((2*ab_data[i].samples_start-total_samples)*1.0)/total_samples;
	b=a*a;
	alpha=ab_data[i].alpha;
	alphabeta=ab_data[i].alphabeta;

	/* compute denominator */
	dr=C_0[n].re*alpha;
	di=C_0[n].im*alpha;
	ds=dr*dr+di*di;
	/* the following is possibly bogus.. better checking for 0 ? 
	   note that we are supposed to be windowing the result anyway, 
	   so, presumably, we don't care about NaNs in odd places.
		*/
	if(ds<=0.0){
		fprintf(stderr,"ds=0.0 for n=%ld\n", n);
		return;
		}

	/* compute numerator */
	nr=H_0[n].re*alphabeta+1;
	ni=H_0[n].im*alphabeta;
	
	/* divide, i.e. multiply by C conjugate and divide by |C|^2 */
	xr=(nr*dr+ni*di)/ds;
	xi=(ni*dr-nr*di)/ds;


	mag=sqrt(xr*xr+xi*xi);

	m[0]=1;
	m[1]=a;
	m[2]=b;
	
	/* compute the estimate using averaging */
	est1_r=phi_r->phi_r_re[0]/ab_free;
	est1_i=phi_r->phi_r_im[0]/ab_free;
	err=fabs(est1_r-xr)/mag;
	if(*max1<err)*max1=err;
	err=fabs(est1_i-xi)/mag;
	if(*max1<err)*max1=err;
	
	
	est3_r=0;
	est3_i=0;
	/* now compute the estimate from the data using quadratic interpolation */
	for(j=0;j<3;j++)
		for(k=0;k<3;k++){
			est3_r+=m[j]*phi->inverse[j][k]*phi_r->phi_r_re[k];
			est3_i+=m[j]*phi->inverse[j][k]*phi_r->phi_r_im[k];
			}
	
	err=fabs(est3_r-xr)/mag;
	if(*max3<err)*max3=err;
	err=fabs(est3_i-xi)/mag;
	if(*max3<err)*max3=err;
	if(nbin==1800*1400){
		/*
		fprintf(stderr,"xr=%g xi=%g\n", xr, xi);
		fprintf(stderr,"est1_r=%g est1_i=%g\n", est1_r, est1_i);
		fprintf(stderr,"est3_r=%g est3_i=%g\n", est3_r, est3_i);
		fprintf(stderr,"     mag=%g max1=%f max3=%f\n", mag,*max1,*max3);
		*/
		
		fprintf(stderr,"%g %g %g %g %g %g %g %g %g\n", xr, xi, est1_r, est1_i, est3_r, est3_i, mag, *max1,*max3);
		}
	}

}

void get_td_calibration(double frequency, long sample, double *re, double *im)
{
float alpha, alphabeta;
double dr,di,ds,nr,ni;
int n;
n=(frequency*60.0+0.5);
if(n>H_free){
	*re=0.0;
	*im=0.0;
	return;
	}

get_alpha_beta(sample, &alpha, &alphabeta);

/* compute denominator */
dr=C_0[n].re*alpha;
di=C_0[n].im*alpha;
ds=dr*dr+di*di;
/* the following is possibly bogus.. better checking for 0 ? 
   note that we are supposed to be windowing the result anyway, 
   so, presumably, we don't care about NaNs in odd places.
	*/
if(ds<=0.0){
	*re=0.0;
	*im=0.0;
	if(!dont_fail_on_zero_ds){
		fprintf(stderr,"ds=0.0 for n=%d\n", n);
		exit(-1);
		}
	return;		
	}

/* compute numerator */
nr=H_0[n].re*alphabeta+1;
ni=H_0[n].im*alphabeta;
	
/* divide, i.e. multiply by C conjugate and divide by |C|^2 */
*re=(nr*dr+ni*di)/ds;
*im=(ni*dr-nr*di)/ds;
}
