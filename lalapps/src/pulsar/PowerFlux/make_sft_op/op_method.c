
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <string.h>

#define __USE_ISOC99 1
#include <math.h>

#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/RealFFT.h>

#include "make_sft_global.h"
#include "op_method.h"

/* set window_sum to 1.0 when no windowing is used 
   this function modifies in place *data_in */
void compute_test_fft(REAL4Vector *data_in, COMPLEX8Vector **phi, int approx, long freq_start, long freq_stop)
{
LALStatus status={level:0, statusPtr: NULL};
RealFFTPlan *fft_plan=NULL;
COMPLEX8Vector *fft=NULL;
REAL4 *p;
int i;
long L,k;

L=(data_in)->length;

LALCCreateVector(&status, &fft, data_in->length/2+1);
TESTSTATUS(&status);

LALCreateForwardRealFFTPlan(&status, &fft_plan, L, 0);
TESTSTATUS(&status);

for(i=0;i<approx;i++){
	fprintf(stderr, "Computing phi[%d]\n", i);
	if(i>0){
		p=(data_in)->data;
		/* multiply input data by \phi_1(s) */		
		for(k=0;k<L;k++){
			*p=((*p)*k)/L;
			p++;
			}
		}
	
	LALForwardRealFFT(&status, fft, data_in, fft_plan);
	TESTSTATUS(&status);
	
	memcpy(phi[i]->data, &(fft->data[freq_start]), (freq_stop-freq_start)*8);
	}
LALCDestroyVector(&status, &fft);

LALDestroyRealFFTPlan(&status, &fft_plan);
TESTSTATUS(&status);
}

/* This is only for approx=3, we cannot do larger because of memory constraints */
void compute_calibrated_periodogram3(COMPLEX8Vector **phi, long freq_start, long freq_stop, COMPLEX8Vector *pgram, double window_sum)
{
long i;
int j,k;
double a,b,c,d;
PHI_DATA3 phi_data;
PHI_RESPONSE3 phi_r;
double max1,max3;

compute_Phi(&phi_data);

for(i=0;i<pgram->length;i++){
	compute_phi_r(&phi_r, i+freq_start);
	if(!(i % (10*1800))){
		max1=0.0;
		max3=0.0;
		}
	compute_estimation_errors(&phi_data,&phi_r,i+freq_start, &a, &b);
	if(max1<a)max1=a;
	if(max3<b)max3=b;
	if(!((i+1) % (10*1800))){
		fprintf(stderr,"max1=%f max3=%f\n", max1, max3);
		}
	a=0.0;
	b=0.0;
	for(j=0;j<3;j++)
		for(k=0;k<3;k++){
			c=phi_data.inverse[j][k]*phi_r.phi_r_re[k];
			d=phi_data.inverse[j][k]*phi_r.phi_r_im[k];
			a+=phi[j]->data[i].re*c-phi[j]->data[i].im*d;
			b+=phi[j]->data[i].re*d+phi[j]->data[i].im*c;
			}
	pgram->data[i].re=a/window_sum;
	pgram->data[i].im=b/window_sum;
	}
}
