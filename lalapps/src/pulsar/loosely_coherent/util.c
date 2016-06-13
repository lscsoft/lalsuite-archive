/*
*  Copyright (C) 2007 Vladimir Dergachev
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/select.h>
#include <unistd.h>
#include <sys/syscall.h>

#include <gsl/gsl_sf_trig.h>

#include "global.h"
#include "util.h"

extern FILE *LOG;

void condor_safe_sleep(int seconds)
{
struct timeval timeout;
timeout.tv_sec=seconds;
timeout.tv_usec=1;
select(0, NULL, NULL, NULL, &timeout);
}


int direct_fcntl(int fd, int cmd, void *arg)
{
    int ret = -1;
    int unmapped_fd;

#ifdef CONDOR
    /* This call would normally go through the file table. However, we want
        to do it locally here, so get the unmapped fd the kernel originally
        gave us for this fd and use that directly. */
    unmapped_fd = _condor_get_unmapped_fd(fd);
#else 
    unmapped_fd = fd;
#endif
    ret = syscall(SYS_fcntl, unmapped_fd, cmd, arg);

    return ret;
}

int fast_fseek(FILE *stream, long offset, int whence)
{
long delta;
char *buf;
if(whence==SEEK_SET) {
	delta=ftell(stream);
	if(delta>offset)
		return(fseek(stream, offset, whence));
	if(delta==offset)return 0;
	delta=offset-delta;
	} else
if(whence==SEEK_CUR) {
	if(offset<0)
		return(fseek(stream, offset, whence));
	if(offset==0)return 0;
	delta=offset;
	} else {
	return(fseek(stream, offset, whence));
	}
if(delta>1000)
	return(fseek(stream, offset, whence));
buf=alloca(delta+1);
return(-(fread(buf, 1, delta, stream)!=delta));
}

void locate_arg(char *line, int length, int arg, int *arg_start, int *arg_stop)
{
int i;
int k;
*arg_start=0;
*arg_stop=0;
for(i=0,k=0;i<length;i++) {
	if((line[i]=='\n') || (line[i]=='\r'))break;
	if((line[i]==' ') || (line[i]=='\t')) {
		k++;
		while((i<length) && ((line[i]==' ') || (line[i]=='\t')))i++;
		}
	if(k==arg)break;
	}
if(k!=arg)return;
*arg_start=i;
while((i<length) && !((line[i]==' ') || (line[i]=='\t') || (line[i]=='\n')))i++;
*arg_stop=i;
if(line[*arg_start]=='"')(*arg_start)++;
if(((*arg_stop)>(*arg_start)) && (line[(*arg_stop)-1]=='"'))(*arg_stop)--;
}

float hann_response(float delta)
{
if(fabs(delta-1.0)<0.5) {
	return(gsl_sf_sinc(delta-1.0)/(delta*(1+delta)));
	}
if(fabs(delta+1.0)<0.5) {
	return(gsl_sf_sinc(delta+1.0)/(delta*(delta-1)));
	}
if(fabs(delta)<0.5) {
	return(gsl_sf_sinc(delta)/(1-delta*delta));
	}
return(sin(M_PI*delta)/(M_PI*delta*(1-delta*delta)));
}


void fill_hann_filterN(float *coeffs, int length, int middle, float mismatch)
{
int i;
for(i=0;i<length;i++) {
	coeffs[i]=hann_response(mismatch-i+middle);
	/* change sign according to e^{pi i delta} factor - we don't care about common phase but we care about -1^n */
	if((i+middle)&1)coeffs[i]=-coeffs[i];
	}
}

/* this is the inverse of covariance matrix of hann windowed gaussian noise  (times 4.0) */

// float R_inv[7][7]={
// 	7,      6,      5,      4,      3,      2,      1,
// 	6,      12,     10,     8,      6,      4,      2,
// 	5,      10,     15,     12,     9,      6,      3,
// 	4,      8,      12,     16,     12,     8,      4,
// 	3,      6,      9,      12,     15,     10,     5,
// 	2,      4,      6,      8,      10,     12,     6,
// 	1,      2,      3,      4,      5,      6,      7
// 	};

float R_inv[7][7]={
	336,    504,    540,    480,    360,    216,    84,
	504,    1071,   1260,   1170,   900,    549,    216,
	540,    1260,   1800,   1800,   1440,   900,    360,
	480,    1170,   1800,   2100,   1800,   1170,   480,
	360,    900,    1440,   1800,   1800,   1260,   540,
	216,    549,    900,    1170,   1260,   1071,   504,
	84,     216,    360,    480,    540,    504,    336
	};

void fill_hann_filter7(float *coeffs, float mismatch)
{
int i, j;
float c[7];
float norm;
for(i=0;i<7;i++) {
	c[i]=hann_response(mismatch-(i-3));
	/* change sign according to e^{pi i delta} factor - we don't care about common phase but we care about -1^n */
	if((i+3)&1)c[i]=-c[i];
	}
norm=0;
for(i=0;i<7;i++) {
	coeffs[i]=0;
	for(j=0;j<7;j++)coeffs[i]+=R_inv[i][j]*c[j];
	
	//coeffs[i]*=1.0/135.0;
	norm+=coeffs[i]*c[i];
	}
//fprintf(stderr, "mismatch=%f norm=%g\n", mismatch, norm);
norm=1.0/norm;
for(i=0;i<7;i++)coeffs[i]*=norm;
}

void fill_diff_hann_filter7(float *coeffs, float mismatch)
{
float c1[7], c2[7];
int i;

fill_hann_filter7(c1, mismatch-0.001);
fill_hann_filter7(c2, mismatch+0.001);
for(i=0;i<7;i++)coeffs[i]=2000.0*(c2[i]-c1[i]);
}

#define N 255
float filter7_table[2*N+1][7];

void tabulate_hann_filter7(void)
{
int i,j;
float f, max;
for(i=0;i<2*N+1;i++) {
	fill_hann_filter7(filter7_table[i], (i-N)/(1.0*N));
	}

max=0;
for(i=0;i<2*N;i++) {
	for(j=0;j<7;j++) {
		f=fabs(filter7_table[i][j]-filter7_table[i+1][j]);
		if(f>max)max=f;
		}
	}
fprintf(stderr, "maximum filter7 sampling error: %g\n", max);
fprintf(LOG, "maximum filter7 sampling error: %g\n", max);
}

void tabulated_fill_hann_filter7(float *coeffs, float mismatch)
{
int k;
k=rintf((mismatch+1.0)*N);
if((k<0)||(k>=(2*N+1))) {
	fprintf(stderr, "mismatch value (%g) outside sampled range requiested.. Bug ?\n", mismatch);
	fill_hann_filter7(coeffs, mismatch);
	return;
	}
memcpy(coeffs, filter7_table[k], 7*sizeof(float));
}

void test_hann(void)
{
int i;
float x;
fprintf(stderr, "hann_response:\n");
for(i=0;i<1000;i++) {
	x=(i-500)*0.01;
	fprintf(stderr, "%f %f\n",x, hann_response(x));
	}
}

void test_filter(void) 
{
float filter[7];
float mismatch;
int i;
for(i=0;i<=100;i++) {
	mismatch=(i-50)*0.02;
	fill_hann_filterN(filter, 7, 3, mismatch);

	fprintf(stderr, "%d %f %f %f %f %f %f %f %f\n", 
		i, mismatch,
		filter[0],
		filter[1],
		filter[2],
		filter[3],
		filter[4],
		filter[5],
		filter[6]
		);
	}
}

VARRAY *new_varray(int item_size)
{
VARRAY *v;
v=do_alloc(1, sizeof(*v));
v->item_size=item_size;
v->free=0;
v->size=100;
v->data=do_alloc(v->size, v->item_size);
return v;
}

void free_varray(VARRAY *v)
{
free(v->data);
v->data=NULL;
free(v);
}

int varray_add(VARRAY *v, void *item)
{
if(v->free>=v->size) {
	void *p;
	v->size=2*v->size+10;
	p=do_alloc(v->size, v->item_size);
	if(v->free>0)memcpy(p, v->data, v->free*v->item_size);
	if(v->data!=NULL)free(v->data);
	v->data=p;
	}
memcpy(& ((((char *)(v->data))[v->free*v->item_size])), item, v->item_size);
v->free++;
return(v->free-1);
}

int is_round235(int n)
{
while(n>1 && (n % 2==0))n=n/2;
while(n>1 && (n % 3==0))n=n/3;
while(n>1 && (n % 5==0))n=n/5;
if(n==1)return 1;
return 0;
}

/* make an integer round by being a product of 2 3 and 5 */
int round235up_int(int n)
{
while(!is_round235(n))n++;
return(n);
}

int round235down_int(int n)
{
while(!is_round235(n))n--;
return(n);
}
