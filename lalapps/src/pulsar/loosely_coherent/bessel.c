#include <stdio.h>
#include <math.h>
#include <string.h>

#include <pmmintrin.h>
#include <xmmintrin.h>

#include <gsl/gsl_sf.h>
#include <lal/LALDatatypes.h>

#include "global.h"
#include "bessel.h"

extern FILE * LOG;

/* gsl functions are not precize for small arguments */
void my_bessel_Jn_array(int count, double rho, double *bessel_values)
{
int order;
int i;
double x, a;

x=rho*100;
order=floor(x*x);
if(order>42)order=42;

//fprintf(stderr, "order=%d rho=%g\n", order, rho);

if(order>count)order=count;

if(order>0)
	gsl_sf_bessel_Jn_array(0, order, rho, bessel_values);

if(order==count)return;

if(order==0)bessel_values[0]=1.0;

a=1.0;
for(i=1;i<order;i++)a=a*(0.5*rho)/i;

for(;i<=count;i++) {
	a=a*(0.5*rho)/i;
	bessel_values[i]=a;
	if(fabs(a)<1e-6) {
		memset(&(bessel_values[i]), 0, (count-i+1)*sizeof(*bessel_values));
		return;
		}
	}
}

void check_my_bessel_Jn(void)
{
int i, j, N;
double a, x, rho;
double val[256];
double err=0;
	
for(j=0;j<4000;j++) {
	rho=(j*1.0)/256.0;
	
	N=127;
	my_bessel_Jn_array(N, rho, val);
	for(i=0;i<=N;i++) {
		//fprintf(stderr, "%d %f\n", i, rho);
		if(fabs(rho)<=0) {
			if(i==0)a=1.0;
				else a=0;
			} else
		if(rho*100<sqrt(i+1) || (i>42)) {
			a=exp(i*log(0.5*rho))/gsl_sf_fact(i);
			} else
			a=gsl_sf_bessel_Jn(i, rho);
		x=fabs(a-val[i]);
		if(x>err)err=x;
		if(x>1e-5) {
			fprintf(stderr, "Bessel approximation error: i=%d rho=%.12g %g for %.12g vs %.12g\n", i, rho, x, a, val[i]);
			}
		}
	}
fprintf(stderr, "Maximum Bessel approximation error: %g\n", err);
fprintf(LOG, "Maximum Bessel approximation error: %g\n", err);
if(err>1e-5) {
	fprintf(stderr, "*** Bessel approximation error too large, exiting\n");
	exit(-1);
	}
}

void make_bessel_filter(COMPLEX8 *filter, int filter_size, COMPLEX8 *coeffs, int coeffs_size, double scale)
{
double *bessel_values;
double *cnorm;
COMPLEX8 *cuniti, *filter_tmp;
int i, j, k, m;
int offset;
COMPLEX8 a, b;
double ad;
double rho;

//coeffs_size=3;

if(!(filter_size & 1)) {
	fprintf(stderr, "*** INTERNAL ERROR: filter size should be odd\n");
	exit(-1);
	}
offset=filter_size>>1;

bessel_values=alloca((offset+1)*sizeof(*bessel_values));
filter_tmp=alloca(filter_size*sizeof(*filter_tmp));
	
cnorm=alloca(coeffs_size*sizeof(*cnorm));
cuniti=alloca(coeffs_size*sizeof(*cuniti));

for(i=0;i<coeffs_size;i++) {
	a=coeffs[i];
	cnorm[i]=sqrt(a.re*a.re+a.im*a.im);
	
	/* normalize coeffs multiplied by 1i*/
	cuniti[i].re=-a.im/cnorm[i];
	cuniti[i].im=a.re/cnorm[i];
	}

//scale=-scale;
	
if(scale<0) {
	scale=-scale;
	
	for(i=0;i<coeffs_size;i++) {
		cuniti[i].re=-cuniti[i].re;
		cuniti[i].im=-cuniti[i].im;
		}
	}
	
// for(i=0;i<coeffs_size;i++) {
// 	fprintf(stderr, "h %d %f %f %f\n", i, cnorm[i], cuniti[i].re, cuniti[i].im);
// 	}
// fprintf(stderr, "\n");

//gsl_sf_bessel_Jn_array(0, offset, scale*2*cnorm[0], bessel_values);
my_bessel_Jn_array(offset, scale*2*cnorm[0], bessel_values);

filter[offset].re=bessel_values[0];
filter[offset].im=0;

// for(i=0;i<=offset;i++) {
// 	fprintf(stderr, " %f", bessel_values[i]);
// 	}
// fprintf(stderr, "\n\n");

a=cuniti[0];
for(i=1;i<=offset;i++) {
	filter[offset+i].re=bessel_values[i]*a.re;
	filter[offset+i].im=bessel_values[i]*a.im;

	filter[offset-i].re=bessel_values[i]*a.re*(1-(i & 1)*2);
	filter[offset-i].im=-bessel_values[i]*a.im*(1-(i & 1)*2);
	
	b.re=a.re*cuniti[0].re-a.im*cuniti[0].im;
	b.im=a.im*cuniti[0].re+a.re*cuniti[0].im;
	a=b;
	}
	
j=1;
while((j<coeffs_size) && (j<offset)) {
	
/*	for(i=0;i<filter_size;i++) {
		fprintf(stderr, "%d %d %f %f\n", j+1, i, filter[i].re, filter[i].im);
		}
	fprintf(stderr, "\n");*/
	rho=scale*2*cnorm[j];
	
	if(fabs(rho)<1e-3){
		j++;
		continue;
		}

	
	memcpy(filter_tmp, filter, filter_size*sizeof(*filter));
	
	//gsl_sf_bessel_Jn_array(0, offset, scale*2*cnorm[j], bessel_values);
	memset(bessel_values, 0, (offset+1)*sizeof(*bessel_values));
	my_bessel_Jn_array(ceil(offset*1.0/(j+1)), rho, bessel_values);
	
// 	for(i=0;i<=offset;i++) {
// 		fprintf(stderr, " %f", bessel_values[i]);
// 		}
// 	fprintf(stderr, "\n\n");

	for(i=0;i<filter_size;i++) {
		filter[i].re=filter_tmp[i].re*bessel_values[0];
		filter[i].im=filter_tmp[i].im*bessel_values[0];
		}

	a=cuniti[j];
	for(i=1;i<=offset;i++) {

		if(fabs(bessel_values[i])>0.0) {
			for(k=0;k<filter_size;k++) {
				m=k+i*(j+1);
				if(m<filter_size) {
					b=filter_tmp[m];
					filter[k].re+=bessel_values[i]*(b.re*a.re+b.im*a.im)*(1-(i & 1)*2);
					filter[k].im+=bessel_values[i]*(-b.re*a.im+b.im*a.re)*(1-(i & 1)*2);
					}

				m=k-i*(j+1);
				if(m>=0) {
					b=filter_tmp[m];
					filter[k].re+=bessel_values[i]*(b.re*a.re-b.im*a.im);
					filter[k].im+=bessel_values[i]*(b.re*a.im+b.im*a.re);
					}
				}
			}
		
		b.re=a.re*cuniti[j].re-a.im*cuniti[j].im;
		b.im=a.im*cuniti[j].re+a.re*cuniti[j].im;
		a=b;
		}
	j++;
	}
	
/* Normalize filter */
ad=0.0;
for(i=0;i<filter_size;i++) {
	ad+=(double)(filter[i].re)*(double)(filter[i].re)+(double)(filter[i].im)*(double)(filter[i].im);
	}
ad=1.0/sqrt(ad);
if(fabs(ad-1.0)>0.01) {
	for(i=0;i<=offset && i<coeffs_size;i++) {
		fprintf(stderr, " %f", cnorm[i]);
		}
	fprintf(stderr, "\nbessel values\n");

	for(i=0;i<=offset;i++) {
		fprintf(stderr, " %f", bessel_values[i]);
		}
	fprintf(stderr, "\n\n");

	for(i=0;i<filter_size;i++) {
		fprintf(stderr, " (%f,%f)", filter[i].re, filter[i].im);
		}
	fprintf(stderr, "\n\n");

	fprintf(stderr, "ad=%f scale=%g coeff[1]=(%g, %g) filter_size=%d\n", ad, scale, coeffs[1].re, coeffs[1].im, filter_size);
	}
for(i=0;i<filter_size;i++) {
	filter[i].re=filter[i].re*ad;
	filter[i].im=filter[i].im*ad;
	}
/* check offset normalization */
#if 0
ad=0.0;
for(i=1;i<filter_size;i++) {
	ad+=(double)(filter[i].re)*(double)(filter[i-1].re)+(double)(filter[i].im)*(double)(filter[i-1].im);
	}
if(fabs(ad)>0.001)fprintf(stderr, "2 ad=%f scale=%g coeff[1]=(%g, %g) filter_size=%d\n", ad, scale, coeffs[1].re, coeffs[1].im, filter_size);
#endif
}

void test_bessel_filter(void)
{
COMPLEX8 harmonic_coeffs[3];
COMPLEX8 filter_coeffs[7];
COMPLEX8 filter_coeffs_gold[7];
int i;

harmonic_coeffs[0].re=-5.32583;
harmonic_coeffs[0].im=-0.772804;
harmonic_coeffs[1].re=-1.31682;
harmonic_coeffs[1].im=-0.102469;
harmonic_coeffs[2].re=-0.583185;
harmonic_coeffs[2].im=-0.0306731;

make_bessel_filter(filter_coeffs, 7, harmonic_coeffs, 3, 0.125);

filter_coeffs_gold[0].re=-0.0630661037184435;
filter_coeffs_gold[1].re=-0.221543839933367;
filter_coeffs_gold[2].re=-0.16305690111018;
filter_coeffs_gold[3].re=0.583937398666762;
filter_coeffs_gold[4].re=-0.00916336621599076;
filter_coeffs_gold[5].re=-0.208374643680681;
filter_coeffs_gold[6].re=-0.0955543870573582; 

filter_coeffs_gold[0].im=0.0233244359530231;
filter_coeffs_gold[1].im=-0.0402574464426947;
filter_coeffs_gold[2].im=-0.504925020787896;
filter_coeffs_gold[3].im=0.0738975042310994;
filter_coeffs_gold[4].im=-0.48804568992295;
filter_coeffs_gold[5].im=-0.140818825977675;
filter_coeffs_gold[6].im=-0.0141210971953729; 

for(i=0;i<7;i++) {
	fprintf(stderr, "test_bessel_filter %d %g %g\n", i-3, filter_coeffs[i].re, filter_coeffs[i].im);
	}

for(i=0;i<7;i++) {
	if(fabs(filter_coeffs[i].re-filter_coeffs_gold[i].re)>1e-6 || fabs(filter_coeffs[i].im-filter_coeffs_gold[i].im)>1e-6) {
		fprintf(stderr, "*** INTERNAL ERROR: mismatch in filter coeffs (%.12f, %.12f) vs (%.12f, %.12f) i=%d\n", 
			filter_coeffs[i].re, filter_coeffs[i].im, filter_coeffs_gold[i].re, filter_coeffs_gold[i].im, i);
		}
	}
}

void shift_fft7(COMPLEX8Vector *fft_out, COMPLEX8Vector *fft_in, COMPLEX8 *filter)
{
int i, j, k;
double a, b;
int nsamples=fft_in->length;
COMPLEX8 *pf, *pd;

if(fft_in->length!=fft_out->length) {
	fprintf(stderr, "*** INTERNAL ERROR: fft lengths do not match %d vs %d\n", fft_in->length, fft_out->length);
	exit(-1);
	}
	
if(nsamples<10) {
	fprintf(stderr, "*** INTERNAL ERROR: cannot filter very small SFTs (%d)\n", nsamples);
	exit(-1);
	}

for(i=0;i<5;i++) {
	a=0.0;
	b=0.0;
	for(j=0;j<7;j++) {
		k=i-3+j;
		if(k<0)k=nsamples+k;
		a+=filter[6-j].re*fft_in->data[k].re-filter[6-j].im*fft_in->data[k].im;
		b+=filter[6-j].re*fft_in->data[k].im+filter[6-j].im*fft_in->data[k].re;
		}
	
	fft_out->data[i].re=a;
	fft_out->data[i].im=b;
	}

for(;i<nsamples-5;i++) {
	a=0.0;
	b=0.0;
	
	#define ADD {\
		a+=pf->re*pd->re-pf->im*pd->im; \
		b+=pf->re*pd->im+pf->im*pd->re; \
		pf++; \
		pd--; \
		}
		
	
	pf=filter;
	pd=&(fft_in->data[i+3]);
/*	for(j=0;j<7;j++) {
		k=i-3+j;
		a+=filter[6-j].re*fft_in->data[k].re-filter[6-j].im*fft_in->data[k].im;
		b+=filter[6-j].re*fft_in->data[k].im+filter[6-j].im*fft_in->data[k].re;
		}*/

	ADD
	ADD
	ADD
	ADD
	ADD
	ADD
	ADD
	
	fft_out->data[i].re=a;
	fft_out->data[i].im=b;
	}

for(;i<nsamples;i++) {
	a=0.0;
	b=0.0;
	for(j=0;j<7;j++) {
		k=i-3+j;
		if(k>=nsamples)k=k-nsamples;
		a+=filter[6-j].re*fft_in->data[k].re-filter[6-j].im*fft_in->data[k].im;
		b+=filter[6-j].re*fft_in->data[k].im+filter[6-j].im*fft_in->data[k].re;
		}
	
	fft_out->data[i].re=a;
	fft_out->data[i].im=b;
	}

}

void shift_fft7_sse(COMPLEX8Vector *fft_out, COMPLEX8Vector *fft_in, COMPLEX8 *filter)
{
int i, j, k;
float a, b;
int nsamples=fft_in->length;
COMPLEX8 *pf, *pd;
float *filter_re, *filter_im, *tmp_re, *tmp_im;
__m128 filter_re1, filter_re2, filter_im1, filter_im2, tmp_re1, tmp_re2, tmp_im1, tmp_im2, a1, a2, b1, b2;

filter_re=aligned_alloca(8*sizeof(*filter_re));
filter_im=aligned_alloca(8*sizeof(*filter_im));
tmp_re=aligned_alloca(8*sizeof(*tmp_re));
tmp_im=aligned_alloca(8*sizeof(*tmp_im));
// tmp_ret=aligned_alloca(8*sizeof(*tmp_re));
// tmp_imt=aligned_alloca(8*sizeof(*tmp_im));

if(fft_in->length!=fft_out->length) {
	fprintf(stderr, "*** INTERNAL ERROR: fft lengths do not match %d vs %d\n", fft_in->length, fft_out->length);
	exit(-1);
	}
	
if(nsamples<10) {
	fprintf(stderr, "*** INTERNAL ERROR: cannot filter very small SFTs (%d)\n", nsamples);
	exit(-1);
	}

for(i=0;i<7;i++) {
	filter_re[i]=filter[i].re;
	filter_im[i]=filter[i].im;
	}
filter_re[7]=0.0;
filter_im[7]=0.0;
tmp_re[7]=0.0;
tmp_im[7]=0.0;

filter_re1=_mm_load_ps(filter_re);
filter_re2=_mm_load_ps(&(filter_re[4]));
filter_im1=_mm_load_ps(filter_im);
filter_im2=_mm_load_ps(&(filter_im[4]));

for(i=0;i<5;i++) {
	a=0.0;
	b=0.0;
	for(j=0;j<7;j++) {
		k=i-3+j;
		if(k<0)k=nsamples+k;
		a+=filter[6-j].re*fft_in->data[k].re-filter[6-j].im*fft_in->data[k].im;
		b+=filter[6-j].re*fft_in->data[k].im+filter[6-j].im*fft_in->data[k].re;
		}
	
	fft_out->data[i].re=a;
	fft_out->data[i].im=b;
	}
	
for(j=0;j<7;j++) {
	tmp_re[j]=fft_in->data[i+3-j].re;
	tmp_im[j]=fft_in->data[i+3-j].im;
	}

tmp_re1=_mm_load_ps(tmp_re);
tmp_im1=_mm_load_ps(tmp_im);

tmp_re2=_mm_load_ps(&(tmp_re[4]));
tmp_im2=_mm_load_ps(&(tmp_im[4]));

for(;i<nsamples-5;i++) {
	#if 0
	a=0.0;
	b=0.0;
	
	#define ADD {\
		a+=(float)pf->re*(float)pd->re-(float)pf->im*(float)pd->im; \
		b+=(float)pf->re*(float)pd->im+(float)pf->im*(float)pd->re; \
		pf++; \
		pd--; \
		}
		
	
	pf=filter;
	pd=&(fft_in->data[i+3]);

	ADD
	ADD
	ADD
	ADD
	ADD
	ADD
	ADD
	
	fft_out->data[i].re=a;
	fft_out->data[i].im=b;
	#endif
	
	a1=_mm_sub_ps(_mm_mul_ps(tmp_re1, filter_re1), _mm_mul_ps(tmp_im1, filter_im1));
	b1=_mm_add_ps(_mm_mul_ps(tmp_re1, filter_im1), _mm_mul_ps(tmp_im1, filter_re1));

	a2=_mm_sub_ps(_mm_mul_ps(tmp_re2, filter_re2), _mm_mul_ps(tmp_im2, filter_im2));
	b2=_mm_add_ps(_mm_mul_ps(tmp_re2, filter_im2), _mm_mul_ps(tmp_im2, filter_re2));	
	
	/* shuffle data and load next elements */
	
	tmp_re1=_mm_shuffle_ps(tmp_re1, tmp_re1, _MM_SHUFFLE(2,1,0,3));
	a1=_mm_add_ps(a1, a2);
	tmp_re2=_mm_move_ss(_mm_shuffle_ps(tmp_re2, tmp_re2, _MM_SHUFFLE(2,1,0,3)), tmp_re1);
	b1=_mm_add_ps(b1, b2);
	tmp_re1=_mm_move_ss(tmp_re1, _mm_set_ss(fft_in->data[i+4].re));
	
	a1=_mm_hadd_ps(a1, b1);
	tmp_im1=_mm_shuffle_ps(tmp_im1, tmp_im1, _MM_SHUFFLE(2,1,0,3));
	a1=_mm_hadd_ps(a1, a1);
	tmp_im2=_mm_move_ss(_mm_shuffle_ps(tmp_im2, tmp_im2, _MM_SHUFFLE(2,1,0,3)), tmp_im1);
	_mm_store_ss(&fft_out->data[i].re, a1);
	tmp_im1=_mm_move_ss(tmp_im1, _mm_set_ss(fft_in->data[i+4].im));	
	_mm_store_ss(&fft_out->data[i].im, _mm_shuffle_ps(a1, a1, _MM_SHUFFLE(3,2,0,1)));
	#if 0
	if(fabs(a-fft_out->data[i].re)>1e-4*fabs(a) || fabs(b-fft_out->data[i].im)>1e-4*fabs(b)) {
		fprintf(stderr, "(%g, %g) vs (%g, %g)\n", fft_out->data[i].re, fft_out->data[i].im, a, b);
		}
	#endif

	}

for(;i<nsamples;i++) {
	a=0.0;
	b=0.0;
	for(j=0;j<7;j++) {
		k=i-3+j;
		if(k>=nsamples)k=k-nsamples;
		a+=filter[6-j].re*fft_in->data[k].re-filter[6-j].im*fft_in->data[k].im;
		b+=filter[6-j].re*fft_in->data[k].im+filter[6-j].im*fft_in->data[k].re;
		}
	
	fft_out->data[i].re=a;
	fft_out->data[i].im=b;
	}

}

void shift_fft9(COMPLEX8Vector *fft_out, COMPLEX8Vector *fft_in, COMPLEX8 *filter)
{
int i, j, k;
double a, b;
int nsamples=fft_in->length;
COMPLEX8 *pf, *pd;

if(fft_in->length!=fft_out->length) {
	fprintf(stderr, "*** INTERNAL ERROR: fft lengths do not match %d vs %d\n", fft_in->length, fft_out->length);
	exit(-1);
	}
	
if(nsamples<10) {
	fprintf(stderr, "*** INTERNAL ERROR: cannot filter very small SFTs (%d)\n", nsamples);
	exit(-1);
	}

for(i=0;i<5;i++) {
	a=0.0;
	b=0.0;
	for(j=0;j<9;j++) {
		k=i-4+j;
		if(k<0)k=nsamples+k;
		a+=filter[8-j].re*fft_in->data[k].re-filter[8-j].im*fft_in->data[k].im;
		b+=filter[8-j].re*fft_in->data[k].im+filter[8-j].im*fft_in->data[k].re;
		}
	
	fft_out->data[i].re=a;
	fft_out->data[i].im=b;
	}

for(;i<nsamples-5;i++) {
	a=0.0;
	b=0.0;
	
	#define ADD {\
		a+=pf->re*pd->re-pf->im*pd->im; \
		b+=pf->re*pd->im+pf->im*pd->re; \
		pf++; \
		pd--; \
		}
		
	
	pf=filter;
	pd=&(fft_in->data[i+4]);
/*	for(j=0;j<7;j++) {
		k=i-3+j;
		a+=filter[6-j].re*fft_in->data[k].re-filter[6-j].im*fft_in->data[k].im;
		b+=filter[6-j].re*fft_in->data[k].im+filter[6-j].im*fft_in->data[k].re;
		}*/

	ADD
	ADD
	ADD
	ADD
	ADD
	ADD
	ADD
	ADD
	ADD
	
	fft_out->data[i].re=a;
	fft_out->data[i].im=b;
	}

for(;i<nsamples;i++) {
	a=0.0;
	b=0.0;
	for(j=0;j<9;j++) {
		k=i-4+j;
		if(k>=nsamples)k=k-nsamples;
		a+=filter[8-j].re*fft_in->data[k].re-filter[8-j].im*fft_in->data[k].im;
		b+=filter[8-j].re*fft_in->data[k].im+filter[8-j].im*fft_in->data[k].re;
		}
	
	fft_out->data[i].re=a;
	fft_out->data[i].im=b;
	}

}

void shift_fft11(COMPLEX8Vector *fft_out, COMPLEX8Vector *fft_in, COMPLEX8 *filter)
{
int i, j, k;
double a, b;
int nsamples=fft_in->length;
COMPLEX8 *pf, *pd;

if(fft_in->length!=fft_out->length) {
	fprintf(stderr, "*** INTERNAL ERROR: fft lengths do not match %d vs %d\n", fft_in->length, fft_out->length);
	exit(-1);
	}
	
if(nsamples<10) {
	fprintf(stderr, "*** INTERNAL ERROR: cannot filter very small SFTs (%d)\n", nsamples);
	exit(-1);
	}

for(i=0;i<7;i++) {
	a=0.0;
	b=0.0;
	for(j=0;j<11;j++) {
		k=i-5+j;
		if(k<0)k=nsamples+k;
		a+=filter[10-j].re*fft_in->data[k].re-filter[10-j].im*fft_in->data[k].im;
		b+=filter[10-j].re*fft_in->data[k].im+filter[10-j].im*fft_in->data[k].re;
		}
	
	fft_out->data[i].re=a;
	fft_out->data[i].im=b;
	}

for(;i<nsamples-7;i++) {
	a=0.0;
	b=0.0;
	
	#define ADD {\
		a+=pf->re*pd->re-pf->im*pd->im; \
		b+=pf->re*pd->im+pf->im*pd->re; \
		pf++; \
		pd--; \
		}
		
	
	pf=filter;
	pd=&(fft_in->data[i+5]);
/*	for(j=0;j<7;j++) {
		k=i-3+j;
		a+=filter[6-j].re*fft_in->data[k].re-filter[6-j].im*fft_in->data[k].im;
		b+=filter[6-j].re*fft_in->data[k].im+filter[6-j].im*fft_in->data[k].re;
		}*/

	ADD
	ADD
	ADD
	ADD
	ADD
	ADD
	ADD
	ADD
	ADD
	ADD
	ADD
	
	fft_out->data[i].re=a;
	fft_out->data[i].im=b;
	}

for(;i<nsamples;i++) {
	a=0.0;
	b=0.0;
	for(j=0;j<11;j++) {
		k=i-5+j;
		if(k>=nsamples)k=k-nsamples;
		a+=filter[10-j].re*fft_in->data[k].re-filter[10-j].im*fft_in->data[k].im;
		b+=filter[10-j].re*fft_in->data[k].im+filter[10-j].im*fft_in->data[k].re;
		}
	
	fft_out->data[i].re=a;
	fft_out->data[i].im=b;
	}

}

void shift_fft11_sse(COMPLEX8Vector *fft_out, COMPLEX8Vector *fft_in, COMPLEX8 *filter)
{
int i, j, k;
float a, b;
int nsamples=fft_in->length;
COMPLEX8 *pf, *pd;
float *filter_re, *filter_im, *tmp_re, *tmp_im;
__m128 filter128_re[3], filter128_im[3], tmp128_re[3], tmp128_im[3], a1, a2, b1, b2, a3, b3;

filter_re=aligned_alloca(12*sizeof(*filter_re));
filter_im=aligned_alloca(12*sizeof(*filter_im));
tmp_re=aligned_alloca(12*sizeof(*tmp_re));
tmp_im=aligned_alloca(12*sizeof(*tmp_im));
// tmp_ret=aligned_alloca(8*sizeof(*tmp_re));
// tmp_imt=aligned_alloca(8*sizeof(*tmp_im));

if(fft_in->length!=fft_out->length) {
	fprintf(stderr, "*** INTERNAL ERROR: fft lengths do not match %d vs %d\n", fft_in->length, fft_out->length);
	exit(-1);
	}
	
if(nsamples<16) {
	fprintf(stderr, "*** INTERNAL ERROR: cannot filter very small SFTs (%d)\n", nsamples);
	exit(-1);
	}

for(i=0;i<11;i++) {
	filter_re[i]=filter[i].re;
	filter_im[i]=filter[i].im;
	}
filter_re[11]=0.0;
filter_im[11]=0.0;
tmp_re[11]=0.0;
tmp_im[11]=0.0;

for(j=0;j<3;j++) {
	filter128_re[j]=_mm_load_ps(&(filter_re[4*j]));
	filter128_im[j]=_mm_load_ps(&(filter_im[4*j]));
	}

for(i=0;i<7;i++) {
	a=0.0;
	b=0.0;
	for(j=0;j<11;j++) {
		k=i-5+j;
		if(k<0)k=nsamples+k;
		a+=filter[10-j].re*fft_in->data[k].re-filter[10-j].im*fft_in->data[k].im;
		b+=filter[10-j].re*fft_in->data[k].im+filter[10-j].im*fft_in->data[k].re;
		}
	
	fft_out->data[i].re=a;
	fft_out->data[i].im=b;
	}
	
for(j=0;j<11;j++) {
	tmp_re[j]=fft_in->data[i+5-j].re;
	tmp_im[j]=fft_in->data[i+5-j].im;
	}

for(j=0;j<3;j++) {
	tmp128_re[j]=_mm_load_ps(&(tmp_re[4*j]));
	tmp128_im[j]=_mm_load_ps(&(tmp_im[4*j]));
	}

for(;i<nsamples-7;i++) {
	#if 0
	a=0.0;
	b=0.0;
	
	#define ADD {\
		a+=(float)pf->re*(float)pd->re-(float)pf->im*(float)pd->im; \
		b+=(float)pf->re*(float)pd->im+(float)pf->im*(float)pd->re; \
		pf++; \
		pd--; \
		}
		
	
	pf=filter;
	pd=&(fft_in->data[i+5]);

	ADD
	ADD
	ADD
	ADD
	ADD
	ADD
	ADD
	ADD
	ADD
	ADD
	ADD
	
	fft_out->data[i].re=a;
	fft_out->data[i].im=b;
	#endif
	
	a1=_mm_sub_ps(_mm_mul_ps(tmp128_re[0], filter128_re[0]), _mm_mul_ps(tmp128_im[0], filter128_im[0]));
	b1=_mm_add_ps(_mm_mul_ps(tmp128_re[0], filter128_im[0]), _mm_mul_ps(tmp128_im[0], filter128_re[0]));

	a2=_mm_sub_ps(_mm_mul_ps(tmp128_re[1], filter128_re[1]), _mm_mul_ps(tmp128_im[1], filter128_im[1]));
	b2=_mm_add_ps(_mm_mul_ps(tmp128_re[1], filter128_im[1]), _mm_mul_ps(tmp128_im[1], filter128_re[1]));

	a3=_mm_sub_ps(_mm_mul_ps(tmp128_re[2], filter128_re[2]), _mm_mul_ps(tmp128_im[2], filter128_im[2]));
	b3=_mm_add_ps(_mm_mul_ps(tmp128_re[2], filter128_im[2]), _mm_mul_ps(tmp128_im[2], filter128_re[2]));

	/* shuffle data and load next elements */
	
	a1=_mm_add_ps(a1, _mm_add_ps(a2, a3));
	b1=_mm_add_ps(b1, _mm_add_ps(b2, b3));

	tmp128_re[0]=_mm_shuffle_ps(tmp128_re[0], tmp128_re[0], _MM_SHUFFLE(2,1,0,3));
	tmp128_re[1]=_mm_shuffle_ps(tmp128_re[1], tmp128_re[1], _MM_SHUFFLE(2,1,0,3));
	tmp128_re[2]=_mm_shuffle_ps(tmp128_re[2], tmp128_re[2], _MM_SHUFFLE(2,1,0,3));
	
	tmp128_re[2]=_mm_move_ss(tmp128_re[2], tmp128_re[1]);
	tmp128_re[1]=_mm_move_ss(tmp128_re[1], tmp128_re[0]);
	tmp128_re[0]=_mm_move_ss(tmp128_re[0], _mm_set_ss(fft_in->data[i+6].re));
	
	a1=_mm_hadd_ps(a1, b1);
	a1=_mm_hadd_ps(a1, a1);
	_mm_store_ss(&fft_out->data[i].re, a1);
	_mm_store_ss(&fft_out->data[i].im, _mm_shuffle_ps(a1, a1, _MM_SHUFFLE(3,2,0,1)));

	tmp128_im[0]=_mm_shuffle_ps(tmp128_im[0], tmp128_im[0], _MM_SHUFFLE(2,1,0,3));
	tmp128_im[1]=_mm_shuffle_ps(tmp128_im[1], tmp128_im[1], _MM_SHUFFLE(2,1,0,3));
	tmp128_im[2]=_mm_shuffle_ps(tmp128_im[2], tmp128_im[2], _MM_SHUFFLE(2,1,0,3));
	
	tmp128_im[2]=_mm_move_ss(tmp128_im[2], tmp128_im[1]);
	tmp128_im[1]=_mm_move_ss(tmp128_im[1], tmp128_im[0]);
	tmp128_im[0]=_mm_move_ss(tmp128_im[0], _mm_set_ss(fft_in->data[i+6].im));

	#if 0
	if(fabs(a-fft_out->data[i].re)>1e-4*fabs(a) || fabs(b-fft_out->data[i].im)>1e-4*fabs(b)) {
		fprintf(stderr, "(%g, %g) vs (%g, %g)\n", fft_out->data[i].re, fft_out->data[i].im, a, b);
		}
	#endif

	}

for(;i<nsamples;i++) {
	a=0.0;
	b=0.0;
	for(j=0;j<11;j++) {
		k=i-5+j;
		if(k>=nsamples)k=k-nsamples;
		a+=filter[10-j].re*fft_in->data[k].re-filter[10-j].im*fft_in->data[k].im;
		b+=filter[10-j].re*fft_in->data[k].im+filter[10-j].im*fft_in->data[k].re;
		}
	
	fft_out->data[i].re=a;
	fft_out->data[i].im=b;
	}

}

void shift_fft(COMPLEX8Vector *fft_out, COMPLEX8Vector *fft_in, COMPLEX8 *filter, int filter_size)
{
int offset, i, j, k, nsamples=fft_out->length;
COMPLEX8 a,b, *p, *o;

if(!(filter_size & 1)) {
	fprintf(stderr, "*** INTERNAL ERROR: filter size should be odd\n");
	exit(-1);
	}

if(fft_in->length!=fft_out->length) {
	fprintf(stderr, "*** INTERNAL ERROR: fft lengths do not match %d vs %d\n", fft_in->length, fft_out->length);
	exit(-1);
	}
	
if(nsamples<=filter_size) {
	fprintf(stderr, "*** INTERNAL ERROR: cannot filter very small SFTs (length=%d filter_length=%d)\n", nsamples, filter_size);
	exit(-1);
	}

#if 0
memcpy(fft_out->data, fft_in->data, nsamples*sizeof(*(fft_out->data)));
return;
#endif

if(filter_size==7) {
	shift_fft7_sse(fft_out, fft_in, filter);
	return;
	}

if(filter_size==9) {
	shift_fft9(fft_out, fft_in, filter);
	return;
	}

if(filter_size==11) {
	shift_fft11_sse(fft_out, fft_in, filter);
	return;
	}

offset=filter_size>>1;

for(i=0;i<nsamples;i++) {
	fft_out->data[i].re=0.0;
	fft_out->data[i].im=0.0;
	}

for(k=-offset;k<=offset;k++) {
	a=filter[offset-k];
	
	o=fft_out->data;
	for(i=0;i<offset;i++) {
		j=i+k;
		if(j<0)j+=nsamples;
		b=fft_in->data[j];
		
		o->re+=a.re*b.re-a.im*b.im;
		o->im+=a.im*b.re+a.re*b.im;
		o++;
		}

	p=&(fft_in->data[i+k]);
	for(;i+offset<nsamples;i++) {
		
		o->re+=a.re*p->re-a.im*p->im;
		o->im+=a.im*p->re+a.re*p->im;
		p++;
		o++;
		}

	for(;i<nsamples;i++) {
		j=i+k;
		if(j>=nsamples)j-=nsamples;
		b=fft_in->data[j];
		
		o->re+=a.re*b.re-a.im*b.im;
		o->im+=a.im*b.re+a.re*b.im;
		o++;
		}

	}

}

void test_bessel(void)
{
check_my_bessel_Jn();
test_bessel_filter();
}