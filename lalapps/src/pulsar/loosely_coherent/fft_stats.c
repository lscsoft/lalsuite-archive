#include <stdio.h>
#include <string.h>
#include <math.h>

#include "global.h"
#include "fft_stats.h"

#define N_IOTA 16
#define N_PSI  32

extern FILE *LOG;

static ALIGNMENT_DATA *acd=NULL;

static void allocate_alignment_coeffs(void)
{
int i, j, k;
ALIGNMENT_COEFFS *ac;

acd=do_alloc(1, sizeof(*acd));

acd->size=N_IOTA*N_PSI+2;
acd->free=0;

acd->coeffs=do_alloc(acd->size, sizeof(*acd->coeffs));

ac=&(acd->coeffs[0]);

ac->iota=0;
ac->psi=0;
acd->free++;

ac=&(acd->coeffs[1]);

ac->iota=M_PI;
ac->psi=0;
acd->free++;

k=2;
for(i=0;i<N_IOTA;i++)
	for(j=0;j<N_PSI;j++) {
		ac=&(acd->coeffs[k]);

		ac->iota=(M_PI*(i+1))/(N_IOTA+2);
		ac->psi=(0.5*M_PI*j)/N_PSI;

		acd->free++;
		k++;
		}
		
for(k=0;k<acd->free;k++) {
	ac=&(acd->coeffs[k]);

	ac->Ax=cos(ac->iota);
	ac->Ap=(1+ac->Ax*ac->Ax)*0.5;
	
	ac->w1_re=ac->Ap*cos(2*ac->psi);
	ac->w1_im=ac->Ax*sin(2*ac->psi);
	ac->w2_re=ac->Ap*sin(2*ac->psi);
	ac->w2_im=-ac->Ax*cos(2*ac->psi);
	
	ac->w11=ac->w1_re*ac->w1_re+ac->w1_im*ac->w1_im;
	ac->w12=ac->w1_re*ac->w2_re+ac->w1_im*ac->w2_im;
	ac->w22=ac->w2_re*ac->w2_re+ac->w2_im*ac->w2_im;
	}

fprintf(LOG, "alignment_coeffs: iota psi Ax Ap w1_re w1_im w2_re w2_im\n");
for(k=0;k<acd->free;k++) {
	ac=&(acd->coeffs[k]);
	fprintf(LOG, "%d %f %f %f %f %f %f %f %f\n", k, ac->iota, ac->psi, ac->Ax, ac->Ap, ac->w1_re, ac->w1_im, ac->w2_re, ac->w2_im);
	}

}

void init_stats(FFT_STATS *st)
{
memset(st, 0, sizeof(*st));
}

void update_stats(FFT_STATS *st_accum, FFT_STATS *st)
{
#define UPDATE_STAT(a)	{\
	if(st->a.value>st_accum->a.value)st_accum->a=st->a; \
	} 
	
UPDATE_STAT(snr);
UPDATE_STAT(ul);
UPDATE_STAT(circ_ul);
}

void log_stats(FILE *f, char *tag, FFT_STATS *st, double ul_adjust)
{
#define LOG(a, adj) \
	fprintf(f, "stats: \"%s\" %s %lg %d %lg %lg %lg %lg %lg %lg\n", tag, #a, st->a.value*adj, st->a.fft_bin, st->a.fft_offset, st->a.iota, st->a.psi, st->a.phi, st->a.z.re, st->a.z.im);
	
LOG(snr, 1)
LOG(ul, ul_adjust)
LOG(circ_ul, ul_adjust)
fprintf(f, "ratio: \"%s\" %g %g %f\n", tag, st->template_count, st->stat_hit_count, st->stat_hit_count/st->template_count);
}

void update_SNR_stats(STAT_INFO *st, COMPLEX8 z1, COMPLEX8 z2, double fpp, double fpc, double fcc, int bin, double fft_offset)
{
int i;
double a, b, x, y, p;
ALIGNMENT_COEFFS *ac;

for(i=0;i<acd->free;i++) {
	ac=&(acd->coeffs[i]);
	x=z1.re*ac->w1_re-z1.im*ac->w1_im+z2.re*ac->w2_re-z2.im*ac->w2_im;
	y=z1.re*ac->w1_im+z1.im*ac->w1_re+z2.re*ac->w2_im+z2.im*ac->w2_re;
	p=x*x+y*y;
	a=fpp*ac->w11+2*fpc*ac->w12+fcc*ac->w22;
	b=p/a;
	
	if(b>st->value) {
		st->value=b;
		st->z.re=x;
		st->z.im=y;
		st->fft_bin=bin;
		st->fft_offset=fft_offset;
		st->alignment_bin=i;
		st->iota=ac->iota;
		st->psi=ac->psi;
		st->phi=atan2(x, y);
		}
		
	}
}

#define UL_CONFIDENCE_LEVEL 1.65

void update_UL_stats(STAT_INFO *st, COMPLEX8 z1, COMPLEX8 z2, double fpp, double fpc, double fcc, int bin, double fft_offset)
{
int i;
double a, b, x, y, p;
ALIGNMENT_COEFFS *ac;

for(i=0;i<acd->free;i++) {
	ac=&(acd->coeffs[i]);
	x=z1.re*ac->w1_re-z1.im*ac->w1_im+z2.re*ac->w2_re-z2.im*ac->w2_im;
	y=z1.re*ac->w1_im+z1.im*ac->w1_re+z2.re*ac->w2_im+z2.im*ac->w2_re;
	p=x*x+y*y;
	a=fpp*ac->w11+2*fpc*ac->w12+fcc*ac->w22;
	b=sqrt(p/(a*a)+UL_CONFIDENCE_LEVEL/a);
	
	if(b>st->value) {
		st->value=b;
		st->z.re=x;
		st->z.im=y;
		st->fft_bin=bin;
		st->fft_offset=fft_offset;
		st->alignment_bin=i;
		st->iota=ac->iota;
		st->psi=ac->psi;
		st->phi=atan2(x, y);
		}
	}
}

void update_circ_UL_stats(STAT_INFO *st, COMPLEX8 z1, COMPLEX8 z2, double fpp, double fpc, double fcc, int bin, double fft_offset)
{
int i;
double a, b, x, y, p;
ALIGNMENT_COEFFS *ac;

for(i=0;i<2;i++) {
	ac=&(acd->coeffs[i]);
	x=z1.re*ac->w1_re-z1.im*ac->w1_im+z2.re*ac->w2_re-z2.im*ac->w2_im;
	y=z1.re*ac->w1_im+z1.im*ac->w1_re+z2.re*ac->w2_im+z2.im*ac->w2_re;
	p=x*x+y*y;
	a=fpp*ac->w11+2*fpc*ac->w12+fcc*ac->w22;
	b=sqrt(p/(a*a)+UL_CONFIDENCE_LEVEL/a);
	
	if(b>st->value) {
		st->value=b;
		st->z.re=x;
		st->z.im=y;
		st->fft_bin=bin;
		st->fft_offset=fft_offset;
		st->alignment_bin=i;
		st->iota=ac->iota;
		st->psi=ac->psi;
		st->phi=atan2(x, y);
		}
	}
}

void compute_fft_stats(FFT_STATS *stats, COMPLEX8Vector *fft1, COMPLEX8Vector *fft2, double fpp, double fpc, double fcc, double fft_offset)
{
float M[4];
float V;
int idx;
int i;

#if 0
V=0;
for(i=0;i<fft1->length;i++) {
	V+=fft1->data[i].re+fft1->data[i].im+fft2->data[i].re+fft2->data[i].im;
	}
stats->ul.value=V;
return;
#endif

for(i=0;i<4;i++)M[i]=0.0;

for(i=0;i<fft1->length;i++) {
	V=fft1->data[i].re*fft1->data[i].re+fft1->data[i].im*fft1->data[i].im+fft2->data[i].re*fft2->data[i].re+fft2->data[i].im*fft2->data[i].im;
	idx=(fft1->data[i].re*fft2->data[i].re+fft1->data[i].im*fft2->data[i].im>=0)*2+(fft1->data[i].im*fft2->data[i].re-fft1->data[i].re*fft2->data[i].im>=0);
	
	if(V>M[idx])M[idx]=V;
	}

for(i=0;i<fft1->length;i++) {
	V=fft1->data[i].re*fft1->data[i].re+fft1->data[i].im*fft1->data[i].im+fft2->data[i].re*fft2->data[i].re+fft2->data[i].im*fft2->data[i].im;
	idx=(fft1->data[i].re*fft2->data[i].re+fft1->data[i].im*fft2->data[i].im>=0)*2+(fft1->data[i].im*fft2->data[i].re-fft1->data[i].re*fft2->data[i].im>=0);
	
	if(V<0.5*M[idx])continue;
	
	update_SNR_stats(&(stats->snr), fft1->data[i], fft2->data[i], fpp, fpc, fcc, i, fft_offset);
	update_UL_stats(&(stats->ul), fft1->data[i], fft2->data[i], fpp, fpc, fcc, i, fft_offset);
	update_circ_UL_stats(&(stats->circ_ul), fft1->data[i], fft2->data[i], fpp, fpc, fcc, i, fft_offset);
	
	stats->stat_hit_count++;
	}
stats->template_count+=fft1->length;
}

void init_fft_stats(void)
{
allocate_alignment_coeffs();	
}