#include <stdio.h>
#include <string.h>
#include <math.h>

#include "global.h"
#include "fft_stats.h"
#include "cmdline.h"

#define N_IOTA 16
#define N_PSI  32

extern FILE *LOG;

extern struct gengetopt_args_info args_info;

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
	ac->w1_im=-ac->Ax*sin(2*ac->psi);
	ac->w2_re=ac->Ap*sin(2*ac->psi);
	ac->w2_im=ac->Ax*cos(2*ac->psi);
	
	ac->w11=ac->w1_re*ac->w1_re+ac->w1_im*ac->w1_im;
	ac->w12=2.0*(ac->w1_re*ac->w2_re+ac->w1_im*ac->w2_im);
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
st->B_stat.value=-1e25;
st->F_stat.value=-1e25;
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

void log_stats(LOOSE_CONTEXT *ctx, FILE *f, char *tag, FFT_STATS *st, double ul_adjust)
{
#define LOG(a, adj) \
	fprintf(f, "stats: \"%s\" \"%s\" %s %lg %d %lg %.12f %lg %.12f %.12f %lg %lg %lg %lg %lg\n", \
		args_info.label_arg, tag, #a, st->a.value*adj, st->a.fft_bin, st->a.fft_offset, \
		st->a.frequency, st->a.spindown, st->a.ra, st->a.dec, st->a.iota, st->a.psi, st->a.phi, st->a.z.re, st->a.z.im);
	
LOG(snr, 1)
LOG(ul, ul_adjust)
LOG(circ_ul, ul_adjust)
LOG(B_stat, 1)
LOG(F_stat, 1)
fprintf(f, "ratio: \"%s\" %g %g %f\n", tag, st->template_count, st->stat_hit_count, st->stat_hit_count/st->template_count);
}

void update_SNR_stats(LOOSE_CONTEXT *ctx, STAT_INFO *st, COMPLEX8 z1, COMPLEX8 z2, int bin, double fft_offset)
{
int i;
double a, b, x, y, p;
double fpp=ctx->weight_pp, fpc=ctx->weight_pc, fcc=ctx->weight_cc;
ALIGNMENT_COEFFS *ac;

for(i=0;i<acd->free;i++) {
	ac=&(acd->coeffs[i]);
	x=z1.re*ac->w1_re-z1.im*ac->w1_im+z2.re*ac->w2_re-z2.im*ac->w2_im;
	y=z1.re*ac->w1_im+z1.im*ac->w1_re+z2.re*ac->w2_im+z2.im*ac->w2_re;
	p=x*x+y*y;
	a=fpp*ac->w11+fpc*ac->w12+fcc*ac->w22;
	b=p/a;
	
	if(b>st->value) {
		st->value=b;
		st->z.re=x;
		st->z.im=y;
		st->fft_bin=bin;
		st->fft_offset=fft_offset;
		st->alignment_bin=i;
		st->frequency=ctx->frequency+st->fft_offset+(1.0-ctx->te_sc->slope)*(2*st->fft_bin>ctx->nsamples ? st->fft_bin-ctx->nsamples : st->fft_bin)/ctx->timebase;
		st->spindown=ctx->spindown;
		st->ra=ctx->ra;
		st->dec=ctx->dec;
		st->iota=ac->iota;
		st->psi=ac->psi;
		st->phi=atan2(x, y);
		}
		
	}
}

#define UL_CONFIDENCE_LEVEL 1.65

void update_UL_stats(LOOSE_CONTEXT *ctx, STAT_INFO *st, COMPLEX8 z1, COMPLEX8 z2, int bin, double fft_offset)
{
int i;
double a, b, x, y, p;
double fpp=ctx->weight_pp, fpc=ctx->weight_pc, fcc=ctx->weight_cc;
ALIGNMENT_COEFFS *ac;

for(i=0;i<acd->free;i++) {
	ac=&(acd->coeffs[i]);
	x=z1.re*ac->w1_re-z1.im*ac->w1_im+z2.re*ac->w2_re-z2.im*ac->w2_im;
	y=z1.re*ac->w1_im+z1.im*ac->w1_re+z2.re*ac->w2_im+z2.im*ac->w2_re;
	p=x*x+y*y;
	a=fpp*ac->w11+fpc*ac->w12+fcc*ac->w22;
	b=sqrt(p/(a*a)+UL_CONFIDENCE_LEVEL/a);
	
	if(b>st->value) {
		st->value=b;
		st->z.re=x;
		st->z.im=y;
		st->fft_bin=bin;
		st->fft_offset=fft_offset;
		st->alignment_bin=i;
		st->frequency=ctx->frequency+st->fft_offset+(1.0-ctx->te_sc->slope)*(2*st->fft_bin>ctx->nsamples ? st->fft_bin-ctx->nsamples : st->fft_bin)/ctx->timebase;
		st->spindown=ctx->spindown;
		st->ra=ctx->ra;
		st->dec=ctx->dec;
		st->iota=ac->iota;
		st->psi=ac->psi;
		st->phi=atan2(x, y);
		}
	}
}

void update_circ_UL_stats(LOOSE_CONTEXT *ctx, STAT_INFO *st, COMPLEX8 z1, COMPLEX8 z2, int bin, double fft_offset)
{
int i;
double a, b, x, y, p;
double fpp=ctx->weight_pp, fpc=ctx->weight_pc, fcc=ctx->weight_cc;
ALIGNMENT_COEFFS *ac;

for(i=0;i<2;i++) {
	ac=&(acd->coeffs[i]);
	x=z1.re*ac->w1_re-z1.im*ac->w1_im+z2.re*ac->w2_re-z2.im*ac->w2_im;
	y=z1.re*ac->w1_im+z1.im*ac->w1_re+z2.re*ac->w2_im+z2.im*ac->w2_re;
	p=x*x+y*y;
	a=fpp*ac->w11+fpc*ac->w12+fcc*ac->w22;
	b=sqrt(p/(a*a)+UL_CONFIDENCE_LEVEL/a);
	
	if(b>st->value) {
		st->value=b;
		st->z.re=x;
		st->z.im=y;
		st->fft_bin=bin;
		st->fft_offset=fft_offset;
		st->alignment_bin=i;
		st->frequency=(ctx->frequency+(1.0-ctx->te_sc->slope)*(2*st->fft_bin>ctx->nsamples ? st->fft_bin-ctx->nsamples : st->fft_bin)/ctx->timebase)*(1+st->fft_offset);
		st->spindown=ctx->spindown;
		st->ra=ctx->ra;
		st->dec=ctx->dec;
		st->iota=ac->iota;
		st->psi=ac->psi;
		st->phi=atan2(x, y);
		}
	}
}

void update_B_stats(LOOSE_CONTEXT *ctx, STAT_INFO *st, COMPLEX8 z1, COMPLEX8 z2, int bin, double fft_offset)
{
int i;
double a, b, x, y, p, v, w, w_total;
double fpp=ctx->weight_pp, fpc=ctx->weight_pc, fcc=ctx->weight_cc;
ALIGNMENT_COEFFS *ac;

v=-1e25;
w_total=0;
for(i=0;i<acd->free;i++) {
	ac=&(acd->coeffs[i]);
	x=z1.re*ac->w1_re-z1.im*ac->w1_im+z2.re*ac->w2_re-z2.im*ac->w2_im;
	y=z1.re*ac->w1_im+z1.im*ac->w1_re+z2.re*ac->w2_im+z2.im*ac->w2_re;
	//p=x*x+y*y;
	//a=fpp*ac->w11+fpc*ac->w12+fcc*ac->w22;
	a=(fcc*ac->w11*x*x-fpc*ac->w12*x*y+fpp*ac->w22*y*y)/(fpp*fcc-0.25*fpc*fpc);
	w=1.0;
	if(w>0) {
		w_total+=w;
		b=a-log(w);
		if(v>b) v+=log1p(exp(b-v));
			else
			v=b+log1p(exp(v-b));
		}
	}

v-=log(w_total);

if(v>st->value) {
	st->value=v;
	st->z.re=-1;
	st->z.im=-1;
	st->fft_bin=bin;
	st->fft_offset=fft_offset;
	st->alignment_bin=i;
	st->frequency=ctx->frequency+st->fft_offset+(1.0-ctx->te_sc->slope)*(2*st->fft_bin>ctx->nsamples ? st->fft_bin-ctx->nsamples : st->fft_bin)/ctx->timebase;
	st->spindown=ctx->spindown;
	st->ra=ctx->ra;
	st->dec=ctx->dec;
	st->iota=-1;
	st->psi=-1;
	st->phi=-1;
	}
}

void update_F_stats_int(LOOSE_CONTEXT *ctx, STAT_INFO *st, COMPLEX8 z1, COMPLEX8 z2, int bin, double fft_offset)
{
int i;
double a, b, x, y, p, v, w, w_total;
double fpp=ctx->weight_pp, fpc=ctx->weight_pc, fcc=ctx->weight_cc;
ALIGNMENT_COEFFS *ac;

v=-1e25;
w_total=0;
for(i=0;i<acd->free;i++) {
	ac=&(acd->coeffs[i]);
	x=z1.re*ac->w1_re-z1.im*ac->w1_im+z2.re*ac->w2_re-z2.im*ac->w2_im;
	y=z1.re*ac->w1_im+z1.im*ac->w1_re+z2.re*ac->w2_im+z2.im*ac->w2_re;
	//p=x*x+y*y;
	//a=fpp*ac->w11+fpc*ac->w12+fcc*ac->w22;
	a=(fcc*ac->w11*x*x-fpc*ac->w12*x*y+fpp*ac->w22*y*y)/(fpp*fcc-0.25*fpc*fpc);
	b=(1.0-ac->Ax*ac->Ax);
	w=b*b*b;
	if(w>0) {
		w_total+=w;
		b=a-log(w);
		if(v>b) v+=log1p(exp(b-v));
			else
			v=b+log1p(exp(v-b));
		}
	//fprintf(stderr, "%g %g %g\n", v, b, w);
	//fprintf(stderr, "%g %g %g\n", adj, sqrt(p)*adj, a*adj);
	}

v-=log(w_total);

if(v>st->value) {
	st->value=v;
	st->z.re=-1;
	st->z.im=-1;
	st->fft_bin=bin;
	st->fft_offset=fft_offset;
	st->alignment_bin=i;
	st->frequency=ctx->frequency+st->fft_offset+(1.0-ctx->te_sc->slope)*(2*st->fft_bin>ctx->nsamples ? st->fft_bin-ctx->nsamples : st->fft_bin)/ctx->timebase;
	st->spindown=ctx->spindown;
	st->ra=ctx->ra;
	st->dec=ctx->dec;
	st->iota=-1;
	st->psi=-1;
	st->phi=-1;
	}
}

void update_F_stats(LOOSE_CONTEXT *ctx, STAT_INFO *st, COMPLEX8 z1, COMPLEX8 z2, int bin, double fft_offset)
{
int i;
double a, b, x, y, p;
double fpp=ctx->weight_pp, fpc=ctx->weight_pc, fcc=ctx->weight_cc;
ALIGNMENT_COEFFS *ac;

for(i=0;i<acd->free;i++) {
	ac=&(acd->coeffs[i]);
	x=z1.re*ac->w1_re-z1.im*ac->w1_im+z2.re*ac->w2_re-z2.im*ac->w2_im;
	y=z1.re*ac->w1_im+z1.im*ac->w1_re+z2.re*ac->w2_im+z2.im*ac->w2_re;
	b=(fcc*ac->w11*x*x-fpc*ac->w12*x*y+fpp*ac->w22*y*y)/(fpp*fcc-0.25*fpc*fpc);
	
	if(b>st->value) {
		st->value=b;
		st->z.re=x;
		st->z.im=y;
		st->fft_bin=bin;
		st->fft_offset=fft_offset;
		st->alignment_bin=i;
		st->frequency=ctx->frequency+st->fft_offset+(1.0-ctx->te_sc->slope)*(2*st->fft_bin>ctx->nsamples ? st->fft_bin-ctx->nsamples : st->fft_bin)/ctx->timebase;
		st->spindown=ctx->spindown;
		st->ra=ctx->ra;
		st->dec=ctx->dec;
		st->iota=ac->iota;
		st->psi=ac->psi;
		st->phi=atan2(x, y);
		}
		
	}
}

void compute_fft_stats(LOOSE_CONTEXT *ctx, FFT_STATS *stats, COMPLEX8Vector *fft1, COMPLEX8Vector *fft2, double fft_offset)
{
float M[4];
float V;
int idx;
int i;
int nsamples=fft1->length;

#if 0
V=0;
for(i=0;i<fft1->length;i++) {
	V+=fft1->data[i].re+fft1->data[i].im+fft2->data[i].re+fft2->data[i].im;
	}
stats->ul.value=V;
return;
#endif

for(i=0;i<4;i++)M[i]=0.0;

for(i=0;i<nsamples;i++) {
	/* crudely skip indices outside Nyquist */
	if(abs(i-(nsamples>>1))<(nsamples>>2))continue;
	
	V=fft1->data[i].re*fft1->data[i].re+fft1->data[i].im*fft1->data[i].im+fft2->data[i].re*fft2->data[i].re+fft2->data[i].im*fft2->data[i].im;
	idx=(fft1->data[i].re*fft2->data[i].re+fft1->data[i].im*fft2->data[i].im>=0)*2+(fft1->data[i].im*fft2->data[i].re-fft1->data[i].re*fft2->data[i].im>=0);
	
	if(V>M[idx])M[idx]=V;
	}

for(i=0;i<nsamples;i++) {
	/* crudely skip indices outside Nyquist */
	if(abs(i-(nsamples>>1))<(nsamples>>2))continue;

	V=fft1->data[i].re*fft1->data[i].re+fft1->data[i].im*fft1->data[i].im+fft2->data[i].re*fft2->data[i].re+fft2->data[i].im*fft2->data[i].im;
	idx=(fft1->data[i].re*fft2->data[i].re+fft1->data[i].im*fft2->data[i].im>=0)*2+(fft1->data[i].im*fft2->data[i].re-fft1->data[i].re*fft2->data[i].im>=0);
	
	if(V<0.5*M[idx])continue;
	
	update_SNR_stats(ctx, &(stats->snr), fft1->data[i], fft2->data[i], (i*2>nsamples ? i-nsamples : i), fft_offset);
	update_UL_stats(ctx, &(stats->ul), fft1->data[i], fft2->data[i], (i*2>nsamples ? i-nsamples : i), fft_offset);
	update_circ_UL_stats(ctx, &(stats->circ_ul), fft1->data[i], fft2->data[i], (i*2>nsamples ? i-nsamples : i), fft_offset);
	update_F_stats(ctx, &(stats->F_stat), fft1->data[i], fft2->data[i], (i*2>nsamples ? i-nsamples : i), fft_offset);
	update_B_stats(ctx, &(stats->B_stat), fft1->data[i], fft2->data[i], (i*2>nsamples ? i-nsamples : i), fft_offset);
	
	stats->stat_hit_count++;
	}
stats->template_count+=fft1->length;
}

void init_fft_stats(void)
{
allocate_alignment_coeffs();	
}