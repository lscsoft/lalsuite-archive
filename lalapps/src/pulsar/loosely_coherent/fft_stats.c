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
st->min_noise_ratio=1e24;
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
#define LOG(a, adj, shift) \
	fprintf(f, "stats: \"%s\" \"%s\" %s %d %lg %d %lg %.12f %lg %.12f %.12f %lg %lg %lg %lg %lg\n", \
		args_info.label_arg, tag, #a, ctx->patch_id, st->a.value*adj+shift, st->a.fft_bin, st->a.fft_offset, \
		st->a.frequency, st->a.spindown, st->a.ra, st->a.dec, st->a.iota, st->a.psi, st->a.phi, st->a.z.re, st->a.z.im);
	
LOG(snr, 1, 0)
LOG(ul, ul_adjust, 0)
LOG(circ_ul, ul_adjust, 0)
LOG(B_stat, 1/log(10.0), log10(ul_adjust))
fprintf(stderr, "B_stat: %g log10(ul_adjust)=%g\n", st->B_stat.value/log(10.0), log10(ul_adjust));
fprintf(stderr, "SNR: %g log10(UL)=%g iota_UL=%f psi_UL=%f\n", st->snr.value, log10(st->ul.value*ul_adjust), st->ul.iota, st->ul.psi);
LOG(F_stat, 1, 0)
fprintf(f, "ratio: \"%s\" \"%s\" %d %g %g %f %f %f %f %f %f %f %f %f %f %g\n", args_info.label_arg, tag, ctx->patch_id, st->template_count, st->stat_hit_count, st->stat_hit_count/st->template_count, st->min_noise_ratio, st->max_noise_ratio, 
	ctx->ratio_SNR, ctx->ratio_UL, ctx->ratio_UL_circ, ctx->ratio_B_stat, ctx->ratio_F_stat, ctx->max_ratio, ctx->noise_adj[0], ctx->noise_adj[1]);
fprintf(f, "weight: \"%s\" \"%s\" %d %g %g %g\n", args_info.label_arg, tag, ctx->patch_id, ctx->weight_pp, ctx->weight_pc, ctx->weight_cc);
}

void compute_SNR_variance1(LOOSE_CONTEXT *ctx)
{
int i,j,k;
double fpp=ctx->weight_pp, fpc=ctx->weight_pc, fcc=ctx->weight_cc;
ALIGNMENT_COEFFS *ac;
double x1_re, x1_im, x2_re, x2_im;
double min_norm, max_norm;
double a, b, c, d, s, aa, bb, cc;
int i_max=100, j_max=100;

double x_min[4], x_max[4];

max_norm=0;
min_norm=1e20;

for(i=0;i<=i_max;i++) {
	a=(1.0*i)/i_max;
	b=sqrt(1.0-a*a);
	for(j=0;j<j_max;j++) {
		c=cos(0.25*2*M_PI*j/j_max);
		s=sin(0.25*2*M_PI*j/j_max);
		
		x1_re=a*c;
		x1_im=-a*s;
		x2_re=b*c;
		x2_im=b*s;
		
		d=0;
		for(k=0;k<acd->free;k++) {
			ac=&(acd->coeffs[k]);
			aa=x1_re*ac->w1_re+x1_im*ac->w1_im+x2_re*ac->w2_re-x2_im*ac->w2_im;
			bb=x1_re*ac->w1_im-x1_im*ac->w1_re+x2_re*ac->w2_im+x2_im*ac->w1_re;
			cc=(aa*aa+bb*bb)/(fpp*ac->w11+fpc*ac->w12+fcc*ac->w22);
			if(cc>d)d=cc;
			}
		fprintf(stderr, "x1=(%f, %f) x2=(%f, %f) d=%g\n", x1_re, x1_im, x2_re, x2_im, d);
		if(d>max_norm) {
			max_norm=d;
			x_max[0]=x1_re;
			x_max[1]=x1_im;
			x_max[2]=x2_re;
			x_max[3]=x2_im;
			}
		if(d<min_norm){
			min_norm=d;
			x_min[0]=x1_re;
			x_min[1]=x1_im;
			x_min[2]=x2_re;
			x_min[3]=x2_im;
			}
		}
	}
fprintf(stderr, "min=%g (%f,%f), (%f, %f)  max=%g (%f,%f), (%f, %f)  ratio=%f\n", min_norm, x_min[0], x_min[1], x_min[2], x_min[3], max_norm, x_max[0], x_max[1], x_max[2], x_max[3],  max_norm/min_norm);
}

void update_SNR_stats(LOOSE_CONTEXT *ctx, STAT_INFO *st, COMPLEX8 z1, COMPLEX8 z2, int bin, double fft_offset)
{
int i;
double a, b, x, y, p;
double fpp, fpc, fcc;
ALIGNMENT_COEFFS *ac;
int nsamples=ctx->nsamples;

/* compute noise adjustment */
a=ctx->noise_adj[0]+(ctx->noise_adj[1]*(fabs(bin)-(nsamples>>3)))/(nsamples>>4);
fpp=a*ctx->weight_pp;
fpc=a*ctx->weight_pc;
fcc=a*ctx->weight_cc;

for(i=0;i<acd->free;i++) {
	ac=&(acd->coeffs[i]);
	x=z1.re*ac->w1_re-z1.im*ac->w1_im+z2.re*ac->w2_re-z2.im*ac->w2_im;
	y=z1.re*ac->w1_im+z1.im*ac->w1_re+z2.re*ac->w2_im+z2.im*ac->w2_re;
	p=x*x+y*y;
	a=fpp*ac->w11+fpc*ac->w12+fcc*ac->w22;
	b=p/a;
	
	if(b>st->value) {
		st->value=b;
		st->z.re=ctx->var_offset[0];
		st->z.im=ctx->var_offset[1];
		st->fft_bin=bin;
		st->fft_offset=fft_offset;
		st->alignment_bin=i;
		st->frequency=(1.0+st->fft_offset/ctx->frequency)*(ctx->frequency+(1.0-ctx->te_sc->slope)*(2*st->fft_bin>ctx->nsamples ? st->fft_bin-ctx->nsamples : st->fft_bin)/ctx->timebase);
		st->spindown=ctx->spindown;
		st->ra=ctx->ra+(ctx->sb_ra[0]-ctx->ra)*ctx->var_offset[0]+(ctx->sb_ra[1]-ctx->ra)*ctx->var_offset[1];
		st->dec=ctx->dec+(ctx->sb_dec[0]-ctx->dec)*ctx->var_offset[0]+(ctx->sb_dec[1]-ctx->dec)*ctx->var_offset[1];
		st->iota=ac->iota;
		st->psi=ac->psi;
		st->phi=atan2(x, y);
		}
		
	}
}


//#define UL_CONFIDENCE_LEVEL 1.65
/* Use 3.0 which is good for both Gaussian and exponential statistic, and all chi2 */
#define UL_CONFIDENCE_LEVEL 3.0

void update_UL_stats(LOOSE_CONTEXT *ctx, STAT_INFO *st, COMPLEX8 z1, COMPLEX8 z2, int bin, double fft_offset)
{
int i;
double a, b, x, y, p, stv;
double fpp, fpc, fcc;
ALIGNMENT_COEFFS *ac;
int nsamples=ctx->nsamples;

/* compute noise adjustment */
a=ctx->noise_adj[0]+(ctx->noise_adj[1]*(fabs(bin)-(nsamples>>3)))/(nsamples>>4);
fpp=a*ctx->weight_pp;
fpc=a*ctx->weight_pc;
fcc=a*ctx->weight_cc;

stv=st->value*st->value;
for(i=0;i<acd->free;i++) {
	ac=&(acd->coeffs[i]);
	x=z1.re*ac->w1_re-z1.im*ac->w1_im+z2.re*ac->w2_re-z2.im*ac->w2_im;
	y=z1.re*ac->w1_im+z1.im*ac->w1_re+z2.re*ac->w2_im+z2.im*ac->w2_re;
	p=x*x+y*y;
	a=1.0/(fpp*ac->w11+fpc*ac->w12+fcc*ac->w22);
	b=p*(a*a)+2*sqrt(p*(a*a*a))+(UL_CONFIDENCE_LEVEL-1)*a;
	
	if(b>stv) {
		stv=b;
		st->value=sqrt(b);
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

/* this function is for ratio estimation using worst case limit of abs(z1)+abs(z2)->inf */
void update_UL_stats_raw(LOOSE_CONTEXT *ctx, STAT_INFO *st, COMPLEX8 z1, COMPLEX8 z2, int bin, double fft_offset)
{
int i;
double a, b, x, y, p, stv;
double fpp, fpc, fcc;
ALIGNMENT_COEFFS *ac;
int nsamples=ctx->nsamples;

/* compute noise adjustment */
a=ctx->noise_adj[0]+(ctx->noise_adj[1]*(fabs(bin)-(nsamples>>3)))/(nsamples>>4);
fpp=a*ctx->weight_pp;
fpc=a*ctx->weight_pc;
fcc=a*ctx->weight_cc;

stv=st->value*st->value;
for(i=0;i<acd->free;i++) {
	ac=&(acd->coeffs[i]);
	x=z1.re*ac->w1_re-z1.im*ac->w1_im+z2.re*ac->w2_re-z2.im*ac->w2_im;
	y=z1.re*ac->w1_im+z1.im*ac->w1_re+z2.re*ac->w2_im+z2.im*ac->w2_re;
	p=x*x+y*y;
	a=1.0/(fpp*ac->w11+fpc*ac->w12+fcc*ac->w22);
	b=p*(a*a);
	
	if(b>stv) {
		stv=b;
		st->value=sqrt(b);
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
double fpp, fpc, fcc;
ALIGNMENT_COEFFS *ac;
int nsamples=ctx->nsamples;

/* compute noise adjustment */
a=ctx->noise_adj[0]+(ctx->noise_adj[1]*(fabs(bin)-(nsamples>>3)))/(nsamples>>4);
fpp=a*ctx->weight_pp;
fpc=a*ctx->weight_pc;
fcc=a*ctx->weight_cc;

for(i=0;i<2;i++) {
	ac=&(acd->coeffs[i]);
	x=z1.re*ac->w1_re-z1.im*ac->w1_im+z2.re*ac->w2_re-z2.im*ac->w2_im;
	y=z1.re*ac->w1_im+z1.im*ac->w1_re+z2.re*ac->w2_im+z2.im*ac->w2_re;
	p=x*x+y*y;
	a=1.0/(fpp*ac->w11+fpc*ac->w12+fcc*ac->w22);
	b=sqrt(p*(a*a)+2*sqrt(p*(a*a*a))+(UL_CONFIDENCE_LEVEL-1)*a);
	//b=sqrt(p/(a*a)+UL_CONFIDENCE_LEVEL/a);
	
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

/* this function is for ratio estimation using worst case limit of abs(z1)+abs(z2)->inf */
void update_circ_UL_stats_raw(LOOSE_CONTEXT *ctx, STAT_INFO *st, COMPLEX8 z1, COMPLEX8 z2, int bin, double fft_offset)
{
int i;
double a, b, x, y, p;
double fpp, fpc, fcc;
ALIGNMENT_COEFFS *ac;
int nsamples=ctx->nsamples;

/* compute noise adjustment */
a=ctx->noise_adj[0]+(ctx->noise_adj[1]*(fabs(bin)-(nsamples>>3)))/(nsamples>>4);
fpp=a*ctx->weight_pp;
fpc=a*ctx->weight_pc;
fcc=a*ctx->weight_cc;

for(i=0;i<2;i++) {
	ac=&(acd->coeffs[i]);
	x=z1.re*ac->w1_re-z1.im*ac->w1_im+z2.re*ac->w2_re-z2.im*ac->w2_im;
	y=z1.re*ac->w1_im+z1.im*ac->w1_re+z2.re*ac->w2_im+z2.im*ac->w2_re;
	p=x*x+y*y;
	a=1.0/(fpp*ac->w11+fpc*ac->w12+fcc*ac->w22);
	b=sqrt(p*(a*a));
	
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

void update_B_stats1(LOOSE_CONTEXT *ctx, STAT_INFO *st, COMPLEX8 z1, COMPLEX8 z2, int bin, double fft_offset)
{
int i;
double a, b, x, y, p, v, w, w_total;
double fpp, fpc, fcc;
ALIGNMENT_COEFFS *ac;
int nsamples=ctx->nsamples;

/* compute noise adjustment */
a=ctx->noise_adj[0]+(ctx->noise_adj[1]*(fabs(bin)-(nsamples>>3)))/(nsamples>>4);
fpp=a*ctx->weight_pp;
fpc=a*ctx->weight_pc;
fcc=a*ctx->weight_cc;

v=-1e25;
w_total=0;
for(i=0;i<acd->free;i++) {
	ac=&(acd->coeffs[i]);
	x=z1.re*ac->w1_re-z1.im*ac->w1_im+z2.re*ac->w2_re-z2.im*ac->w2_im;
	y=z1.re*ac->w1_im+z1.im*ac->w1_re+z2.re*ac->w2_im+z2.im*ac->w2_re;
	//p=x*x+y*y;
	//a=fpp*ac->w11+fpc*ac->w12+fcc*ac->w22;
	a=(fcc*ac->w11*x*x-2*fpc*ac->w12*x*y+fpp*ac->w22*y*y)/(fpp*fcc-fpc*fpc);
	w=1.0;
	if(w>0) {
		w_total+=w;
		b=a-log(w);
		if(v>b+2) v+=0;
			else
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
	st->frequency=(1.0+st->fft_offset/ctx->frequency)*(ctx->frequency+(1.0-ctx->te_sc->slope)*(2*st->fft_bin>ctx->nsamples ? st->fft_bin-ctx->nsamples : st->fft_bin)/ctx->timebase);
	st->spindown=ctx->spindown;
	st->ra=ctx->ra;
	st->dec=ctx->dec;
	st->iota=-1;
	st->psi=-1;
	st->phi=-1;
	}
}


void update_B_stats2(LOOSE_CONTEXT *ctx, STAT_INFO *st, COMPLEX8 z1, COMPLEX8 z2, int bin, double fft_offset)
{
int i;
double a, b, x, y, p, v, w, w_total;
double fpp, fpc, fcc;
ALIGNMENT_COEFFS *ac;
int nsamples=ctx->nsamples;

/* compute noise adjustment */
a=ctx->noise_adj[0]+(ctx->noise_adj[1]*(fabs(bin)-(nsamples>>3)))/(nsamples>>4);
fpp=a*ctx->weight_pp;
fpc=a*ctx->weight_pc;
fcc=a*ctx->weight_cc;

v=-1e55;
w_total=0;
for(i=0;i<acd->free;i++) {
	ac=&(acd->coeffs[i]);
	x=z1.re*ac->w1_re-z1.im*ac->w1_im+z2.re*ac->w2_re-z2.im*ac->w2_im;
	y=z1.re*ac->w1_im+z1.im*ac->w1_re+z2.re*ac->w2_im+z2.im*ac->w2_re;
	//p=x*x+y*y;
	//a=fpp*ac->w11+fpc*ac->w12+fcc*ac->w22;
	b=fpp*fcc-fpc*fpc;
	a=0.5*(fcc*ac->w11*x*x-2*fpc*ac->w12*x*y+fpp*ac->w22*y*y)/(b);
	a=a-0.5*log(b)-log1p(exp(-a*(1.6+0.07*a*a-0.0005*a*a*a*a)))-0.5*log(2*M_PI);
	w=1.0;
	if(w>0) {
		w_total+=w;
		b=a-log(w);
		if(v>b+2) v+=0;
			else
		if(v>b) v=v+log1p(exp(b-v));
			else
			v=b+log1p(exp(v-b));
		}
	}

//v-=log(w_total);

if(v>st->value) {
	st->value=v;
	st->z.re=-1;
	st->z.im=-1;
	st->fft_bin=bin;
	st->fft_offset=fft_offset;
	st->alignment_bin=i;
	st->frequency=(1.0+st->fft_offset/ctx->frequency)*(ctx->frequency+(1.0-ctx->te_sc->slope)*(2*st->fft_bin>ctx->nsamples ? st->fft_bin-ctx->nsamples : st->fft_bin)/ctx->timebase);
	st->spindown=ctx->spindown;
	st->ra=ctx->ra;
	st->dec=ctx->dec;
	st->iota=-1;
	st->psi=-1;
	st->phi=-1;
	}
}

void update_B_stats3(LOOSE_CONTEXT *ctx, STAT_INFO *st, COMPLEX8 z1, COMPLEX8 z2, int bin, double fft_offset)
{
int i;
double a, b, x, y, p, v, w, w_total;
double fpp, fpc, fcc;
ALIGNMENT_COEFFS *ac;
int nsamples=ctx->nsamples;

/* compute noise adjustment */
a=ctx->noise_adj[0]+(ctx->noise_adj[1]*(fabs(bin)-(nsamples>>3)))/(nsamples>>4);
fpp=a*ctx->weight_pp;
fpc=a*ctx->weight_pc;
fcc=a*ctx->weight_cc;

v=-1e55;
w_total=0;
for(i=0;i<acd->free;i++) {
	ac=&(acd->coeffs[i]);
	x=z1.re*ac->w1_re-z1.im*ac->w1_im+z2.re*ac->w2_re-z2.im*ac->w2_im;
	y=z1.re*ac->w1_im+z1.im*ac->w1_re+z2.re*ac->w2_im+z2.im*ac->w2_re;
	p=x*x+y*y;
	b=fpp*ac->w11+fpc*ac->w12+fcc*ac->w22;
	//b=fpp*fcc-fpc*fpc;
	//a=0.5*(fcc*ac->w11*x*x-2*fpc*ac->w12*x*y+fpp*ac->w22*y*y)/(b);
	a=p/b;
	p=0.5*a-0.5*log(b)-log1p(exp(-a*(1.6+0.07*a*a-0.0005*a*a*a*a)))-0.5*log(2*M_PI);
	w=1.0;
	if(w>0) {
		w_total+=w;
		b=p-log(w);
		if(v>b+2) v+=0;
			else
		if(v>b) v=v+log1p(exp(b-v));
			else
			v=b+log1p(exp(v-b));
		}
	}

//v-=log(w_total);

if(v>st->value) {
	st->value=v;
	st->z.re=-1;
	st->z.im=-1;
	st->fft_bin=bin;
	st->fft_offset=fft_offset;
	st->alignment_bin=i;
	st->frequency=(1.0+st->fft_offset/ctx->frequency)*(ctx->frequency+(1.0-ctx->te_sc->slope)*(2*st->fft_bin>ctx->nsamples ? st->fft_bin-ctx->nsamples : st->fft_bin)/ctx->timebase);
	st->spindown=ctx->spindown;
	st->ra=ctx->ra;
	st->dec=ctx->dec;
	st->iota=-1;
	st->psi=-1;
	st->phi=-1;
	}
}

double theta_appr(double x)
{
double a0, a2, a4, b0, b2, b4, y;
a0=7.71990148904874;
a2=19.0337266315871;
a4=5.20172247607555;
b0=7.7201854234519;
b2=21.1533518190664;
b4=4.28188537828521;

y=x*x;
return(exp(sqrt(0.25+y))/sqrt(sqrt(4*M_PI*M_PI*y+16*exp(2)))*(a0+a2*y+a4*y*y+y*y*y)/(b0+b2*y+b4*y*y+y*y*y));
}

/* fast, moderate precision log that makes use of large dynamic range of double x */
static inline double fast_log(double x)
{
//fprintf(stderr, "%g\n", x);
return(log(x));
}

double log_theta_appr(double x)
{
double a0, a2, a4, b0, b2, b4, y, y2;
a0=7.71990148904874;
a2=19.0337266315871;
a4=5.20172247607555;
b0=7.7201854234519;
b2=21.1533518190664;
b4=4.28188537828521;
if(fabs(x)<1.0) {
	y=x*x;
	return((sqrt(0.25+y))-0.25*fast_log((4*M_PI*M_PI*y+16*expf(2)))+fast_log((a0+a2*y+a4*y*y+y*y*y)/(b0+b2*y+b4*y*y+y*y*y)));
	} else {
	y=x*x;
	y2=1.0/y;
	return((sqrt(0.25+y))-0.25*fast_log((4*M_PI*M_PI*y+16*expf(2)))+fast_log((a0*y2*y2*y2+a2*y2*y2+a4*y2+1.0)/(b0*y2*y2*y2+b2*y2*y2+b4*y2+1.0)));
	}
}

void update_B_stats(LOOSE_CONTEXT *ctx, STAT_INFO *st, COMPLEX8 z1, COMPLEX8 z2, int bin, double fft_offset)
{
int i;
double a, b, x, y, p, v, w, w_total;
double fpp, fpc, fcc;
ALIGNMENT_COEFFS *ac;
int nsamples=ctx->nsamples;

/* compute noise adjustment */
a=ctx->noise_adj[0]+(ctx->noise_adj[1]*(fabs(bin)-(nsamples>>3)))/(nsamples>>4);
fpp=a*ctx->weight_pp;
fpc=a*ctx->weight_pc;
fcc=a*ctx->weight_cc;

v=-1e55;
w_total=0;
for(i=0;i<acd->free;i++) {
	ac=&(acd->coeffs[i]);
	x=z1.re*ac->w1_re-z1.im*ac->w1_im+z2.re*ac->w2_re-z2.im*ac->w2_im;
	y=z1.re*ac->w1_im+z1.im*ac->w1_re+z2.re*ac->w2_im+z2.im*ac->w2_re;
	p=x*x+y*y;
	b=fpp*ac->w11+fpc*ac->w12+fcc*ac->w22;
	a=p/b;
	
/*	if(a>5)fprintf(stderr, "a=%g b=%g p=%g", a, b, p);*/
	p=log_theta_appr(0.5*a)-0.5*log(b);
/*	if(a>5)fprintf(stderr, " log(B-stat)=%g\n", p);*/
	w=1.0;
	if(w>0) {
		w_total+=w;
		b=p-log(w);
		if(v>b+2) v+=0;
			else
		if(v>b) v=v+log1p(exp(b-v));
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
	st->frequency=(1.0+st->fft_offset/ctx->frequency)*(ctx->frequency+(1.0-ctx->te_sc->slope)*(2*st->fft_bin>ctx->nsamples ? st->fft_bin-ctx->nsamples : st->fft_bin)/ctx->timebase);
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
double a, b, x, y, p;
double fpp, fpc, fcc;
ALIGNMENT_COEFFS *ac;
int nsamples=ctx->nsamples;
COMPLEX16 w1, w2;

/* compute noise adjustment */
a=ctx->noise_adj[0]+(ctx->noise_adj[1]*(fabs(bin)-(nsamples>>3)))/(nsamples>>4);
fpp=a*ctx->weight_pp;
fpc=a*ctx->weight_pc;
fcc=a*ctx->weight_cc;

b=(fcc*(z1.re*z1.re+z1.im*z1.im)-2*fpc*(z1.re*z2.re+z1.im*z2.im)+fpp*(z2.re*z2.re+z2.im*z2.im))/(fpp*fcc-fpc*fpc);

if(b>st->value) {
	w1.re=a*(fcc*z1.re-fpc*z2.re);
	w1.im=a*(fcc*z1.im-fpc*z2.im);
	
	w2.re=a*(-fpc*z1.re+fpp*z2.re);
	w2.im=a*(-fpc*z1.im+fpp*z2.im);

	st->value=b;
	st->z.re=w1.re;
	st->z.im=w1.im;
	st->fft_bin=bin;
	st->fft_offset=fft_offset;
	st->alignment_bin=-1;
	st->frequency=ctx->frequency+st->fft_offset+(1.0-ctx->te_sc->slope)*(2*st->fft_bin>ctx->nsamples ? st->fft_bin-ctx->nsamples : st->fft_bin)/ctx->timebase;
	st->spindown=ctx->spindown;
	st->ra=ctx->ra;
	st->dec=ctx->dec;
	
	a=1.0/(fpp*fcc-fpc*fpc);
	
	st->iota=-1;
	st->psi=w2.re;
	st->phi=w2.im;
	}
}

double compute_stats_func_ratio(LOOSE_CONTEXT *ctx, void (*stats_func)(LOOSE_CONTEXT *ctx, STAT_INFO *st, COMPLEX8 z1, COMPLEX8 z2, int bin, double fft_offset))
{
int i,j,k;
double fpp=ctx->weight_pp, fpc=ctx->weight_pc, fcc=ctx->weight_cc;
ALIGNMENT_COEFFS *ac;
double min_norm, max_norm;
double a, b, c, d, s, aa, bb, cc;
int i_max=100, j_max=100;
STAT_INFO st;
COMPLEX8 z1, z2;
double norm1, norm2;

double x_min[4], x_max[4];

max_norm=-1e25;
min_norm=1e25;

//norm=sqrt(fpp+fcc);

/* make sure we use both z1 and z2 */
if(fpp<0.2*fcc)fpp=0.2*fcc;
if(fcc<0.2*fpp)fcc=0.2*fpp;

norm1=sqrt(fpp);
norm2=sqrt(fcc);

for(i=0;i<=i_max;i++) {
	a=(1.0*i)/i_max;
	b=sqrt(1.0-a*a);
	for(j=0;j<j_max;j++) {
		c=cos(0.25*2*M_PI*j/j_max);
		s=sin(0.25*2*M_PI*j/j_max);
		
// 		z1.re=norm*a*c;
// 		z1.im=-norm*a*s;
// 		z2.re=norm*b*c;
// 		z2.im=norm*b*s;

		z1.re=norm1*a*c;
		z1.im=-norm1*a*s;
		z2.re=norm2*b*c;
		z2.im=norm2*b*s;

		d=0;
		st.value=0;
		stats_func(ctx, &st, z1, z2, 0, 0);

/*		fprintf(stderr, "x1=(%f, %f) x2=(%f, %f) d=%g\n", z1.re, z1.im, z2.re, z2.im, st.value);*/
		if(st.value>max_norm) {
			max_norm=st.value;
			x_max[0]=z1.re;
			x_max[1]=z1.im;
			x_max[2]=z2.re;
			x_max[3]=z2.im;
			}
		if(st.value<min_norm){
			min_norm=st.value;
			x_min[0]=z1.re;
			x_min[1]=z1.im;
			x_min[2]=z2.re;
			x_min[3]=z2.im;
			}
		}
	}
/*fprintf(stderr, "min=%g (%f,%f), (%f, %f)  max=%g (%f,%f), (%f, %f)  ratio=%f\n", min_norm, x_min[0], x_min[1], x_min[2], x_min[3], max_norm, x_max[0], x_max[1], x_max[2], x_max[3],  max_norm/min_norm);*/
return (max_norm/min_norm);
}

void compute_stats_variance(LOOSE_CONTEXT *ctx)
{
ctx->noise_adj[0]=1.0;
ctx->noise_adj[1]=0.0;
ctx->ratio_SNR=compute_stats_func_ratio(ctx, update_SNR_stats);
ctx->ratio_UL=compute_stats_func_ratio(ctx, update_UL_stats_raw);
ctx->ratio_UL=ctx->ratio_UL*ctx->ratio_UL;
ctx->ratio_UL_circ=compute_stats_func_ratio(ctx, update_circ_UL_stats_raw);
ctx->ratio_UL_circ=ctx->ratio_UL_circ*ctx->ratio_UL_circ;
ctx->ratio_B_stat=compute_stats_func_ratio(ctx, update_B_stats);
ctx->ratio_F_stat=compute_stats_func_ratio(ctx, update_F_stats);

ctx->max_ratio=ctx->ratio_SNR;
if(ctx->ratio_UL>ctx->max_ratio)ctx->max_ratio=ctx->ratio_UL;
if(ctx->ratio_UL_circ>ctx->max_ratio)ctx->max_ratio=ctx->ratio_UL_circ;
//if(ctx->ratio_B_stat>ctx->max_ratio)ctx->max_ratio=ctx->ratio_B_stat;
if(ctx->ratio_F_stat>ctx->max_ratio)ctx->max_ratio=ctx->ratio_F_stat;
fprintf(stderr, "ratio: SNR=%f UL=%f UL_circ=%f B_stat=%f F_stat=%f fpp=%g fpc=%g fcc=%g\n", 
		ctx->ratio_SNR,
		ctx->ratio_UL,
		ctx->ratio_UL_circ,
		ctx->ratio_B_stat,
		ctx->ratio_F_stat,
		ctx->weight_pp,
		ctx->weight_pc,
		ctx->weight_cc
       );
}

void compute_fft_stats(LOOSE_CONTEXT *ctx, FFT_STATS *stats, COMPLEX8Vector *fft1, COMPLEX8Vector *fft2, double fft_offset)
{
float M[4];
float V;
int idx;
int i;
int nsamples=fft1->length;
int template_count=0, noise_count;
double sum, sum2, noise_level;
double fpp=ctx->weight_pp, fpc=ctx->weight_pc, fcc=ctx->weight_cc;
FILE *f;
int skip_band;
int center=nsamples>>1;

#if 0
V=0;
for(i=0;i<fft1->length;i++) {
	V+=fft1->data[i].re+fft1->data[i].im+fft2->data[i].re+fft2->data[i].im;
	}
stats->ul.value=V;
return;
#endif

/* WARNING !!!! this function must match compute ratio */

/* make sure we use both z1 and z2 */
if(fpp<0.2*fcc)fpp=0.2*fcc;
if(fcc<0.2*fpp)fcc=0.2*fpp;

for(i=0;i<4;i++)M[i]=0.0;

sum=0;
sum2=0;
template_count=0;
noise_level=(fcc*ctx->weight_pp+fpp*ctx->weight_cc)*20;
noise_count=0;
// f=fopen("sum.txt", "w");

/* Precompute maximum norm fcc*|X1|+fpp*|X2| as well as small linear deviation in noise level with index */

//skip_band=(0.5*(nsamples/ctx->timebase)-args_info.focus_f0_delta_arg)*ctx->timebase;
skip_band=floor((0.5+0.5*(1.0-args_info.focus_f0_useful_fraction_arg))*(nsamples>>1));

//fprintf(stderr, "nsamples=%d skip_band=%d fraction=%f\n", nsamples, skip_band, skip_band*1.0/(nsamples>>2));

for(i=0;i<nsamples;i++) {
	/* crudely skip indices outside Nyquist */
	//if(abs(i-(nsamples>>1))<(nsamples>>2))continue;
	if(abs(i-center)<skip_band)continue;
	
	V=fcc*(fft1->data[i].re*fft1->data[i].re+fft1->data[i].im*fft1->data[i].im)+fpp*(fft2->data[i].re*fft2->data[i].re+fft2->data[i].im*fft2->data[i].im);
	idx=0;
	//idx=(fft1->data[i].re*fft2->data[i].re+fft1->data[i].im*fft2->data[i].im>=0)*2+(fft1->data[i].im*fft2->data[i].re-fft1->data[i].re*fft2->data[i].im>=0);
	
	if(V>M[idx])M[idx]=V;
	template_count++;

	if(V<noise_level) {
		noise_count++;
		sum+=V;
		sum2+=(V*((i>(nsamples>>1) ? nsamples-i : i)-(((nsamples>>1)-skip_band)>>1)))/(((nsamples>>1)-skip_band)>>2);
/*		fprintf(f, "%d %d %g %g %g %g %g %g %g\n", i, ((i>(nsamples>>1) ? nsamples-i : i)-(skip_band>>1)), V/(ctx->weight_pp+ctx->weight_cc), (sum/noise_count)/(ctx->weight_pp+ctx->weight_cc), (sum2/noise_count)/(ctx->weight_pp+ctx->weight_cc), fft1->data[i].re, fft1->data[i].im, fft2->data[i].re, fft2->data[i].im);*/
		}
	//fprintf(stderr, "V=%g V_norm=%g sum=%g sum_norm=%g\n", V, V/(ctx->weight_pp+ctx->weight_cc), sum/template_count, (sum/template_count)/(ctx->weight_pp+ctx->weight_cc));
	}
	
ctx->noise_adj[0]=(sum/noise_count)/(fcc*ctx->weight_pp+fpp*ctx->weight_cc);
ctx->noise_adj[1]=(sum2/noise_count)/(fcc*ctx->weight_pp+fpp*ctx->weight_cc);
/*fprintf(stderr, "%g %g skip_band=%d nsamples>>1=%d noise_count=%d\n", ctx->noise_adj[0], ctx->noise_adj[1], skip_band, nsamples>>1, noise_count);*/
// fclose(f);
// exit(-1);

for(i=0;i<nsamples;i++) {
	/* crudely skip indices outside Nyquist */
	//if(abs(i-(nsamples>>1))<(nsamples>>2))continue;
	if(abs(i-center)<skip_band)continue;

	V=fcc*(fft1->data[i].re*fft1->data[i].re+fft1->data[i].im*fft1->data[i].im)+fpp*(fft2->data[i].re*fft2->data[i].re+fft2->data[i].im*fft2->data[i].im);
/*	V=fft1->data[i].re*fft1->data[i].re+fft1->data[i].im*fft1->data[i].im+fft2->data[i].re*fft2->data[i].re+fft2->data[i].im*fft2->data[i].im;*/
	idx=0;
	//idx=(fft1->data[i].re*fft2->data[i].re+fft1->data[i].im*fft2->data[i].im>=0)*2+(fft1->data[i].im*fft2->data[i].re-fft1->data[i].re*fft2->data[i].im>=0);
	
	if(V*ctx->max_ratio<M[idx])continue;

	if(V*ctx->ratio_SNR>M[idx])
		update_SNR_stats(ctx, &(stats->snr), fft1->data[i], fft2->data[i], (i*2>nsamples ? i-nsamples : i), fft_offset);
	if(V*ctx->ratio_UL>M[idx])
		update_UL_stats(ctx, &(stats->ul), fft1->data[i], fft2->data[i], (i*2>nsamples ? i-nsamples : i), fft_offset);
	if(V*ctx->ratio_UL_circ>M[idx])
		update_circ_UL_stats(ctx, &(stats->circ_ul), fft1->data[i], fft2->data[i], (i*2>nsamples ? i-nsamples : i), fft_offset);
	if(V*ctx->ratio_F_stat>M[idx])
		update_F_stats(ctx, &(stats->F_stat), fft1->data[i], fft2->data[i], (i*2>nsamples ? i-nsamples : i), fft_offset);
	
	if(V*ctx->ratio_F_stat>M[idx])
		update_B_stats(ctx, &(stats->B_stat), fft1->data[i], fft2->data[i], (i*2>nsamples ? i-nsamples : i), fft_offset);
	
	stats->stat_hit_count++;
	}
stats->template_count+=template_count;
sum/=noise_count*(ctx->weight_pp+ctx->weight_cc);
if(sum>stats->max_noise_ratio)stats->max_noise_ratio=sum;
if(sum<stats->min_noise_ratio)stats->min_noise_ratio=sum;
//fprintf(stderr, "%g %g\n", stats->stat_hit_count*100.0/stats->template_count, ctx->max_ratio);
//fprintf(stderr, "sum=%g template_count=%d\n", sum, template_count);
}

void init_fft_stats(void)
{
allocate_alignment_coeffs();	
}