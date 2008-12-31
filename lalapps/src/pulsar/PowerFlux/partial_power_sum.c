SUFFIX(PARTIAL_POWER_SUM) *SUFFIX(allocate_partial_power_sum)(int pps_bins)
{
SUFFIX(PARTIAL_POWER_SUM) *r;
int i;

r=do_alloc(1, sizeof(*r));

r->type=sizeof(REAL);
r->nbins=pps_bins;
r->offset=(pps_bins-useful_bins)>>1;

r->weight_pppp=do_alloc(pps_bins, sizeof(*(r->weight_pppp)));
r->weight_pppc=do_alloc(pps_bins, sizeof(*(r->weight_pppc)));
r->weight_ppcc=do_alloc(pps_bins, sizeof(*(r->weight_ppcc)));
r->weight_pccc=do_alloc(pps_bins, sizeof(*(r->weight_pccc)));
r->weight_cccc=do_alloc(pps_bins, sizeof(*(r->weight_cccc)));

r->power_pp=do_alloc(pps_bins, sizeof(*(r->power_pp)));
r->power_pc=do_alloc(pps_bins, sizeof(*(r->power_pc)));
r->power_cc=do_alloc(pps_bins, sizeof(*(r->power_cc)));

for(i=0;i<pps_bins;i++) {
	r->weight_pppp[i]=NAN;
	r->weight_pppc[i]=NAN;
	r->weight_ppcc[i]=NAN;
	r->weight_pccc[i]=NAN;
	r->weight_cccc[i]=NAN;

	r->power_pp[i]=NAN;
	r->power_pc[i]=NAN;
	r->power_cc[i]=NAN;
	}

r->c_weight_pppp=NAN;
r->c_weight_pppc=NAN;
r->c_weight_ppcc=NAN;
r->c_weight_pccc=NAN;
r->c_weight_cccc=NAN;

r->weight_arrays_non_zero=1;
r->collapsed_weight_arrays=0;
return r;
}

void SUFFIX(zero_partial_power_sum)(SUFFIX(PARTIAL_POWER_SUM) *pps)
{
int i;
int pps_bins=pps->nbins;
if(pps->type!=sizeof(REAL)) {
	fprintf(stderr, "*** INTERNAL ERROR: %s power sum real type mismatch type=%d sizeof(REAL)=%ld\n",
		__FUNCTION__,
		pps->type, sizeof(REAL));
	exit(-1);
	}

for(i=0;i<pps_bins;i++) {
	pps->weight_pppp[i]=0;
	pps->weight_pppc[i]=0;
	pps->weight_ppcc[i]=0;
	pps->weight_pccc[i]=0;
	pps->weight_cccc[i]=0;

	pps->power_pp[i]=0;
	pps->power_pc[i]=0;
	pps->power_cc[i]=0;
	}

pps->c_weight_pppp=0;
pps->c_weight_pppc=0;
pps->c_weight_ppcc=0;
pps->c_weight_pccc=0;
pps->c_weight_cccc=0;

pps->weight_arrays_non_zero=0;
pps->collapsed_weight_arrays=0;
}

void SUFFIX(accumulate_partial_power_sum)(SUFFIX(PARTIAL_POWER_SUM) *accum, SUFFIX(PARTIAL_POWER_SUM) *partial)
{
int i;
REAL *a1, *a2, *a3, *a4, *a5, *p1, *p2, *p3, *p4, *p5;
int pps_bins=accum->nbins;
int shift;
if(accum->type!=sizeof(REAL)) {
	fprintf(stderr, "*** INTERNAL ERROR: %s power sum real type mismatch type=%d sizeof(REAL)=%ld\n",
		__FUNCTION__,
		accum->type, sizeof(REAL));
	exit(-1);
	}
if(partial->type!=sizeof(REAL)) {
	fprintf(stderr, "*** INTERNAL ERROR: %s power sum real type mismatch type=%d sizeof(REAL)=%ld\n",
		__FUNCTION__,
		partial->type, sizeof(REAL));
	exit(-1);
	}
shift=partial->offset-accum->offset;
if( (shift<0) || (shift+pps_bins>partial->nbins)) {
	fprintf(stderr, "*** INTERNAL ERROR: %s power sums do not overlap shift=%d pps_bins=%d partial->offset=%d accum->offset=%d partial->nbins=%d\n",
		__FUNCTION__,
		shift, pps_bins, partial->offset, accum->offset, partial->nbins);
	exit(-1);
	}

a1=accum->power_pp;
a2=accum->power_pc;
a3=accum->power_cc;

p1=&(partial->power_pp[shift]);
p2=&(partial->power_pc[shift]);
p3=&(partial->power_cc[shift]);

for(i=0;i<pps_bins;i++) {
	(*a1)+=*p1;
	(*a2)+=*p2;
	(*a3)+=*p3;

	a1++;
	p1++;
	a2++;
	p2++;
	a3++;
	p3++;
	}

if(partial->weight_arrays_non_zero) {
	a1=accum->weight_pppp;
	a2=accum->weight_pppc;
	a3=accum->weight_ppcc;
	a4=accum->weight_pccc;
	a5=accum->weight_cccc;

	p1=&(partial->weight_pppp[shift]);
	p2=&(partial->weight_pppc[shift]);
	p3=&(partial->weight_ppcc[shift]);
	p4=&(partial->weight_pccc[shift]);
	p5=&(partial->weight_cccc[shift]);

	for(i=0;i<pps_bins;i++) {
		(*a1)+=*p1;
		(*a2)+=*p2;
		(*a3)+=*p3;
		(*a4)+=*p4;
		(*a5)+=*p5;
	
		a1++;
		p1++;
		a2++;
		p2++;
		a3++;
		p3++;
		a4++;
		p4++;
		a5++;
		p5++;
		}
	accum->weight_arrays_non_zero=1;
	}

accum->c_weight_pppp+=partial->c_weight_pppp;
accum->c_weight_pppc+=partial->c_weight_pppc;
accum->c_weight_ppcc+=partial->c_weight_ppcc;
accum->c_weight_pccc+=partial->c_weight_pccc;
accum->c_weight_cccc+=partial->c_weight_cccc;
}

void SUFFIX(dump_partial_power_sum)(FILE *out, SUFFIX(PARTIAL_POWER_SUM) *pps)
{
int i;
int pps_bins=pps->nbins;
if(pps->type!=sizeof(REAL)) {
	fprintf(stderr, "*** INTERNAL ERROR: %s power sum real type mismatch type=%d sizeof(REAL)=%ld\n",
		__FUNCTION__,
		pps->type, sizeof(REAL));
	exit(-1);
	}

fprintf(out, "%g %g %g %g %g", 
	pps->c_weight_pppp,
	pps->c_weight_pppc,
	pps->c_weight_ppcc,
	pps->c_weight_pccc,
	pps->c_weight_cccc );

for(i=0;i<pps_bins;i++) fprintf(out, " %g", pps->weight_pppp[i]);
for(i=0;i<pps_bins;i++) fprintf(out, " %g", pps->weight_pppc[i]);
for(i=0;i<pps_bins;i++) fprintf(out, " %g", pps->weight_ppcc[i]);
for(i=0;i<pps_bins;i++) fprintf(out, " %g", pps->weight_pccc[i]);
for(i=0;i<pps_bins;i++) fprintf(out, " %g", pps->weight_cccc[i]);

for(i=0;i<pps_bins;i++) fprintf(out, " %g", pps->power_pp[i]);
for(i=0;i<pps_bins;i++) fprintf(out, " %g", pps->power_pc[i]);
for(i=0;i<pps_bins;i++) fprintf(out, " %g", pps->power_cc[i]);
}

void SUFFIX(free_partial_power_sum)(SUFFIX(PARTIAL_POWER_SUM) *pps)
{
if(pps->type!=sizeof(REAL)) {
	fprintf(stderr, "*** INTERNAL ERROR: %s power sum real type mismatch type=%d sizeof(REAL)=%ld\n",
		__FUNCTION__,
		pps->type, sizeof(REAL));
	exit(-1);
	}

free(pps->weight_pppp);
free(pps->weight_pppc);
free(pps->weight_ppcc);
free(pps->weight_pccc);
free(pps->weight_cccc);

free(pps->power_pp);
free(pps->power_pc);
free(pps->power_cc);

memset(pps, 0, sizeof(*pps));
free(pps);
}
