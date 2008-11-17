SUFFIX(PARTIAL_POWER_SUM) *SUFFIX(allocate_partial_power_sum)(void)
{
SUFFIX(PARTIAL_POWER_SUM) *r;
int i;

r=do_alloc(1, sizeof(*r));

r->type=sizeof(REAL);

r->weight_pppp=do_alloc(useful_bins, sizeof(*(r->weight_pppp)));
r->weight_pppc=do_alloc(useful_bins, sizeof(*(r->weight_pppc)));
r->weight_ppcc=do_alloc(useful_bins, sizeof(*(r->weight_ppcc)));
r->weight_pccc=do_alloc(useful_bins, sizeof(*(r->weight_pccc)));
r->weight_cccc=do_alloc(useful_bins, sizeof(*(r->weight_cccc)));

r->power_pp=do_alloc(useful_bins, sizeof(*(r->power_pp)));
r->power_pc=do_alloc(useful_bins, sizeof(*(r->power_pc)));
r->power_cc=do_alloc(useful_bins, sizeof(*(r->power_cc)));

for(i=0;i<useful_bins;i++) {
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
if(pps->type!=sizeof(REAL)) {
	fprintf(stderr, "*** INTERNAL ERROR: %s power sum real type mismatch type=%d sizeof(REAL)=%ld\n",
		__FUNCTION__,
		pps->type, sizeof(REAL));
	exit(-1);
	}

for(i=0;i<useful_bins;i++) {
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

for(i=0;i<useful_bins;i++) {
	accum->power_pp[i]+=partial->power_pp[i];
	accum->power_pc[i]+=partial->power_pc[i];
	accum->power_cc[i]+=partial->power_cc[i];
	}

if(partial->weight_arrays_non_zero) {
	for(i=0;i<useful_bins;i++) {
		accum->weight_pppp[i]+=partial->weight_pppp[i];
		accum->weight_pppc[i]+=partial->weight_pppc[i];
		accum->weight_ppcc[i]+=partial->weight_ppcc[i];
		accum->weight_pccc[i]+=partial->weight_pccc[i];
		accum->weight_cccc[i]+=partial->weight_cccc[i];
		}
	accum->weight_arrays_non_zero=1;
	}

accum->c_weight_pppp+=partial->c_weight_pppp;
accum->c_weight_pppc+=partial->c_weight_pppc;
accum->c_weight_ppcc+=partial->c_weight_ppcc;
accum->c_weight_pccc+=partial->c_weight_pccc;
accum->c_weight_cccc+=partial->c_weight_cccc;
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
