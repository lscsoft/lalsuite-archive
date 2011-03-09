SUFFIX(PARTIAL_POWER_SUM) *SUFFIX(allocate_partial_power_sum)(int pps_bins, int cross_terms_present)
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
	
if(cross_terms_present) {
	r->power_im_pc=do_alloc(pps_bins, sizeof(*(r->power_pc)));
	for(i=0;i<pps_bins;i++) {
		r->power_im_pc[i]=NAN;
		}
	} else 
	r->power_im_pc=NULL;

r->c_weight_pppp=NAN;
r->c_weight_pppc=NAN;
r->c_weight_ppcc=NAN;
r->c_weight_pccc=NAN;
r->c_weight_cccc=NAN;
r->c_weight_im_ppcc=NAN;

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
	
if(pps->power_im_pc!=NULL) {
	for(i=0;i<pps_bins;i++) {
		pps->power_im_pc[i]=0;
		}
	}

pps->c_weight_pppp=0;
pps->c_weight_pppc=0;
pps->c_weight_ppcc=0;
pps->c_weight_pccc=0;
pps->c_weight_cccc=0;
pps->c_weight_im_ppcc=0;

pps->weight_arrays_non_zero=0;
pps->collapsed_weight_arrays=0;
}

void SUFFIX(randomize_partial_power_sum)(SUFFIX(PARTIAL_POWER_SUM) *pps)
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
	pps->weight_pppp[i]=(rand()*100.0)/RAND_MAX;
	pps->weight_pppc[i]=(rand()*100.0)/RAND_MAX;
	pps->weight_ppcc[i]=(rand()*100.0)/RAND_MAX;
	pps->weight_pccc[i]=(rand()*100.0)/RAND_MAX;
	pps->weight_cccc[i]=(rand()*100.0)/RAND_MAX;

	pps->power_pp[i]=(rand()*100.0)/RAND_MAX;
	pps->power_pc[i]=(rand()*100.0)/RAND_MAX;
	pps->power_cc[i]=(rand()*100.0)/RAND_MAX;
	}

if(pps->power_im_pc!=NULL) {
	for(i=0;i<pps_bins;i++) {
		pps->power_im_pc[i]=(rand()*100.0)/RAND_MAX;
		}
	}

pps->c_weight_pppp=(rand()*100.0)/RAND_MAX;
pps->c_weight_pppc=(rand()*100.0)/RAND_MAX;
pps->c_weight_ppcc=(rand()*100.0)/RAND_MAX;
pps->c_weight_pccc=(rand()*100.0)/RAND_MAX;
pps->c_weight_cccc=(rand()*100.0)/RAND_MAX;
pps->c_weight_im_ppcc= (pps->power_im_pc!=NULL) ? (rand()*100.0)/RAND_MAX : 0.0;

pps->weight_arrays_non_zero=1;
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

if(partial->power_im_pc!=NULL) {
	if(accum->power_im_pc==NULL) {
		fprintf(stderr, "*** INTERNAL ERROR: attempt to accumulate partial power sum with cross terms into accumulator without\n");
		exit(-1);
		}
	a1=accum->power_im_pc;
	p1=&(partial->power_im_pc[shift]);
	for(i=0;i<pps_bins;i++) {
		(*a1)+=*p1;
		
		a1++;
		p1++;
		}
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
accum->c_weight_im_ppcc+=partial->c_weight_im_ppcc;
}

int SUFFIX(compare_partial_power_sums)(char *prefix, SUFFIX(PARTIAL_POWER_SUM) *ref, SUFFIX(PARTIAL_POWER_SUM) *test, REAL rel_tolerance, REAL rel_abs_tolerance)
{
int i;
int pps_bins=ref->nbins;

if(ref->type!=sizeof(REAL)) {
	fprintf(stderr, "%sreference power sum real type mismatch type=%d sizeof(REAL)=%ld\n",
		prefix,
		ref->type, sizeof(REAL));
	return -1;
	}
if(test->type!=sizeof(REAL)) {
	fprintf(stderr, "%stest power sum real type mismatch type=%d sizeof(REAL)=%ld\n",
		prefix,
		test->type, sizeof(REAL));
	return -2;
	}

#define TEST(field, format, tolerance) \
	if( (ref->field!=test->field) && !(fabs(test->field-ref->field)<tolerance*(fabs(ref->field)+fabs(test->field)))) { \
		fprintf(stderr, "%s" #field " fields do not match ref=" format " test=" format " test-ref=" format "\n", \
			prefix, \
			ref->field, test->field, test->field-ref->field); \
		return -3; \
		}

TEST(offset, "%d", -1)
TEST(nbins, "%d", -1)
TEST(weight_arrays_non_zero, "%d", -1)
TEST(collapsed_weight_arrays, "%d", -1)

TEST(c_weight_pppp, "%g", rel_tolerance)
TEST(c_weight_pppc, "%g", rel_tolerance)
TEST(c_weight_ppcc, "%g", rel_tolerance)
TEST(c_weight_pccc, "%g", rel_tolerance)
TEST(c_weight_cccc, "%g", rel_tolerance)
TEST(c_weight_im_ppcc, "%g", rel_tolerance)

#undef TEST

#define TEST(field) \
	{ \
	REAL max_field=0.0; \
	for(i=0;i<pps_bins;i++) { \
		if((ref->field[i]!=test->field[i]) && !(fabs(test->field[i]-ref->field[i])<rel_tolerance*(fabs(test->field[i])+fabs(ref->field[i])))) { \
			fprintf(stderr, "%s" #field "[%d] fields mismatch ref=%g test=%g test-ref=%g\n", prefix, i, ref->field[i], test->field[i], test->field[i]-ref->field[i]); \
			return -4; \
			} \
		if(fabs(ref->field[i])>max_field)max_field=fabs(ref->field[i]); \
		if(fabs(test->field[i])>max_field)max_field=fabs(test->field[i]); \
		} \
	for(i=0;i<pps_bins;i++) { \
		if((ref->field[i]!=test->field[i]) && !(fabs(test->field[i]-ref->field[i])<rel_abs_tolerance*max_field)) { \
			fprintf(stderr, "%s" #field "[%d] fields mismatch ref=%g test=%g test-ref=%g max_field=%g\n", prefix, i, ref->field[i], test->field[i], test->field[i]-ref->field[i], max_field); \
			return -5; \
			} \
		} \
	}

TEST(power_pp)
TEST(power_pc)
TEST(power_cc)
	
if(ref->power_im_pc!=NULL || test->power_im_pc!=NULL) {
	if(ref->power_im_pc==NULL || test->power_im_pc==NULL) {
		fprintf(stderr, "mismatch: test->power_im_pc=%p ref->power_im_pc=%p\n", test->power_im_pc, ref->power_im_pc);
		} else {
		TEST(power_im_pc)
		}
	}
	

if(ref->weight_arrays_non_zero) {
	TEST(weight_pppp)
	TEST(weight_pppc)
	TEST(weight_ppcc)
	TEST(weight_pccc)
	TEST(weight_cccc)
	}

#undef TEST
return 0;
}

void SUFFIX(cblas_accumulate_partial_power_sum)(SUFFIX(PARTIAL_POWER_SUM) *accum, SUFFIX(PARTIAL_POWER_SUM) *partial)
{
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

CBLAS_AXPY(pps_bins, 1.0, &(partial->power_pp[shift]), 1, accum->power_pp, 1);
CBLAS_AXPY(pps_bins, 1.0, &(partial->power_pc[shift]), 1, accum->power_pc, 1);
CBLAS_AXPY(pps_bins, 1.0, &(partial->power_cc[shift]), 1, accum->power_cc, 1);

if(partial->power_im_pc!=NULL) {
	if(accum->power_im_pc==NULL) {
		fprintf(stderr, "*** INTERNAL ERROR: attempt to accumulate partial power sum with cross terms into accumulator without\n");
		exit(-1);
		}

	CBLAS_AXPY(pps_bins, 1.0, &(partial->power_im_pc[shift]), 1, accum->power_im_pc, 1);
	}


if(partial->weight_arrays_non_zero) {
	CBLAS_AXPY(pps_bins, 1.0, &(partial->weight_pppp[shift]), 1, accum->weight_pppp, 1);
	CBLAS_AXPY(pps_bins, 1.0, &(partial->weight_pppc[shift]), 1, accum->weight_pppc, 1);
	CBLAS_AXPY(pps_bins, 1.0, &(partial->weight_ppcc[shift]), 1, accum->weight_ppcc, 1);
	CBLAS_AXPY(pps_bins, 1.0, &(partial->weight_pccc[shift]), 1, accum->weight_pccc, 1);
	CBLAS_AXPY(pps_bins, 1.0, &(partial->weight_cccc[shift]), 1, accum->weight_cccc, 1);

	accum->weight_arrays_non_zero=1;
	}

accum->c_weight_pppp+=partial->c_weight_pppp;
accum->c_weight_pppc+=partial->c_weight_pppc;
accum->c_weight_ppcc+=partial->c_weight_ppcc;
accum->c_weight_pccc+=partial->c_weight_pccc;
accum->c_weight_cccc+=partial->c_weight_cccc;
accum->c_weight_im_ppcc+=partial->c_weight_im_ppcc;
}

void SUFFIX(sse_accumulate_partial_power_sum)(SUFFIX(PARTIAL_POWER_SUM) *accum, SUFFIX(PARTIAL_POWER_SUM) *partial)
{
int pps_bins=accum->nbins;
int shift;
if(accum->type!=sizeof(float)) {
	fprintf(stderr, "*** INTERNAL ERROR: %s power sum real type mismatch type=%d but must be 4 for sse mode\n",
		__FUNCTION__,
		accum->type);
	exit(-1);
	}
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

SSE_SUM(pps_bins, &(partial->power_pp[shift]), accum->power_pp);
SSE_SUM(pps_bins, &(partial->power_pc[shift]), accum->power_pc);
SSE_SUM(pps_bins, &(partial->power_cc[shift]), accum->power_cc);

if(partial->power_im_pc!=NULL) {
	if(accum->power_im_pc==NULL) {
		fprintf(stderr, "*** INTERNAL ERROR: attempt to accumulate partial power sum with cross terms into accumulator without\n");
		exit(-1);
		}

	SSE_SUM(pps_bins, &(partial->power_im_pc[shift]), accum->power_im_pc);
	}

if(partial->weight_arrays_non_zero) {
	SSE_SUM(pps_bins, &(partial->weight_pppp[shift]), accum->weight_pppp);
	SSE_SUM(pps_bins, &(partial->weight_pppc[shift]), accum->weight_pppc);
	SSE_SUM(pps_bins, &(partial->weight_ppcc[shift]), accum->weight_ppcc);
	SSE_SUM(pps_bins, &(partial->weight_pccc[shift]), accum->weight_pccc);
	SSE_SUM(pps_bins, &(partial->weight_cccc[shift]), accum->weight_cccc);

	accum->weight_arrays_non_zero=1;
	}

accum->c_weight_pppp+=partial->c_weight_pppp;
accum->c_weight_pppc+=partial->c_weight_pppc;
accum->c_weight_ppcc+=partial->c_weight_ppcc;
accum->c_weight_pccc+=partial->c_weight_pccc;
accum->c_weight_cccc+=partial->c_weight_cccc;
accum->c_weight_im_ppcc+=partial->c_weight_im_ppcc;
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

fprintf(out, "%g %g %g %g %g %g", 
	pps->c_weight_pppp,
	pps->c_weight_pppc,
	pps->c_weight_ppcc,
	pps->c_weight_pccc,
	pps->c_weight_cccc,
	pps->c_weight_im_ppcc
       );

for(i=0;i<pps_bins;i++) fprintf(out, " %g", pps->weight_pppp[i]);
for(i=0;i<pps_bins;i++) fprintf(out, " %g", pps->weight_pppc[i]);
for(i=0;i<pps_bins;i++) fprintf(out, " %g", pps->weight_ppcc[i]);
for(i=0;i<pps_bins;i++) fprintf(out, " %g", pps->weight_pccc[i]);
for(i=0;i<pps_bins;i++) fprintf(out, " %g", pps->weight_cccc[i]);

for(i=0;i<pps_bins;i++) fprintf(out, " %g", pps->power_pp[i]);
for(i=0;i<pps_bins;i++) fprintf(out, " %g", pps->power_pc[i]);
for(i=0;i<pps_bins;i++) fprintf(out, " %g", pps->power_cc[i]);

if(pps->power_im_pc!=NULL) {
	for(i=0;i<pps_bins;i++) fprintf(out, " %g", pps->power_im_pc[i]);
	} else {
	fprintf(out, " power_im_pc=NULL");
	}
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

if(pps->power_im_pc!=NULL)free(pps->power_im_pc);

memset(pps, 0, sizeof(*pps));
free(pps);
}
