
typedef struct {
	int type;
	int nbins;

	REAL c_weight_pppp;
	REAL c_weight_pppc;
	REAL c_weight_ppcc;
	REAL c_weight_pccc;
	REAL c_weight_cccc;
	REAL c_weight_im_ppcc;

	REAL *weight_pppp;
	REAL *weight_pppc;
	REAL *weight_ppcc;
	REAL *weight_pccc;
	REAL *weight_cccc;
	/* REAL *weight_im_ppcc; commented out as loosely coherent search does not use bin-avoidance and does not use these arrays */	

	/* power sums - plus^2, plus*cross and cross^2*/
	REAL *power_pp;
	REAL *power_pc;
	REAL *power_cc;
	REAL *power_im_pc; /* imaginary part for loosely coherent statistic */

	int offset; /* this is the index of the bin with index 0 in the output (i.e. firstbin) */

	int weight_arrays_non_zero;
	int collapsed_weight_arrays;
	} SUFFIX(PARTIAL_POWER_SUM);

SUFFIX(PARTIAL_POWER_SUM) *SUFFIX(allocate_partial_power_sum)(int pps_bins, int cross_terms_present);
void SUFFIX(zero_partial_power_sum)(SUFFIX(PARTIAL_POWER_SUM) *pps);
void SUFFIX(randomize_partial_power_sum)(SUFFIX(PARTIAL_POWER_SUM) *pps);
void SUFFIX(accumulate_partial_power_sum)(SUFFIX(PARTIAL_POWER_SUM) *accum, SUFFIX(PARTIAL_POWER_SUM) *partial);
int SUFFIX(compare_partial_power_sums)(char *prefix, SUFFIX(PARTIAL_POWER_SUM) *ref, SUFFIX(PARTIAL_POWER_SUM) *test, REAL rel_tolerance, REAL rel_abs_tolerance);
void SUFFIX(cblas_accumulate_partial_power_sum)(SUFFIX(PARTIAL_POWER_SUM) *accum, SUFFIX(PARTIAL_POWER_SUM) *partial);
void SUFFIX(sse_accumulate_partial_power_sum)(SUFFIX(PARTIAL_POWER_SUM) *accum, SUFFIX(PARTIAL_POWER_SUM) *partial);
void SUFFIX(free_partial_power_sum)(SUFFIX(PARTIAL_POWER_SUM) *pps);
void SUFFIX(dump_partial_power_sum)(FILE *out, SUFFIX(PARTIAL_POWER_SUM) *pps);

