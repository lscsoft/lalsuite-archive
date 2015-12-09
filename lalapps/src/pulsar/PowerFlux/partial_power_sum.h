
typedef REAL * __restrict__  SUFFIX(REALPTR);

typedef struct {
	int type;
	int nbins;

	REAL c_weight_pppp;
	REAL c_weight_pppc;
	REAL c_weight_ppcc;
	REAL c_weight_pccc;
	REAL c_weight_cccc;
	REAL c_weight_im_ppcc;

	REAL * __restrict__ weight_pppp;
	REAL * __restrict__ weight_pppc;
	REAL * __restrict__ weight_ppcc;
	REAL * __restrict__ weight_pccc;
	REAL * __restrict__ weight_cccc;
	/* REAL *weight_im_ppcc; commented out as loosely coherent search does not use bin-avoidance and does not use these arrays */	

	/* power sums - plus^2, plus*cross and cross^2*/
	REAL * __restrict__ power_pp;
	REAL * __restrict__ power_pc;
	REAL * __restrict__ power_cc;
	REAL * __restrict__ power_im_pc; /* imaginary part for loosely coherent statistic */

	REAL * __restrict__ memory_pool;
	
	int memory_pool_size;
	
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

