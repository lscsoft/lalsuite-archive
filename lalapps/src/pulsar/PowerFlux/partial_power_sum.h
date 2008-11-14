
typedef struct {
	REAL c_weight_pppp;
	REAL c_weight_pppc;
	REAL c_weight_ppcc;
	REAL c_weight_pccc;
	REAL c_weight_cccc;

	REAL *weight_pppp;
	REAL *weight_pppc;
	REAL *weight_ppcc;
	REAL *weight_pccc;
	REAL *weight_cccc;

	/* power sums - plus^2, plus*cross and cross^2*/
	REAL *power_pp;
	REAL *power_pc;
	REAL *power_cc;

	int weight_arrays_non_zero;
	int collapsed_weight_arrays;
	} SUFFIX(PARTIAL_POWER_SUM);

SUFFIX(PARTIAL_POWER_SUM) *SUFFIX(allocate_partial_power_sum)(void);
void SUFFIX(zero_partial_power_sum)(SUFFIX(PARTIAL_POWER_SUM) *pps);
void SUFFIX(accumulate_partial_power_sum)(SUFFIX(PARTIAL_POWER_SUM) *accum, SUFFIX(PARTIAL_POWER_SUM) *partial);
void SUFFIX(free_partial_power_sum)(SUFFIX(PARTIAL_POWER_SUM) *pps);

