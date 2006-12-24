#ifndef __CANDIDATES_H__
#define __CANDIDATES_H__

typedef struct {
		
	int polarization_index;
	
	int rank;
	int opt_rank;
	float score;
	int point_index;
	int domain_size;
	float ul;
	float S_band;
	float M_band;
	float max_dx;
	float frequency;
	float ra;
	float dec;
	float spindown;
	float psi;
	float iota;
	float weight_ratio;
	int skyband;
	float coherence_score;
	float chi_sq;
	float power_cor;
	float snr;
	float strain;
	float strain_err;
	float total_weight;
	float f_max;
	float ifo_freq;
	float ifo_freq_sd;
	} CANDIDATE;

void identify_candidates(void);

#endif
