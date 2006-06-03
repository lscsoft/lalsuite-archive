#ifndef __CANDIDATES_H__
#define __CANDIDATES_H__

typedef struct {
		
	int point_index;
	int polarization_index;
	
	int rank;
	int domain_size;
	float max_dx;
	float frequency;
	float ra;
	float dec;
	float spindown;
	float a_plus;
	float a_cross;
	float weight_ratio;
	int skyband;
	float coherence_score;
	} CANDIDATE;

void identify_candidates(void);

#endif
