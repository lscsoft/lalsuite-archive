/*
*  Copyright (C) 2007 Vladimir Dergachev
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

#ifndef __CANDIDATES_H__
#define __CANDIDATES_H__

typedef struct {
	int better_candidate;
		
	int max_dx_polarization_index;
	
	int rank;
	int opt_rank;
	float score;
	int point_index;
	int domain_size;
	float ul;
	float S_band;
	float M_band;
	float max_dx;
	double frequency;
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
	double f_max;
	double ifo_freq;
	float ifo_freq_sd;
	} CANDIDATE;

void init_candidates(void);
void identify_candidates(void);
void output_candidates(FILE *fout);
void compute_scores(CANDIDATE *cand, int debug);
void single_pass_compute_matched_snr(CANDIDATE *cand);
void single_pass_compute_simple_snr(CANDIDATE *cand);

#endif
