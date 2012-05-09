#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <lal/FrequencySeries.h>

struct twod_waveform_interpolant {

	gsl_vector *svd_basis;
 
	/* See http://arxiv.org/pdf/1108.5618v1.pdf  This represents the C
 	 * matrix of formula (8) without mu.  Note that you specify a separate waveform
	 * interpolant object for each mu 
	 */

	gsl_matrix_complex *C_KL;

		
};
	
struct twod_waveform_interpolant_array {
	struct twod_waveform_interpolant *interp;
	unsigned int size;
	double param1_min;
	double param1_max;
	double param2_min;
	double param2_max;
	double inner_param1_min;
	double inner_param2_min;
	double inner_param1_max;
	double inner_param2_max;
};

struct twod_waveform_interpolant_manifold {
	REAL8FrequencySeries *psd;
	unsigned int patches_in_eta;
	unsigned int patches_in_mc;
	unsigned int number_templates_along_eta;	
	unsigned int number_templates_along_mc;
	struct twod_waveform_interpolant_array *interp_arrays;
	double eta_padding;
	double mc_padding;
	double inner_param1_min;
	double inner_param1_max;
	double inner_param2_min;
	double inner_param2_max;
	double outer_param1_min;
	double outer_param1_max;
	double outer_param2_min;
	double outer_param2_max;
};

int free_waveform_interp_objects(struct twod_waveform_interpolant_array *);

struct twod_waveform_interpolant_array* new_waveform_interpolant_array_from_svd_bank(gsl_matrix *svd_bank);

struct twod_waveform_interpolant_patches* interpolants_over_patches(REAL8FrequencySeries *psd, int N_patches);
