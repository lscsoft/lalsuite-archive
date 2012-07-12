#include <lal/LALInference.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <lal/FrequencySeries.h>



int free_waveform_interp_objects(struct twod_waveform_interpolant_array *);

struct twod_waveform_interpolant_array* new_waveform_interpolant_array_from_svd_bank(gsl_matrix *svd_bank);

struct twod_waveform_interpolant_patches* interpolants_over_patches(REAL8FrequencySeries *psd, int N_patches);

int XLALInferenceDestroyInterpManifold(struct twod_waveform_interpolant_manifold*);

struct twod_waveform_interpolant* new_waveform_interpolant_from_svd_bank(gsl_matrix*);

struct twod_waveform_interpolant_manifold* interpolants_manifold_init(REAL8FrequencySeries*, unsigned int, unsigned int , int , int , double , double , double , double , double , double , double , double , double , double, unsigned int, double );

int interpolate_waveform_from_mchirp_and_eta(struct twod_waveform_interpolant_array*, gsl_vector_complex*, double, double);

struct twod_waveform_interpolant_manifold *XLALInferenceCreateInterpManifold(REAL8FrequencySeries* , double, double, double, double, double, double, double, double);

int index_into_patch(struct twod_waveform_interpolant_manifold*, double, double);

double compute_chirp_time(double, double, double, int, double);

int dewhiten_template_wave(gsl_vector_complex*, COMPLEX16TimeSeries*, COMPLEX16FrequencySeries*, COMPLEX16FrequencySeries* 
,COMPLEX16FFTPlan*, REAL8FrequencySeries*, double, double);

double ffinal(double);
