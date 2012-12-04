/* Use old lal complex structures for now */
#define LAL_USE_OLD_COMPLEX_STRUCTS

#include <LALInferenceInterps.h>
#include <math.h>
#include <complex.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_gamma.h>
/* LAL includes */
#include <lal/LALDatatypes.h>
#include <lal/Units.h>
#include <lal/TimeFreqFFT.h>
#include <lal/ComplexFFT.h>
#include <lal/RealFFT.h>
#include <lal/TimeSeries.h>
#include <lal/FrequencySeries.h>
#include <lal/LALConstants.h>
#include <lal/Sequence.h>
#include <assert.h>
#include <lal/LALComplex.h>
#include <lal/VectorOps.h>

/*
 * Data structure methods
 */



int free_waveform_interp_objects(struct twod_waveform_interpolant_array * interps) {
	unsigned int i;
	struct twod_waveform_interpolant *interp = interps->interp;
	for (i = 0; i < interps->size; i++, interp++) {
		gsl_matrix_complex_free(interp->C_KL);
		gsl_vector_free(interp->svd_basis);	
		//free(&interps->interp[i]);
	}
	//free(&interps->interp[i]);
	//free(interps);
	return 0;

}


int XLALInferenceDestroyInterpManifold(struct twod_waveform_interpolant_manifold *manifold){

	unsigned int i;
		
	for(i=0; i < manifold->patches_in_eta*manifold->patches_in_mc; i++){

		free_waveform_interp_objects(&manifold->interp_arrays[i]);
	}	
	free(manifold->interp_arrays);	
	free(manifold);

	return 0;
}
struct twod_waveform_interpolant * new_waveform_interpolant_from_svd_bank(gsl_matrix *svd_bank){

	unsigned int i,j;

	struct twod_waveform_interpolant *interp = (struct twod_waveform_interpolant *) calloc(svd_bank->size2, sizeof(struct twod_waveform_interpolant));

	for (i = 0; i < svd_bank->size2; i++) {
		interp[i].svd_basis = gsl_vector_calloc(svd_bank->size1);
		for (j = 0; j < svd_bank->size1; j++){
			gsl_vector_set(interp[i].svd_basis, j, gsl_matrix_get(svd_bank, j, i));

		}
	}
	return interp;
}

struct twod_waveform_interpolant_manifold* interpolants_manifold_init(REAL8FrequencySeries *psd, unsigned int patches_in_eta, unsigned int patches_in_mc, int number_templates_along_mc, int number_templates_along_eta, double mc_min, double mc_max, double eta_min, double eta_max, double outer_mc_min, double outer_mc_max, double outer_eta_min, double outer_eta_max, double mc_padding, double eta_padding){

	struct twod_waveform_interpolant_manifold * output = (struct twod_waveform_interpolant_manifold *) calloc( 1, sizeof(struct twod_waveform_interpolant_manifold));

	output->interp_arrays = (struct twod_waveform_interpolant_array *) calloc(patches_in_eta*patches_in_mc, sizeof(struct twod_waveform_interpolant_array));
	output->psd = psd;
	output->number_templates_along_eta = number_templates_along_eta; 
	output->number_templates_along_mc = number_templates_along_mc;
	output->patches_in_eta = patches_in_eta;
	output->patches_in_mc = patches_in_mc;
	output->mc_padding = mc_padding;
	output->eta_padding = eta_padding;
	output->inner_param1_min = mc_min;
	output->inner_param1_max = mc_max;
	output->inner_param2_min = eta_min;
	output->inner_param2_max = eta_max;
	output->outer_param1_min = outer_mc_min;
	output->outer_param1_max = outer_mc_max;
	output->outer_param2_min = outer_eta_min;
	output->outer_param2_max = outer_eta_max;
	
	return output;
}


/*
 * Formula (3) in http://arxiv.org/pdf/1108.5618v1.pdf
 */ 



static int projection_coefficient(gsl_vector *svd_basis, gsl_matrix *template_bank, gsl_matrix_complex *M_xy, unsigned int N_mc, unsigned int M_eta){
	unsigned int i,j,k;
	/*project svd basis onto SPA's to get coefficients*/
	gsl_complex M;
	double M_real, M_imag;
	/* compute M_xy at fixed mu */

	assert(N_mc * M_eta == template_bank->size2 / 2);

	GSL_SET_COMPLEX(&M, 0, 0);

	for (k =0; k < template_bank->size2 / 2; k++){

		j = k % M_eta;	

		if ( !(k % M_eta) ){

			i = k/M_eta;
		}
		gsl_vector_view spa_waveform_real = gsl_matrix_column(template_bank, 2*k);
		gsl_vector_view spa_waveform_imag = gsl_matrix_column(template_bank, 2*k+1);
		gsl_blas_ddot(&spa_waveform_real.vector, svd_basis, &M_real);
	        gsl_blas_ddot(&spa_waveform_imag.vector, svd_basis, &M_imag);
	
		GSL_SET_COMPLEX(&M, M_real, M_imag);
		gsl_matrix_complex_set(M_xy, i, j, M);
	}


	return 0;
}

static gsl_matrix_complex *measure_m0_phase(gsl_matrix_complex* M_xy, gsl_matrix_complex* phase_M0_xy){

	double phi_tmp;
	gsl_complex z_tmp;

	for(unsigned int i=0; i < M_xy->size1; i++){
		for(unsigned int j=0; j < M_xy->size2; j++){

			phi_tmp = atan2( GSL_IMAG(gsl_matrix_complex_get(M_xy, i, j) ), GSL_REAL(gsl_matrix_complex_get(M_xy, i, j)) ) ;	
			GSL_SET_COMPLEX(&z_tmp, cos(phi_tmp), -sin(phi_tmp));
			gsl_matrix_complex_set(phase_M0_xy, i, j, z_tmp);

		}
	}

	return phase_M0_xy;
}

static gsl_matrix_complex *rotate_M_xy(gsl_matrix_complex* M_xy, gsl_matrix_complex* phase_M0_xy){

			
	gsl_matrix_complex_mul_elements(M_xy, phase_M0_xy);	

	
	return M_xy;
}

/*
 * Formula (4) in http://arxiv.org/pdf/1108.5618v1.pdf
 */

static double chebyshev_node(int j, int J_max) {
	return cos( M_PI * (j + 0.5) / ( J_max )  );
}

/* 
 * Formula (5) in http://arxiv.org/pdf/1108.5618v1.pdf
 */



static double onedCheby(double x, int J, int J_max) {
	double w=0;
	double T_J=0;
	int count;
	int i;
	count =  J/2 + 1;

	for(i=0; i < count; i++){
	
		T_J+=pow( (x*x - 1.), i )*pow( x,(J-2*i) ) * (double) gsl_sf_choose(J, 2*i);
	
		
	}

	if (J !=0 ){
                w = 1. / sqrt((J_max)/2.);
        }
        else{
                w = 1. / sqrt(J_max);
        }

	T_J*=w;

	return T_J;

}

static double twodCheby(double x, int K, int K_max, double y, int L, int L_max) {
	return onedCheby(x, K, K_max) * onedCheby(y, L, L_max);
}


/* 
 * Formula (7) in http://arxiv.org/pdf/1108.5618v1.pdf
 */


/* you do this for every mu */
static gsl_matrix_complex  *compute_C_KL(gsl_vector *x_k, gsl_vector *y_l, gsl_matrix_complex *M) {
	int K, L, k, l;
	gsl_complex out;
	int k_max = x_k->size;
	int l_max = y_l->size;
	gsl_matrix_complex *C_KL = gsl_matrix_complex_calloc(k_max, l_max);
	gsl_complex tmp;

	GSL_SET_COMPLEX(&out, 0, 0);

	for (K = 0; K < k_max; K++) {
		for (L = 0; L < l_max; L++) {
			for (k = 0; k < k_max; k++) {
				for (l = 0; l < l_max; l++) {
					tmp = gsl_complex_mul_real(gsl_matrix_complex_get(M, k, l), twodCheby(gsl_vector_get(x_k, k), K, k_max, gsl_vector_get(y_l, l), L, l_max));
					out = gsl_complex_add(out, tmp);
				}
			}

			gsl_matrix_complex_set(C_KL, K, L, out);
			GSL_SET_COMPLEX(&out, 0, 0);
		}
	}

	return C_KL;
}


/* 
 * Formula (8) in http://arxiv.org/pdf/1108.5618v1.pdf
 */

static gsl_complex compute_M_xy(gsl_matrix_complex *C_KL, double x, double y) {
	int K_max, L_max;
	int K, L;
	gsl_complex M; 
	gsl_complex tmp;
	
	K_max = C_KL->size1;
	L_max = C_KL->size2;


	GSL_SET_COMPLEX(&M, 0, 0);
	GSL_SET_COMPLEX(&tmp, 0, 0);

	for (K = 0; K < K_max; K++) {
		for (L = 0; L < L_max; L++) {
			/* FIXME, is this correct?? */
			tmp =gsl_complex_mul_real(gsl_matrix_complex_get(C_KL, K, L),twodCheby(x, K, K_max, y, L, L_max));
			M = gsl_complex_add(M, tmp);
		}

	}	

	return M;
}

/*
 * generic cheby utilities
 */

static double map_coordinate_to_cheby(double c_min, double c_max, double c) {

	return  2. * ( c - c_min ) / (c_max - c_min) - 1.;
}

static double map_node_to_coordinate(double c_min, double c_max, double x_node) {


	return  c_min + ( ( x_node + 1.) / 2. ) * (c_max - c_min);

}

/* waveform template stuff */

static double mc2mass1(double mc, double eta)
/* mass 1 (the smaller one) for given mass ratio & chirp mass */
{
 double root = sqrt(0.25-eta);
 double fraction = (0.5+root) / (0.5-root);
 return mc * (pow(1+fraction,0.2) / pow(fraction,0.6));
}


static double mc2mass2(double mc, double eta)
/* mass 2 (the greater one) for given mass ratio & chirp mass */
{
 double root = sqrt(0.25-eta);
 double inversefraction = (0.5-root) / (0.5+root);
 return mc * (pow(1+inversefraction,0.2) / pow(inversefraction,0.6));
}





static double compute_chirp_time (double m1, double m2, double fLower, int order, double chi)
	{

	/* variables used to compute chirp time */
	double c0T, c2T, c3T, c4T, c5T, c6T, c6LogT, c7T;
	double xT, x2T, x3T, x4T, x5T, x6T, x7T, x8T;
	double m = m1 + m2;
	double eta = m1 * m2 / m / m;

	c0T = c2T = c3T = c4T = c5T = c6T = c6LogT = c7T = 0.;

	/* Switch on PN order, set the chirp time coeffs for that order */
	switch (order)
		{
		case 8:
		case 7:
			c7T = LAL_PI * (14809.0 * eta * eta - 75703.0 * eta / 756.0 - 15419335.0 / 127008.0);
		case 6:
			c6T = LAL_GAMMA * 6848.0 / 105.0 - 10052469856691.0 / 23471078400.0 + LAL_PI * LAL_PI * 128.0 / 3.0 + eta * (3147553127.0 / 3048192.0 - LAL_PI * LAL_PI * 451.0 / 12.0) - eta * eta * 15211.0 / 1728.0 + eta * eta * eta * 25565.0 / 1296.0 + log (4.0) * 6848.0 / 105.0;
     			c6LogT = 6848.0 / 105.0;
		case 5:
			c5T = 13.0 * LAL_PI * eta / 3.0 - 7729.0 / 252.0 - (0.4*565.*(-146597. + 135856.*eta + 17136.*eta*eta)*chi/(2268.*(-113. + 76.*eta))); /* last term is 0 if chi is 0*/
		case 4:
			c4T = 3058673.0 / 508032.0 + eta * (5429.0 / 504.0 + eta * 617.0 / 72.0) + (0.4*63845.*(-81. + 4.*eta)*chi*chi/(8.*pow(-113. + 76.*eta, 2.))); /* last term is 0 if chi is 0*/
			c3T = -32.0 * LAL_PI / 5.0 + (0.4*113.*chi/3.); /* last term is 0 if chi is 0*/
			c2T = 743.0 / 252.0 + eta * 11.0 / 3.0;
			c0T = 5.0 * m * LAL_MTSUN_SI / (256.0 * eta);
			break;
		default:
			fprintf (stderr, "ERROR!!!\n");
			break;
		}

	/* This is the PN parameter v evaluated at the lower freq. cutoff */
	xT = pow (LAL_PI * m * LAL_MTSUN_SI * fLower, 1.0 / 3.0);
	x2T = xT * xT;
	x3T = xT * x2T;
	x4T = x2T * x2T;
	x5T = x2T * x3T;
	x6T = x3T * x3T;
	x7T = x3T * x4T;
	x8T = x4T * x4T;

	/* Computes the chirp time as tC = t(v_low)    */
	/* tC = t(v_low) - t(v_upper) would be more    */
	/* correct, but the difference is negligble.   */

	/* This formula works for any PN order, because */
	/* higher order coeffs will be set to zero.     */

	return c0T * (1 + c2T * x2T + c3T * x3T + c4T * x4T + c5T * x5T + (c6T + c6LogT * log (xT)) * x6T + c7T * x7T) / x8T;
}

static double ffinal(double m_total){
	
	// Compute frequency at Schwarzschild ISCO 

	double f_isco;
	
	f_isco = pow(2., ceil( log( (1./LAL_PI)*( pow(6.,-3./2.) )*( pow((m_total)*LAL_MTSUN_SI,-1.) ) )  / log(2.) ) );  //Next highest power of 2 of f_isco 
	
	return f_isco;
}


static int compute_max_chirp_time_and_max_frequency(double mc_min, double mc_max, double eta_min, double eta_max, double f_min, double *f_max, double *t_max) {

	double m1, m2, mc, eta, t_max_tmp;
	int N_mc = 10;
	int M_eta = 10;
	for(int i = 0; i < N_mc; i++){
		for(int j = 0; j < M_eta; j++){

			mc = mc_min + i*(mc_max - mc_min)/(N_mc - 1);
			eta = eta_min + j*(eta_max - eta_min)/(M_eta - 1);
			m1 = mc2mass1(mc, eta);
			m2 = mc2mass2(mc, eta);
			t_max_tmp = compute_chirp_time(m1, m2, f_min, 4, 0);
			if (t_max_tmp > *t_max ) *t_max = t_max_tmp;
			

		}
	} 
	eta_min +=0.;	
	double mt_min = mc_min/(pow(eta_max, 3./5.));
	*f_max = ffinal(mt_min);
	*f_max = 2048.;
	return 0;
}

/*static int compute_min_chirp_time(double mc_min, double mc_max, double eta_min, double eta_max, double f_min){
	double m1, m2, mc, eta, t_min_tmp;
	int N_mc = 10;
	int M_eta = 10;

	eta = eta_min + 1.*(eta_max - eta_min)/(M_eta - 1);
	mc = mc_min + 1.*(mc_max - mc_min)/(N_mc - 1);
	m1 = mc2mass1(mc, eta);
	m2 = mc2mass2(mc, eta);
	double t_min = compute_chirp_time(m1, m2, f_min, 4, 0);
	
	for(int i = 0; i < N_mc; i++){
		for(int j = 0; j < M_eta; j++){

			mc = mc_min + i*(mc_max - mc_min)/(N_mc - 1);
                        eta = eta_min + j*(eta_max - eta_min)/(M_eta - 1);
                        m1 = mc2mass1(mc, eta);
                        m2 = mc2mass2(mc, eta);
                        t_min_tmp = compute_chirp_time(m1, m2, f_min, 4, 0);
                        if (t_min_tmp < t_min ) t_min = t_min_tmp;
			}
		}
	return t_min;
}
*/


static int generate_template(double m1, double m2, double f_low, double order, COMPLEX16FrequencySeries *hOfF, COMPLEX16TimeSeries *tmp_tseries, COMPLEX16FFTPlan *fwdplan, int shift_flag){

	double templateReal, templateImag, twopit, f, re, im, t_shift, f_ref;
	
	REAL8TimeSeries *hplus = NULL;
	REAL8TimeSeries *hcross = NULL;

	//LIGOTimeGPS epoch = LIGOTIMEGPSZERO;

	memset(tmp_tseries->data->data, 0, tmp_tseries->data->length * sizeof (COMPLEX16));

        unsigned int n = hOfF->data->length;


	double m1SI = m1*LAL_MSUN_SI;
	double m2SI = m2*LAL_MSUN_SI;

	double d = 1e6 * LAL_PC_SI;

	XLALSimInspiralTaylorT4PNGenerator( &hplus, &hcross, 0, 1, tmp_tseries->deltaT, m1SI, m2SI, f_low - 10, d, 0, order, order);
	//XLALSimIMREOBNRv2DominantMode(&hplus, &hcross, 0, tmp_tseries->deltaT, m1SI, m2SI, f_low - 15, d, 0);

	memset (hOfF->data->data, 0, hOfF->data->length * sizeof (COMPLEX16));


	if (shift_flag){

		f_ref = 157.;//152.84375;
        	t_shift = -compute_chirp_time(m1, m2, f_ref, order, 0);
	        //int number_of_bins_to_shift = ceil(t_shift/tmp_tseries->deltaT);
		for( unsigned int i =0; i < hplus->data->length; i++){

			tmp_tseries->data->data[tmp_tseries->data->length - 1 - i].re = hplus->data->data[hplus->data->length - 1 - i]; 
			tmp_tseries->data->data[tmp_tseries->data->length - 1 - i].im = hcross->data->data[hplus->data->length - 1 - i];
		}
		/*unsigned int i_max_minus_one = 99;
		for ( unsigned int i =0; i < i_max_minus_one + 1 ; i++){
			tmp_tseries->data->data[tmp_tseries->data->length - 1 - i_max_minus_one - number_of_bins_to_shift + i].re *= (1 - tanh(i));
			tmp_tseries->data->data[tmp_tseries->data->length - 1 - i_max_minus_one - number_of_bins_to_shift + i].im *= (1 - tanh(i));
			tmp_tseries->data->data[tmp_tseries->data->length - 1 - number_of_bins_to_shift - hplus->data->length - 1 + i].re *= tanh(i);
			tmp_tseries->data->data[tmp_tseries->data->length - 1 - number_of_bins_to_shift - hplus->data->length - 1 + i].im *= tanh(i);	
		
		}*/

	}
	else{
		f_ref = 0;
		t_shift = 0;
		for( unsigned int i =0; i < hplus->data->length; i++){
			tmp_tseries->data->data[tmp_tseries->data->length - 1 - i].re = hplus->data->data[hplus->data->length - 1 - i];
			tmp_tseries->data->data[tmp_tseries->data->length - 1 - i].im = hcross->data->data[hplus->data->length - 1 - i];
		}
	}

	XLALCOMPLEX16TimeFreqFFT (hOfF, tmp_tseries, fwdplan);	

	for (unsigned int l = 0; l < hOfF->data->length/2 + 1; l++){
		templateReal = hOfF->data->data[hOfF->data->length - 1 - (n/2 - 1) + l].re;

      		templateImag = hOfF->data->data[hOfF->data->length - 1 - (n/2 - 1) + l].im;
      		twopit = LAL_TWOPI * (t_shift);
      		f = ((double) l) * hOfF->deltaF;
      		re = cos(twopit * f );
      		im =  -sin(twopit * f );
      		hOfF->data->data[hOfF->data->length - 1 - (n/2 - 1) + l].re = -templateReal*re + templateImag*im;
     	 	hOfF->data->data[hOfF->data->length - 1 - (n/2 - 1) + l].im = -templateImag*re - templateReal*im;
		if (l < hOfF->data->length/2){
			hOfF->data->data[l].re = 0;
			hOfF->data->data[l].im = 0;
		}
		/*
                if(l < ceil(f_low/hOfF->deltaF)){
                        hOfF->data->data[l].re = 0;
                        hOfF->data->data[l].im = 0;

                }*/
	}
	XLALDestroyREAL8TimeSeries(hplus);
	XLALDestroyREAL8TimeSeries(hcross);		

	return 0;
}

static int get_psd_from_file(REAL8FrequencySeries *series,  const char *fname){
	
	FILE *fp = fopen(fname, "r");
	double deltaF = series->deltaF;
	double f0 = series->f0;
	double f;
	const gsl_interp_type *t = gsl_interp_linear;
	gsl_interp *g_interp = gsl_interp_alloc (t, 65537);
	gsl_interp_accel *acc = gsl_interp_accel_alloc();
	double freq[65537];
	double psd[65537];
	
	unsigned int i = 0;
	while (!feof(fp)) {
		if (i >= 65537) break;
		fscanf(fp, "%lf %le\n", &freq[i], &psd[i]);
		i++;
	}

	gsl_interp_init(g_interp, freq, psd, 65537);
	for (i = 0; i < series->data->length; i++) {
		f = f0 + i * deltaF;
		series->data->data[i] = gsl_interp_eval(g_interp, freq, psd, f, acc);
	}
	gsl_interp_free(g_interp);
	gsl_interp_accel_free (acc);
	fclose(fp);
	return 0;
}

static int freq_to_time_fft(COMPLEX16FrequencySeries *fseries, COMPLEX16TimeSeries *tseries, COMPLEX16FFTPlan *revplan){
	
	XLALCOMPLEX16FreqTimeFFT(tseries, fseries, revplan);
	return 0;
}

/*
 * High level functions
 */

static int compute_working_length_and_sample_rate(double chirp_time, double f_max, unsigned int *working_length, double *sample_rate, unsigned int *length_max) {
	
	double duration = pow(2., ceil(log(chirp_time) / log(2.))); /* see SPADocstring in _spawaveform.c */
	
	*sample_rate = pow(2., ceil(log(2.* f_max) / log(2.)));
	*length_max =  round(*sample_rate * duration);
	*working_length = 32.0*(*sample_rate);//(unsigned int) round(pow(2., ceil(log(*length_max + round(32.0 * *(sample_rate))) / log(2.))));
	return 0;
}

static int initialize_time_and_freq_series(REAL8FrequencySeries **psd_ptr, COMPLEX16FrequencySeries **fseries_ptr, COMPLEX16FrequencySeries **fseries_for_ifft_ptr, COMPLEX16TimeSeries **tseries_ptr, COMPLEX16FFTPlan **revplan_ptr, double mc_min, double mc_max, double eta_min, double eta_max, double f_min, unsigned int *length_max, int tukey_flag){

	double f_max=0;
	double t_max=0;
	double deltaT=0;
	double sample_rate=0;
	double deltaF=0;
	double working_duration=0;
	unsigned int working_length=0;
	LIGOTimeGPS epoch = LIGOTIMEGPSZERO;

	REAL8FrequencySeries *psd = NULL;
	COMPLEX16TimeSeries *tseries = NULL;
	COMPLEX16FrequencySeries *fseries = NULL;
	COMPLEX16FrequencySeries *fseries_for_ifft = NULL;
	COMPLEX16FFTPlan *revplan = NULL;
	
	compute_max_chirp_time_and_max_frequency(mc_min, mc_max, eta_min, eta_max, f_min, &f_max, &t_max);
	//t_max*=2;
	fprintf(stderr, "f_max %e t_max %e\n", f_max, t_max);

	compute_working_length_and_sample_rate(t_max, f_max, &working_length, &sample_rate, length_max);

	deltaT = 1. / sample_rate;
        working_duration = (working_length / sample_rate);
	deltaF = 1. / working_duration;
	
	/*deltaF = 0.031250;
	deltaT = 1. / sample_rate;
        working_duration = 1./deltaF;
        working_length = round(working_duration*sample_rate);
	*/
	fprintf(stderr, "working_length %d sample_rate %e length_max: %d deltaT: %f deltaF: %f working_duration: %f\n", working_length, sample_rate, *length_max, deltaT, deltaF, working_duration);
	/* set up time series */	
	tseries = XLALCreateCOMPLEX16TimeSeries(NULL, &epoch, 0., deltaT, &lalDimensionlessUnit, working_length);
	memset (tseries->data->data, 0, tseries->data->length * sizeof (COMPLEX16));
	
	// 0 and positive frequencies only 
	fseries = XLALCreateCOMPLEX16FrequencySeries(NULL, &epoch, 0, deltaF, &lalDimensionlessUnit, working_length);
	memset (fseries->data->data, 0, fseries->data->length * sizeof (COMPLEX16));	
	psd = XLALCreateREAL8FrequencySeries(NULL, &epoch, 0, deltaF, &lalDimensionlessUnit, working_length / 2 + 1);
	memset (psd->data->data, 0, psd->data->length * sizeof (REAL8));
	
	/* FIXME: get_psd_from_file needs to change so that it takes in a real8freqseries and returns a new real8freqseries interpolated from the input */
	get_psd_from_file(psd, "H1-32sPSD.dat");
	//interpolate_psd_from_REAL8FrequencySeries(psd, psd_to_interpolate);

	if (tukey_flag == 1){
	
		REAL8TimeSeries* t_dom_psd;
		REAL8TimeSeries* tmp_storage;
	
		REAL8Window* window;

		COMPLEX16FrequencySeries* complex_psd;
	
		window = XLALCreateTukeyREAL8Window(working_length, 0.5);

		t_dom_psd = XLALCreateREAL8TimeSeries(NULL, &epoch, 0, deltaT, &lalDimensionlessUnit, working_length);
		
		tmp_storage = XLALCreateREAL8TimeSeries(NULL, &epoch, 0, deltaT, &lalDimensionlessUnit, working_length);
	
		complex_psd = XLALCreateCOMPLEX16FrequencySeries(NULL, &epoch, 0, deltaF, &lalDimensionlessUnit, working_length / 2 + 1);
	
		for (unsigned int i = 0; i < psd->data->length; i++){
			complex_psd->data->data[i].re = psd->data->data[i];
			complex_psd->data->data[i].im = 0;
		}
	
		REAL8FFTPlan *psd_plan;
		psd_plan = XLALCreateReverseREAL8FFTPlan(working_length, 1);	
		
		XLALREAL8FreqTimeFFT(t_dom_psd, complex_psd, psd_plan);			

		for (unsigned int i = 0; i < t_dom_psd->data->length; i++){
			tmp_storage->data->data[i] = t_dom_psd->data->data[i];
			
		}


		for (unsigned int i = 0; i < t_dom_psd->data->length; i++){ 
	
			if (i >= working_length/2){
				t_dom_psd->data->data[i] = tmp_storage->data->data[ i - working_length/2 ];
			}
			
			if (i < working_length/2){
				t_dom_psd->data->data[i] = tmp_storage->data->data[ i + working_length/2 ];
			
			}		

		}		
    
		/* Tukey PSD
		 */
 
		for (unsigned int i = 0; i < t_dom_psd->data->length; i++){
			t_dom_psd->data->data[i] *= window->data->data[i];

		}

		for (unsigned int i = 0; i < tmp_storage->data->length; i++){

			tmp_storage->data->data[i] = t_dom_psd->data->data[i];
		}
		for (unsigned int i = 0; i < t_dom_psd->data->length; i++){
	
			if (i >= working_length/2){
	                        t_dom_psd->data->data[i] = tmp_storage->data->data[ i - working_length/2];
	                }
			
			if (i < working_length/2){
	                        t_dom_psd->data->data[i] = tmp_storage->data->data[ i + working_length/2];
	
	                }	
	
	        }
	
		REAL8FFTPlan *fwd_psd_plan;
		fwd_psd_plan = XLALCreateForwardREAL8FFTPlan(working_length, 1);
       
		XLALREAL8TimeFreqFFT(complex_psd, t_dom_psd, fwd_psd_plan);
		for (unsigned int i = 0; i < complex_psd->data->length; i++){
			psd->data->data[i] = complex_psd->data->data[i].re;
		}

	   	XLALDestroyCOMPLEX16FrequencySeries(complex_psd);
        	XLALDestroyREAL8TimeSeries(t_dom_psd);
        	XLALDestroyREAL8TimeSeries(tmp_storage);
        	XLALDestroyREAL8FFTPlan(fwd_psd_plan);
        	XLALDestroyREAL8FFTPlan(psd_plan);
        	XLALDestroyREAL8Window(window);


	}
	
	// full frequency series	
	fseries_for_ifft = XLALCreateCOMPLEX16FrequencySeries(NULL, &epoch, 0, deltaF, &lalDimensionlessUnit, working_length);
	memset (fseries_for_ifft->data->data, 0, fseries_for_ifft->data->length * sizeof (COMPLEX16));	
	revplan = XLALCreateReverseCOMPLEX16FFTPlan(working_length, 1);
	
	*psd_ptr = psd;
	*fseries_ptr = fseries;
	*fseries_for_ifft_ptr = fseries_for_ifft;
	*tseries_ptr = tseries;
	*revplan_ptr = revplan;	
	
	return 0;

}

/*static int add_quadrature_phase(COMPLEX16FrequencySeries* fseries, COMPLEX16FrequencySeries* fseries_for_ifft){
	
	unsigned int n = fseries_for_ifft->data->length;	

	fseries->data->data[0].re = 0;
	fseries->data->data[0].im = 0;
	
	if( ! (n % 2) ){
		for (unsigned int i=1; i < (n/2); i++){		
			fseries_for_ifft->data->data[ fseries_for_ifft->data->length - 1 - (n/2 - 1) + i ].re = fseries->data->data[i].re / sqrt(2.*fseries->deltaF);
			fseries_for_ifft->data->data[ fseries_for_ifft->data->length - 1 - (n/2 - 1) + i ].im = fseries->data->data[i].im / sqrt(2.*fseries->deltaF);
		}
	}

	return 0;

}*/

static gsl_vector *even_param_spacing(double min, double max, int count) {
	gsl_vector *out = gsl_vector_calloc(count);
	int i;
	for (i = 0; i < count; i++) {
		gsl_vector_set(out, i, min + (i/(count-1.))*(max - min));
	}
	return out;
}

static gsl_vector *node_param_spacing(double min, double max, gsl_vector* x_nodes) {

	gsl_vector *out = gsl_vector_calloc(x_nodes->size);
	unsigned int i;

	for (i = 0; i < out->size; i++) {

		gsl_vector_set( out, i, map_node_to_coordinate( min, max, gsl_vector_get(x_nodes, i) ) );

	}
	return out;
}

static gsl_vector *raw_nodes(int count) {
	int i;
	gsl_vector *out = gsl_vector_calloc(count);
	for (i = 0; i < count; i++){
		gsl_vector_set(out, i, chebyshev_node(i, count));
	}
	return out;
}

static gsl_vector *linear_space_for_interpolation(double min, double max, unsigned int count){

        gsl_vector *out = gsl_vector_calloc(count);
        double x;
        for(unsigned int i=0; i < count; i++){

                x = -1 + (i/(count-1.))*( 1 - (-1) );
                gsl_vector_set( out, i, map_node_to_coordinate( min, max, x) );

        }

        return out;
}



static int generate_whitened_template(	double m1, double m2, 
					double f_min, unsigned int length_max,
					int order, REAL8FrequencySeries* psd, gsl_vector* template_real,
					gsl_vector* template_imag, COMPLEX16TimeSeries* tseries, COMPLEX16FrequencySeries* fseries,
					COMPLEX16FrequencySeries* fseries_for_ifft, COMPLEX16FFTPlan* revplan, COMPLEX16TimeSeries *tmp_tseries, COMPLEX16FFTPlan *fwdplan_for_tdom_wave, int shift_flag) {

	//gsl_vector_complex *z_t = gsl_vector_complex_calloc(length_max);
	//gsl_complex complex_norm, z_tmp_element;
	//double norm;

	generate_template(m1, m2, f_min, order, fseries, tmp_tseries, fwdplan_for_tdom_wave, shift_flag);
	//XLALWhitenCOMPLEX16FrequencySeries(fseries, psd);
	//add_quadrature_phase(fseries, fseries_for_ifft);

	int n = fseries->data->length;

	for ( unsigned int l = 0; l < psd->data->length; l++){

		fseries->data->data[fseries->data->length - 1 - (n/2 - 1) + l].re /= sqrt(psd->data->data[l]);
		fseries->data->data[fseries->data->length - 1 - (n/2 - 1) + l].im /= sqrt(psd->data->data[l]);

	}
	freq_to_time_fft(fseries, tseries, revplan);
        for(unsigned int l = 0 ; l < length_max; l++){

                gsl_vector_set(template_real, l, tseries->data->data[tseries->data->length - 1 - (length_max - 1) + l].re);
                gsl_vector_set(template_imag, l, tseries->data->data[tseries->data->length - 1 - (length_max- 1) + l].im);		
		//gsl_vector_complex_set(z_t, l, z_tmp_element);
	}
	//gsl_blas_zdotc(z_t, z_t, &complex_norm);

        //norm = 1./sqrt( gsl_complex_abs( complex_norm ) );

  	//gsl_vector_scale (template_real, norm);
        
  	//gsl_vector_scale (template_imag, norm);
        
        //gsl_vector_complex_free(z_t);
	fseries_for_ifft->data->data[0].re = 0;
	return 0;
} 

static gsl_matrix *create_templates_from_mc_and_eta(gsl_vector *mcvec, gsl_vector *etavec, double f_min, int length_max, REAL8FrequencySeries* psd, COMPLEX16TimeSeries* tseries, COMPLEX16FrequencySeries* fseries, COMPLEX16FrequencySeries* fseries_for_ifft, COMPLEX16FFTPlan *revplan){
       /*
 	* N_mc is number of points on M_c grid
 	* viceversa for M_eta
 	*/ 
	unsigned int i,j;
	unsigned int k =0;
	double sample_rate;
	double working_duration;
        int working_length;
	double eta, mc, m1, m2;
	LIGOTimeGPS epoch = LIGOTIMEGPSZERO;

	working_length = fseries_for_ifft->data->length;
	working_duration = working_length*tseries->deltaT;
	sample_rate = round(working_length / working_duration);

	gsl_vector *template_real = gsl_vector_calloc(length_max);
	gsl_vector *template_imag = gsl_vector_calloc(length_max);

	COMPLEX16TimeSeries *tmp_tseries; 
	COMPLEX16FFTPlan *fwdplan_for_tseries;

	tmp_tseries = XLALCreateCOMPLEX16TimeSeries(NULL, &epoch, 0, 1/sample_rate, &lalDimensionlessUnit, tseries->data->length);
        memset (tmp_tseries->data->data, 0, tmp_tseries->data->length * sizeof (COMPLEX16));
        fwdplan_for_tseries = XLALCreateForwardCOMPLEX16FFTPlan(tmp_tseries->data->length, 1);

	/* gsl_matrix *A will contain template bank */
	gsl_matrix *A = gsl_matrix_calloc(length_max, 2 * mcvec->size * etavec->size);

	for ( i = 0; i < mcvec->size ; i++){
		for ( j = 0; j < etavec->size ; j++){
			
			eta = gsl_vector_get(etavec, j);
			mc = gsl_vector_get(mcvec, i);

                        m1 = mc2mass1(mc, eta);
                        m2 = mc2mass2(mc, eta);

			generate_whitened_template(m1, m2, f_min, length_max, 4, psd, template_real, template_imag, tseries, fseries, fseries_for_ifft, revplan, tmp_tseries, fwdplan_for_tseries, 1);
	
			gsl_matrix_set_col(A, 2*k,  template_real);
			gsl_matrix_set_col(A, 2*k+1, template_imag);

			k+=1;
		}
	}
	
	gsl_vector_free(template_real);
	gsl_vector_free(template_imag);

	XLALDestroyCOMPLEX16TimeSeries(tmp_tseries);
	XLALDestroyCOMPLEX16FFTPlan(fwdplan_for_tseries);

	return A;
}

static gsl_matrix *create_svd_basis_from_template_bank(gsl_matrix* template_bank){
	
	double tolerance;
	double norm_s;
	double sum_s = 0.;
	unsigned int n;

	gsl_matrix *output;
	gsl_matrix_view template_view;
	gsl_matrix *V;
	gsl_vector *S; 

	/* Work space matrices */
	gsl_matrix *gX;
	gsl_vector *gW;

	gX = gsl_matrix_calloc(template_bank->size2, template_bank->size2);
 	gW = gsl_vector_calloc(template_bank->size2);	

	V = gsl_matrix_calloc(template_bank->size2, template_bank->size2);
	S = gsl_vector_calloc(template_bank->size2);

	
	gsl_linalg_SV_decomp_mod(template_bank, gX, V, S, gW);



	fprintf(stderr, "SVD completed.\n");


	tolerance = 0.999999;
	norm_s = pow(gsl_blas_dnrm2(S), 2.);

	if (isnan(norm_s) ) exit(1);

	gsl_matrix_free(gX);
        gsl_vector_free(gW);
        gsl_matrix_free(V);	

	/*FIXME make this more sophisticated if you care */

	sum_s = 0;
	for (n = 0; n < S->size; n++) {
		sum_s += gsl_vector_get(S, n) * gsl_vector_get(S, n);
		if (sqrt(sum_s / norm_s) >= tolerance) break;
		}

	if (n % 2)
		n += 1;


	fprintf(stderr,"SVD: using %d basis templates\n:", n);
	fprintf(stderr, "norm: %f", norm_s);	
	template_view = gsl_matrix_submatrix(template_bank, 0, 0, template_bank->size1, n);
	output = gsl_matrix_calloc(template_bank->size1, n);
	gsl_matrix_memcpy(output, &template_view.matrix);
	
	gsl_vector_free(S);

	return output;
}

/* FIXME use a better name */
int interpolate_waveform_from_mchirp_and_eta(struct twod_waveform_interpolant_array *interps, gsl_vector_complex *h_t, double mchirp, double eta) { 

	unsigned int i;
	gsl_complex M;
	double x, y;
	//struct twod_waveform_interpolant *interp = interps->interp;

	gsl_vector_view h_t_real = gsl_vector_complex_real(h_t); 
	gsl_vector_view h_t_imag = gsl_vector_complex_imag(h_t);


	x = map_coordinate_to_cheby(interps->param1_min, interps->param1_max, mchirp);
	y = map_coordinate_to_cheby(interps->param2_min, interps->param2_max, eta);
	/* this is the loop over mu */
	for (i = 0; i < interps->size; i++) {
		M = compute_M_xy(interps->interp[i].C_KL, x, y);
		gsl_blas_daxpy (GSL_REAL(M), interps->interp[i].svd_basis, &h_t_real.vector);
		gsl_blas_daxpy (GSL_IMAG(M), interps->interp[i].svd_basis, &h_t_imag.vector);
	}
	
	return 0;
	
}

static int pad_parameter_bounds(double mc_min, double mc_max, double eta_min, double eta_max, double *outer_eta_min, double *outer_eta_max, double *outer_mc_min, double *outer_mc_max, int number_templates_along_eta, int number_templates_along_mc, int number_of_templates_to_pad, double *mc_padding, double *eta_padding){


	/* pad parameter by number_of_templates_to_pad templates */
	*eta_padding = number_of_templates_to_pad*(eta_max-eta_min)/number_templates_along_eta;
	*mc_padding = number_of_templates_to_pad*(mc_max-mc_min)/number_templates_along_mc; 

	/* check that padded eta_min and eta_max are physical */
	if(0. <= eta_min - *eta_padding){
		
		*outer_eta_min = eta_min - *eta_padding;
	}
	
	else *outer_eta_min = eta_min;

	if (0.25 >= eta_max + *eta_padding){

		*outer_eta_max = eta_max + *eta_padding;
	}

	else *outer_eta_max = eta_max;
	
	/* check that padded mc_min is physical */
	if (0. <= mc_min - *mc_padding){	

		*outer_mc_min = mc_min - *mc_padding;

	}

	else *outer_mc_min = mc_min;

	*outer_mc_max = mc_max + *mc_padding;

	return 0;

}

static int make_patch_from_manifold(struct twod_waveform_interpolant_manifold *manifold){

	unsigned int patches_in_eta, patches_in_mc;
	unsigned int i,j,k;
	double pad_mc, pad_eta;
	double mc_min, mc_max, eta_min, eta_max;
	double eta_width, mc_width;

	patches_in_mc = manifold->patches_in_mc;
	patches_in_eta = manifold->patches_in_eta;

	/* find index (j,k) of the required patch on the manifold from i */
	
	for(i=0; i < patches_in_mc*patches_in_eta; i++){

	
		k = i % patches_in_eta;
	
		if ( !(i % patches_in_eta) ){
	
			j= i/patches_in_eta;
		}	
	
		
		pad_mc = manifold->mc_padding;
		pad_eta = manifold->eta_padding;

		mc_width =  (manifold->inner_param1_max - manifold->inner_param1_min)/patches_in_mc;
		eta_width =  (manifold->inner_param2_max - manifold->inner_param2_min)/patches_in_eta;

		if ( 0. <= manifold->inner_param1_min + j*mc_width - pad_mc){

			mc_min = manifold->inner_param1_min + j*mc_width - pad_mc;	
			manifold->interp_arrays[i].param1_min = mc_min;		
			manifold->interp_arrays[i].inner_param1_min = mc_min + pad_mc;
		}

		else{
			mc_min = manifold->inner_param1_min + j*mc_width;
			manifold->interp_arrays[i].param1_min = mc_min;
			manifold->interp_arrays[i].inner_param1_min = mc_min; 
		}


		mc_max = manifold->inner_param1_min + mc_width + j*mc_width + pad_mc;		
		manifold->interp_arrays[i].param1_max = mc_max;
		manifold->interp_arrays[i].inner_param1_max = mc_max - pad_mc;

		if (0. <= manifold->inner_param2_min + j*eta_width - pad_eta){

			eta_min = manifold->inner_param2_min + k*eta_width - pad_eta;
			manifold->interp_arrays[i].param2_min = eta_min;
			manifold->interp_arrays[i].inner_param2_min = eta_min + pad_eta;
		}	

		else{
		
			eta_min = manifold->inner_param2_min + k*eta_width;
			manifold->interp_arrays[i].param2_min = eta_min;
			manifold->interp_arrays[i].inner_param2_min = eta_min;	
	
		}
		if (0.25 >=  manifold->inner_param2_min + eta_width + k*eta_width + pad_eta){

			eta_max = manifold->inner_param2_min + eta_width + k*eta_width + pad_eta;	
			manifold->interp_arrays[i].param2_max = eta_max;
			manifold->interp_arrays[i].inner_param2_max = eta_max - pad_eta;
		}		
		else{

			eta_max = manifold->inner_param2_min + eta_width + k*eta_width;
			manifold->interp_arrays[i].param2_max = eta_max;
			manifold->interp_arrays[i].inner_param2_max = eta_max;

		}

	}



	return 0;
}

static int populate_interpolants_on_patches(struct twod_waveform_interpolant_manifold *manifold, COMPLEX16FrequencySeries *fseries, COMPLEX16FrequencySeries *fseries_for_ifft, COMPLEX16TimeSeries *tseries, COMPLEX16FFTPlan *revplan, int length_max, double f_min){

	int i;
	gsl_vector *x_nodes;
	gsl_vector *y_nodes;
	gsl_vector *mchirps_even;
	gsl_vector *etas_even;
	gsl_vector *mchirps_nodes;
	gsl_vector *etas_nodes;
	gsl_matrix_complex *M_xy;
	gsl_matrix *templates;
	gsl_matrix *svd_basis;
	gsl_matrix *templates_at_nodes;
	gsl_matrix_complex *phase_M0_xy;
	unsigned int N_mc, M_eta;
	int number_of_patches;
	double mc_min, mc_max, eta_min, eta_max;
	
	number_of_patches = manifold->patches_in_mc * manifold->patches_in_eta;
	
	/* loop over each patch */
	
	

	for(i = 0; i < number_of_patches; i++){

		N_mc = manifold->number_templates_along_mc;
		M_eta = manifold->number_templates_along_eta;
		mc_min = manifold->interp_arrays[i].param1_min;
		mc_max = manifold->interp_arrays[i].param1_max;
		eta_min = manifold->interp_arrays[i].param2_min;
		eta_max = manifold->interp_arrays[i].param2_max;	

		mchirps_even = even_param_spacing(mc_min, mc_max, N_mc);
		etas_even = even_param_spacing(eta_min, eta_max, M_eta);

		templates = create_templates_from_mc_and_eta(mchirps_even, etas_even, f_min, length_max, manifold->psd, tseries, fseries, fseries_for_ifft, revplan);

		gsl_vector_free(etas_even);
		gsl_vector_free(mchirps_even);
	
		svd_basis = create_svd_basis_from_template_bank(templates);	
	
		manifold->interp_arrays[i].size = svd_basis->size2;	
		manifold->interp_arrays[i].interp = new_waveform_interpolant_from_svd_bank(svd_basis);	


		/* Compute new template bank at colocation points-> project onto basis vectors
 		 * to get matrix of coefficients for C_KL computation */
	
	        x_nodes = raw_nodes(N_mc);
	        y_nodes = raw_nodes(M_eta);
	
		mchirps_nodes = node_param_spacing(mc_min, mc_max, x_nodes);
		etas_nodes = node_param_spacing(eta_min, eta_max, y_nodes);
	
		templates_at_nodes = create_templates_from_mc_and_eta(mchirps_nodes, etas_nodes, f_min, length_max, manifold->psd, tseries, fseries, fseries_for_ifft, revplan);
		
		phase_M0_xy = gsl_matrix_complex_calloc(mchirps_nodes->size, etas_nodes->size);
		M_xy = gsl_matrix_complex_calloc(mchirps_nodes->size, etas_nodes->size);
	
		for (unsigned int j = 0; j < manifold->interp_arrays[i].size; j++) {       
	
			projection_coefficient(manifold->interp_arrays[i].interp[j].svd_basis, templates_at_nodes, M_xy, N_mc, M_eta);
	
			if(j==0){
	
				measure_m0_phase(M_xy, phase_M0_xy);
			}
	
			rotate_M_xy(M_xy, phase_M0_xy);
	
			manifold->interp_arrays[i].interp[j].C_KL = compute_C_KL(x_nodes, y_nodes, M_xy);           
	
			
		}

		gsl_vector_free(etas_nodes);
                gsl_vector_free(mchirps_nodes);
                gsl_vector_free(x_nodes);
                gsl_vector_free(y_nodes);
                gsl_matrix_complex_free(M_xy);
                gsl_matrix_free(templates_at_nodes);
                gsl_matrix_complex_free(phase_M0_xy);
                gsl_matrix_free(templates);	
		gsl_matrix_free(svd_basis);
	
	} 
	return 0;

}

int index_into_patch(struct twod_waveform_interpolant_manifold *manifold, double p1, double p2) {
        unsigned int i;
        //struct twod_waveform_interpolant_array *out = manifold->interp_arrays;
        for (i = 0; i < manifold->patches_in_eta*manifold->patches_in_mc; i++){
                if ((p1 >= manifold->interp_arrays[i].inner_param1_min) && (p2 >= manifold->interp_arrays[i].inner_param2_min) && (p1 <= manifold->interp_arrays[i].inner_param1_max) && (p2 <= manifold->interp_arrays[i].inner_param2_max)){
		
			 break;
		}
    	}

        return i;
}






struct twod_waveform_interpolant_manifold *XLALInferenceCreateInterpManifold(double mc_min, double mc_max, double eta_min, double eta_max, double f_min){

	/* Hard code for now. FIXME: figure out way to optimize and automate patching, given parameter bounds */

        unsigned int patches_in_eta = 2;
        unsigned int patches_in_mc = 2;
        unsigned int number_templates_along_eta = 15;
        unsigned int number_templates_along_mc = 15;
	unsigned int number_of_templates_to_pad = 1;
        unsigned int number_of_patches;

        unsigned int length_max = 0;

	double outer_eta_min, outer_eta_max, outer_mc_min, outer_mc_max, mc_padding, eta_padding;

	struct twod_waveform_interpolant_manifold *manifold = NULL;
	
        COMPLEX16FrequencySeries *fseries = NULL;
	REAL8FrequencySeries *psd_for_template_bank = NULL;
	COMPLEX16FrequencySeries *fseries_for_ifft = NULL;
        COMPLEX16TimeSeries *tseries = NULL;
	COMPLEX16FFTPlan *revplan = NULL;

	pad_parameter_bounds(mc_min, mc_max, eta_min, eta_max, &outer_eta_min, &outer_eta_max, &outer_mc_min, &outer_mc_max, number_templates_along_mc, number_templates_along_eta, number_of_templates_to_pad, &mc_padding, &eta_padding);

	/* FIXME: initialize_time_and_freq_series needs to take in psd_to_interpolate and return psd for use in template bank */	
	initialize_time_and_freq_series(&psd_for_template_bank, &fseries, &fseries_for_ifft, &tseries, &revplan, outer_mc_min, outer_mc_max, outer_eta_min, outer_eta_max, f_min, &length_max, 0);

	manifold = interpolants_manifold_init(psd_for_template_bank, patches_in_eta, patches_in_mc, number_templates_along_mc, number_templates_along_eta, mc_min, mc_max, eta_min, eta_max, outer_mc_min, outer_mc_max, outer_eta_min, outer_eta_max, mc_padding, eta_padding);


	number_of_patches = patches_in_eta*patches_in_mc;

	for(unsigned int i=0; i < number_of_patches; i++){
		
		make_patch_from_manifold(manifold);

	}	

	populate_interpolants_on_patches(manifold, fseries, fseries_for_ifft, tseries, revplan, length_max, f_min);

	XLALDestroyCOMPLEX16TimeSeries(tseries);
        XLALDestroyCOMPLEX16FrequencySeries(fseries_for_ifft);
        XLALDestroyCOMPLEX16FrequencySeries(fseries);
        XLALDestroyCOMPLEX16FFTPlan(revplan);

	return manifold;

}


/*static int dewhiten_template_wave(gsl_vector_complex* template, COMPLEX16TimeSeries *dewhitened_tseries, COMPLEX16FrequencySeries *dewhitened_fseries, COMPLEX16FrequencySeries *fseries_for_dewhitening, COMPLEX16FFTPlan *fwdplan_for_dewhitening, REAL8FrequencySeries* psd, double f_min, double f_isco){

	unsigned int k, l;
	double deltaF;
	deltaF = fseries_for_dewhitening->deltaF;	

               	for (k = 0; k < template->size; k++){
			dewhitened_tseries->data->data[dewhitened_tseries->data->length - 1 - (template->size -1) + k].re = GSL_REAL(gsl_vector_complex_get(template, k));
			dewhitened_tseries->data->data[dewhitened_tseries->data->length - 1 - (template->size -1) + k].im = GSL_IMAG(gsl_vector_complex_get(template, k));
               	}
		
		XLALCOMPLEX16TimeFreqFFT(fseries_for_dewhitening, dewhitened_tseries, fwdplan_for_dewhitening);			


		for (l = 0; l < dewhitened_fseries->data->length; l++){

				if(l < ceil(f_min/deltaF) || l > f_isco/deltaF){
					dewhitened_fseries->data->data[l].re = 0;
					dewhitened_fseries->data->data[l].im = 0;
				

				}

				else{
					dewhitened_fseries->data->data[l].re = fseries_for_dewhitening->data->data[l + dewhitened_tseries->data->length/2 - 1 ].re*sqrt(psd->data->data[ l ]/(2.*deltaF)) ;
					dewhitened_fseries->data->data[l].im = fseries_for_dewhitening->data->data[l + dewhitened_tseries->data->length/2 - 1].im*sqrt(psd->data->data[ l ]/(2.*deltaF)) ;
				}

		   }			



	return 0;

}
*/

static int compute_overlap_whitened_waveform(struct twod_waveform_interpolant_manifold *manifold, COMPLEX16FrequencySeries *fseries_for_ifft, COMPLEX16TimeSeries *tseries,  double length_max, unsigned int New_N_mc, unsigned int New_M_eta, double f_min){

	FILE *list_of_overlaps;
	double Overlap=0.;
	double m1, m2, eta, mc, mc_min, mc_max, eta_min, eta_max;
	int working_length, patch_index;
	double working_duration, sample_rate;

	gsl_complex dotc1, dotc2, dotc3, complex_zero, z_tmp_element;

	gsl_vector *template_real;
	gsl_vector *template_imag;
	gsl_vector *mchirps_interps;
	gsl_vector *etas_interps;
	gsl_vector_complex *z_tmp;
	gsl_vector_complex *h_t;
	
	working_length = fseries_for_ifft->data->length;
	working_duration = working_length*tseries->deltaT;
	sample_rate = round(working_length / working_duration);
	
	z_tmp = gsl_vector_complex_calloc(length_max);

	mc_min = manifold->inner_param1_min;
	mc_max = manifold->inner_param1_max;
	eta_min = manifold->inner_param2_min;
	eta_max = manifold->inner_param2_max;

	mchirps_interps = linear_space_for_interpolation(mc_min, mc_max, New_N_mc);
        etas_interps = linear_space_for_interpolation(eta_min, eta_max, New_M_eta);

	GSL_SET_COMPLEX(&dotc1, 0, 0);
	GSL_SET_COMPLEX(&dotc2, 0 ,0);
	GSL_SET_COMPLEX(&dotc3, 0 ,0);
	GSL_SET_COMPLEX(&complex_zero, 0 ,0);

	template_real = gsl_vector_calloc(length_max);
	template_imag = gsl_vector_calloc(length_max);	
	h_t = gsl_vector_complex_calloc(length_max);
	LIGOTimeGPS epoch = LIGOTIMEGPSZERO;
	
	COMPLEX16TimeSeries *tmp_tseries;	
	COMPLEX16FFTPlan *fwdplan_for_tseries;
	COMPLEX16FrequencySeries *htilde;
	COMPLEX16FFTPlan *revplan;

	htilde = XLALCreateCOMPLEX16FrequencySeries(NULL, &epoch, 0, fseries_for_ifft->deltaF, &lalDimensionlessUnit, working_length);
        memset (htilde->data->data, 0, htilde->data->length * sizeof (COMPLEX16));	

        tmp_tseries = XLALCreateCOMPLEX16TimeSeries(NULL, &epoch, 0, 1/sample_rate, &lalDimensionlessUnit, tseries->data->length);
	memset (tmp_tseries->data->data, 0, tmp_tseries->data->length * sizeof (COMPLEX16));

	fwdplan_for_tseries = XLALCreateForwardCOMPLEX16FFTPlan(tmp_tseries->data->length, 1);

	revplan = XLALCreateReverseCOMPLEX16FFTPlan(working_length, 1);

	list_of_overlaps = fopen("overlaps_whitened.txt","w");		

	for (unsigned int i =0; i <  mchirps_interps->size; i++){
		for (unsigned int j =0; j <  etas_interps->size; j++){
			gsl_vector_complex_set_all(h_t, complex_zero);
			gsl_vector_complex_set_all(z_tmp, complex_zero);

                        eta = gsl_vector_get(etas_interps, j);
                        mc = gsl_vector_get(mchirps_interps, i);						
			patch_index = index_into_patch(manifold, mc, eta);
                        interpolate_waveform_from_mchirp_and_eta(&manifold->interp_arrays[patch_index], h_t, mc, eta);
			m1 = mc2mass1(mc, eta);
                        m2 = mc2mass2(mc, eta);
			generate_whitened_template(m1, m2, f_min, length_max, 4, manifold->psd, template_real, template_imag, tseries, htilde, fseries_for_ifft, revplan, tmp_tseries, fwdplan_for_tseries, 1);
			for(unsigned int l = 0; l < length_max; l++){
				GSL_SET_COMPLEX(&z_tmp_element, gsl_vector_get(template_real, l), gsl_vector_get(template_imag, l) );
				gsl_vector_complex_set(z_tmp, l, z_tmp_element);

			}

			gsl_blas_zdotc(z_tmp, h_t, &dotc1);
			gsl_blas_zdotc(h_t, h_t, &dotc2);
			gsl_blas_zdotc(z_tmp, z_tmp, &dotc3);	
			
		 	Overlap = ( gsl_complex_abs( dotc1 ) / sqrt( gsl_complex_abs( dotc2 ) ) / sqrt( gsl_complex_abs( dotc3 ) ) );
	
			fprintf(list_of_overlaps,"mc = %f, eta=%f, overlap=%e\n", mc, eta, Overlap);
			fprintf(stderr, "mc = %f, eta=%f, overlap=%e\n", mc, eta, Overlap);	
			GSL_SET_COMPLEX(&dotc1, 0, 0);
		        GSL_SET_COMPLEX(&dotc2, 0 ,0);
		        GSL_SET_COMPLEX(&dotc3, 0 ,0);

			//XLALDestroyCOMPLEX16FrequencySeries(htilde);	

			for(unsigned int m=0; m < h_t->size; m++){

				gsl_vector_complex_set(h_t, m, dotc1);
			}	


			}
		}		

	fclose(list_of_overlaps);

	gsl_vector_complex_free(h_t);
	gsl_vector_complex_free(z_tmp);
	gsl_vector_free(template_imag);
	gsl_vector_free(template_real);
	gsl_vector_free(mchirps_interps);
	gsl_vector_free(etas_interps);

	return 0;
}

int main(void){


	unsigned int length_max = 0;
	unsigned int New_N_mc = 50, New_M_eta = 50;
	double mc_min = 7.2;
	double eta_min = 0.175;
	double mc_max = 7.6;
	double eta_max = 0.25;
	double f_min = 40.0;

	struct twod_waveform_interpolant_manifold *manifold = NULL;

	REAL8FrequencySeries *psd = NULL;
        COMPLEX16FrequencySeries *fseries = NULL;
	COMPLEX16FrequencySeries *fseries_for_ifft = NULL;
        COMPLEX16TimeSeries *tseries = NULL;
	COMPLEX16FFTPlan *revplan = NULL;

	manifold = XLALInferenceCreateInterpManifold(mc_min, mc_max, eta_min, eta_max, f_min);

	initialize_time_and_freq_series(&psd, &fseries, &fseries_for_ifft, &tseries, &revplan, mc_min, mc_max, eta_min, eta_max, f_min, &length_max, 0);

	compute_overlap_whitened_waveform(manifold, fseries_for_ifft, tseries, length_max, New_N_mc, New_M_eta, f_min);

	XLALDestroyCOMPLEX16TimeSeries(tseries);
        XLALDestroyCOMPLEX16FrequencySeries(fseries_for_ifft);
        XLALDestroyCOMPLEX16FrequencySeries(fseries);
        XLALDestroyCOMPLEX16FFTPlan(revplan);
        XLALDestroyREAL8FrequencySeries(psd);
	
	XLALInferenceDestroyInterpManifold(manifold);	

	return 0;

}
