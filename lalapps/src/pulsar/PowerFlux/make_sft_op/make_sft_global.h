#ifndef __MAKE_SFT_GLOBAL_H__
#define __MAKE_SFT_GLOBAL_H__

#include <lal/LALStdlib.h>

#define MAKE_SFT_VERSION	"make_sft 0.1.1"

/* 64 bit integer datatype */

#ifndef INT64
#define INT64 long long
#endif


#define OUTPUT_MODE_TEXT 	0
#define OUTPUT_MODE_BINARY 	1
#define OUTPUT_MODE_GEO		2

#define DONT_FAKE_DATA			0
#define FAKE_DATA_CONSTANT		1
#define FAKE_DATA_RANDOM_GAUSSIAN	2
#define FAKE_DATA_SINE_WAVE		3
#define FAKE_DATA_CONSTANT_FILL		4
#define FAKE_DATA_CALIBRATED_SINE_WAVE	5

typedef struct {
	int command;
	double phase;
	double frequency;
	double amplitude;
	} FAKE_DATA_COMMAND;

typedef struct {
	int type;
	long tail_size;
	char *debug;
	} MAKE_SFT_WINDOW;


/* make_sft.c */

void set_channel(char *p);
void set_debug_output(char *p);

void add_suppress_region(double line1, double line2, int N);

void print_settings(FILE *f);

/* convenience routines */
int translate_window_type(char *text);


/* make_sft_conf.c */
int yylex(void);

/* make_sft_calibration.c */
void load_R(char *name);
void load_C(char *name);
 /* call this before using H or C data - i.e. before using calibrate_sft() */
void post_init_response_files(void);


void add_alpha_beta(INT64 gps, REAL4 alpha, REAL4 alphabeta);
void load_alpha_beta(char *name);
/* call this before using alpha/beta constants - i.e. before using calibrate_sft() -
   this sorts the array by gps time, among all other things */
void post_init_alpha_beta(void);

void calibrate_fft(COMPLEX8Vector *fft, long first_sample);

#include "op_method.h"

void compute_Phi(PHI_DATA3 *phi_data);
void compute_phi_r(PHI_RESPONSE3 *phi_response, long n);
void compute_estimation_errors(PHI_DATA3 *phi, PHI_RESPONSE3 *phi_r, long n, double *max1, double *max3);

void get_td_calibration(double frequency, long sample, double *re, double *im);

/* make_sft_utils.c */

void fault_file(char *filename);
char *dup_quoted_name(char *p);

void *do_alloc(long a, long b);

REAL4 sum_r4_squares(REAL4 *array, long count);

REAL4 sum_positive_r4(REAL4 *array, long count);

long count_bitmap_bits(unsigned char *bitmap, long length);

void fault_file(char *filename);

void print_REAL4Vector(REAL4Vector *vec, char *file, char *tag, int mode, long bin_start, long bin_stop);
void print_COMPLEX8Vector(COMPLEX8Vector *vec, char *file, char *tag, int mode, long bin_start, long bin_stop);

void print_window_settings(FILE *f, char *window_name, MAKE_SFT_WINDOW *window);
int translate_window_type(char *text);

/* make_sft_td_data */

void add_frame_file(char *filename, double TStart, double TEnd);
void add_start_stop_frame_file(char *text);

void allocate_data(void);

void add_fake_data_command(int command, double phase, double frequency, double amplitude);

void assimilate_segment(char *line);
void assimilate_subsecond_segment(char *line);
void assimilate_sample_segment(char *line);
void print_data_stats(void);
void verify_loaded_data(void);
void generate_fake_data(void);
void print_data(char *file, char *tag);
void print_fake_data_settings(FILE *f);
void linear_interpolate_gaps(void);


#define TESTSTATUS( status ) \
  { if ( (status)->statusCode ) { \
  	fprintf(stderr,"** LAL status encountered in file \"%s\" function \"%s\" line %d\n", __FILE__, __FUNCTION__, __LINE__);\
  	REPORTSTATUS((status)); exit(-1); \
	}  }

#endif
