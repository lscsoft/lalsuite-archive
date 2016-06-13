%option noyywrap

WHITESPACE	[ \t]+
DIGIT	[0-9]
NUMBER	{DIGIT}+
REAL	[0-9eE+\-.]+

%{

#include "make_sft_global.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

extern INT64 gps_start, gps_end;
extern long samples_per_second;

int mode;

extern int output_mode;
extern int precision;
extern long freq_start, freq_stop;
extern long overlap;
extern char *output_calibrated_data;
extern char *output_power;
extern char *output_sft;
extern char bypass_calibration, bypass_freq_window, bypass_lines, 
	bypass_highpass_filter, bypass_lowpass_filter, bypass_first_window, bypass_first_fft;
	
extern REAL4 highpass_filter_f, highpass_filter_a, lowpass_filter_f, lowpass_filter_a;
extern int highpass_filter_order, lowpass_filter_order;

extern MAKE_SFT_WINDOW window1, freq_window, window2;

extern char trace_power, trace_calibration_power;

extern char *channel;

#define MODE_NONE			0
#define MODE_FRAME_FILES		1
#define MODE_SEGMENTS			2
#define MODE_SUBSECOND_SEGMENTS 	3
#define MODE_SAMPLE_SEGMENTS		4
#define MODE_START_STOP_FRAME_FILES	5

double conf_a, conf_b, conf_c;
long conf_i;
unsigned char *conf_p=NULL;

int dont_fail_on_missing_data=0, dont_fail_on_missing_calibration=0, dont_fail_on_zero_ds=0;


%}

%%

#.*\n		{ /* ignore comments */ }

FRAME_FILES{WHITESPACE}?\n	{ mode=MODE_FRAME_FILES; }

START_STOP_FRAME_FILES{WHITESPACE}?\n	{ mode=MODE_START_STOP_FRAME_FILES; }

SEGMENTS{WHITESPACE}?\n	{ mode=MODE_SEGMENTS; }

SUBSECOND_SEGMENTS{WHITESPACE}?\n	{ mode=MODE_SUBSECOND_SEGMENTS; }

SAMPLE_SEGMENTS{WHITESPACE}?\n		{ mode=MODE_SAMPLE_SEGMENTS; }

OUTPUT_MODE_TEXT{WHITESPACE}?\n	{
	output_mode=OUTPUT_MODE_TEXT;
	}

OUTPUT_MODE_BINARY{WHITESPACE}?\n	{
	output_mode=OUTPUT_MODE_BINARY;
	}
	
OUTPUT_MODE_GEO{WHITESPACE}?\n	{
	output_mode=OUTPUT_MODE_GEO;
	}

CHANNEL{WHITESPACE}\".*\"{WHITESPACE}?\n	{
	if(channel!=NULL)free(channel);
	channel=dup_quoted_name(yytext);
	}

DURATION{WHITESPACE}{NUMBER}{WHITESPACE}{NUMBER}{WHITESPACE}?\n	{
	sscanf(yytext+8,"%lld %lld", &gps_start, &gps_end);
	allocate_data();
	}
	
HIGHPASS_FILTER{WHITESPACE}{NUMBER}{WHITESPACE}{REAL}{WHITESPACE}{REAL}{WHITESPACE}?\n	{
	sscanf(yytext+15,"%d %g %g", &highpass_filter_order, &highpass_filter_f, &highpass_filter_a);
	}
	
LOWPASS_FILTER{WHITESPACE}{NUMBER}{WHITESPACE}{REAL}{WHITESPACE}{REAL}{WHITESPACE}?\n	{
	sscanf(yytext+14,"%d %g %g", &lowpass_filter_order, &lowpass_filter_f, &lowpass_filter_a);
	}

OVERLAP{WHITESPACE}{NUMBER}{WHITESPACE}?\n	{
	sscanf(yytext+7, "%ld", &overlap);
	}
	
FREQ_WINDOW_TAIL_SIZE{WHITESPACE}{NUMBER}{WHITESPACE}?\n	{
	sscanf(yytext+21, "%ld", &conf_i);
	if(conf_i<0){
		fprintf(stderr, "** Window tail size must be non-negative\n");
		exit(-1);
		}
	if(samples_per_second<0){
		fprintf(stderr, "** FREQ_WINDOW_TAIL_SIZE must follow valid SAMPLES_PER_SECOND command\n");
		exit(-1);
		}
	if(conf_i>30*samples_per_second){
		fprintf(stderr, "** Window tail size is too big\n");
		exit(-1);
		}
	freq_window.tail_size=conf_i;	
	}

WINDOW(1|2)_TAIL_SIZE{WHITESPACE}{NUMBER}{WHITESPACE}?\n	{
	sscanf(yytext+17, "%ld", &conf_i);
	if(conf_i<0){
		fprintf(stderr, "** Window %c tail size must be non-negative\n", yytext[6]);
		exit(-1);
		}
	if(samples_per_second<0){
		fprintf(stderr, "** WINDOW%c_TAIL_SIZE must follow valid SAMPLES_PER_SECOND command\n", yytext[6]);
		exit(-1);
		}
	if(conf_i>30*samples_per_second){
		fprintf(stderr, "** Window tail size is too big\n");
		exit(-1);
		}
	switch(yytext[6]){
		case '1':
			window1.tail_size=conf_i;	
			break;
		case '2':
			window2.tail_size=conf_i;	
			break;
		default:
			fprintf(stderr, "** INTERNAL ERROR - window tail size specification (%c)\n", yytext[6]);
			exit(-1);
		}
	}

WINDOW(1|2)_TYPE{WHITESPACE}\".*\"{WHITESPACE}?\n	{
	switch(yytext[6]){
		case '1':
			window1.type=translate_window_type(yytext+12);	
			break;
		case '2':
			window2.type=translate_window_type(yytext+12);	
			break;
		default:
			fprintf(stderr, "** INTERNAL ERROR - window type specification (%c)\n", yytext[6]);
			exit(-1);
		}
	}

FREQ_WINDOW_TYPE{WHITESPACE}\".*\"{WHITESPACE}?\n	{
	freq_window.type=translate_window_type(yytext+16);
	}


SAMPLES_PER_SECOND{WHITESPACE}{NUMBER}\n	{
	sscanf(yytext+18, "%ld", &samples_per_second);
	allocate_data();
	}

DEBUG{NUMBER}{WHITESPACE}\".*\"{WHITESPACE}?\n	{
	set_debug_output(yytext);
	}

DEBUG_WINDOW(1|2){WHITESPACE}\".*\"{WHITESPACE}?\n	{
	{ MAKE_SFT_WINDOW *w;
		switch(yytext[12]){
			case '1':
				w=&window1;
				break;
			case '2':
				w=&window2;
				break;
			default:
				fprintf(stderr, "** INTERNAL ERROR - window debug specification (%c)\n", yytext[12]);
				exit(-1);
			}
		if(w->debug!=NULL)free(w->debug);
		w->debug=dup_quoted_name(yytext+13);
		}
	}

DEBUG_FREQ_WINDOW{WHITESPACE}\".*\"{WHITESPACE}?\n	{
	if(freq_window.debug!=NULL)free(freq_window.debug);
	freq_window.debug=dup_quoted_name(yytext+17);
	}

OUTPUT_CALIBRATED_DATA{WHITESPACE}\".*\"{WHITESPACE}?\n	{
	if(output_calibrated_data!=NULL)free(output_calibrated_data);
	output_calibrated_data=dup_quoted_name(yytext+23);
	}

OUTPUT_POWER{WHITESPACE}\".*\"{WHITESPACE}?\n	{
	if(output_power!=NULL)free(output_power);
	output_power=dup_quoted_name(yytext+12);
	}
OUTPUT_SFT{WHITESPACE}\".*\"{WHITESPACE}?\n	{
	if(output_sft!=NULL)free(output_sft);
	output_sft=dup_quoted_name(yytext+10);
	}

OUTPUT_BAND{WHITESPACE}{REAL}{WHITESPACE}{REAL}{WHITESPACE}?\n	{
	if((gps_start<0)||(gps_end<0)){
		fprintf(stderr,"** OUTPUT_BAND command must follow DURATION command\n");
		exit(-1);
		}
	sscanf(yytext+11,"%lg %lg", &conf_a, &conf_b);
	freq_start=lrint(conf_a*(gps_end-gps_start));
	freq_stop=lrint(conf_b*(gps_end-gps_start));
	}
	
OUTPUT_BAND_BINS{WHITESPACE}{NUMBER}{WHITESPACE}{NUMBER}{WHITESPACE}?\n	{
	sscanf(yytext+16, "%ld %ld", &freq_start, &freq_stop);
	}

PRECISION{WHITESPACE}{NUMBER}{WHITESPACE}?\n	{
	sscanf(yytext+9, "%d", &precision);
	}
	
ALPHA_BETA_FILE{WHITESPACE}\".*\"{WHITESPACE}?\n	{
	conf_p=dup_quoted_name(yytext+15);
	load_alpha_beta(conf_p);
	free(conf_p);
	conf_p=NULL;
	}

R_FILE{WHITESPACE}\".*\"{WHITESPACE}?\n	{
	conf_p=dup_quoted_name(yytext+6);
	load_R(conf_p);
	free(conf_p);
	conf_p=NULL;
	}

C_FILE{WHITESPACE}\".*\"{WHITESPACE}?\n	{
	conf_p=dup_quoted_name(yytext+6);
	load_C(conf_p);
	free(conf_p);
	conf_p=NULL;
	}

SUPPRESS{WHITESPACE}{REAL}{WHITESPACE}{REAL}{WHITESPACE}{NUMBER}{WHITESPACE}?\n	{
	sscanf(yytext+8, "%lg %lg %ld", &conf_a, &conf_b, &conf_i);
	add_suppress_region(conf_a, conf_b, conf_i);
	}

DONT_FAIL_ON_MISSING_DATA{WHITESPACE}?\n	{
	dont_fail_on_missing_data=1;
	}	
	
DONT_FAIL_ON_MISSING_CALIBRATION{WHITESPACE}?\n	{
	dont_fail_on_missing_calibration=1;
	}

DONT_FAIL_ON_ZERO_DS{WHITESPACE}?\n	{
	dont_fail_on_zero_ds=1;
	}

BYPASS_FREQ_WINDOW{WHITESPACE}?\n	{
	bypass_freq_window=1;
	}
	
BYPASS_CALIBRATION{WHITESPACE}?\n	{
	bypass_calibration=1;
	}
	
BYPASS_LINE_SUPPRESSION{WHITESPACE}?\n	{
	bypass_lines=1;
	}
	
BYPASS_HIGHPASS_FILTER{WHITESPACE}?\n	{
	bypass_highpass_filter=1;
	}
	
BYPASS_LOWPASS_FILTER{WHITESPACE}?\n	{
	bypass_lowpass_filter=1;
	}

BYPASS_FIRST_WINDOW{WHITESPACE}?\n	{
	bypass_first_window=1;
	}
	
BYPASS_FIRST_FFT{WHITESPACE}?\n	{
	bypass_first_fft=1;
	bypass_calibration=1;
	bypass_lines=1;
	bypass_freq_window=1;
	}

FAKE_DATA_CONSTANT{WHITESPACE}{REAL}{WHITESPACE}?\n	{
	sscanf(yytext+18,"%lg", &conf_a);
	add_fake_data_command(FAKE_DATA_CONSTANT, 0.0, 0.0, conf_a);
	}

FAKE_DATA_CONSTANT_FILL{WHITESPACE}{REAL}{WHITESPACE}?\n	{
	sscanf(yytext+23,"%lg", &conf_a);
	add_fake_data_command(FAKE_DATA_CONSTANT_FILL, 0.0, 0.0, conf_a);
	}

FAKE_DATA_RANDOM_GAUSSIAN{WHITESPACE}{REAL}{WHITESPACE}?\n	{
	sscanf(yytext+18,"%lg", &conf_a);
	add_fake_data_command(FAKE_DATA_RANDOM_GAUSSIAN, 0.0, 0.0, conf_a);
	}

FAKE_DATA_SINE_WAVE{WHITESPACE}{REAL}{WHITESPACE}{REAL}{WHITESPACE}{REAL}{WHITESPACE}?\n	{
	sscanf(yytext+19,"%lg %lg %lg", &conf_a, &conf_b, &conf_c);
	add_fake_data_command(FAKE_DATA_SINE_WAVE, (conf_a*M_PI)/180.0, conf_b, conf_c);
	}

FAKE_DATA_CALIBRATED_SINE_WAVE{WHITESPACE}{REAL}{WHITESPACE}{REAL}{WHITESPACE}{REAL}{WHITESPACE}?\n	{
	sscanf(yytext+30,"%lg %lg %lg", &conf_a, &conf_b, &conf_c);
	add_fake_data_command(FAKE_DATA_SINE_WAVE, (conf_a*M_PI)/180.0, conf_b, conf_c);
	}

TRACE_POWER{WHITESPACE}?\n	{
	trace_power=1;
	}

TRACE_CALIBRATION_POWER{WHITESPACE}?\n	{
	trace_calibration_power=1;
	}

.*	{
	switch(mode){
		case MODE_NONE:
			fprintf(stderr,"Could not parse line \"%s\"\n", yytext);
			break;
		case MODE_FRAME_FILES:
			add_frame_file(yytext, -1, -1);
			break;
		case MODE_START_STOP_FRAME_FILES:
			add_start_stop_frame_file(yytext);
			break;
		case MODE_SEGMENTS:
			assimilate_segment(yytext);
			break;
		case MODE_SUBSECOND_SEGMENTS:
			assimilate_subsecond_segment(yytext);
			break;
		case MODE_SAMPLE_SEGMENTS:
			assimilate_sample_segment(yytext);
			break;
		default:
			fprintf(stderr,"Unimplemented mode %d while parsing line \"%s\"\n", mode, yytext);
		}			
	}


{WHITESPACE}?\n	{
	/* skip whitespace and new lines */
	}
%%
