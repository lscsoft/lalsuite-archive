
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <string.h>

#define __USE_ISOC99 1
#include <math.h>

#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/RealFFT.h>
#include <lal/Window.h>
#include <lal/BandPassTimeSeries.h>

#include "make_sft_global.h"

#include <FrameL.h>

extern char * channel;
extern INT64 gps_start, gps_end;
extern long samples_per_second;
extern int precision;
extern long total_samples;
extern int precision;
extern float *data;

long freq_start=-1, freq_stop=-1;

int output_mode=OUTPUT_MODE_TEXT;


/* Overlap specifies by how many samples nearby 60 second segments should overlap */
long overlap=0;

/* filenames to write debug data. NULL if no debugging has been requested */
char *file_debug1=NULL, *file_debug2=NULL;

/* options to bypass parts of processing */
char bypass_calibration=0, bypass_freq_window=0, bypass_lines=0, 
	bypass_highpass_filter=0, bypass_lowpass_filter=0, bypass_first_window=0, bypass_first_fft=0;


/* should we trace total power during computation ? */
char trace_power=0, trace_calibration_power=0;

/* filename to output results */
char *output_power=NULL, *output_calibrated_data=NULL, *output_sft=NULL;

/* Constants for Butterworth filter */
REAL4 highpass_filter_f=100.0, highpass_filter_a=0.5, lowpass_filter_f=-1, lowpass_filter_a=-1.0;
int highpass_filter_order=5, lowpass_filter_order=-1;

/* large refers to the large SFT that produces power statistic */
REAL4Vector *sft_in=NULL;
REAL4Vector *sft_out=NULL;
COMPLEX8Vector *sft_out2=NULL;
RealFFTPlan *large_forward=NULL;

/* After calibration SFT needs to be windowed to get rid of bogus data outside of calibration range 
   we need to cut off at least some of low frequencies as DC value is always bogus,
   as documented in LAL software manual. */
long freq_window_start=30*60, freq_window_stop=2000*60;
char *debug_freq_window=NULL;

MAKE_SFT_WINDOW freq_window={
	type: Hann, 
	tail_size: 512,
	debug: NULL
	};
	
MAKE_SFT_WINDOW window1={
	type: Hann, 
	tail_size: 512,
	debug: NULL
	};

MAKE_SFT_WINDOW window2={
	type: Hann, 
	tail_size: 512,
	debug: NULL
	};


typedef struct {
	double line1;
	double line2;
	int N;
	} LINE_SUPPRESS_REGION;
	
LINE_SUPPRESS_REGION *suppress_lines=NULL;
int sl_size=0;
int sl_free=0;
	

void add_suppress_region(double line1, double line2, int N)
{
LINE_SUPPRESS_REGION *p;
if(sl_free>=sl_size){
	sl_size=2*sl_size+10;
	p=do_alloc(sl_size, sizeof(*p));
	if(sl_free>0)memcpy(p, suppress_lines, sl_free*sizeof(*suppress_lines));
	if(suppress_lines!=NULL)free(suppress_lines);
	suppress_lines=p;
	}
suppress_lines[sl_free].line1=line1;
suppress_lines[sl_free].line2=line2;
suppress_lines[sl_free].N=N;
sl_free++;
}

void print_suppress_lines(FILE *f)
{
int i;
for(i=0;i<sl_free;i++){
	fprintf(f,"#    %s suppressed: %lg through %lg including %d harmonics\n",
		bypass_lines ? "NOT" : "", 
		suppress_lines[i].line1, suppress_lines[i].line2, suppress_lines[i].N);
	}
}

void print_bypass_options(FILE *f)
{
int none=1;
fprintf(f,"#    bypass options:");
if(bypass_highpass_filter){
	fprintf(f, " HIGHPASS_FILTER");
	none=0;
	}
if(bypass_first_window){
	fprintf(f, " FIRST_WINDOW");
	none=0;
	}
if(bypass_calibration){
	fprintf(f, " CALIBRATION");
	none=0;
	}
if(bypass_lines){
	fprintf(f, " LINE_SUPPRESSION");
	none=0;
	}

if(bypass_freq_window){
	fprintf(f, " FREQ_WINDOW");
	none=0;
	}

if(bypass_first_fft){
	fprintf(f, " FIRST_FFT");
	none=0;
	}

if(none){
	fprintf(f, " none");
	none=0;
	}
fprintf(f, "\n");
}

void print_filter_settings(FILE *f)
{
fprintf(f, "#    highpass filter N: %d\n", highpass_filter_order);
fprintf(f, "#    highpass filter F: %.*g\n", precision, highpass_filter_f);
fprintf(f, "#    highpass filter A: %.*g\n", precision, highpass_filter_a);
fprintf(f, "#    lowpass filter N: %d\n", lowpass_filter_order);
fprintf(f, "#    lowpass filter F: %.*g\n", precision, lowpass_filter_f);
fprintf(f, "#    lowpass filter A: %.*g\n", precision, lowpass_filter_a);
}



void print_settings(FILE *f)
{
fprintf(f,"#    starting bin:	    %ld\n", freq_start);
fprintf(f,"#    stopping bin:       %ld\n", freq_stop);
print_window_settings(f, "first", &window1);
print_window_settings(f, "frequency", &freq_window);
print_window_settings(f, "second", &window2);
if(overlap!=0)fprintf(f,"#    overlap: %ld\n", overlap);
print_suppress_lines(f);
print_bypass_options(f);
print_filter_settings(f);
print_fake_data_settings(f);
}

void calibrate_data(void)
{
LALStatus status={level:0, statusPtr:NULL};
REAL4Vector *series=NULL, *window1_vec=NULL, *window2_vec=NULL, *freq_window_vec=NULL;
RealFFTPlan *forward=NULL, *reverse=NULL;
COMPLEX8Vector *fft=NULL;
long i,j,k,m;
REAL4 a;
long nsamples;
LALWindowParams window_parameters;
PassBandParamStruc filterpar;
REAL4TimeSeries chan; /* must zero the f0 field */
long sl_start, sl_end;
REAL4 sl_left, sl_right;

nsamples=60*samples_per_second;

if(overlap>=nsamples){
	fprintf(stderr, "** Overlap (currently %ld) must be smaller than nsamples (currently %ld)\n", overlap, nsamples);
	exit(-1);
	}
	
if(overlap<0){
	fprintf(stderr, "** Overlap (currently %ld) must be positive\n", overlap);
	exit(-1);
	}

fprintf(stderr,"Starting data calibration, total_samples=%ld window nsamples=%ld\n", total_samples, nsamples);

/* allocate variables */

LALSCreateVector(&status, &series, nsamples);
TESTSTATUS(&status);

   
LALSCreateVector(&status, &window1_vec, window1.tail_size*2);
TESTSTATUS(&status);

LALSCreateVector(&status, &freq_window_vec, freq_window.tail_size*2);
TESTSTATUS(&status);

LALSCreateVector(&status, &window2_vec, window2.tail_size*2);
TESTSTATUS(&status);

LALCCreateVector(&status, &fft, nsamples/2+1);
TESTSTATUS(&status);

LALCreateForwardRealFFTPlan(&status, &forward, nsamples, 0);
TESTSTATUS(&status);

LALCreateReverseRealFFTPlan(&status, &reverse, nsamples, 0);
TESTSTATUS(&status);

/* prepare window */
window_parameters.length=window1.tail_size*2;
window_parameters.type=window1.type; 
for(i=0;i<window1.tail_size*2;i++)window1_vec->data[i]=1.0;
LALWindow(&status, window1_vec, &window_parameters);
TESTSTATUS(&status);

print_REAL4Vector(window1_vec, window1.debug, "WINDOW1, DEBUG", output_mode, 0, window1_vec->length);

if(window2.tail_size>0){
	/* prepare window */
	window_parameters.length=window2.tail_size*2;
	window_parameters.type=window2.type; 
	for(i=0;i<window2.tail_size*2;i++)window2_vec->data[i]=1.0;
	LALWindow(&status, window2_vec, &window_parameters);
	TESTSTATUS(&status);
	}

print_REAL4Vector(window2_vec, window1.debug, "WINDOW2, DEBUG", output_mode, 0, window2_vec->length);


/* prepare frequency window */
window_parameters.length=freq_window.tail_size*2;
window_parameters.type=freq_window.type; 
for(i=0;i<freq_window.tail_size*2;i++)freq_window_vec->data[i]=1.0;
LALWindow(&status, freq_window_vec, &window_parameters);
TESTSTATUS(&status);

print_REAL4Vector(freq_window_vec, debug_freq_window, "FREQUENCY WINDOW, DEBUG", output_mode, 0, freq_window_vec->length);


/* setup structure to hold input for Butterworth filter */
chan.data=NULL;
strncpy(chan.name,channel,LALNameLength);
chan.f0=0;
chan.name[LALNameLength-1]=0; /* make sure it is null-terminated */
chan.deltaT=1.0/samples_per_second;
chan.epoch.gpsSeconds=0; /* no need */
chan.epoch.gpsNanoSeconds=0;

/* initialize output */
for(i=0;i<overlap;i++){
	sft_in->data[i]=0.0;
	}

/* perform calibration */
for(i=0;i+nsamples<=total_samples;i+=(nsamples-overlap)){
	fprintf(stderr, "\ti=%ld nsamples=%ld\n", i, nsamples);
	if(trace_calibration_power){
		fprintf(stderr, "\tcal_start power=%.*g\n",
			precision, sum_r4_squares(&(data[i]), nsamples)/samples_per_second);
		}
		
	/* copy data */
	for(k=0;k<nsamples;k++){		
		series->data[k]=data[i+k];
		}

	if(!bypass_highpass_filter && (highpass_filter_f>0) && (highpass_filter_a>0)){
		/* Setup Butterworth filter */
		filterpar.name="Butterworth High Pass";
		filterpar.nMax=highpass_filter_order;
		filterpar.f2=highpass_filter_f;
		filterpar.a2=highpass_filter_a;
		/* values that are 'not given' = out of range */
		filterpar.f1=-1.0;
		filterpar.a1=-1.0;
		
		/* REAL4Sequence is the same as REAL4Vector - according to lal/Datatypes.h */
		chan.data=(REAL4Sequence *)series;
		LALDButterworthREAL4TimeSeries(&status, &chan, &filterpar);
		TESTSTATUS(&status);

		if(trace_calibration_power){
			fprintf(stderr, "\thigh_pass power=%.*g\n",
				precision, sum_r4_squares(series->data, nsamples)/samples_per_second);
			}
		}
		
	if(!bypass_lowpass_filter && (lowpass_filter_f>0) && (lowpass_filter_a > 0)){
		/* Setup Butterworth filter */
		filterpar.name="Butterworth Low Pass";
		filterpar.nMax=lowpass_filter_order;
		filterpar.f1=lowpass_filter_f;
		filterpar.a1=lowpass_filter_a;
		/* values that are 'not given' = out of range */
		/* why 2*nsamples ? LAL used to have a bug in it where
		   the pairs were sorted before deciding whether the filter 
		   is lowpass or highpass. Therefore, we need to specify a 
		   large, out of range, frequency f so that we get a low-pass filter.
		   The maximum frequency is Nyquist, hence 2*nsamples is definitely out of range */
		filterpar.f2=2*nsamples;
		filterpar.a2=-1.0;

		/* REAL4Sequence is the same as REAL4Vector - according to lal/Datatypes.h */
		chan.data=(REAL4Sequence *)series;
		LALDButterworthREAL4TimeSeries(&status, &chan, &filterpar);
		TESTSTATUS(&status);

		if(trace_calibration_power){
			fprintf(stderr, "\tlow_pass power=%.*g\n",
				precision, sum_r4_squares(series->data, nsamples)/samples_per_second);
			}
		}

	if(!bypass_first_window){
		/* window the result */
		for(k=0;k<window1.tail_size;k++){		
			series->data[k]*=window1_vec->data[k];
			}
		for(k=nsamples-window1.tail_size;k<nsamples;k++){		
			series->data[k]*=window1_vec->data[nsamples-k-1];
			}

		if(trace_calibration_power){
			fprintf(stderr, "\twindow1 power=%.*g\n",
				precision, sum_r4_squares(series->data, nsamples)/samples_per_second);
			}
		}
		
	if(!bypass_first_fft){
		LALForwardRealFFT(&status, fft, series, forward);
		TESTSTATUS(&status);

		if(trace_calibration_power){
			fprintf(stderr, "\tfft_power=%.*g\n",
				precision, (2.0*sum_r4_squares((REAL4 *)fft->data, fft->length*2)-fft->data[0].re*fft->data[0].re)/(1.0*nsamples*samples_per_second));
			}
		}

	if(!bypass_calibration){
		calibrate_fft(fft, i);	

		if(trace_calibration_power){
			fprintf(stderr, "\tcalibration power=%.*g\n",
				precision, (2.0*sum_r4_squares((REAL4 *)fft->data, fft->length*2)-fft->data[0].re*fft->data[0].re)/(1.0*nsamples*samples_per_second));
			}
		}
	
	if(!bypass_lines){
		/* suppress lines */
		for(k=0;k<sl_free;k++){
			for(j=1;j<=suppress_lines[k].N;j++){
				sl_start=lrint(60.0*j*suppress_lines[k].line1);
				sl_end=lrint(60.0*j*suppress_lines[k].line2);
				if(sl_end<1)continue;
				if(sl_start>=(nsamples/2-1))break;
				if(sl_start<1)sl_start=1;
				if(sl_end>=(nsamples/2-1))sl_end=nsamples/2-2;
				if(sl_start>sl_end)continue;
				sl_left=fft->data[sl_start-1].re;
				sl_right=fft->data[sl_end+1].re;
				for(m=sl_start;m<=sl_end;m++){
					fft->data[m].re=(sl_left*(sl_end-m+1)+sl_right*(m-sl_start+1))/(sl_end-sl_start+2);
					}			
				sl_left=fft->data[sl_start-1].im;
				sl_right=fft->data[sl_end+1].im;
				for(m=sl_start;m<=sl_end;m++){
					fft->data[m].im=(sl_left*(sl_end-m+1)+sl_right*(m-sl_start+1))/(sl_end-sl_start+2);
					}			
				}
			}

		if(trace_calibration_power){
			fprintf(stderr, "\tline_sup power=%.*g\n",
				precision, (2.0*sum_r4_squares((REAL4 *)fft->data, fft->length*2)-fft->data[0].re*fft->data[0].re)/(1.0*nsamples*samples_per_second));
			}
		}
		
	if(!bypass_freq_window){
		/* window the result */
		for(k=0;k<fft->length;k++){
			if((k<freq_window_start)||(k>=freq_window_stop)){
				fft->data[k].re=0.0;
				fft->data[k].im=0.0;
				continue;
				}
			if(k<(freq_window_start+freq_window.tail_size)){
				a=freq_window_vec->data[k-freq_window_start];
				fft->data[k].re*=a;
				fft->data[k].im*=a;
				}
			if(k>=(freq_window_stop-freq_window.tail_size)){
				a=freq_window_vec->data[freq_window_stop-k-1];
				fft->data[k].re*=a;
				fft->data[k].im*=a;
				}
			}

		if(trace_calibration_power){
			fprintf(stderr, "\tfreq_window power=%.*g\n",
				precision, (2.0*sum_r4_squares((REAL4 *)fft->data, fft->length*2)-fft->data[0].re*fft->data[0].re)/(1.0*nsamples*samples_per_second));
			}
		}
	
	if(!bypass_first_fft){
		LALReverseRealFFT(&status, series, fft, reverse);
		TESTSTATUS(&status);

		if(trace_calibration_power){
			fprintf(stderr, "\tinv_fft power=%.*g\n",
				precision, sum_r4_squares(series->data, nsamples)/(1.0*nsamples*nsamples*samples_per_second));
			}
		}

	if(!bypass_first_fft){
		/* copy data and window it again */
		for(k=0;k<overlap;k++){
			/* we need to divide by nsamples to obtain true
			   FFT inverse */
			   
			if(k<window2.tail_size)sft_in->data[i+k]+=series->data[k]*window2_vec->data[k]/nsamples;
				else
			if(k>=(nsamples-window2.tail_size))sft_in->data[i+k]+=series->data[k]*window2_vec->data[nsamples-k-1]/nsamples;
				else 
				sft_in->data[i+k]+=series->data[k]/nsamples;
			}
		for(k=overlap;k<nsamples;k++){
			/* we need to divide by nsamples to obtain true
			   FFT inverse */
			   
			if(k<window2.tail_size)sft_in->data[i+k]=series->data[k]*window2_vec->data[k]/nsamples;
				else
			if(k>=(nsamples-window2.tail_size))sft_in->data[i+k]=series->data[k]*window2_vec->data[nsamples-k-1]/nsamples;
				else 
				sft_in->data[i+k]=series->data[k]/nsamples;
			}

		if(trace_calibration_power){
			fprintf(stderr, "\twindow2 power=%.*g\n",
				precision, sum_r4_squares(&(sft_in->data[i]), nsamples)/samples_per_second);
			}
		} else {
		/* copy data and window it again */
		for(k=0;k<overlap;k++){
			/* we need to divide by nsamples to obtain true
			   FFT inverse */
			   
			if(k<window2.tail_size)sft_in->data[i+k]+=series->data[k]*window2_vec->data[k];
				else
			if(k>=(nsamples-window2.tail_size))sft_in->data[i+k]+=series->data[k]*window2_vec->data[nsamples-k-1];
				else 
				sft_in->data[i+k]+=series->data[k];
			}
		for(k=overlap;k<nsamples;k++){
			/* we need to divide by nsamples to obtain true
			   FFT inverse */
			   
			if(k<window2.tail_size)sft_in->data[i+k]=series->data[k]*window2_vec->data[k];
				else
			if(k>=(nsamples-window2.tail_size))sft_in->data[i+k]=series->data[k]*window2_vec->data[nsamples-k-1];
				else 
				sft_in->data[i+k]=series->data[k];
			}
		if(trace_calibration_power){
			fprintf(stderr, "\twindow2 power=%.*g\n",
				precision, sum_r4_squares(&(sft_in->data[i]), nsamples)/samples_per_second);
			}
		}
	}

/* deallocate variables */

LALDestroyRealFFTPlan(&status, &reverse);
TESTSTATUS(&status);

LALDestroyRealFFTPlan(&status, &forward);
TESTSTATUS(&status);

LALCDestroyVector(&status, &fft);
TESTSTATUS(&status);

LALSDestroyVector(&status, &window1_vec);
TESTSTATUS(&status);

LALSDestroyVector(&status, &freq_window_vec);
TESTSTATUS(&status);

LALSDestroyVector(&status, &window2_vec);
TESTSTATUS(&status);

LALSDestroyVector(&status, &series);
TESTSTATUS(&status);
}

void set_debug_output(char *p)
{
switch(p[5]){
	case '1':
		if(file_debug1!=NULL)free(file_debug1);
		file_debug1=dup_quoted_name(p);
		break;
	case '2':
		if(file_debug2!=NULL)free(file_debug2);
		file_debug2=dup_quoted_name(p);
		break;
	default:
		fprintf(stderr,"** Could not parse statement \"%s\"\n", p);
		exit(-1);
	}
}

#define BUF_SIZE 200000
char buffer[BUF_SIZE];

int main(int argc, char *argv[])
{
LALStatus status={level:0, statusPtr: NULL};

fprintf(stderr,"make_sft version %s\n", MAKE_SFT_VERSION);
fprintf(stderr,"Using frame library %s\n", FrLibVersionF());
fprintf(stderr,"Using LAL version %s, CVS tag %s\n", LAL_VERSION, LAL_CVS_TAG);

yylex();

if(freq_start<0)freq_start=0;
if(freq_stop<0)freq_stop=total_samples;

/* post_init various subsystems */

post_init_response_files();
post_init_alpha_beta();

/* print settings on stderr */
print_settings(stderr);

/* start processing data */
print_data_stats();
print_data(file_debug1, "debug 1");

verify_loaded_data();

generate_fake_data();

linear_interpolate_gaps();
print_data(file_debug2, "debug 2");

if(overlap!=0){
	/* allocate space for large vectors */
	LALSCreateVector(&status, &(sft_in), total_samples);
	TESTSTATUS(&status);
	} else {
	/* this is a hack, but it significantly reduces memory footprint */
	sft_in=do_alloc(1,sizeof(*sft_in));
	sft_in->length=total_samples;
	sft_in->data=data;
	}

LALCreateForwardRealFFTPlan(&status, &large_forward, total_samples, 0);
TESTSTATUS(&status);

if(trace_power){
	fprintf(stderr, "Before calibration total power=%.*g\n",
		precision, sum_r4_squares(data, total_samples)/samples_per_second);
	}

calibrate_data();

if(trace_power){
	fprintf(stderr, "After calibration total power=%.*g\n",
		precision, sum_r4_squares(sft_in->data, sft_in->length)/samples_per_second);
	}

/* free *data which is huge and not needed anymore */
if(overlap!=0){
	free(data);
	data=NULL;
	}

print_REAL4Vector(sft_in, output_calibrated_data, "CALIBRATED DATA", output_mode, 0, sft_in->length);

if(output_power!=NULL){
	LALSCreateVector(&status, &(sft_out), total_samples/2+1);
	TESTSTATUS(&status);

	LALRealPowerSpectrum(&status, sft_out, sft_in, large_forward);
	TESTSTATUS(&status);

	/* GEO format is used only for SFTs */
	print_REAL4Vector(sft_out, output_power, "POWER", output_mode, freq_start, freq_stop);

	if(trace_power){
		fprintf(stderr, "Power SFT total power=%.*g\n",
			precision, (2.0*sum_positive_r4(sft_out->data, sft_out->length)-sft_out->data[0])/(1.0*total_samples*samples_per_second));
		}

	LALSDestroyVector(&status, &sft_out);
	TESTSTATUS(&status);
	}

if(output_sft!=NULL){
	LALCCreateVector(&status, &(sft_out2), total_samples/2+1);
	TESTSTATUS(&status);

	LALForwardRealFFT(&status, sft_out2, sft_in, large_forward);
	TESTSTATUS(&status);

	/* GEO format is used only for SFTs */
	print_COMPLEX8Vector(sft_out2, output_sft, "SFT", output_mode, freq_start, freq_stop);

	if(trace_power){
		fprintf(stderr, "SFT total power=%.*g\n",
			precision, (2.0*sum_r4_squares((REAL4 *)sft_out2->data, sft_out2->length*2)-sft_out2->data[0].re*sft_out2->data[0].re)/(1.0*total_samples*samples_per_second));
		}

	LALCDestroyVector(&status, &sft_out2);
	TESTSTATUS(&status);
	}

/* deallocate space for large vectors */
LALDestroyRealFFTPlan(&status, &large_forward);
TESTSTATUS(&status);

if(overlap!=0){
	LALSDestroyVector(&status, &sft_in);
	TESTSTATUS(&status);
	} else {
	/* why bother freeing a huge vector when we are exiting anyway ?
	   this is important as freeing causes memory allocator to 
	   touch a lot of RAM (possibly having it swapped in) for no
	   reason whatsoever */
	}


return 0;
}
