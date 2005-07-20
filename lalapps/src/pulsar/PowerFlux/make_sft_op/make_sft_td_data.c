/*

  Time domain functions and data:
  
  data loading, simulation and manipulation.

*/
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

#include <FrIO.h>
#include <FrameL.h>

#include "make_sft_global.h"

extern int precision;
extern int dont_fail_on_missing_data;

typedef struct {
	char *name;
	INT64 gps_start;
	INT64 gps_end;
	} frame_file;

frame_file *files=NULL;
int files_size=0;
int files_free=0;
char *channel=NULL;
INT64 gps_start=-1, gps_end=-1;
long samples_per_second=-1;

/* data we read from the frames */
#ifdef DOUBLE_DATA
double *data=NULL;
#else
float *data=NULL;
#endif

/* bitmap that indicates which data is present */
unsigned char *data_present=NULL;
/* bitmap that indicates valid data segments */
unsigned char *segments=NULL;
long total_samples=-1;

#define DATA_PRESENT(i)		(data_present[(i)>>3] & (1<<((i)&7)))


/* generate test data */
FAKE_DATA_COMMAND *fake_data;
int fake_data_size=0;
int fake_data_free=0;

void add_fake_data_command(int command, double phase, double frequency, double amplitude)
{
FAKE_DATA_COMMAND *p;
int i;
if(fake_data_free>=fake_data_size){
	fake_data_size=2*fake_data_size+10;
	p=do_alloc(fake_data_size, sizeof(FAKE_DATA_COMMAND));
	if(fake_data_free>0)memcpy(p, fake_data, fake_data_free*sizeof(FAKE_DATA_COMMAND));
	if(fake_data!=NULL)free(fake_data);
	fake_data=p;
	}
i=fake_data_free;
fake_data_free++;
fake_data[i].command=command;
fake_data[i].phase=phase;
fake_data[i].frequency=frequency;
fake_data[i].amplitude=amplitude;
}

void print_fake_data_settings(FILE *f)
{
int i;
for(i=0;i<fake_data_free;i++)
	switch(fake_data[i].command){
		case DONT_FAKE_DATA:
			fprintf(f, "#    fake data: not enabled\n");
			break;
		case FAKE_DATA_CONSTANT:
			fprintf(f, "#    fake data: CONSTANT\n");
			fprintf(f, "#    \tfd amplitude: %lg\n", fake_data[i].amplitude);
			break;
		case FAKE_DATA_CONSTANT_FILL:
			fprintf(f, "#    fake data: CONSTANT_FILL\n");
			fprintf(f, "#    \tfd amplitude: %lg\n", fake_data[i].amplitude);
			break;
		case FAKE_DATA_RANDOM_GAUSSIAN:
			fprintf(f, "#    fake data: RANDOM_GAUSSIAN\n");
			fprintf(f, "#    \tfd amplitude: %lg\n", fake_data[i].amplitude);
			break;
		case FAKE_DATA_SINE_WAVE:
			fprintf(f, "#    fake data: SINE_WAVE\n");
			fprintf(f, "#    \tfd frequency: %lg\n", fake_data[i].frequency);
			fprintf(f, "#    \tfd amplitude: %lg\n", fake_data[i].amplitude);
			fprintf(f, "#    \tfd phase: %lg\n", fake_data[i].phase);
			break;
		case FAKE_DATA_CALIBRATED_SINE_WAVE:
			fprintf(f, "#    fake data: CALIBRATED_SINE_WAVE\n");
			fprintf(f, "#    \tfd frequency: %lg\n", fake_data[i].frequency);
			fprintf(f, "#    \tfd amplitude: %lg\n", fake_data[i].amplitude);
			fprintf(f, "#    \tfd phase: %lg\n", fake_data[i].phase);
			break;
		default:
			fprintf(f, "#    fake data: **UNKNOWN**\n");
			break;
		}

}

void generate_fake_data(void)
{
long j,i;
double a,b;
if(fake_data_free>0){
	fprintf(stderr, "Generating fake data\n");
	}
for(i=0;i<total_samples;i++)
	for(j=0;j<fake_data_free;j++)
	   switch(fake_data[j].command){
		case FAKE_DATA_CONSTANT_FILL:
			data[i]=fake_data[j].amplitude;
			break;
		case FAKE_DATA_CONSTANT:
			data[i]+=fake_data[j].amplitude;
			break;
		case FAKE_DATA_RANDOM_GAUSSIAN:
			do{
				a=(1.0*random())/RAND_MAX;
				b=(6.0*random())/RAND_MAX-3.0;
				} while(sqrt(2*M_PI)*a> exp(-b*b/2.0));
			data[i]+=fake_data[j].amplitude*b;
			break;
		case FAKE_DATA_SINE_WAVE:
			data[i]+=fake_data[j].amplitude*sin((2.0*M_PI*i*fake_data[j].frequency)/samples_per_second + fake_data[j].phase);
			break;		
		case FAKE_DATA_CALIBRATED_SINE_WAVE:
			get_td_calibration(fake_data[j].frequency, i, &a, &b);
			data[i]+=fake_data[j].amplitude*sin((2.0*M_PI*i*fake_data[j].frequency)/samples_per_second + fake_data[j].phase+atan2(a,b))/sqrt(a*a+b*b);
			break;		
		}

/* check whether we need to mess with data_present bits */
/* for later, if needed: simply do data_present[i]=segments[i] */
for(j=0;j<fake_data_free;j++)
	if(fake_data[j].command==FAKE_DATA_CONSTANT_FILL){
		for(i=0;i<(total_samples+7)>>3;i++)data_present[i]=0xff;
		break;
		}
}


void allocate_data(void)
{
if((gps_start<0)||(gps_end<0)||(samples_per_second<0))return;
total_samples=(gps_end-gps_start)*samples_per_second;
data=calloc(total_samples, sizeof(*data));
if(data==NULL){
	fprintf(stderr,"** Could not allocate memory for sample storage. gps_start=%lld gps_end=%lld samples_per_second=%ld\n",
		gps_end, gps_start, samples_per_second);
	exit(-1);
	}
fprintf(stderr, "Allocated storage for %ld samples (%.03f MB)\n", total_samples,(total_samples*sizeof(*data))/(1024.0*1024.0));
data_present=do_alloc((total_samples+7)>>3, 1);
segments=do_alloc((total_samples+7)>>3, 1);
/* set all bits to 0 - no data */
memset(data_present, 0, (total_samples+7)>>3);
memset(segments, 0, (total_samples+7)>>3);
}

void verify_loaded_data(void)
{
long j,i;
unsigned char c;
/* mark as not present samples outside valid segments */
j=-1;
for(i=0;i<(total_samples+7)>>3;i++){
	c=data_present[i]&segments[i];
	data_present[i]=c;
	if(c!=segments[i]){
		for(j=i+1;j<(total_samples+7)>>3;j++){
			c=data_present[j]&segments[j];
			if(c==segments[j])break;
			}
		fprintf(stderr,"** Not all data specified in segments has been loaded, in particular see gps time from %lld through %lld\n",
			gps_start+((i<<3)/samples_per_second), gps_start+((j<<3)+7)/samples_per_second);
		i=j;
		}
	}
if((j>=0) && !dont_fail_on_missing_data){
	fprintf(stderr,"** Not all requested data was found, exit as requested\n");
	exit(-1);
	}
}
/* just for testing. Compute and print trends of FrVect 
   slightly buggy - does not handle dimensions */
void print_vect_trends(FrVect *vect, long offset, long step)
{
int max,min,a;
double mean, sqr;
float maxf,minf, af;
long i,k;
switch(vect->type){
	case FR_VECT_2S:
		for(i=offset;i+step<=vect->nData;i+=step){
			max=vect->dataS[i];
			min=vect->dataS[i];
			mean=0.0;
			sqr=0.0;
			for(k=0;k<step;k++){
				a=vect->dataS[i+k];
				if(a<min)min=a;
				if(a>max)max=a;
				mean+=a;
				sqr+=a*a;
				}
			fprintf(stderr,"vect=\"%s\" min=%d max=%d mean=%g rms=%g\n",vect->name, min, max, mean/step, sqrt(sqr/step));
			}
		break;
	case FR_VECT_4R:
		for(i=offset;i+step<=vect->nData;i+=step){
			maxf=vect->dataF[i];
			minf=vect->dataF[i];
			mean=0.0;
			sqr=0.0;
			for(k=0;k<step;k++){
				af=vect->dataF[i+k];
				if(af<minf)minf=af;
				if(af>maxf)maxf=af;
				mean+=af;
				sqr+=af*af;
				}
			fprintf(stderr,"vect=\"%s\" min=%g max=%g mean=%g rms=%g\n",vect->name, minf, maxf, mean/step, sqrt(sqr/step));
			}
		break;
	case FR_VECT_8R:
		break;
	}

}

#if 0
/*
Slower version that does handle vect->next
*/
void assimilate_data_next(FrVect *vect)
{
long i;
float *p1;
float *p2;
unsigned char *p3;
float *p2_flag;
long offset;
long nsamples,bit;
unsigned char mask;
offset=(vect->GTime-gps_start)*samples_per_second;
nsamples=vect->nData;
fprintf(stderr,"vect: gps_start %f  length %ld, data length %ld\n", vect->GTime, vect->nData/samples_per_second, vect->nData);
i=0;
if(offset<0){
        nsamples+=offset;
        i=-offset;
        offset=0;
        }
if(offset+nsamples>total_samples)nsamples=total_samples-offset;
if(nsamples<=0)return; /* no data to assimilate */

switch(vect->type){
        case FR_VECT_4R:
		p2=&(vect->dataF[i]);
		p2_flag=&(vect->next->dataF[i]);
		for(;i<nsamples;i++){
			fprintf(stderr,"*** ERROR INCOMPLETE CODE\n");
			}
		break;
	default:
	fprintf(stderr,"** format %d for FrVect data is not supported yet\n", vect->type);
	}
}

#endif

/* 
   does not handle dimensions or vect->next 
*/

void assimilate_data(FrVect *vect)
{
long i;
#ifdef DOUBLE_DATA
	double *p1;
#else
	float *p1;
#endif
float *p2;
unsigned char *p3;
double *p2d;
long offset;
long nsamples,bit;
if(vect->next!=NULL){
	fprintf(stderr,"** Cannot handle vect->next!=NULL\n");
	fprintf(stderr,"** vect->GTime=%g\n", vect->GTime);
	exit(-1);
	}
offset=(vect->GTime-gps_start)*samples_per_second;
nsamples=vect->nData;
fprintf(stderr,"vect: gps_start %f  length %ld, data length %ld\n", vect->GTime, vect->nData/samples_per_second, vect->nData);
i=0;
if(offset<0){
	nsamples+=offset;
	i=-offset;
	offset=0;
	}
if(offset+nsamples>total_samples)nsamples=total_samples-offset;
if(nsamples<=0)return; /* no data to assimilate */

switch(vect->type){
	case FR_VECT_4R:
		p2=&(vect->dataF[i]);
		p1=&(data[offset]);
		p3=&(data_present[offset>>3]);
			#define COPY_REAL *p1=*p2;p1++;p2++;
		/* start copying - we might not have offset to be a multiple of 8 */
		bit=offset&7;
		i=0;
		if(bit){
			*p3|=((1<<(8-bit))-1)<<bit;
			p3++;
			for(;(bit)&&(i<nsamples);i++,bit--){ COPY_REAL }
			}
		/* continue copying in chunks of 8 */
		for(;i+8<=nsamples;i+=8){
			/* copy 8 reals to their new places */
			COPY_REAL
			COPY_REAL
			COPY_REAL
			COPY_REAL
			COPY_REAL
			COPY_REAL
			COPY_REAL
			COPY_REAL
			/* mark all 8 bits at once */
			*p3=0xff;p3++;
			}			
		/* copy and mark the remaining bits */
		*p3|=(1<<(nsamples-i))-1;
		for(;i<nsamples;i++){ COPY_REAL }
		#undef COPY_REAL
		break;
	case FR_VECT_8R:
		p2d=&(vect->dataD[i]);
		p1=&(data[offset]);
		p3=&(data_present[offset>>3]);
			#define COPY_REAL *p1=*p2d;p1++;p2d++;
		/* start copying - we might not have offset to be a multiple of 8 */
		bit=offset&7;
		i=0;
		if(bit){
			*p3|=((1<<(8-bit))-1)<<bit;
			p3++;
			for(;(bit)&&(i<nsamples);i++,bit--){ COPY_REAL }
			}
		/* continue copying in chunks of 8 */
		for(;i+8<=nsamples;i+=8){
			/* copy 8 reals to their new places */
			COPY_REAL
			COPY_REAL
			COPY_REAL
			COPY_REAL
			COPY_REAL
			COPY_REAL
			COPY_REAL
			COPY_REAL
			/* mark all 8 bits at once */
			*p3=0xff;p3++;
			}			
		/* copy and mark the remaining bits */
		*p3|=(1<<(nsamples-i))-1;
		for(;i<nsamples;i++){ COPY_REAL }
		#undef COPY_REAL
		break;
	case FR_VECT_2S:
	default:
		fprintf(stderr,"** format %d for FrVect data is not supported yet\n", vect->type);
		break;
	}
}

/* this function assumes that each file contains one or more frames
   and that the data spans contiguous time segment (or, at least,
   the segments from different files do not overlap */

void add_frame_file(char *filename, double TStart, double TEnd)
{
char *p;
FrFile *ff;
FrVect *vect;
frame_file *fdata;
/* check that sample storage has been allocated */
if(data==NULL){
	fprintf(stderr, "** You must use DURATION and SAMPLES_PER_SECOND commands before FRAME_FILES\n");
	exit(-1);
	}

/* remove leading whitespace */
while((*filename==' ')||(*filename=='\t'))filename++;
/* break at first newline */
for(p=filename;*p;p++)if(*p=='\n'){*p=0; break; }
/* remove trailing whitespace */
while((*p==' ')||(*p=='\t')){
	*p=0;
	if(p==filename)break;
	p--;
	}

fprintf(stderr,"Adding frame file \"%s\"\n", filename);
if(files_free>=files_size){
	files_size=2*files_size+10;
	/* reuse p.. we do not care it is the wrong type, void * is sufficient */
	p=(char *)do_alloc(files_size, sizeof(frame_file));
	if(files_free>0)memcpy(p, files, files_free*sizeof(frame_file));
	if(files!=NULL)free(files);
	files=(frame_file *)p;
	}
fdata=&(files[files_free]);
fdata->name=strdup(filename);
fault_file(filename);
ff=FrFileINew(filename);
if(ff==NULL){
	fprintf(stderr,"** Error opening \"%s\" for reading frame fdata\n", filename);
	return;
	}
if(TStart<0)TStart=FrFileITStart(ff);
if(TEnd<0)TEnd=FrFileITEnd(ff);
fprintf(stderr,"Opened frame file \"%s\" start=%f end=%f\n", filename, TStart, TEnd);
if((TStart > TEnd) || (TStart < 0)){
	fprintf(stderr, "** TStart > TEnd || TStart < 0 on file %s\n", filename);
	FrFileIEnd(ff);
	return;
	}
/*  This does not appear to work - vect returns bogus values for start time.
if(TStart < gps_start)TStart=gps_start;
if(TEnd > gps_end)TEnd=gps_end;
*/
/* ok, file is opened and we _are_ able to read from it.. 
   validate the entry and read the fdata */
files_free++;
fdata->gps_start=llrint(TStart);
fdata->gps_end=llrint(TEnd);
errno=0;
vect=FrFileIGetV(ff, channel, TStart, TEnd-TStart);
if((vect==NULL)||(errno!=0)){
	fprintf(stderr,"** Could not find data for channel \"%s\" [%.10lg,%.10lg] in file \"%s\" (errno=%s)\n", channel, TStart, TEnd, filename, strerror(errno));
	} else {
	/* print_vect_trends(vect, 0, 16384); */
	assimilate_data(vect);
	FrVectFree(vect);
	}
FrFileIEnd(ff);
}

void add_start_stop_frame_file(char *text)
{
double TStart=-1.0, TEnd=-1.0;
int name_start;
if(sscanf(text, "%lg %lg%n", &TStart, &TEnd, &name_start)!=2){
	fprintf(stderr,"** Error parsing start_stop_framefile: \"%s\"=(%.10lg %.10lg %d)\n", text, TStart, TEnd, name_start);
	return;
	}
add_frame_file(&(text[name_start]), TStart, TEnd);
}

void mark_segment(long start, long end)
{
long i;
fprintf(stderr,"marking segment %ld %ld\n", start, end);
if(start<0)start=0;
if(end>total_samples)end=total_samples;
for(i=start;i<end;i++){
	segments[i>>3]|=1<<(i & 7);
	if(!(i&7))while((i+8<=end)){
		segments[i>>3]=0xff;
		i+=8;
		}
	}
}

void assimilate_segment(char *line)
{
INT64 start,end;
if(segments==NULL){
	fprintf(stderr,"** SEGMENTS *must* follow DURATION and SAMPLES_PER_SECOND commands\n");
	exit(-1);
	}
sscanf(line,"%lld %lld", &start, &end);
mark_segment((start-gps_start)*samples_per_second, (end-gps_start)*samples_per_second);
}

void assimilate_subsecond_segment(char *line)
{
double start,end;
if(segments==NULL){
	fprintf(stderr,"** SUBSECOND_SEGMENTS *must* follow DURATION and SAMPLES_PER_SECOND commands\n");
	exit(-1);
	}
sscanf(line,"%lg %lg", &start, &end);
mark_segment(lrint((start-gps_start)*samples_per_second), lrint((end-gps_start)*samples_per_second));
}

void assimilate_sample_segment(char *line)
{
long start,end;
if(segments==NULL){
	fprintf(stderr,"** SAMPLE_SEGMENTS *must* follow DURATION and SAMPLES_PER_SECOND commands\n");
	exit(-1);
	}
sscanf(line,"%ld %ld", &start, &end);
mark_segment(start, end);
}

void print_data_stats(void)
{
fprintf(stderr,"Total samples %ld , samples loaded %ld, segment samples %ld\n", 
	total_samples, 
	count_bitmap_bits(data_present, total_samples),
	count_bitmap_bits(segments, total_samples));
}

void linear_interpolate_gaps(void)
{
unsigned char c;
float a,b;
long i,start,end,gap_start, gap_end;
/* fill in the start with constant value */
for(i=0;i<total_samples;i+=8){
	c=data_present[i>>3];
	if(c){
		while(!((c>>(i&7))&1))i++;
		break;
		}
	}
if(i>=total_samples){
	fprintf(stderr,"** linear_interpolate_gaps: no data found\n");
	exit(-1);
	}
a=data[i];
start=i;
while(i>0){
	i--;
	data[i]=a;
	}
/* fill in the ending with constant value */
for(i=(total_samples-1)&~7;i>0;i-=8){
	c=data_present[i>>3];
	if(c){
		i+=7;
		while(!((c>>(i&7))&1))i--;
		break;
		}
	}
a=data[i];
end=i;
i++;
while(i<total_samples){
	data[i]=a;
	i++;
	}
/* now linear interpolate gaps */

for(i=start;i<end;i++){
	/* skip in bounds of 8 */
	while(!(i&7) && (data_present[i>>3]==0xff) && (i<end))i+=8;
	/* check whether we found a gap */
	if((i<end) && !DATA_PRESENT(i)){
		gap_start=i-1;
		a=data[gap_start];
		while(!DATA_PRESENT(i))i++;
		gap_end=i;
		b=data[gap_end];
		for(i=gap_start+1;i<gap_end;i++)data[i]=(a*(gap_end-i)+b*(i-gap_start))/(gap_end-gap_start);
		i=gap_end;
		}
	}
}

void print_data(char *file, char *tag)
{
FILE *f;
long i;
int count=0;
/* Do we really want debug data ? */
if(file==NULL)return;
fprintf(stderr,"Writing debugging info (tag \"%s\") into file \"%s\"\n",
	tag, file);
f=fopen(file, "w");
while(f==NULL){
	if(count<60){
		sleep(1);
		f=fopen(file, "w");
		count++;
		continue;
		}
	perror(file);
	return;
	}
fprintf(f,"#\n#      Tag: %s\n#\n", tag);
fprintf(f,"#    Channel:  \"%s\"\n#\n", channel);
fprintf(f,"#    GPS start: %lld\n", gps_start);
fprintf(f,"#    GPS end:   %lld\n", gps_end);
fprintf(f,"#    samples_per_second: %ld\n", samples_per_second);
fprintf(f,"#    total_samples:      %ld\n", total_samples);
fprintf(f,"#    version: %s\n", MAKE_SFT_VERSION);
print_settings(f);
fprintf(f,"#\n#\n");
for(i=0;i<total_samples;i++){
	fprintf(f, "%.*g\n", precision, data[i]);
	}
fclose(f);
}

