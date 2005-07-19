/*

   Miscellaneious handy functions

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

#include "make_sft_global.h"

extern char * channel;
extern INT64 gps_start, gps_end;
extern long samples_per_second;
extern long total_samples;

int precision=6;

struct {
	int type;
	char *name; 
	} windows[]={
	{Hann, "Hann"},
	{Welch, "Welch"},
	{Bartlett, "Bartlett"},
	/* Default is Hann */
	{Hann, NULL}
	};


void *do_alloc(long a, long b)
{
void *p;
int count=0;
p=calloc(a,b);
while(p==NULL){
	fprintf(stderr,"** Failed to allocated %ld chunks of %ld bytes (%ld bytes total)\n", a, b, a*b);
	#ifdef __EXIT_ON_LOW_MEMORY__
	exit(-1);
	#endif
	sleep(10);
	p=calloc(a,b);
	count++;
	if(count>100)exit(-1);
	}
return p;
}

int c_r8(REAL8 *a, REAL8 *b)
{
if(*a<*b)return -1;
if(*a>*b)return 1;
return 0;
}

REAL8 sum_r8_squares(REAL8 *array, long count)
{
int i,c;
REAL8 *squares, a;
squares=do_alloc(count, sizeof(*squares));
for(i=0;i<count;i++)squares[i]=array[i]*array[i];
qsort(squares, count, sizeof(*squares), c_r8);
c=count;
while(c>1){
	for(i=0;(i+2)<=c;i+=2){
		squares[i>>1]=squares[i]+squares[i+1];
		}
	if(c & 1){
		squares[(c>>1)]=squares[c-1];
		c=(c>>1)+1;
		} else {
		c=(c>>1);
		}
	}
a=squares[0];
free(squares);
return a;
}

int c_r4(REAL4 *a, REAL4 *b)
{
if(*a<*b)return -1;
if(*a>*b)return 1;
return 0;
}


REAL4 sum_r4_squares(REAL4 *array, long count)
{
int i,c;
REAL4 *squares, a;
squares=do_alloc(count, sizeof(REAL4));
for(i=0;i<count;i++)squares[i]=array[i]*array[i];
qsort(squares, count, sizeof(REAL4), c_r4);
c=count;
while(c>1){
	for(i=0;(i+2)<=c;i+=2){
		squares[i>>1]=squares[i]+squares[i+1];
		}
	if(c & 1){
		squares[(c>>1)]=squares[c-1];
		c=(c>>1)+1;
		} else {
		c=(c>>1);
		}
	}
a=squares[0];
free(squares);
return a;
}

REAL4 sum_positive_r4(REAL4 *array, long count)
{
int i,c;
REAL4 *backup, a;
backup=do_alloc(count, sizeof(REAL4));
memcpy(backup, array, count*sizeof(REAL4));
qsort(backup, count, sizeof(REAL4), c_r4);
c=count;
while(c>1){
	for(i=0;(i+2)<=c;i+=2){
		backup[i>>1]=backup[i]+backup[i+1];
		}
	if(c & 1){
		backup[(c>>1)]=backup[c-1];
		c=(c>>1)+1;
		} else {
		c=(c>>1);
		}
	}
a=backup[0];
free(backup);
return a;
}

long count_bitmap_bits(unsigned char *bitmap, long length)
{
long i;
long bit_count=0;
unsigned char c;
int k;
for(i=0;i<((length>>3));i++){
	c=bitmap[i];
	if(!c)continue;
	if(c==0xff){
		bit_count+=8;
		continue;
		}
	for(k=0;k<8;k++)if(c & (1<<k))bit_count++;	
	}
return bit_count;
}


/* When we try to access a file it may not be present - automounter
   might not be ready right away.

   Thus attempt to access a file until it is there.
   
   400=5*60+100.
*/

void fault_file(char *filename)
{
struct stat buf;
int count=0;
while(stat(filename, &buf)<0){
	count++;
	if(count==1){
		fprintf(stderr,"Did not find file \"%s\", will continue trying for 400 seconds\n", filename);
		}
	if(count>400){
		fprintf(stderr,"** Could not establish existence of file \"%s\"\n", filename);
		exit(-1);	
		}
	sleep(1);
	}
if(count>0){
	fprintf(stderr, "Found file \"%s\" after %d attempts\n", filename, count);
	}
}

/* insert comments so that the header ends on a 512-12 byte boundary (12 bytes are using for binary+EGIB signature) */
void align_file(FILE *file)
{
long i;
char padding[512];
i=ftell(file);
if(i<0){
	perror("align_file");
	return;
	}
i=i % 500;
i=499-i;
if(!i)return;
memset(padding,'-',512);
memcpy(padding, "# PADDING", 9);
padding[i-1]='\n';
padding[i]=0;
fprintf(file,"%s",padding);
}

void print_REAL4Vector(REAL4Vector *vec, char *file, char *tag, int mode, long bin_start, long bin_stop)
{
FILE *f;
long i;
long x;
unsigned char *s=(unsigned char *)&x;
int count=0;
/* Do we really want debug data ? */
if(file==NULL)return;
fprintf(stderr,"Writing REAL4Vector data (tag \"%s\") into file \"%s\"\n",
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

if(bin_start<0)bin_start=0;
if(bin_stop>vec->length)bin_stop=vec->length;

fprintf(f,"#\n#      Tag: %s\n#\n", tag);
fprintf(f,"#    Channel:  \"%s\"\n#\n", channel);
fprintf(f,"#    GPS start: %lld\n", gps_start);
fprintf(f,"#    GPS end:   %lld\n", gps_end);
fprintf(f,"#    samples_per_second: %ld\n", samples_per_second);
fprintf(f,"#    total_samples:      %ld\n", total_samples);
fprintf(f,"#    version: %s\n", MAKE_SFT_VERSION);
print_settings(f);
fprintf(f,"#\n#\n");
switch(mode){
	case OUTPUT_MODE_TEXT:
		for(i=bin_start;i+4<=bin_stop;i+=4){
			fprintf(f, "%.*g\n%.*1$g\n%.*1$g\n%.*1$g\n", 
				precision, vec->data[i],
				vec->data[i+1],
				vec->data[i+2],
				vec->data[i+3]);
			}
		for(;i<bin_stop;i++){
			fprintf(f, "%.*g\n", precision, vec->data[i]);
			}
		break;
	case OUTPUT_MODE_GEO: /* GEO format exists only for SFTs */
	case OUTPUT_MODE_BINARY:
		align_file(f);
		fprintf(f, "binary\n");
		x=('B'<<24)|('i'<<16)|('g'<<8)|('E');
		fwrite(s, 1, 4, f);
		fprintf(f,"\n");
		fwrite(&(vec->data[bin_start]), sizeof(*(vec->data)), bin_stop-bin_start, f);
		break;
	default:
		fprintf(stderr, "Don't know how to write using mode %d, exiting.\n", mode);
		exit(-1);
	}
fclose(f);
}

void print_COMPLEX8Vector(COMPLEX8Vector *vec, char *file, char *tag, int mode, long bin_start, long bin_stop)
{
FILE *f;
long i;
long x;
unsigned char *s=(unsigned char *)&x;
int count=0;
INT4 int4;
REAL8 real8;
REAL4 a,factor;
/* Do we really want debug data ? */
if(file==NULL)return;
fprintf(stderr,"Writing COMPLEX8Vector data (tag \"%s\") into file \"%s\"\n",
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

if(bin_start<0)bin_start=0;
if(bin_stop>vec->length)bin_stop=vec->length;

/* header first */

switch(mode){
	case OUTPUT_MODE_TEXT:
	case OUTPUT_MODE_BINARY:
		fprintf(f,"#\n#      Tag: %s\n#\n", tag);
		fprintf(f,"#    Channel:  \"%s\"\n#\n", channel);
		fprintf(f,"#    GPS start: %lld\n", gps_start);
		fprintf(f,"#    GPS end:   %lld\n", gps_end);
		fprintf(f,"#    samples_per_second: %ld\n", samples_per_second);
		fprintf(f,"#    total_samples:      %ld\n", total_samples);
		fprintf(f,"#    version: %s\n", MAKE_SFT_VERSION);
		print_settings(f);
		fprintf(f,"#\n#\n");
		break;
	case OUTPUT_MODE_GEO:
		/* endian */
		real8=1;
		fwrite(&real8, sizeof(real8), 1, f);
		/* gps_sec */
		int4=gps_start;
		fwrite(&int4, sizeof(int4), 1, f);
		/* gps_nsec - always 0 in this program */
		int4=0;
		fwrite(&int4, sizeof(int4), 1, f);
		/* timebase */
		real8=gps_end-gps_start;
		fwrite(&real8, sizeof(real8), 1, f);
		/* firstfreqindex */
		int4=bin_start;
		fwrite(&int4, sizeof(int4), 1, f);
		/* nsamples */
		int4=bin_stop-bin_start;
		fwrite(&int4, sizeof(int4), 1, f);
		break;
	default:
		fprintf(stderr, "Don't know how to write using mode %d, exiting.\n", mode);
		exit(-1);
	}
		
switch(mode){
	case OUTPUT_MODE_TEXT:
		for(i=bin_start;i+2<=bin_stop;i+=2){
			fprintf(f, "%.*g %.*1$g\n%.*1$g %.*1$g\n", 
				precision, vec->data[i].re,
				vec->data[i].im,
				vec->data[i+1].re,
				vec->data[i+1].im);
			}
		for(;i<bin_stop;i++){
			fprintf(f, "%.*g %.*1g\n", precision, vec->data[i].re, vec->data[i].im);
			}
		break;
	case OUTPUT_MODE_BINARY:
		align_file(f);
		fprintf(f, "binary\n");
		x=('B'<<24)|('i'<<16)|('g'<<8)|('E');
		fwrite(s, 1, 4, f);
		fprintf(f,"\n");
		fwrite(&(vec->data[bin_start]), sizeof(*(vec->data)), bin_stop-bin_start, f);
		break;
	case OUTPUT_MODE_GEO: /* GEO format exists only for SFTs */
		/* we need to renormalize data properly */
		factor=(1.0*(bin_stop-bin_start))/(total_samples/2+1);
		for(i=bin_start;i<bin_stop;i++){
			a=vec->data[i].re*factor;
			fwrite(&(a), sizeof(a), 1, f);
			a=vec->data[i].im*factor;
			fwrite(&(a), sizeof(a), 1, f);
			}
		break;
	default:
		fprintf(stderr, "Don't know how to write using mode %d, exiting.\n", mode);
		exit(-1);
	}
fclose(f);
}


char *dup_quoted_name(char *p)
{
char *r;
/* find opening doublequote and copy name starting after it */
while(*p && *p!='"')p++;
if(*p)p++;
r=strdup(p);
p=r;
/* find ending doublequote and zero it, use backslash for escaping */
while(*p && *p!='"'){
	if(*p=='\\')p++;
	if(*p)p++;
	}
*p=0;
return r;
}

void print_window_settings(FILE *f, char *window_name, MAKE_SFT_WINDOW *window)
{
int i;
fprintf(f,"#   %s window type: ", window_name);
for(i=0;windows[i].name!=NULL;i++)
	if(windows[i].type==window->type){
		fprintf(f, "%s\n", windows[i].name);
		break;
		}
if(windows[i].name==NULL)fprintf(f, "**UNKNOWN**\n");
fprintf(f,"#   %s window tail size: %ld\n", window_name, window->tail_size);
}

int translate_window_type(char *text)
{
int i;
char *type;
type=dup_quoted_name(text);
for(i=0;windows[i].name!=NULL;i++){
	if(!strcasecmp(type, windows[i].name))break;
	}
free(type);
return windows[i].type;
}
