#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>

long starting_bin, stopping_bin;

unsigned char s[2000];

long long gps_start, gps_stop;
long file_starting_bin, file_stopping_bin;

typedef float REAL4;
REAL4 *data=NULL;
long data_size=0;

void * do_alloc(long a, long b)
{
void *r;
r=calloc(a,b);
while(r==NULL){
	fprintf(stderr,"Could not allocate %ld units of %ld bytes each (%ld bytes total)\n", a, b, a*b);
	sleep(1);
	r=calloc(a,b);
	}
return r;
}

int c_r4(REAL4 *a, REAL4 *b)
{
if(*a<*b)return -1;
if(*a>*b)return 1;
return 0;
}

/* this function modifies numbers in the array */
REAL4 sum_positive_r4(REAL4 *array, long count)
{
int i,c;
qsort(array, count, sizeof(REAL4), c_r4);
c=count;
while(c>1){
	for(i=0;i<c;i+=2){
		array[i>>1]=array[i]+array[i+1];
		}
	if(c & 1){
		array[(c>>1)]=array[c-1];
		c=(c>>1)+1;
		} else {
		c=(c>>1);
		}
	}
return array[0];
}

void compute_statistics(REAL4 *array, long count, REAL4 *average, REAL4 *median)
{

int i,c;
qsort(array, count, sizeof(REAL4), c_r4);
if(count & 1)*median=array[count>>1];
	else *median=(array[count>>1]+array[(count>>1)-1])/2.0;
c=count;
while(c>1){
	for(i=0;i<c;i+=2){
		array[i>>1]=array[i]+array[i+1];
		}
	if(c & 1){
		array[(c>>1)]=array[c-1];
		c=(c>>1)+1;
		} else {
		c=(c>>1);
		}
	}
*average=array[0]/count;
}

void process_tag(unsigned char *s)
{
unsigned char *ptag, *pval;
ptag=s;
while((*ptag==' ')||(*ptag=='\t'))ptag++;
pval=ptag;
while(*pval && (*pval!=':'))pval++;
if(*pval==':')pval++;
if(!strncmp(ptag, "GPS start", 9)){
	sscanf(pval, "%lld", &gps_start);	
	} else
if(!strncmp(ptag, "GPS end", 7)){
	sscanf(pval, "%lld", &gps_stop);
	} else
if(!strncmp(ptag, "starting bin", 12)){
	sscanf(pval, "%ld", &file_starting_bin);
	} else
if(!strncmp(ptag, "stopping bin", 12)){
	sscanf(pval, "%ld", &file_stopping_bin);
	}
}

int main(int argc, char *argv[])
{
int i;
FILE *f=NULL;
double start_freq, stop_freq;
long BigE;
unsigned char *endianness=(unsigned char *)&BigE;
REAL4 average, median;

BigE=('B'<<24)|('i'<<16)|('g'<<8)|('E');
if(argc<4){
	fprintf(stderr, "Usage: %s start_freq stop_freq file1 file2 ...\n", argv[0]);
	exit(-1);
	}
sscanf(argv[1], "%lg", &start_freq);
sscanf(argv[2], "%lg", &stop_freq);
if(stop_freq<start_freq){
	fprintf(stderr, "Stop frequency must be greater or equal to the start frequency\n");
	exit(-1);
	}
printf("# %d files\n", argc-3);
printf("# Start frequency %lg\n", start_freq);
printf("# Stop frequency %lg\n", stop_freq);
printf("#Gps_time average median\n");
for(i=3;i<argc;i++){
	gps_start=-1;
	gps_stop=-1;
	file_starting_bin=-1;
	file_stopping_bin=-1;
	f=fopen(argv[i], "r");
	if(f==NULL){
		fprintf(stderr, "Could not open file \"%s\" for reading:", argv[i]);
		perror("");
		continue;
		}
	while(!feof(f)){
		fgets(s, 2000, f);
		if(!s[0])continue;
		if(s[0]=='#'){
			process_tag(s+1);
			continue;
			} else
		if(!strncmp(s, "binary", 6)){
			fgets(s, 2000, f);
			if(strncmp(s, endianness, 4)){
				fprintf(stderr, "Could not process file \"%s\" because of different endianness.\n", argv[i]);
				break;
				}
			if((gps_start<0)||(gps_stop<=0)||(file_starting_bin<0)||(file_stopping_bin<0)){
				fprintf(stderr, "Could not process file \"%s\" because of some of the required tags are missing.\n", argv[i]);
				break;
				}
			starting_bin=lrint(start_freq*(gps_stop-gps_start));
			stopping_bin=lrint(stop_freq*(gps_stop-gps_start));
			if(stopping_bin==starting_bin)stopping_bin=starting_bin+1;
			if((starting_bin < file_starting_bin)||(starting_bin>=file_stopping_bin)||
			   (stopping_bin > file_stopping_bin)||(stopping_bin<=starting_bin)){
				fprintf(stderr, "Could not process file \"%s\" because requested data segment is not present.\n", argv[i]);
				fprintf(stderr, "\tstarting_bin=%ld stopping_bin=%ld file_starting_bin=%ld file_stopping_bin=%ld\n", 
					starting_bin, stopping_bin, file_starting_bin, file_stopping_bin);
				break;
			   	}
			if(starting_bin != file_starting_bin)
				if(fseek(f, (starting_bin-file_starting_bin)*sizeof(REAL4), SEEK_CUR)<0){
					fprintf(stderr, "Unable to seek to data location in file \"%s\":", argv[i]);
					perror("");
					break;
					}
			if((stopping_bin-starting_bin) > data_size){
				if(data!=NULL)free(data);
				data_size=stopping_bin-starting_bin;
				data=do_alloc(data_size, sizeof(*data));
				}
			if(fread(data, sizeof(REAL4), stopping_bin-starting_bin, f)<(stopping_bin-starting_bin)){
				fprintf(stderr, "Unable to read data from file \"%s\":", argv[i]);
				perror("");
				break;
				}
			compute_statistics(data, stopping_bin-starting_bin, &average, &median);
			printf("%lld %g %g\n",gps_start, average, median);
			break;
			} else {
			fprintf(stderr, "Could not process file \"%s\" because it is not in binary format.\n", argv[i]);
			break;
			}
		}
	fclose(f);
	}
fflush(stdout);
return 0;
}
