#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>

#define TESTIO(a)   { \
	if((a)<0){ \
		perror(""); \
		exit(-1); \
		} \
	}

typedef double REAL8;
typedef float REAL4;

typedef struct {
	REAL8 key;
	int gps;
	int nsec;
	REAL8 timebase;
	int bin_start;
	int nbins;
	} HEADER;

int main(int argc, char *argv[])
{
int fd;
unsigned char *p;
REAL4 *data;
struct stat buf;
HEADER *header;
long i;

if(argc<2){
	fprintf(stderr,"Usage: %s filename\n", argv[0]);
	exit(-1);
	}

TESTIO(fd=open(argv[1], O_RDWR));
TESTIO(fstat(fd, &buf));
p=mmap(NULL, buf.st_size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
TESTIO(p==NULL?-1:0);
close(fd);

header=p;
data=p+sizeof(*header);

fprintf(stderr,"key=%g\n", header->key);
fprintf(stderr,"gps=%d\n", header->gps);
fprintf(stderr,"nsec=%d\n", header->nsec);
fprintf(stderr,"timebase=%g\n", header->timebase);
fprintf(stderr,"bin_start=%d\n", header->nsec);
fprintf(stderr,"nbins=%d\n", header->nbins);

/* Modify file here : */

if(header->key!=1.0){
	fprintf(stderr,"file %s has unknown key\n", argv[1]);
	exit(-3);
	}

if(header->timebase >0) {
	fprintf(stderr,"file %s is correct - leaving alone\n", argv[1]);
	exit(-2);
	}
fprintf(stderr,"Adjusting file %s\n", argv[1]);

header->timebase=-header->timebase;
return 0;
(volatile)header->key=1.5;
for(i=0;i<header->nbins*2;i++)
	data[i]=data[i]*(header->nbins*1.0)/(0.5*1800.0*16384.0);
(volatile)header->key=1.0;
return 0;
}
