/*
*  Copyright (C) 2007 Vladimir Dergachev
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

/**
 * \file
 * \ingroup pulsarApps
 * \author Vladimir Dergachev
 */

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
