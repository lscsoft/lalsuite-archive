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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "global.h"
#include "intervals.h"


INTERVAL_SET *new_interval_set(void)
{
INTERVAL_SET *is;
is=do_alloc(1, sizeof(*is));
is->free=0;
is->size=10;
is->i=do_alloc(is->size, sizeof(*is->i));
return is;
}

void free_interval_set(INTERVAL_SET *is)
{
free(is->i);
free(is);
}

void add_interval(INTERVAL_SET *is, INT64 start, INT64 stop)
{
INTERVAL *p;
if(is->free>=is->size){
	is->size+=is->size+10;
	p=do_alloc(is->size, sizeof(*p));
	if(is->free>0)memcpy(p, is->i, is->free*sizeof(*(is->i)));
	free(is->i);
	is->i=p;
	}
is->i[is->free].gps_start=start;
is->i[is->free].gps_stop=stop;
is->free++;
}

void add_intervals_from_file(INTERVAL_SET *is, char *filename)
{
FILE *fin;
char s[2000];
INT64 start, stop;

fprintf(stderr, "Loading segments from file \"%s\"\n", filename);
fin=fopen(filename, "r");
if(fin==NULL) {
	perror("Could not open file");
	return;
	}

while(!feof(fin)){
	fgets(s, 2000, fin);
	sscanf(s, "%lld %lld\n", &start, &stop);
	/* fprintf(stderr, "%lld %lld\n", start, stop); */
	add_interval(is, start, stop);
	}
fclose(fin);
}

/* an easy check, uncomplicated by technicalities (like hash table) */
int check_intervals(INTERVAL_SET *is, INT64 gps)
{
long i;
if(is==NULL)return -1; /* we are not really checking anything */
for(i=0;i<is->free;i++){
	if((gps>=is->i[i].gps_start)&&(gps<is->i[i].gps_stop))return 1;
	}
return 0;
}

