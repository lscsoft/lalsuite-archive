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

fin=fopen(filename, "r");
if(fin==NULL)return;

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

