#ifndef __INTERVALS_H__
#define __INTERVALS_H__

typedef struct {
	INT64 gps_start;
	INT64 gps_stop;
	} INTERVAL;
	
typedef struct {
	INTERVAL *i;
	long free;
	long size;
	} INTERVAL_SET;

INTERVAL_SET *new_interval_set(void);
void add_interval(INTERVAL_SET *is, INT64 start, INT64 stop);
void add_intervals_from_file(INTERVAL_SET *is, char *filename);
int check_intervals(INTERVAL_SET *is, INT64 gps);


#endif
