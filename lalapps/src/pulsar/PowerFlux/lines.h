#ifndef __LINES_H__
#define __LINES_H__

/* constants for lines bitmap */
#define LINE_CANDIDATE 	(1<<0)
#define LINE_HIGH	(1<<1)
#define LINE_VERY_HIGH 	(1<<2)
#define LINE_CLUSTERED	(1<<3)
#define LINE_ISOLATED	(1<<4)

#define LINE_YES	(1<<7)

typedef struct {
	int x0;
	int x1;  /* starting and ending bin indices - i.e. window we are analyzing */
	unsigned char *lines; /* bitmap: bit 0 indicates we have a line */
	
	int nlines; /* maximum number of lines to keep */
	int *lines_list; /* list of line numbers */
	
	int nlittle; /* number of high rank values we can safely ignore and still 
		          get good estimate of variance */
	/* statistics, just for debugging and reporting */
	double median;
	double qmost;
	double qlines;
	} LINES_REPORT;

LINES_REPORT *make_lines_report(int x0, int x1, int nlines);
void free_lines_report(LINES_REPORT *lr);
void detect_lines_d(double *z, LINES_REPORT *lr);
void detect_lines_f(float *z, LINES_REPORT *lr);
void detect_background_lines(double *mean);
void print_lines_report(FILE *f,LINES_REPORT *lr,char *tag);

#endif
