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
void detect_lines_d(double *z, LINES_REPORT *lr, char *tag);
void detect_lines_f(float *z, LINES_REPORT *lr, char *tag);
void print_lines_report(FILE *f,LINES_REPORT *lr, char *tag);

#endif
