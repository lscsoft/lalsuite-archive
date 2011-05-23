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

#ifndef __UTIL_H__
#define __UTIL_H__

void condor_safe_sleep(int seconds);
int direct_fcntl(int fd, int cmd, void *arg);
int fast_fseek(FILE *stream, long offset, int whence);

/* locate argument number arg in the character line of length length */
void locate_arg(char *line, int length, int arg, int *arg_start, int *arg_stop);
float hann_response(float delta);
void fill_hann_filterN(float *coeffs, int length, int middle, float mismatch);
void fill_hann_filter7(float *coeffs, float mismatch);
void fill_diff_hann_filter7(float *coeffs, float mismatch);
void tabulate_hann_filter7(void);
void tabulated_fill_hann_filter7(float *coeffs, float mismatch);

typedef struct {
	void *data;
	int item_size;
	int size;
	int free;
	} VARRAY;

#define VELT(v, type, i)  (*((type *) &( ( (char *) (v->data) )[v->item_size*(i)])))

//#define VELT(v, type, i)  ( (((type) *)v->data) [(i)] )

VARRAY *new_varray(int item_size);
void free_varray(VARRAY *v);
int varray_add(VARRAY *v, void *item);

int is_round235(int n);
int round235up_int(int n);
int round235down_int(int n);

#endif
