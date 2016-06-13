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

#ifndef __RASTERMAGIC_H__
#define __RASTERMAGIC_H__

#include <stdarg.h>

typedef struct {
	int width;
	int height;
	int stride;   /* offset, in bytes, between lines, 
	                  assumed to be the same for all arrays;
			  */
	int step;     /* offset, in bytes, between pixels */
	int flag;    /* flag */
	unsigned char *red; /* pointer to red pixels */
	unsigned char *green; /* etc */
	unsigned char *blue; /* etc */
	} RGBPic;


#define RGBPIC_FLAG_FREE_RED	1
#define RGBPIC_FLAG_FREE_GREEN	1
#define RGBPIC_FLAG_FREE_BLUE	1

#define COLOR(R,G,B)	(((R)<<16)|((G)<<8)|((B)))
#define RED_COLOR(a)	(((a)>>16)&0xff)
#define GREEN_COLOR(a)	(((a)>>8)&0xff)
#define BLUE_COLOR(a)	(((a))&0xff)

#define DRAW_POINT(pic,x,y,color)	{\
		int c0=(color); \
		RGBPic *p0=(pic);\
		int offset; \
		offset=(y)*((p0)->stride)+(x)*((p0)->step); \
		(p0)->red[offset]=RED_COLOR(c0); \
		(p0)->green[offset]=GREEN_COLOR(c0); \
		(p0)->blue[offset]=BLUE_COLOR(c0); \
		}
		
		

RGBPic *make_RGBPic(int width, int height);

void free_RGBPic(RGBPic *p);


void RGBPic_vprintf(RGBPic *p, int x, int y, 
	int fg_color, int bg_color,
	const char *format, va_list ap);
	
void RGBPic_printf(RGBPic *p, int x, int y, 
	int fg_color, int bg_color,
	const char *format, ...);

/* same as above, but prints text vertically, bottom to top */
void RGBPic_vprintf_v(RGBPic *p, int x, int y, 
	int fg_color, int bg_color,
	const char *format, va_list ap);
	
void RGBPic_printf_v(RGBPic *p, int x, int y, 
	int fg_color, int bg_color,
	const char *format, ...);

/* all coordinates are inclusive, i.e. both pixels (x1,y1) and (x2,y2) are painted */
void RGBPic_clear_area(RGBPic *p, int color, 
	int x1, int y1, int x2, int y2);

void RGBPic_draw_line(RGBPic *p, int color, int x1, int y1, int x2, int y2);

void RGBPic_dump_ppm(char *filename, RGBPic *p);
void RGBPic_dump_png(char *filename, RGBPic *p);

typedef struct {
	/* user-settable part */
	  /* display limits */
	double lower_x;
	double upper_x;
	double lower_y;
	double upper_y;
	
	  /* size in pixels, colors, flags */
	int width, height;
	int grid_color;
	int bg_color;
	int fg_color;
	int  logscale_x;
	int  logscale_y;
	
	  /* the following values are computed automatically by draw_grid, 
	     do not touch ! */
	  
	  /* actual upper limits - after rationalization and grid fitting */
	  /* they are never smaller than the ones requested */
	double lx0;
	double lx1;
	double ly0;
	double ly1;
	  /* location of plot */
	int px0;
	int px1;
	int py0;
	int py1;
	} PLOT;

  /* create new plot structure with reasonable defaults */
PLOT *make_plot(int width, int height);
void free_plot(PLOT *plot);
void adjust_plot_limits_f(PLOT *plot, float *x, float *y, int count, int step_x, int step_y, int replace);
void adjust_plot_limits_d(PLOT *plot, double *x, double *y, int count, int step_x, int step_y, int replace);
void draw_grid(RGBPic *p, PLOT *plot, int x, int y);
void draw_points_f(RGBPic *p, PLOT *plot, int color, float *xx, float *yy, int count, int step_x, int step_y);
void draw_points_d(RGBPic *p, PLOT *plot, int color, double *xx, double *yy, int count, int step, int step_y);

typedef struct {
	int ncolors;
	int *color;
	} PALETTE;

PALETTE * make_hue_palette(int ncolors);
void free_palette(PALETTE *p);

typedef struct {
	int x_pixels_per_point;
	int y_pixels_per_point;
	
	double lower_z;
	double upper_z;
	
	int bg_color;
	int fg_color;
	
	/* flags.. */
	int flip_x;   /* flip x axis */
	int flip_y;   /* flip y axis */
	int swap_xy;  /* swap x and y axis */
	int logscale_z;
	
	/* nbands - assign value less than 2 to disable */
	int nbands;
	int dec_band_color;
	
	PALETTE *palette; /* precomputed palette */
	
	/* the following values are computed by the code, 
	   do not touch ! */
	int actual_width;
	int actual_height;
	int key_width;
	int key_height;
	} DENSITY_MAP;

DENSITY_MAP *make_density_map(int x_ppp, int y_ppp, int nbands);
void free_density_map(DENSITY_MAP *dm);
void adjust_density_map_limits_f(DENSITY_MAP *plot, float *z, int count, int step, int replace);
void adjust_density_map_limits_d(DENSITY_MAP *dm, double *z, int count, int step, int replace);
void draw_density_map_key(RGBPic *p, DENSITY_MAP *dm, int x, int y);
void draw_density_map_f(RGBPic *p, int x, int y, DENSITY_MAP *dm, float *z, int x_count, int y_count, int step_x, int step_y);
void draw_density_map_d(RGBPic *p, int x, int y, DENSITY_MAP *dm, double *z, int x_count, int y_count, int step_x, int step_y);

void plot_single_density_map_f(RGBPic *p, DENSITY_MAP *dm, 
	float *z, int x_count, int y_count, int step_x, int step_y);
void draw_density_map_d(RGBPic *p, int x, int y, 
	DENSITY_MAP *dm, 
	double *z, int x_count, int y_count, int step_x, int step_y);
void plot_single_density_map_d(RGBPic *p, DENSITY_MAP *dm, 
	double *z, int x_count, int y_count, int step_x, int step_y);

#include "grid.h"
void plot_grid_f(RGBPic *p, SKY_GRID *grid, float *z, int step);
void plot_grid_d(RGBPic *p, SKY_GRID *grid, double *z, int step);

#endif
