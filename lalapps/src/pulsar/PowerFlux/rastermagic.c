#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include "rastermagic.h"
#include "global.h"
#include "cmdline.h"

extern char *output_dir;
extern struct gengetopt_args_info args_info;

void *do_alloc(long a, long b);

/* include bitmaps directly - no reason to pollute namespace */
#include "font_8x16.c"

RGBPic *make_RGBPic(long width, long height)
{
RGBPic *p;

p=do_alloc(1, sizeof(*p));
p->width=width;
p->height=height;
p->stride=width;
p->step=1;
p->red=do_alloc(width*height,1);
p->green=do_alloc(width*height,1);
p->blue=do_alloc(width*height,1);
p->flag=RGBPIC_FLAG_FREE_RED|RGBPIC_FLAG_FREE_GREEN|RGBPIC_FLAG_FREE_BLUE;
return p;
}


void free_RGBPic(RGBPic *p)
{
if(p->flag & RGBPIC_FLAG_FREE_RED){
	free(p->red);
	p->red=NULL;
	}
if(p->flag & RGBPIC_FLAG_FREE_GREEN){
	free(p->green);
	p->green=NULL;
	}
if(p->flag & RGBPIC_FLAG_FREE_BLUE){
	free(p->blue);
	p->blue=NULL;
	}
free(p);
}

void RGBPic_vprintf(RGBPic *p, long x, long y, 
	long fg_color, long bg_color,
	const char *format, va_list ap)
{
unsigned char *s,*q,*dst;
long size=1024,count;
long z,i,j,step;
unsigned char c1,c2,line,line2;

/* easy sanity checking first */
if(y<0)return;
if((y+8)>=p->height)return;

s=do_alloc(size,1);
while((count=vsnprintf(s, size, format, ap))>size){
	free(s);
	size=count+1;
	s=do_alloc(size,1);
	}
/* go through characters one by one */
z=x;
step=p->step;
for(q=s;*q;q++){
	/* handle newlines */
	if(*q=='\n'){
		z=x;
		y+=16;
		if((y+16)>=p->height){
			free(s);
			return;
			}
		continue;
		}
	/* don't print if the character cannot be completely fitted */
	if((z+8)>=p->width)continue; /* we continue since we might have a newline
	                                 which resets z */
	for(i=0;i<16;i++){
		/* get font bitmap */
		line=fontdata_8x16[(*q)*16+i];
		#define PAINT_LINE(a,fc,bc)	{\
			line2=line; \
			/* get first pixel to paint */ \
			dst=p->a+(y+i)*p->stride+z*p->step; \
			/* get background and foreground colors */ \
			c1=fc; \
			c2=bc; \
			for(j=0;j<8;j++){ \
				if(line2 & 0x80){ \
					*dst=c1; \
					} else { \
					*dst=c2; \
					} \
				dst+=step; \
				line2<<=1; \
				} \
			}
		PAINT_LINE(red, RED_COLOR(fg_color), RED_COLOR(bg_color));
		PAINT_LINE(green, GREEN_COLOR(fg_color), GREEN_COLOR(bg_color));
		PAINT_LINE(blue, BLUE_COLOR(fg_color), BLUE_COLOR(bg_color));
		#undef PAINT_LINE
		}
	z+=8;
	}
free(s);
}

void RGBPic_printf(RGBPic *p, long x, long y, 
	long fg_color, long bg_color,
	const char *format, ...)
{
va_list ap;
va_start(ap, format);
RGBPic_vprintf(p, x, y, 
	fg_color, bg_color,
	format, ap);
va_end(ap);
}

void RGBPic_vprintf_v(RGBPic *p, long x, long y, 
	long fg_color, long bg_color,
	const char *format, va_list ap)
{
unsigned char *s,*q,*dst;
long size=1024,count;
long z,i,j;
unsigned char c1,c2,line,line2;

/* easy sanity checking first */
if(y<16)return;
if((y+8)>=p->height)return;

s=do_alloc(size,1);
while((count=vsnprintf(s, size, format, ap))>size){
	free(s);
	size=count+1;
	s=do_alloc(size,1);
	}
/* go through characters one by one */
z=y;
for(q=s;*q;q++){
	/* handle newlines */
	if(*q=='\n'){
		z=y;
		x+=16;
		if((x+16)>=p->height){
			free(s);
			return;
			}
		continue;
		}
	/* don't print if the character cannot be complitely fitted */
	if(z<8)continue; /* we continue since we might have a newline
	                                 which resets z */
	for(i=0;i<16;i++){
		/* get font bitmap */
		line=fontdata_8x16[(*q)*16+i];
		#define PAINT_LINE(a,fc,bc)	{\
			line2=line; \
			/* get first pixel to paint */		 \
			dst=p->a+(z)*p->stride+(x+i)*p->step; \
			/* get background and foreground colors */ \
			c1=fc; \
			c2=bc; \
			for(j=0;j<8;j++){ \
				if(line2 & 0x80){ \
					*dst=c1; \
					} else { \
					*dst=c2; \
					} \
				dst-=p->stride; \
				line2<<=1; \
				} \
			}
		PAINT_LINE(red, RED_COLOR(fg_color), RED_COLOR(bg_color));
		PAINT_LINE(green, GREEN_COLOR(fg_color), GREEN_COLOR(bg_color));
		PAINT_LINE(blue, BLUE_COLOR(fg_color), BLUE_COLOR(bg_color));
		#undef PAINT_LINE
		}
	z-=8;
	}
free(s);
}

void RGBPic_printf_v(RGBPic *p, long x, long y, 
	long fg_color, long bg_color,
	const char *format, ...)
{
va_list ap;
va_start(ap, format);
RGBPic_vprintf_v(p, x, y, 
	fg_color, bg_color,
	format, ap);
va_end(ap);
}

void RGBPic_dump_ppm(char *filename, RGBPic *p)
{
FILE *fout;
long i,j;
char s[20000];

snprintf(s,20000,"%s%s", output_dir, filename);
fout=fopen(s, "wb");
if(fout==NULL){
	fprintf(stderr,"Error dumping %ldx%ld picture to \"%s\" in PPM format:",
		p->width, p->height, filename);
	perror("");
	return;
	}
fprintf(fout,"P6\n%ld\n%ld\n255\n", p->width, p->height);
for(i=0;i<p->height;i++)
	for(j=0;j<p->width;j++){
		fputc(p->red[i*p->stride+j*p->step],fout);
		fputc(p->green[i*p->stride+j*p->step],fout);
		fputc(p->blue[i*p->stride+j*p->step],fout);
		}
fclose(fout);
}

#if 0
void RGBPic_dump_png(char *filename, RGBPic *p)
{
char *tmpfile;
char s[20000];
char *old_output_dir=output_dir;
tmpfile=tmpnam(NULL);
output_dir="";
RGBPic_dump_ppm(tmpfile, p);
output_dir=old_output_dir;
memset(s,0,20000);
snprintf(s, 19999, "%s \"%s\" > \"%s%s\" ; rm -f \"%s\"",
	args_info.pnmtopng_arg, tmpfile, output_dir, filename, tmpfile);
fprintf(stderr,"Executing: %s\n", s);
system(s);
}
#endif

void RGBPic_clear_area(RGBPic *p, long color, long x1, long y1, long x2, long y2)
{
unsigned char c;
long i;long j;

/* check dimensions */
if(x1<0)x1=0;
if(y1<0)y1=0;
if(x2>=p->width)x2=p->width-1;
if(y2>=p->height)y2=p->height-1;
if(x1>x2)return;
if(y1>y2)return;

c=RED_COLOR(color);
for(i=y1;i<=y2;i++)
	for(j=x1;j<=x2;j++){
		p->red[i*p->stride+j*p->step]=c;
		}
c=GREEN_COLOR(color);
for(i=y1;i<=y2;i++)
	for(j=x1;j<=x2;j++){
		p->green[i*p->stride+j*p->step]=c;
		}
c=BLUE_COLOR(color);
for(i=y1;i<=y2;i++)
	for(j=x1;j<=x2;j++){
		p->blue[i*p->stride+j*p->step]=c;
		}
}

/* Bresenham algorithm */
void RGBPic_draw_line(RGBPic *p, long color, long x1, long y1, long x2, long y2)
{
long x,y,dx,dy,adx,ady,eps;
/* special case vertical lines */
if(x1==x2){
	if(y2<y1){
		y=y2;
		y2=y1;
		y1=y;
		}
	for(y=y1;y<=y2;y++)DRAW_POINT(p,x1,y,color);
	return;
	}
if(x1>x2){
	x=x2;
	x2=x1;
	x1=x;
	
	y=y2;
	y2=y1;
	y1=y;
	}
/* special case horizontal lines */
if(y1==y2){
	for(x=x1;x<=x2;x++)DRAW_POINT(p,x,y1,color);
	return;
	}
dx=x2-x1;
dy=y2-y1;
if(dy<0)ady=-dy;
	else ady=dy;
eps=0;
if(dx>ady){
	y=y1;
	for(x=x1;x<=x2;x++){
		DRAW_POINT(p,x,y,color);
		eps+=ady;
		if((eps<<1)>=dx){
			if(dy<0)y--;
				else y++;
			eps-=dx;
			}
		}
	} else {
	if(dy<0){
		x=x2;
		x2=x1;
		x1=x;
		
		y=y2;
		y2=y1;
		y1=y;		
		}
	x=x1;
	dx=x2-x1;
	dy=y2-y1;
	if(dx<0)adx=-dx;
		else adx=dx;
	for(y=y1;y<=y2;y++){
		DRAW_POINT(p,x,y,color);
		eps+=adx;
		if((eps<<1)>=dy){
			if(dx<0)x--;
				else x++;
			eps-=dy;
			}
		}
	}
}

PLOT *make_plot(long width, long height)
{
PLOT *r;
r=do_alloc(1,sizeof(*r));
r->width=width;
r->height=height;
r->grid_color=COLOR(200,200,200);
r->bg_color=COLOR(255,255,255);
r->fg_color=COLOR(0,0,0);
r->logscale_x=0;
r->logscale_y=0;
return r;
}

void free_plot(PLOT *plot)
{
free(plot);
}

void adjust_plot_limits_f(PLOT *plot, float *x, float *y, int count, int step_x, int step_y, int replace)
{
float x0,x1,y0,y1;
int i;
if(count<1)return;
x0=*x;
x1=*x;
y0=*y;
y1=*y;

for(i=1;i<count;i++){
	if(x0>x[i*step_x])x0=x[i*step_x];
	if(x1<x[i*step_x])x1=x[i*step_x];
	if(y0>y[i*step_y])y0=y[i*step_y];
	if(y1<y[i*step_y])y1=y[i*step_y];
	}
if(replace){
	plot->lower_x=x0;
	plot->upper_x=x1;
	plot->lower_y=y0;
	plot->upper_y=y1;
	} else {
	if(plot->lower_x>x0)plot->lower_x=x0;
	if(plot->upper_x<x1)plot->upper_x=x1;
	if(plot->lower_y>y0)plot->lower_y=y0;
	if(plot->upper_y>y1)plot->upper_y=y1;
	}
}

void adjust_plot_limits_d(PLOT *plot, double *x, double *y, int count, int step_x, int step_y, int replace)
{
double x0,x1,y0,y1;
int i;
if(count<1)return;
x0=*x;
x1=*x;
y0=*y;
y1=*y;

for(i=1;i<count;i++){
	if(x0>x[i*step_x])x0=x[i*step_x];
	if(x1<x[i*step_x])x1=x[i*step_x];
	if(y0>y[i*step_y])y0=y[i*step_y];
	if(y1<y[i*step_y])y1=y[i*step_y];
	}
if(replace){
	plot->lower_x=x0;
	plot->upper_x=x1;
	plot->lower_y=y0;
	plot->upper_y=y1;
	} else {
	if(plot->lower_x>x0)plot->lower_x=x0;
	if(plot->upper_x<x1)plot->upper_x=x1;
	if(plot->lower_y>y0)plot->lower_y=y0;
	if(plot->upper_y>y1)plot->upper_y=y1;
	}
}

#define exp10(a)	exp(M_LN10*(a))

void draw_grid(RGBPic *p, PLOT *plot, long x, long y)
{
double a,dx,dy,dsx,dsy,e10x,e10y;
long i,j,k;
char s[200];

if(plot->lower_x > plot->upper_x){
	a=plot->upper_x;
	plot->upper_x=plot->lower_x;
	plot->lower_x=a;
	}
if(plot->lower_y > plot->upper_y){
	a=plot->upper_y;
	plot->upper_y=plot->lower_y;
	plot->lower_y=a;
	}

/* the height of the letters is 16 pixels */
plot->px0=x+24;
plot->py0=x+3;
plot->px1=y+plot->width-4;
plot->py1=y+plot->height-25;

dx=plot->upper_x-plot->lower_x;
dy=plot->upper_y-plot->lower_y;

e10x=exp10(floor(log10(dx)));
e10y=exp10(floor(log10(dy)));

plot->lx0=floor(plot->lower_x*2/e10x)*0.5*e10x;
plot->lx1=ceil(plot->upper_x*2/e10x)*0.5*e10x;

plot->ly0=floor(plot->lower_y*2/e10y)*0.5*e10y;
plot->ly1=ceil(plot->upper_y*2/e10y)*0.5*e10y;

dx=plot->lx1-plot->lx0;
dy=plot->ly1-plot->ly0;

dsx=(48*dx)/(plot->px1-plot->px0+1);
a=exp10(floor(log10(fabs(dsx))));
dsx=ceil(dsx/(5*a))*(5*a);
if(e10x<a)e10x=a;

dsy=(48*dy)/(plot->py1-plot->py0+1);
a=exp10(floor(log10(fabs(dsy))));
dsy=ceil(dsy/(5*a))*(5*a);
if(e10y<a)e10y=a;

if((e10x>0.05)&&(e10x<20))e10x=1;
if((e10y>0.05)&&(e10y<20))e10y=1;

RGBPic_clear_area(p,plot->bg_color,x,y,x+plot->width-1,y+plot->height-1);
RGBPic_draw_line(p,plot->fg_color,plot->px0-1,plot->py0-1,plot->px1+1,plot->py0-1);
RGBPic_draw_line(p,plot->fg_color,plot->px0-1,plot->py1+1,plot->px1+1,plot->py1+1);
RGBPic_draw_line(p,plot->fg_color,plot->px0-1,plot->py0-1,plot->px0-1,plot->py1+1);
RGBPic_draw_line(p,plot->fg_color,plot->px1+1,plot->py0-1,plot->px1+1,plot->py1+1);
/* draw x marks */
//TRACE("x marks")
for(a=rint(plot->lx0/dsx)*dsx; a<=plot->lx1;a+=dsx){
	/* special case a=0.0 or we run into precision problems */
	if(fabs(a)<0.5*dsx)a=0.0;
	j=rint(((a-plot->lx0)*(plot->px1-plot->px0+1))/dx);
	k=sprintf(s,"%g",a/e10x);
	if((j>=32)&&(j<(plot->px1-plot->px0-64))){
		RGBPic_draw_line(p,plot->grid_color,plot->px0+j,plot->py0,plot->px0+j,plot->py1);
		RGBPic_draw_line(p,plot->fg_color,plot->px0+j,plot->py1+1,plot->px0+j,plot->py1+4);
		RGBPic_printf(p,plot->px0+j-4*k,plot->py1+5,plot->fg_color, plot->bg_color, s);
		}
	}
if(e10x!=1){
	k=sprintf(s,"*%4.0g", e10x);
	RGBPic_printf(p,x+plot->width-1-8*k,plot->py1+5,plot->fg_color, plot->bg_color, s);
	}
/* draw y marks */
//TRACE("y marks")
for(a=rint(plot->ly0/dsy)*dsy; a<=plot->ly1;a+=dsy){
	/* special case a=0.0 or we run into precision problems */
	if(fabs(a)<0.5*dsy)a=0.0;
	j=rint(((plot->ly1-a)*(plot->py1-plot->py0+1))/dy);	
	k=sprintf(s,"%g",a/e10y);
	
	if((j>=64)&&(j<(plot->py1-plot->py0-32))){
		RGBPic_draw_line(p,plot->grid_color,plot->px0,plot->py0+j,plot->px1,plot->py0+j);
		RGBPic_draw_line(p,plot->fg_color,x+19,plot->py0+j,x+23,plot->py0+j);
		RGBPic_printf_v(p,x,plot->py0+j+4*k,plot->fg_color, plot->bg_color, s);
		}
	}
if(e10y!=1){
	k=sprintf(s,"*%4.0g", e10y);
	RGBPic_printf_v(p,x,y+8*k,plot->fg_color, plot->bg_color, s);
	}
if(plot->logscale_x || plot->logscale_y)
	TRACE("Not implemented yet")

}

void draw_points_f(RGBPic *p, PLOT *plot, long color, float *xx, float *yy, long count, long step_x, long step_y)
{
long i;
double x,y,dx,dy;
long px,py;
dx=plot->px1-plot->px0;
if(plot->lx1>plot->lx0)dx/=plot->lx1-plot->lx0;
	else dx=0;
dy=plot->py1-plot->py0;
if(plot->ly1>plot->ly0)dy/=plot->ly1-plot->ly0;
	else dy=0;
for(i=0;i<count;i++){
	x=xx[i*step_x];
	if(x<plot->lx0)continue;
	if(x>plot->lx1)continue;
	y=yy[i*step_y];
	if(y<plot->ly0)continue;
	if(y>plot->ly1)continue;
	px=plot->px0+rint((x-plot->lx0)*dx);
	if(px<plot->px0)continue;
	if(px>plot->px1)continue;
	py=plot->py1-rint((y-plot->ly0)*dy);
	if(py<plot->py0)continue;
	if(py>plot->py1)continue;
	DRAW_POINT(p,px,py,color)
	DRAW_POINT(p,px-1,py,color)
	DRAW_POINT(p,px,py-1,color)
	DRAW_POINT(p,px+1,py,color)
	DRAW_POINT(p,px,py+1,color)
	}
}

void draw_points_d(RGBPic *p, PLOT *plot, long color, double *xx, double *yy, long count, long step_x, long step_y)
{
long i;
double x,y,dx,dy;
long px,py;
dx=plot->px1-plot->px0;
if(plot->lx1>plot->lx0)dx/=plot->lx1-plot->lx0;
	else dx=0;
dy=plot->py1-plot->py0;
if(plot->ly1>plot->ly0)dy/=plot->ly1-plot->ly0;
	else dy=0;
for(i=0;i<count;i++){
	x=xx[i*step_x];
	if(x<plot->lx0)continue;
	if(x>plot->lx1)continue;
	y=yy[i*step_y];
	if(y<plot->ly0)continue;
	if(y>plot->ly1)continue;
	px=plot->px0+rint((x-plot->lx0)*dx);
	if(px<plot->px0)continue;
	if(px>plot->px1)continue;
	py=plot->py1-rint((y-plot->ly0)*dy);
	if(py<plot->py0)continue;
	if(py>plot->py1)continue;
	DRAW_POINT(p,px,py,color)
	DRAW_POINT(p,px-1,py,color)
	DRAW_POINT(p,px,py-1,color)
	DRAW_POINT(p,px+1,py,color)
	DRAW_POINT(p,px,py+1,color)
	}
}

static inline long hue_z_to_color(float z0)
{
long color;
float H,S,B; /* HSB coordinates */
int r,g,b; /* R, G, B coordinates */
float c,s,r1;

B=0.7;
S=1.0;
H=(4.5*z0*M_PI-M_PI)/3.0;
		
c=cos(H)*0.5;
s=sin(H)*0.5;		
r1=sqrt(3)/2;
		
r=255*B*((1-S)+S*(0.5-0.5*s-r1*c));
g=255*B*((1-S)+S*(0.5+s));
b=255*B*((1-S)+S*(0.5-0.5*s+r1*c));
		
		
#define CLAMP(a)	{ if((a)<0)a=0; if ((a)>255) a=255; }
CLAMP(r)
CLAMP(g)
CLAMP(b)
		
color=COLOR(r,g,b);
return color;
}

PALETTE * make_hue_palette(long ncolors)
{
PALETTE *p;
long i;
p=do_alloc(1,sizeof(*p));
p->ncolors=ncolors;
p->color=do_alloc(p->ncolors, sizeof(*(p->color)));
for(i=0;i<p->ncolors;i++){
	p->color[i]=hue_z_to_color((1.0*i)/(p->ncolors-1));
	}

return p;
}

void free_palette(PALETTE *p)
{
free(p->color);
free(p);
}

DENSITY_MAP *make_density_map(int x_ppp, int y_ppp)
{
DENSITY_MAP *dm;
dm=do_alloc(1, sizeof(*dm));
dm->x_pixels_per_point=x_ppp;
dm->y_pixels_per_point=y_ppp;
dm->bg_color=COLOR(255,255,255);
dm->fg_color=COLOR(0,0,0);
dm->flip_x=1;
dm->flip_y=1;
dm->swap_xy=0;
dm->logscale_z=0;
dm->nbands=args_info.nbands_arg;
dm->dec_band_color=COLOR(127,127,127);
dm->palette=make_hue_palette(230);
return dm;
}

void free_density_map(DENSITY_MAP *dm)
{
free_palette(dm->palette);
free(dm);
}

void adjust_density_map_limits_f(DENSITY_MAP *dm, float *z, int count, int step, int replace)
{
float z0,z1;
int i;
if(count<1)return;
z0=*z;
z1=*z;

for(i=1;i<count;i++){
	if(z0>z[i*step])z0=z[i*step];
	if(z1<z[i*step])z1=z[i*step];
	}
if(replace){
	dm->lower_z=z0;
	dm->upper_z=z1;
	} else {
	if(dm->lower_z>z0)dm->lower_z=z0;
	if(dm->upper_z<z1)dm->upper_z=z1;
	}
}

void adjust_density_map_limits_d(DENSITY_MAP *dm, double *z, int count, int step, int replace)
{
double z0,z1;
int i;
if(count<1)return;
z0=*z;
z1=*z;

for(i=1;i<count;i++){
	if(z0>z[i*step])z0=z[i*step];
	if(z1<z[i*step])z1=z[i*step];
	}
if(replace){
	dm->lower_z=z0;
	dm->upper_z=z1;
	} else {
	if(dm->lower_z>z0)dm->lower_z=z0;
	if(dm->upper_z<z1)dm->upper_z=z1;
	}
}

void adjust_masked_density_map_limits_f(DENSITY_MAP *dm, float *z, int *mask, int count, int step, int replace)
{
float z0,z1;
int i, set;
set=0;

for(i=0;i<count;i++){
	if(mask[i*step]<0)continue;
	if(!set){
		z0=z[i*step];
		z1=z0;
		set=1;
		}
	if(z0>z[i*step])z0=z[i*step];
	if(z1<z[i*step])z1=z[i*step];
	}
if(!set)return;
if(replace){
	dm->lower_z=z0;
	dm->upper_z=z1;
	} else {
	if(dm->lower_z>z0)dm->lower_z=z0;
	if(dm->upper_z<z1)dm->upper_z=z1;
	}
}

void adjust_masked_density_map_limits_d(DENSITY_MAP *dm, double *z, int *mask, int count, int step, int replace)
{
double z0,z1;
int i, set;
set=0;

for(i=1;i<count;i++){
	if(mask[i*step]<0)continue;
	if(!set){
		z0=z[i*step];
		z1=z0;
		set=1;
		}
	if(z0>z[i*step])z0=z[i*step];
	if(z1<z[i*step])z1=z[i*step];
	}
if(!set)return;
if(replace){
	dm->lower_z=z0;
	dm->upper_z=z1;
	} else {
	if(dm->lower_z>z0)dm->lower_z=z0;
	if(dm->upper_z<z1)dm->upper_z=z1;
	}
}

static inline long z_to_color(PALETTE *p,float z0)
{
long ii;

#if 0   /* testing */
return hue_z_to_color(z0);
#endif

ii=floor(z0*p->ncolors);
if(ii>=p->ncolors)ii=p->ncolors-1;
if(ii<0)ii=0;

return p->color[ii];
}

void draw_density_map_key(RGBPic *p, DENSITY_MAP *dm, long x, long y)
{
int i,j;
float z;
long color;
dm->key_width=32+8*11;
dm->key_height=196;

RGBPic_clear_area(p,dm->bg_color,x,y,x+dm->key_width-1,y+dm->key_height-1);

for(i=0;i<dm->key_height;i++){
	z=1.0-(1.0*i)/(dm->key_height-1);
	color=z_to_color(dm->palette,z);
	for(j=0;j<20;j++)DRAW_POINT(p, x+j,y+i,color);
	}

/* print upper and lower values */
RGBPic_printf(p,x+23,y+1,dm->fg_color, dm->bg_color, "%g", dm->upper_z);
RGBPic_printf(p,x+23,y+dm->key_height-1-16,dm->fg_color, dm->bg_color, "%g", dm->lower_z);
}


void draw_density_map_f(RGBPic *p, long x, long y, 
	DENSITY_MAP *dm, 
	float *z, int x_count, int y_count, int step_x, int step_y)
{
int i,j,k,m;
float z0,dz;
long color,x0,y0;
int lz;
char tmp[10];

if(dm->flip_x){
	z=z+(x_count-1)*step_x;
	step_x=-step_x;
	}

if(dm->flip_y){
	z=z+(y_count-1)*step_y;
	step_y=-step_y;
	}

if(dm->swap_xy){
	i=x_count;
	j=step_x;
	x_count=y_count;
	step_x=step_y;
	y_count=i;
	step_y=j;
	}

dm->actual_width=dm->x_pixels_per_point*x_count;
dm->actual_height=dm->y_pixels_per_point*y_count;

lz=dm->logscale_z;


if(lz && (dm->lower_z<0)){
	fprintf(stderr,"numbers must be positive for logscale_z mode\n");
	lz=0;
	}

if(lz){
	dz=log10(dm->upper_z)-log10(dm->lower_z);
	if(dz<=0)dz=1;
	} else {
	dz=dm->upper_z-dm->lower_z;
	if(dz<=0)dz=1;
	}

RGBPic_clear_area(p,dm->bg_color,x,y,x+dm->actual_width-1,y+dm->actual_height-1);

#if 0
if(dm->nbands>=2){
	for(j=0;j<=dm->nbands;j++){
		RGBPic_draw_line(p, dm->dec_band_color, 
			x, y+(j*dm->actual_height-1)/dm->nbands,
			x+dm->actual_width-1, y+(j*dm->actual_height-1)/dm->nbands);
		}
	for(j=0;j<dm->nbands;j++){
		RGBPic_printf(p,x+2,y+((2*j+1)*(dm->actual_height-1))/(2*dm->nbands),
			dm->fg_color, dm->bg_color, "%d", j);
		}
	}
#endif

for(j=0;j<y_count;j++)
	for(i=0;i<x_count;i++){
		z0=z[i*step_x+j*step_y]-dm->lower_z;
		if(z0<-0.0001)continue;
		if(lz){
			z0=log10(z0)/dz; /* normalize so it is between 0 and 1 */
			} else {
			z0=z0/dz; /* normalize so it is between 0 and 1 */
			}
		if(z0>1.0001)continue;

		color=z_to_color(dm->palette,z0);

		x0=x+i*dm->x_pixels_per_point;
		y0=y+j*dm->y_pixels_per_point;
		for(m=0;m<dm->y_pixels_per_point;m++)
			for(k=0;k<dm->x_pixels_per_point;k++)
				DRAW_POINT(p,x0+k,y0+m,color);
		}
}

void layout_density_map_plot(RGBPic *p, DENSITY_MAP *dm, int max_x_count, int max_y_count)
{
long x;
RGBPic_clear_area(p,dm->bg_color,0,0,p->width-1,p->height-1);
x=p->width-120;
dm->x_pixels_per_point=(x-10)/max_x_count;
if(dm->x_pixels_per_point<1){
	fprintf(stderr,"Picture is too small to hold %d x points\n", max_x_count);
	return;
	}
dm->y_pixels_per_point=p->height/max_y_count;
if(dm->y_pixels_per_point<1){
	fprintf(stderr,"Picture is too small to hold %d y points\n", max_y_count);
	return;
	}
/* Make points square */
if(dm->x_pixels_per_point<dm->y_pixels_per_point){
	dm->y_pixels_per_point=dm->x_pixels_per_point;
	} else
if(dm->x_pixels_per_point>dm->y_pixels_per_point){
	dm->x_pixels_per_point=dm->y_pixels_per_point;
	}
draw_density_map_key(p, dm, x, 0);
}

void plot_single_density_map_f(RGBPic *p, DENSITY_MAP *dm, 
	float *z, int x_count, int y_count, int step_x, int step_y)
{
layout_density_map_plot(p, dm, x_count, y_count);
draw_density_map_f(p, 0, 0, dm, z, x_count, y_count, step_x, step_y);
}

void draw_density_map_d(RGBPic *p, long x, long y, 
	DENSITY_MAP *dm, 
	double *z, int x_count, int y_count, int step_x, int step_y)
{
int i,j,k,m;
double z0,dz;
long color,x0,y0;

if(dm->flip_x){
	z=z+(x_count-1)*step_x;
	step_x=-step_x;
	}

if(dm->flip_y){
	z=z+(y_count-1)*step_y;
	step_y=-step_y;
	}

if(dm->swap_xy){
	i=x_count;
	j=step_x;
	x_count=y_count;
	step_x=step_y;
	y_count=i;
	step_y=j;
	}

dm->actual_width=dm->x_pixels_per_point*x_count;
dm->actual_height=dm->y_pixels_per_point*y_count;

dz=dm->upper_z-dm->lower_z;
if(fabsf(dz)<=0)dz=1;

RGBPic_clear_area(p,dm->bg_color,x,y,x+dm->actual_width-1,y+dm->actual_height-1);

#if 0
if(dm->nbands>=2){
	for(j=0;j<=dm->nbands;j++){
		RGBPic_draw_line(p, dm->dec_band_color, 
			x, y+(j*dm->actual_height-1)/dm->nbands,
			x+dm->actual_width-1, y+(j*dm->actual_height-1)/dm->nbands);
		}
	for(j=0;j<dm->nbands;j++){
		RGBPic_printf(p,x+2,y+((2*j+1)*(dm->actual_height-1))/(2*dm->nbands),
			dm->fg_color, dm->bg_color, "%d", j);
		}
	}
#endif

if(dm->logscale_z)fprintf(stderr,"logscale_z is not implemented yet\n");

for(j=0;j<y_count;j++)
	for(i=0;i<x_count;i++){
		z0=z[i*step_x+j*step_y]-dm->lower_z;
		if(z0<-0.0001)continue;
		z0=z0/dz; /* normalize so it is between 0 and 1 */
		if(z0>1.0001)continue;

		color=z_to_color(dm->palette,z0);

		x0=x+i*dm->x_pixels_per_point;
		y0=y+j*dm->y_pixels_per_point;
		for(m=0;m<dm->y_pixels_per_point;m++)
			for(k=0;k<dm->x_pixels_per_point;k++)
				DRAW_POINT(p,x0+k,y0+m,color);
		}

}

void plot_single_density_map_d(RGBPic *p, DENSITY_MAP *dm, 
	double *z, int x_count, int y_count, int step_x, int step_y)
{
layout_density_map_plot(p, dm, x_count, y_count);
draw_density_map_d(p, 0, 0, dm, z, x_count, y_count, step_x, step_y);
}

void plot_sin_theta_f(RGBPic *p, long x, long y, DENSITY_MAP *dm, SKY_GRID *grid, float *z, int step)
{
int i,j,k,kk,m, shift;
float z0,dz;
long color,x0,y0;
int lz;
SIN_THETA_SKY_GRID_PRIV *priv;

if(strcmp(grid->name,"sin theta")){
	fprintf(stderr,"Cannot apply function %s to grid %s\n", __FUNCTION__,
		grid->name);
	exit(-1);
	}

priv=grid->grid_priv;

dm->actual_width=grid->max_n_ra*dm->x_pixels_per_point;
dm->actual_height=grid->max_n_dec*dm->y_pixels_per_point;

lz=dm->logscale_z;


if(lz && (dm->lower_z<0)){
	fprintf(stderr,"numbers must be positive for logscale_z mode\n");
	lz=0;
	}

if(lz){
	dz=log10(dm->upper_z)-log10(dm->lower_z);
	if(dz<=0)dz=1;
	} else {
	dz=dm->upper_z-dm->lower_z;
	if(dz<=0)dz=1;
	}
RGBPic_clear_area(p,dm->bg_color,x,y,x+dm->actual_width-1,y+dm->actual_height-1);

#if 0
if(dm->nbands>=2){
	for(j=0;j<=dm->nbands;j++){
		RGBPic_draw_line(p, dm->dec_band_color, 
			x, y+(j*dm->actual_height-1)/dm->nbands,
			x+dm->actual_width-1, y+(j*dm->actual_height-1)/dm->nbands);
		}
	for(j=0;j<dm->nbands;j++){
		RGBPic_printf(p,x+2,y+((2*j+1)*(dm->actual_height-1))/(2*dm->nbands),
			dm->fg_color, dm->bg_color, "%d", j);
		}
	}
#endif

kk=0;
for(j=0;j<priv->num_dec;j++){
	shift=(grid->max_n_ra-priv->num_ra[j])>>1;
	for(i=0;i<priv->num_ra[j];i++){
		z0=z[kk]-dm->lower_z;
		kk+=step;
		if(z0<-0.0001)continue;
		if(lz){
			z0=log10(z0)/dz; /* normalize so it is between 0 and 1 */
			} else {
			z0=z0/dz; /* normalize so it is between 0 and 1 */
			}
		if(z0>1.0001)continue;

		color=z_to_color(dm->palette,z0);

		if(dm->flip_x){
			x0=x+(shift+priv->num_ra[j]-i-1)*dm->x_pixels_per_point;
			} else {
			x0=x+(shift+i)*dm->x_pixels_per_point;
			}
		if(dm->flip_y){
			y0=y+(priv->num_dec-j-1)*dm->y_pixels_per_point;
			} else {
			y0=y+j*dm->y_pixels_per_point;
			}
		for(m=0;m<dm->y_pixels_per_point;m++)
			for(k=0;k<dm->x_pixels_per_point;k++)
				DRAW_POINT(p,x0+k,y0+m,color);
		}
	}
}

void plot_grid_f(RGBPic *p, SKY_GRID *grid, float *z, int step)
{
if(!strcmp(grid->name,"arcsin") || !strcmp(grid->name,"plain rectangular")){
	DENSITY_MAP *dm;
	RECT_SKY_GRID_PRIV *priv=grid->grid_priv;
	dm=make_density_map(1,1);
	adjust_masked_density_map_limits_f(dm, z, grid->band, grid->npoints, step, 1);
	plot_single_density_map_f(p, dm, z, priv->num_ra, priv->num_dec, step*priv->num_dec, step);
	free_density_map(dm);
	return;
	}
if(!strcmp(grid->name,"sin theta")){
	DENSITY_MAP *dm;
	dm=make_density_map(1,1);
	adjust_masked_density_map_limits_f(dm, z, grid->band, grid->npoints, step, 1);
	layout_density_map_plot(p, dm, grid->max_n_ra, grid->max_n_dec);
	plot_sin_theta_f(p, 0, 0, dm, grid, z, step);
	return;
	}

fprintf(stderr,"**ERROR: do not know how to plot sky grid \"%s\"\n",grid->name);
exit(-1);
}

void plot_sin_theta_d(RGBPic *p, long x, long y, DENSITY_MAP *dm, SKY_GRID *grid, double *z, int step)
{
int i,j,k,kk,m, shift;
double z0,dz;
long color,x0,y0;
int lz;
SIN_THETA_SKY_GRID_PRIV *priv;

if(strcmp(grid->name,"sin theta")){
	fprintf(stderr,"Cannot apply function %s to grid %s\n", __FUNCTION__,
		grid->name);
	exit(-1);
	}

priv=grid->grid_priv;

dm->actual_width=grid->max_n_ra*dm->x_pixels_per_point;
dm->actual_height=grid->max_n_dec*dm->y_pixels_per_point;


lz=dm->logscale_z;


if(lz && (dm->lower_z<0)){
	fprintf(stderr,"numbers must be positive for logscale_z mode\n");
	lz=0;
	}

if(lz){
	dz=log10(dm->upper_z)-log10(dm->lower_z);
	if(dz<=0)dz=1;
	} else {
	dz=dm->upper_z-dm->lower_z;
	if(dz<=0)dz=1;
	}
RGBPic_clear_area(p,dm->bg_color,x,y,x+dm->actual_width-1,y+dm->actual_height-1);

kk=0;
for(j=0;j<priv->num_dec;j++){
	shift=(grid->max_n_ra-priv->num_ra[j])>>1;
	for(i=0;i<priv->num_ra[j];i++){
		z0=z[kk]-dm->lower_z;
		kk++;
		if(z0<-0.0001)continue;
		if(lz){
			z0=log10(z0)/dz; /* normalize so it is between 0 and 1 */
			} else {
			z0=z0/dz; /* normalize so it is between 0 and 1 */
			}
		if(z0>1.0001)continue;

		color=z_to_color(dm->palette,z0);

		if(dm->flip_x){
			x0=x+(shift+priv->num_ra[j]-i-1)*dm->x_pixels_per_point;
			} else {
			x0=x+(shift+i)*dm->x_pixels_per_point;
			}
		if(dm->flip_y){
			y0=y+(priv->num_dec-j-1)*dm->y_pixels_per_point;
			} else {
			y0=y+j*dm->y_pixels_per_point;
			}
		for(m=0;m<dm->y_pixels_per_point;m++)
			for(k=0;k<dm->x_pixels_per_point;k++)
				DRAW_POINT(p,x0+k,y0+m,color);
		}
	}
}

void plot_grid_d(RGBPic *p, SKY_GRID *grid, double *z, int step)
{
if(!strcmp(grid->name,"arcsin") || !strcmp(grid->name,"plain rectangular")){
	DENSITY_MAP *dm;
	RECT_SKY_GRID_PRIV *priv=grid->grid_priv;
	dm=make_density_map(1,1);
	adjust_masked_density_map_limits_d(dm, z, grid->band, grid->npoints, step, 1);
	plot_single_density_map_d(p, dm, z, priv->num_ra, priv->num_dec, step*priv->num_dec, step);
	free_density_map(dm);
	return;
	}
if(!strcmp(grid->name,"sin theta")){
	DENSITY_MAP *dm;
	dm=make_density_map(1,1);
	adjust_masked_density_map_limits_d(dm, z, grid->band, grid->npoints, step, 1);
	layout_density_map_plot(p, dm, grid->max_n_ra, grid->max_n_dec);
	plot_sin_theta_d(p, 0, 0, dm, grid, z, step);
	return;
	}
fprintf(stderr,"**ERROR: do not know how to plot sky grid \"%s\"\n",grid->name);
exit(-1);
}
