#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "rastermagic.h"

void *do_alloc(long a, long b)
{
void *r;
r=calloc(a,b);
while(r==NULL){
	fprintf(stderr,"Could not allocate %ld chunks of %ld bytes each (%ld bytes total)\n",a,b,a*b);
	sleep(10);
	r=calloc(a,b);
	}
return r;
}

char *output_dir="";

int main(int argc, char *argv[])
{
RGBPic *p;
PLOT *plot;

p=make_RGBPic(800,600);
RGBPic_clear_area(p,COLOR(60,80,50),0,0,p->width,p->height);

plot=make_plot(640,480);
plot->lower_x=3.1415;
plot->upper_x=12423;
plot->lower_y=2.34E-10;
plot->upper_y=43.5E-8;
draw_grid(p,plot,50,50);
free_plot(plot);

//RGBPic_printf(p,30,30,COLOR(255,240,210),COLOR(1,2,3),"Hello all !\n line 2 %d %p",argc,p);
//RGBPic_printf_v(p,300,300,COLOR(255,240,210),COLOR(1,2,3),"Hello all 2 !\n line 2 %d %p",argc,p);
RGBPic_draw_line(p,COLOR(40,255,80), 40,50,60,70);
RGBPic_draw_line(p,COLOR(250,255,80), 30,50,60,30);

RGBPic_dump_ppm("test1.ppm",p);
free_RGBPic(p);
return 0;
}
