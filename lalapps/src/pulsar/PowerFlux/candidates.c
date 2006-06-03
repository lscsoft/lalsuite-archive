#include <stdio.h>
#include <stdlib.h>
/* We need this define to get NAN values */
#define __USE_ISOC99
#include <math.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <string.h>

#include "global.h"
#include "rastermagic.h"
#include "cmdline.h"
#include "hookup.h"
#include "grid.h"
#include "polarization.h"
#include "statistics.h"
#include "dataset.h"
#include "candidates.h"

extern DATASET *datasets;
extern int d_free;


extern POLARIZATION_RESULTS *polarization_results;
extern int nlinear_polarizations, ntotal_polarizations;

extern SKY_GRID *fine_grid, *patch_grid;
extern FILE *LOG;

extern double spindown;
extern double resolution;

extern int subinstance;
extern char *subinstance_name;

extern int nbins, first_bin, side_cut, useful_bins;

INT64 spindown_start;

CANDIDATE *candidate=NULL;
int candidate_free=0;
int candidate_size=-1;

extern struct gengetopt_args_info args_info;

char s[20000];

SUM_TYPE *max_dx=NULL;
short *polarization_index=NULL;
int *max_dx_order=NULL;
int *dec_order=NULL;
int *inverse_dec_order=NULL;
int *ra_order=NULL;
int *inverse_ra_order=NULL;
int *max_dx_local_map=NULL;

float noise_floor;

int max_dx_compare(int *a1, int *a2) 
{
if(max_dx[*a1]<max_dx[*a2])return 1;
if(max_dx[*a1]>max_dx[*a2])return -1;
return 0;
}

int dec_compare(int *a1, int *a2) 
{
float ra_dist, dec_dist;
dec_dist=fine_grid->latitude[*a1]-fine_grid->latitude[*a2];
if(dec_dist>0)return 1;
if(dec_dist<0)return -1;

ra_dist=cos(fine_grid->latitude[*a1]+fine_grid->latitude[*a2])*(fine_grid->longitude[*a1]-fine_grid->longitude[*a2]);
if(ra_dist<0)return -1;
if(ra_dist>0)return 1;

return 0;
}

int ra_compare(int *a1, int *a2) 
{
float ra_dist, dec_dist;

ra_dist=cos(0.5*(fine_grid->latitude[*a1]+fine_grid->latitude[*a2]))*(fine_grid->longitude[*a1]-fine_grid->longitude[*a2]);
if(ra_dist<-1.5*resolution)return -1;
if(ra_dist>1.5*resolution)return 1;

dec_dist=fine_grid->latitude[*a1]-fine_grid->latitude[*a2];
if(dec_dist>0)return 1;
if(dec_dist<0)return -1;

return 0;
}

int min_i=0, max_i=0;

static int sweep_points_ra(int index, int mark) 
{
int j,m, count=0;
for(j=inverse_dec_order[index]+1;j<fine_grid->npoints;j++) {
	m=dec_order[j];
	if(max_dx_local_map[m]>=0)break;
	if(max_dx[m]<noise_floor)break;
	if(max_dx[m]>max_dx[dec_order[j-1]]+args_info.snr_precision_arg)break;
	max_dx_local_map[m]=mark;
	if(m<min_i)min_i=m;
	if(m>max_i)max_i=m;
	count++;
	}
for(j=inverse_dec_order[index]-1;j>=0;j--) {
	m=dec_order[j];
	if(max_dx_local_map[m]>=0)break;
	if(max_dx[m]<noise_floor)break;
	if(max_dx[m]>max_dx[dec_order[j+1]]+args_info.snr_precision_arg)break;
	max_dx_local_map[m]=mark;
	if(m<min_i)min_i=m;
	if(m>max_i)max_i=m;
	count++;
	}
return count;
}

static int sweep_points_dec(int index, int mark) 
{
int j,m, count=0;
for(j=inverse_ra_order[index]+1;j<fine_grid->npoints;j++) {
	m=ra_order[j];
	if(max_dx_local_map[m]>=0)break;
	if(max_dx[m]<noise_floor)break;
	if(abs(fine_grid->latitude[m]-fine_grid->latitude[ra_order[j-1]])>0.1)break;
	if(max_dx[m]>max_dx[ra_order[j-1]]+args_info.snr_precision_arg)break;
	max_dx_local_map[m]=mark;
	if(m<min_i)min_i=m;
	if(m>max_i)max_i=m;
	count++;
	}
for(j=inverse_ra_order[index]-1;j>=0;j--) {
	m=ra_order[j];
	if(max_dx_local_map[m]>=0)break;
	if(max_dx[m]<noise_floor)break;
	if(abs(fine_grid->latitude[m]-fine_grid->latitude[ra_order[j-1]])>0.1)break;
	if(max_dx[m]>max_dx[ra_order[j+1]]+args_info.snr_precision_arg)break;
	max_dx_local_map[m]=mark;
	if(m<min_i)min_i=m;
	if(m>max_i)max_i=m;
	count++;
	}
return count;
}

/* slow but sure.. */
static int sweep_neighbourhood(int index, int mark)
{
int j,m, count=0;
float fudge=resolution*3.1;
float cosfudge=1.0-cos(fudge);

for(j=inverse_dec_order[index]+1;j<fine_grid->npoints;j++) {
	m=dec_order[j];

	if(abs(fine_grid->latitude[m]-fine_grid->latitude[index])>fudge)break;

	if(max_dx_local_map[m]>=0)continue;
	
	if(max_dx[m]>max_dx[index])continue;
	

	if(fast_spherical_distance(fine_grid->longitude[m], fine_grid->latitude[m],
				fine_grid->longitude[index], fine_grid->latitude[index])<cosfudge) {
		
		max_dx_local_map[m]=mark;
		count++;
		}
	}
	
for(j=inverse_dec_order[index]-1;j>=0;j--) {
	m=dec_order[j];

	if(abs(fine_grid->latitude[m]-fine_grid->latitude[index])>fudge)break;

	if(max_dx_local_map[m]>=0)continue;

	if(max_dx[m]>max_dx[index])continue;
	

	if(fast_spherical_distance(fine_grid->longitude[m], fine_grid->latitude[m],
				fine_grid->longitude[index], fine_grid->latitude[index])<cosfudge) {
	
		max_dx_local_map[m]=mark;
		count++;
		}
	}
return count;
}

void sweep_points(int index, int mark)
{
int i, count;
unsigned char *p;

max_i=index;
min_i=index;
max_dx_local_map[index]=mark;

p=do_alloc(fine_grid->npoints, sizeof(*p));
for(i=0;i<fine_grid->npoints;i++)*p=0;

while(1) {
	count=0;
	for(i=min_i;i<=max_i;i++) {
		if(p[i])continue;
		if(max_dx_local_map[i]==mark){
			count+=sweep_neighbourhood(i, mark);
			p[i]=1;
			}
		}
	fprintf(stderr, "mark=%d count=%d\n", mark, count);
	if(count==0){ 
		free(p);
		return;
		}
	}
}

void expand_candidates(void)
{
CANDIDATE *p;
candidate_size=candidate_size*2+10;
p=do_alloc(candidate_size, sizeof(*p));
if(candidate_free>0)memcpy(p, candidate, candidate_free*sizeof(*p));
free(candidate);
candidate=p;
}

void dump_candidate(int i)
{
float doppler;
FILE *fout;
char s[20];
DATASET *d;
int index=candidate[i].point_index;
float spindown=candidate[i].spindown;
int b0, b1, bin_shift, j, k, b;
float a, a_plus, a_cross;
POLARIZATION *pl;

snprintf(s, 20, "cand_%d.dat", i);
fprintf(stderr, "%s\n", s);
fout=fopen(s, "w");
fprintf(fout, "dataset\tk\tgps\tdoppler\tspindown\ta_plus\ta_cross\texpTMedian");
for(b=0;b<nbins-2*side_cut;b++)
	fprintf(fout, "\tb%d", b);
fprintf(fout, "\n");
/* loop over datasets */
for(j=0;j<d_free;j++) {
	d=&(datasets[j]);
	pl=&(d->polarizations[candidate[i].polarization_index]);	
	
	/* process single SFT */
	for(k=0;k<d->free;k++){
		/* Get amplitude response */
		a_plus=F_plus(k, fine_grid, index, pl->AM_coeffs);
		a_cross=F_plus(k, fine_grid, index, pl->conjugate->AM_coeffs);

		doppler=fine_grid->e[0][index]*d->detector_velocity[3*k+0]+
			fine_grid->e[1][index]*d->detector_velocity[3*k+1]+
			fine_grid->e[2][index]*d->detector_velocity[3*k+2];

		bin_shift=-rint((first_bin+nbins*0.5)*doppler+1800.0*spindown*(d->gps[k]-spindown_start));

		b0=side_cut-bin_shift;
		b1=(nbins-side_cut)-bin_shift;

		if((b0<0)||(b1>nbins)){
			fprintf(stderr,"Working frequency range obscured by bin_shift shift: bin_shift=%d index=%d i=%d\n",
				bin_shift, index, i);
			exit(-1);
			}
		fprintf(fout, "%s\t%d\t%lld\t%g\t%g\t%f\t%f\t%g", 
			d->name, 
			k, 
			d->gps[k], 
			doppler, 
			spindown*(d->gps[k]-spindown_start),
			a_plus,
			a_cross,
			d->expTMedians[k]);
		for(b=b0; b<b1;b++) {
			a=d->bin[k*nbins+b].im;
			if(a>0)
				fprintf(fout, "\t%g+%gi", d->bin[k*nbins+b].re, a);
				else 
				fprintf(fout, "\t%g%gi", d->bin[k*nbins+b].re, a);
			}

		fprintf(fout, "\n");
		}
	}
fclose(fout);
}

void coherence_score(int i)
{
float doppler;
FILE *fout;
DATASET *d;
int index=candidate[i].point_index;
float spindown=candidate[i].spindown;
float frequency=candidate[i].frequency;
int b0, b1, j, k, b;
float a, score, a_plus, a_cross;
POLARIZATION *pl;
double weight, total_weight, *response, *mismatch, f, x, y;
int window=5;
int *signal_bin;
COMPLEX *p0, *p1, *p2, w, z;
struct {
	double re;
	double im;
	} *phase_sum;

total_weight=0;
phase_sum=do_alloc(2*window+1, sizeof(*phase_sum));
for(b=0;b<2*window+1;b++) {
	phase_sum[b].re=0.0;
	phase_sum[b].im=0.0;
	}

/* loop over datasets */
for(j=0;j<d_free;j++) {
	d=&(datasets[j]);
	pl=&(d->polarizations[candidate[i].polarization_index]);	
	
	response=do_alloc(d->free, sizeof(*response));
	signal_bin=do_alloc(d->free, sizeof(*signal_bin));
	mismatch=do_alloc(d->free, sizeof(*mismatch));
	
	
	/* process single SFT */
	for(k=0;k<d->free;k++){
		/* Get amplitude response */
		a_plus=F_plus(k, fine_grid, index, pl->AM_coeffs);
		a_cross=F_plus(k, fine_grid, index, pl->conjugate->AM_coeffs);
		
		response[k]=pl->plus_factor*a_plus*a_plus+pl->cross_factor*a_cross*a_cross;

		doppler=fine_grid->e[0][index]*d->detector_velocity[3*k+0]+
			fine_grid->e[1][index]*d->detector_velocity[3*k+1]+
			fine_grid->e[2][index]*d->detector_velocity[3*k+2];

		f=frequency+frequency*doppler+spindown*(d->gps[k]-spindown_start);

		signal_bin[k]=rint(1800.0*f-first_bin);
		mismatch[k]=1800.0*f-first_bin-signal_bin[k];

		if(!k || (signal_bin[k]<b0))b0=signal_bin[k];
		if(!k || (signal_bin[k]>b1))b1=signal_bin[k];

		}
	fprintf(stderr, "b0=%d b1=%d\n", b0, b1);

	for(k=1;k< (d->free-1);k++){
		/* skip SFT with asymmetric gaps */
		if(d->gps[k]-d->gps[k-1]!=d->gps[k+1]-d->gps[k])continue;
	
		weight=((response[k-1]*d->expTMedians[k-1]+2*response[k]*d->expTMedians[k]+response[k+1]*d->expTMedians[k+1])*d->weight);
		total_weight+=weight;
	
		p0=&(d->bin[(k-1)*nbins+signal_bin[k-1]-window]);
		p1=&(d->bin[k*nbins+signal_bin[k]-window]);
		p2=&(d->bin[(k+1)*nbins+signal_bin[k+1]-window]);

		for(b=0; b< (2*window+1); b++) {
			x=sqrt(p0[b].re*p0[b].re+p0[b].im*p0[b].im);
			y=sqrt(p1[b].re*p1[b].re+p1[b].im*p1[b].im);

			w.re=(p0[b].re*p1[b].re+p0[b].im*p1[b].im)/(x*y);
			w.im=(p0[b].re*p1[b].im-p0[b].im*p1[b].re)/(x*y);

			x=sqrt(p2[b].re*p2[b].re+p2[b].im*p2[b].im);
			z.re=(p2[b].re*p1[b].re+p2[b].im*p1[b].im)/(x*y);
			z.im=(p2[b].re*p1[b].im-p2[b].im*p1[b].re)/(x*y);
			
			phase_sum[b].re+=weight*(w.re*z.re-w.im*z.im);
			phase_sum[b].im+=weight*(w.im*z.re+w.re*z.im);
			}
		}
		
	free(response);
	free(signal_bin);
	free(mismatch);
	response=NULL;
	signal_bin=NULL;
	mismatch=NULL;
	}

fprintf(stderr, "total_weight=%g\n", total_weight);

score=0;
for(b=0;b<2*window+1;b++) {
	phase_sum[b].re/=total_weight;
	phase_sum[b].im/=total_weight;
	a=sqrt(phase_sum[b].re*phase_sum[b].re+phase_sum[b].im*phase_sum[b].im);
	if(i==0)fprintf(stderr, " %1.2f", a);
	if(a>score)score=a;
	}
if(!i)fprintf(stderr, "\n");

candidate[i].coherence_score=score;
free(phase_sum);
}

void identify_candidates(void)
{
RGBPic *p;
PLOT *plot;
SUM_TYPE a;
int i,k,m;
float *a_f;


candidate_free=0;
candidate_size=args_info.max_candidates_arg;
candidate=do_alloc(candidate_size, sizeof(*candidate));

max_dx=do_alloc(fine_grid->npoints, sizeof(*max_dx));
max_dx_local_map=do_alloc(fine_grid->npoints, sizeof(*max_dx_local_map));
max_dx_order=do_alloc(fine_grid->npoints, sizeof(*max_dx_order));

dec_order=do_alloc(fine_grid->npoints, sizeof(*dec_order));
inverse_dec_order=do_alloc(fine_grid->npoints, sizeof(*inverse_dec_order));

ra_order=do_alloc(fine_grid->npoints, sizeof(*ra_order));
inverse_ra_order=do_alloc(fine_grid->npoints, sizeof(*inverse_ra_order));

a_f=do_alloc(fine_grid->npoints, sizeof(*a_f));
polarization_index=do_alloc(fine_grid->npoints, sizeof(*polarization_index));

if(fine_grid->max_n_dec<800){
	p=make_RGBPic(fine_grid->max_n_ra*(800/fine_grid->max_n_dec)+140, fine_grid->max_n_dec*(800/fine_grid->max_n_dec));
	} else 
	p=make_RGBPic(fine_grid->max_n_ra+140, fine_grid->max_n_dec);	

plot=make_plot(p->width, p->height);

for(i=0;i<fine_grid->npoints;i++) {
	max_dx[i]=0;
	polarization_index[i]=-1;
	max_dx_order[i]=i;
	dec_order[i]=i;
	ra_order[i]=i;
	max_dx_local_map[i]=-1;
	}

for(i=0;i<fine_grid->npoints;i++) {
	for(k=0;k<ntotal_polarizations;k++){
		a=polarization_results[k].skymap.max_dx[i];
		if(a>max_dx[i]) {
			max_dx[i]=a;
			polarization_index[i]=k;
			}
		}
	}
qsort(max_dx_order, fine_grid->npoints, sizeof(*max_dx_order), max_dx_compare);
qsort(dec_order, fine_grid->npoints, sizeof(*dec_order), dec_compare);
qsort(ra_order, fine_grid->npoints, sizeof(*ra_order), ra_compare);

for(i=0;i<fine_grid->npoints;i++) {
	inverse_dec_order[dec_order[i]]=i;
	inverse_ra_order[ra_order[i]]=i;
	}

/* Take noise floor to be bottom 1/3 of max_dx */
noise_floor=max_dx[max_dx_order[(fine_grid->npoints*2)/3]];
fprintf(LOG, "noise_floor: %f\n", noise_floor);
fprintf(stderr, "noise_floor: %f\n", noise_floor);
	
/* Find all local maxima, starting with largest.. */
for(i=0;i<fine_grid->npoints;i++) {
	/* if(i>50)break; /* !!!!!!!!!!!!!!!! TODO - REMOVE !!!!!!!!!!!!!!!!!! */
	k=max_dx_order[i];
	/* skip masked points.. */
	if(polarization_index[k]<0)continue;
	
	/* is the point marked already ? */
	if(max_dx_local_map[k]>=0){
		sweep_neighbourhood(k, max_dx_local_map[k]);
		continue;
		}
	
	/* record found maximum.. */
	candidate[candidate_free].point_index=k;
	
	/* mark and sweep nearby points */
	max_dx_local_map[k]=candidate_free;
	sweep_neighbourhood(k, candidate_free);

	candidate_free++;	
	if(candidate_free>=candidate_size)expand_candidates();
	}

fprintf(LOG, "candidates: polarization_index rank point_index domain_size max_dx frequency ra dec spindown a_plus a_cross weight_ratio skyband coherence_score total\n");
fprintf(stderr, "candidates: polarization_index rank point_index domain_size max_dx frequency ra dec spindown a_plus a_cross weight_ratio skyband coherence_score total\n");


for(i=0;i<candidate_free;i++) {
	candidate[i].domain_size=0;
	}
	
for(i=0;i<fine_grid->npoints;i++) {
	k=max_dx_local_map[i];
	if(k>=0)candidate[k].domain_size++;
	}

m=0;	
for(i=0;i<candidate_free;i++) {
	k=candidate[i].point_index;
	fprintf(stderr, "k=%d\n", k);
	candidate[i].polarization_index=polarization_index[k];
	fprintf(stderr, "polarization_index=%d\n", polarization_index[k]);
	candidate[i].rank=i;
	candidate[i].max_dx=max_dx[k];
	if(candidate[i].max_dx>noise_floor)m++;
	candidate[i].ra=fine_grid->longitude[k];
	candidate[i].dec=fine_grid->latitude[k];
	candidate[i].spindown=spindown;
	candidate[i].a_plus=polarization_results[polarization_index[k]].plus_proj;
	candidate[i].a_cross=polarization_results[polarization_index[k]].cross_proj;
	candidate[i].frequency=polarization_results[polarization_index[k]].skymap.freq_map[k];
	candidate[i].weight_ratio=1.0-polarization_results[polarization_index[k]].skymap.max_sub_weight[k]/polarization_results[polarization_index[k]].skymap.total_weight[k];
	candidate[i].skyband=fine_grid->band[k];
	
	coherence_score(i);

	fprintf(LOG, "candidate: %d %d %d %d %f %f %f %f %g %f %f %f %d %f %d\n",
		candidate[i].polarization_index,
		candidate[i].rank,
		candidate[i].point_index,
		candidate[i].domain_size,
		candidate[i].max_dx,
		candidate[i].frequency,
		candidate[i].ra,
		candidate[i].dec,
		candidate[i].spindown,
		candidate[i].a_plus,
		candidate[i].a_cross,
		candidate[i].weight_ratio,
		candidate[i].skyband,
		candidate[i].coherence_score,
		candidate_free);
	fprintf(stderr, "candidate: %d %d %d %d %f %f %f %f %g %f %f %f %d %f %d\n",
		candidate[i].polarization_index,
		candidate[i].rank,
		candidate[i].point_index,
		candidate[i].domain_size,
		candidate[i].max_dx,
		candidate[i].frequency,
		candidate[i].ra,
		candidate[i].dec,
		candidate[i].spindown,
		candidate[i].a_plus,
		candidate[i].a_cross,
		candidate[i].weight_ratio,
		candidate[i].skyband,
		candidate[i].coherence_score,
		candidate_free);


	//if(i<5)dump_candidate(i);
	}

fprintf(LOG, "high_candidates_count: %d\n", m);
fprintf(stderr, "high_candidates_count: %d\n", m);

fprintf(LOG, "candidates_count: %d\n", candidate_free);
fprintf(stderr, "candidates_count: %d\n", candidate_free);

snprintf(s, 19999, "%smax_dx.png", subinstance_name);
if(clear_name_png(s)){
	plot_grid_f(p, fine_grid, max_dx, 1);
	RGBPic_dump_png(s, p);
	}
snprintf(s, 19999, "%smax_dx.dat", subinstance_name);
dump_floats(s, max_dx, fine_grid->npoints, 1);

for(i=0;i<fine_grid->npoints;i++) {
	a_f[max_dx_order[i]]=fine_grid->npoints-i;
	}
snprintf(s, 19999, "%smax_dx_order.png", subinstance_name);
if(clear_name_png(s)){
	plot_grid_f(p, fine_grid, a_f, 1);
	RGBPic_dump_png(s, p);
	}
snprintf(s, 19999, "%smax_dx_order.dat", subinstance_name);
dump_ints(s, max_dx_order, fine_grid->npoints, 1);

for(i=0;i<fine_grid->npoints;i++) {
	a_f[dec_order[i]]=i;
	}
snprintf(s, 19999, "%sdec_order.png", subinstance_name);
if(clear_name_png(s)){
	plot_grid_f(p, fine_grid, a_f, 1);
	RGBPic_dump_png(s, p);
	}

for(i=0;i<fine_grid->npoints;i++) {
	a_f[ra_order[i]]=i;
	}
snprintf(s, 19999, "%sra_order.png", subinstance_name);
if(clear_name_png(s)){
	plot_grid_f(p, fine_grid, a_f, 1);
	RGBPic_dump_png(s, p);
	}

for(i=0;i<fine_grid->npoints;i++) {
	a_f[i]=max_dx_local_map[i]<0 ? -1 : candidate_free-max_dx_local_map[i];
	}
snprintf(s, 19999, "%smax_dx_local_map.png", subinstance_name);
if(clear_name_png(s)){
	plot_grid_f(p, fine_grid, a_f, 1);
	RGBPic_dump_png(s, p);
	}
snprintf(s, 19999, "%smax_dx_local_map.dat", subinstance_name);
dump_ints(s, max_dx_local_map, fine_grid->npoints, 1);

for(i=0;i<fine_grid->npoints;i++) {
	a_f[i]=polarization_index[i];
	}
snprintf(s, 19999, "%smax_dx_polarization_index.png", subinstance_name);
if(clear_name_png(s)){
	plot_grid_f(p, fine_grid, a_f, 1);
	RGBPic_dump_png(s, p);
	}
snprintf(s, 19999, "%smax_dx_polarization_index.dat", subinstance_name);
dump_shorts(s, polarization_index, fine_grid->npoints, 1);


free(a_f);
a_f=NULL;
free_plot(plot);
free_RGBPic(p);

free(max_dx);
free(max_dx_local_map);
free(max_dx_order);
max_dx=NULL;
max_dx_local_map=NULL;
max_dx_order=NULL;

free(dec_order);
free(inverse_dec_order);
dec_order=NULL;
inverse_dec_order=NULL;

free(ra_order);
free(inverse_ra_order);
ra_order=NULL;
inverse_ra_order=NULL;

free(candidate);
free(polarization_index);
candidate=NULL;
polarization_index=NULL;
}
