#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>
#include <signal.h>
#include <wait.h>
#include <string.h>
#include <math.h>

typedef float REAL4;

void * do_alloc(long a, long b)
{
void *r;
r=calloc(a,b);
while(r==NULL){
	fprintf(stderr,"Could not allocate %ld units of %ld bytes each (%ld bytes total)\n", a, b, a*b);
	sleep(1);
	r=calloc(a,b);
	}
return r;
}

/* This function returns the value of exponential distribution function
   with parameter lambda. Needed for Kolmogorov-Smirnov test 
*/
REAL4 exponential_distribution(REAL4 x, REAL4 lambda)
{
if(x<=0.0)return 0.0;
return (1-exp(-lambda*x));
}

int c_r4(REAL4 *a, REAL4 *b)
{
if(*a<*b)return -1;
if(*a>*b)return 1;
return 0;
}

/* this function modifies numbers in the array */
REAL4 sum_positive_r4(REAL4 *array, long count)
{
int i,c;
qsort(array, count, sizeof(REAL4), c_r4);
c=count;
while(c>1){
	for(i=0;i<c;i+=2){
		array[i>>1]=array[i]+array[i+1];
		}
	if(c & 1){
		array[(c>>1)]=array[c-1];
		c=(c>>1)+1;
		} else {
		c=(c>>1);
		}
	}
return array[0];
}

/* this function modifies numbers in the array 
   It assumes that input numbers are positive - as far as finding
   average is concerned 
   
   The parameter lambda is estimated from the same data - so strictly speaking this is
   a modified KS test. 
   
   Need to find out how this affects critical values for Dn.
   Similar situation for Normal distribution results in much lower critical values.
   (See Bickel and Doksum)
   
   */
   
/* Critical values of KS test (from Bickel and Doksum). Does not apply directly 

alpha=0.01
n	10	20	30	40	50	60	80	n>80
	.489	.352	.290	.252	.226	.207	.179	1.628/(sqrt(n)+0.12+0.11/sqrt(n))
	
alpha=0.05
n	10	20	30	40	50	60	80	n>80
	.409	.294	.242	.210	.188	.172	.150	1.358/(sqrt(n)+0.12+0.11/sqrt(n))

*/
void Kolmogorov_Smirnov_test(REAL4 *array, long count, REAL4 *average, REAL4 *median, REAL4* lambda, REAL4 *Dn, REAL4 outlier_level, long *outliers)
{
int i,c;
REAL4 *backup, a,b;
qsort(array, count, sizeof(REAL4), c_r4);
backup=alloca(sizeof(REAL4)* count);
memcpy(backup, array, sizeof(REAL4)*count);
if(count & 1){
	*median=array[count>>1];
	} else {
	*median=(array[count>>1]+array[(count>>1)-1])/2.0;
	}
c=count;
while(c>1){
	for(i=0;i<c;i+=2){
		array[i>>1]=array[i]+array[i+1];
		}
	if(c & 1){
		array[(c>>1)]=array[c-1];
		c=(c>>1)+1;
		} else {
		c=(c>>1);
		}
	}
*average=array[0]/count;
if(*average==0.0){
	/* this is definitely not an exponential distribution 
	   However 0.0 average would mean - no noise and no signal 
	   so we assign Dn and outliers a special value of -1
	*/
	*Dn=-1.0;
	*outliers=-1;
	return;
	}
a=-2.0;
#if 0
*lambda=1/(*average);
#else
*lambda=log(2.0)/(*median);
#endif
*outliers=0;
for(i=0;i<count;i++){
	b=(((i+1)*1.0)/count)-exponential_distribution(backup[i], *lambda);
	if(b>outlier_level)(*outliers)++;
	if(a<b)a=b;
	b=exponential_distribution(backup[i], *lambda)-((i*1.0)/count);
	if(b>outlier_level)(*outliers)++;
	if(a<b)a=b;
	}
*Dn=a;
}

#define BUFFER_SIZE (256*1024)

#define MODE_TEXT		1
#define MODE_BINARY		2
#define MODE_BINARY_SWAPPED	3

struct {
	unsigned char * filename;
	FILE *fin;
	unsigned char *buffer;
	int fd[2];
	REAL4 weight;
	pid_t pid;
	int mode;
	} * input;

long num_inputs=0;
int precision=6;

REAL4 *data=NULL, *data_squared=NULL, *ks_data=NULL;

void open_pipe(char *command, char *filename, int num)
{
unsigned char *pipe_cmd;
fclose(input[num].fin);
input[num].fin=NULL;
if(pipe(input[num].fd)<0){
	fprintf(stderr,"Could not open pipe while starting \"%s %s\":", command, filename);
	perror("");
	exit(-1);
	}
input[num].pid=fork();
if(input[num].pid==-1){
	fprintf(stderr,"Could not fork while starting \"%s %s\":", command, filename);
	perror("");
	exit(-1);
	}
if(input[num].pid==0){
	/* child */
	/* close input */
	close(0);
	/* redirect output into pipe */
	dup2(input[num].fd[1], 1);
	/* close old descriptor */
	close(input[num].fd[1]);
	/* make command */
	pipe_cmd=alloca(strlen(filename)+strlen(command)+10);
	sprintf(pipe_cmd, "%s %s", command, filename);
	/* start command */
	execl("/bin/sh", "sh", "-c", pipe_cmd, NULL);
	/* we should never reach this */
	fprintf(stderr,"Error starting \"sh -c %s\":", pipe_cmd);
	perror("");
	close(1);
	exit(-1);
	}
/* close descriptor for writing - we do not need it */
close(input[num].fd[1]);
input[num].fin=fdopen(input[num].fd[0], "r");
}

void attach_input_file(char *filename, int num)
{
unsigned char s[20];
int count=0;
input[num].fin=fopen(filename, "r");
while(input[num].fin==NULL){
	if(count<15){
		sleep(1);
		input[num].fin=fopen(filename, "r");
		count++;
		continue;
		}		
	fprintf(stderr,"Could not open file \"%s\":", filename);
	perror("");
	exit(-1);
	}
fgets( s, 20, input[num].fin);
fseek(input[num].fin, 0, SEEK_SET);

/* use magic to tell compressed files */
    /* gzip */
if((s[0]==31)&&(s[1]==139)){
	open_pipe("gzip -c -d", filename, num);
	} else
    /* bzip2 */
if((s[0]=='B')&&(s[1]=='Z')&&(s[2]=='h')){
	open_pipe("bzip2 -c -d", filename, num);
	}
}

void kill_decompressors(void)
{
int i;
for(i=0;i<num_inputs;i++){
	if(input[i].pid>0){
		kill(input[i].pid, SIGTERM);
		}
	}
sleep(1);
for(i=0;i<num_inputs;i++){
	if(input[i].pid>0){
		kill(input[i].pid, SIGKILL);
		}
	}
}

void sig_handler(int a)
{
fprintf(stderr,"Received signal %d\n", a);
fclose(stdout);
kill_decompressors();
exit(0);
}

void sig_child_handler(int a)
{
int i, status;
fprintf(stderr,"Decompressor exited:");
for(i=0;i<num_inputs;i++){
	if((input[i].pid>0)){
		waitpid(input[i].pid, &status, WNOHANG);
		if(WIFEXITED(status)){
			fprintf(stderr,"\n\t[%d pid=%d file=%s]", i, input[i].pid, input[i].filename);
			input[i].pid=0;
			}
		}
	}
fprintf(stderr,"\n");
}

inline void swap_real4(REAL4 *a)
{
unsigned char *c;
unsigned char b,d;
c=(unsigned char *)a;
b=c[0];
d=c[1];
c[0]=c[3];
c[1]=c[2];
c[2]=d;
c[3]=b;
}

unsigned char *output_buffer;
REAL4 total_weight=0.0;

#define LINE_SIZE 200000

int main(int argc, char *argv[])
{
int i;
unsigned char *s;
REAL4 sum;
REAL4 sum_squared, average, rayleigh, variance, ks_average, ks_Dn, outlier_level, ks_median, ks_lambda;
long outliers;
long r;
unsigned char *endianness=(unsigned char *)&r;

/* determine our endianness */
r=('B'<<24)|('i'<<16)|('g'<<8)|('E');

output_buffer=do_alloc(BUFFER_SIZE, 1);
setvbuf(stdout, output_buffer, _IOFBF, BUFFER_SIZE);

total_weight=0.0;

num_inputs=(argc-1)/2;
input=do_alloc(num_inputs, sizeof(*input));
for(i=0;i<num_inputs;i++){
	input[i].filename=argv[2*i+2];
	input[i].fin=NULL;
	input[i].buffer=do_alloc(BUFFER_SIZE,1);
	input[i].fd[0]=-1;
	input[i].fd[1]=-1;
	input[i].weight=atof(argv[2*i+1]);
	input[i].pid=0;
	input[i].mode=MODE_TEXT;
	attach_input_file(argv[2*i+2],i);
	setvbuf(input[i].fin, input[i].buffer, _IOFBF, BUFFER_SIZE);
	}

data=do_alloc(num_inputs, sizeof(*data));
data_squared=do_alloc(num_inputs, sizeof(*data_squared));
ks_data=do_alloc(num_inputs, sizeof(*ks_data));
/* compute total weight */
for(i=0;i<num_inputs;i++){
	data[i]=input[i].weight;
	}
total_weight=sum_positive_r4(data, num_inputs);
		
s=do_alloc(LINE_SIZE, 1);

printf("#\n");
printf("#         Tag: WEIGHTED SUM\n");
for(i=0;i<num_inputs;i++){
	printf("#      file: %s  weight %g\n", input[i].filename, input[i].weight);
	}
printf("#\n");
printf("# Average Median Lambda KS KS_count\n#\n");

signal(SIGINT, sig_handler);
signal(SIGTERM, sig_handler);
signal(SIGPIPE, sig_handler);
/* just for debugging: */
/* signal(SIGCHLD, sig_child_handler); */

outlier_level=1.358/(sqrt(num_inputs)+0.12+0.11/sqrt(num_inputs));

while(1){
	for(i=0;i<num_inputs;i++){
		data[i]=0.0;
		if(feof(input[i].fin) || ferror(input[i].fin)){
			fprintf(stderr,"End of input on file \"%s\"\n", input[i].filename);
			fclose(stdout);
			kill_decompressors();
			exit(0);
			}
		switch(input[i].mode){
			case MODE_TEXT:
				fgets(s, LINE_SIZE, input[i].fin);
				while((s[0]=='#')||!s[0]){
					printf("# %s: %s", input[i].filename, s+1);
					if(feof(input[i].fin) || ferror(input[i].fin)){
						fprintf(stderr,"End of input on file \"%s\"\n", input[i].filename);
						fclose(stdout);
						kill_decompressors();
						exit(0);
						}
					fgets(s, LINE_SIZE, input[i].fin);
					}
				if(s[0]!='b'){
					data[i]=strtod(s,NULL);
					break;
					}
				fgets(s, LINE_SIZE, input[i].fin);
				if(!strncmp(endianness, s, 4))
					input[i].mode=MODE_BINARY;
					else 
					input[i].mode=MODE_BINARY_SWAPPED;
			case MODE_BINARY:
			case MODE_BINARY_SWAPPED:
				fread(&(data[i]), sizeof(data[i]), 1, input[i].fin);
				if(input[i].mode==MODE_BINARY_SWAPPED)swap_real4(&(data[i]));
				break;
			default:
				fprintf(stderr,"Internal error, input[%d=%s].mode=%d\n", i, input[i].filename, input[i].mode);
				kill_decompressors();
				exit(-1);
				break;
			}
		ks_data[i]=data[i];
		data_squared[i]=data[i]*data[i]*input[i].weight;
		data[i]*=input[i].weight;
		}
	Kolmogorov_Smirnov_test(ks_data, num_inputs, &ks_average, &ks_median, &ks_lambda, &ks_Dn, outlier_level, &outliers);
	printf("%.*g %.*1$g %.*1$g %.*1$g %ld\n", precision, ks_average, ks_median, ks_lambda, ks_Dn, outliers);
	#if 0
	sum=sum_positive_r4(data, num_inputs);
	sum_squared=sum_positive_r4(data_squared, num_inputs);
	average=sum/total_weight;
	if(num_inputs>1)
		variance=((sum_squared*total_weight-sum*sum)*num_inputs)/(total_weight*total_weight*(num_inputs-1));
		else variance=0;
	if(fabs(sum)>0)
		rayleigh=sqrt(variance)/average;
		else
		rayleigh=-1;
	printf("%.*g %.*1$g %.*1$g %.*1$g\n", precision, sum, average, variance, rayleigh);
	#endif
	}
return 0;
}
