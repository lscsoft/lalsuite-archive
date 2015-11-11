#include <stdio.h>
#include <stdlib.h>
#define _GNU_SOURCE
#include <pthread.h>
#include <string.h>
#include "global.h"
#include "jobs.h"
#include "cmdline.h"

extern FILE *LOG;
extern struct gengetopt_args_info args_info;

/* Low level interface */

#ifdef LINUX_THREADS

pthread_t  * thread_id;
int max_threads=-1;
int num_threads=0;
int threads_started=0;
MUTEX thread_num_mutex;
CONDITION thread_not_needed;
cpu_set_t cpuset;

int cpu_count(void)
{
FILE *cpuinfo;
int count;
char s[1000];
cpuinfo=fopen("/proc/cpuinfo", "r");
if(cpuinfo==NULL) {
	perror("Opening /proc/cpuinfo");
	fprintf(stderr, "Could not determine CPU count, defaulting to 2\n");
	fprintf(LOG, "Could not determine CPU count, defaulting to 2\n");
	return 2;
	}

count=0;
while(!feof(cpuinfo)) {
	fgets(s, 999, cpuinfo);
	if(!strncmp(s, "processor", 9))count++;
	}
fclose(cpuinfo);
if(count<1) {
	fprintf(stderr, "Could not determine CPU count, defaulting to 2\n");
	fprintf(LOG, "Could not determine CPU count, defaulting to 2\n");
	return 2;
	}
return count;
}

void *thread_cruncher(long i)
{

while(1) {
	//fprintf(stderr, "Thread %d waiting\n", i);
	wait_for_more_jobs();
	thread_mutex_lock(&thread_num_mutex);
	while(i>=num_threads)thread_cond_wait(&thread_not_needed, &thread_num_mutex);
	thread_mutex_unlock(&thread_num_mutex);
	//fprintf(stderr, "Thread %d running\n", i);
	pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset);
	while(do_single_job(i));
	}

return NULL; /* make compiler happy */
}

void init_threads(int mt)
{
max_threads=mt-1;
if(max_threads<0) {
	max_threads=cpu_count()-1;
	}
if(max_threads<0) max_threads=0;

/* Store initial cpu set as Intel OpenMP compiler changes affinities for non-owned threads */
pthread_getaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset);

thread_id=do_alloc(max_threads, sizeof(*thread_id));
num_threads=0;
threads_started=0;
thread_cond_init(&thread_not_needed);
thread_mutex_init(&thread_num_mutex);
set_concurrency(1);
fprintf(stderr, "maximum threads: %d\n", max_threads+1);
fprintf(LOG, "maximum threads: %d\n", max_threads+1);
}

void thread_mutex_init(MUTEX *mutex)
{
pthread_mutex_init(mutex, NULL);
}

void thread_cond_init(CONDITION *cond)
{
pthread_cond_init(cond, NULL);
}

int get_max_threads(void)
{
return(max_threads+1);
}

int get_concurrency(void)
{
int a;
thread_mutex_lock(&thread_num_mutex);
a=num_threads;
thread_mutex_unlock(&thread_num_mutex);
/* return one more because we have num_threads running in addition
   to the main one */
return (a+1);
}

void set_concurrency(int num)
{
long i;
thread_mutex_lock(&thread_num_mutex);
//fprintf(stderr, "set_concurrency(%d) max_threads=%d\n", num, max_threads);
num_threads=num-1;
if(num_threads>max_threads)num_threads=max_threads;
for(i=threads_started;i<num_threads;i++)
	if(pthread_create(&(thread_id[i]), NULL,(void * (*)(void *))  thread_cruncher, (void *)i)<0){
		fprintf(stderr,"Could not spawn thread %ld:", i);
		perror("");
		} // else fprintf(stderr, "Spawned thread %ld\n", i);
if(num_threads>threads_started)
	threads_started=num_threads;
thread_cond_broadcast(&thread_not_needed);
thread_mutex_unlock(&thread_num_mutex);
}

#else

void init_threads(int mt)
{
fprintf(stderr, "maximum threads: 1 (no_thread)\n");
fprintf(LOG, "maximum threads: 1 (no_thread)\n");
}

int get_max_threads(void)
{
return 1;
}

int get_concurrency(void)
{
return 1;
}

void set_concurrency(int num)
{
}

#endif

/* High level interface - jobs */

typedef struct {
	void *data;
	JOB_CRUNCHER job;
	} JOB;

MUTEX jobs_mutex;
CONDITION wait_for_more_jobs_condition;
CONDITION all_done_condition;
JOB *job=NULL;
long jobs_free=0;
long jobs_size=1000;
long next_job_to_do=0;
long jobs_submitted=0;
long jobs_done=0;
long jobs_reset_count=0;

void init_jobs(void)
{
jobs_free=0;
jobs_size=100;
job=do_alloc(jobs_size, sizeof(JOB));
jobs_submitted=0;
jobs_done=0;
next_job_to_do=0;
thread_mutex_init(&jobs_mutex);
thread_cond_init(&wait_for_more_jobs_condition);
thread_cond_init(&all_done_condition);
}

void expand_jobs(void)
{
JOB *p;
jobs_size=2*jobs_size+10;
p=do_alloc(jobs_size, sizeof(JOB));
if(jobs_free>0)memcpy(p, job, jobs_free*sizeof(JOB));
free(job);
job=p;
}

void submit_job(JOB_CRUNCHER f, void *data)
{
set_concurrency(get_max_threads());

thread_mutex_lock(&jobs_mutex);
if(jobs_free>=jobs_size)expand_jobs();
jobs_submitted++;
job[jobs_free].data=data;
job[jobs_free].job=f;
jobs_free++;
//thread_cond_broadcast(&wait_for_more_jobs_condition);
thread_cond_signal(&wait_for_more_jobs_condition);
thread_mutex_unlock(&jobs_mutex);
}

int do_single_job(int thread_id)
{
long i, count, stop;
thread_mutex_lock(&jobs_mutex);
i=next_job_to_do;
if(i>=jobs_free){
	thread_mutex_unlock(&jobs_mutex);
	return 0;
	}
if(jobs_free>i+1000)count=5;
	else count=1;
next_job_to_do+=count;
thread_mutex_unlock(&jobs_mutex);
//fprintf(stderr, "Running job %ld on thread %d\n", i, thread_id);
stop=count+i;
for(;i<stop;i++)
	job[i].job(thread_id, job[i].data);
thread_mutex_lock(&jobs_mutex);
jobs_done+=count;
if(jobs_done==jobs_submitted) {
	//fprintf(stderr, "jobs_done=%ld tid=%d i=%ld\n", jobs_done, thread_id, i);
	thread_cond_broadcast(&all_done_condition);
	set_concurrency(1);
	thread_mutex_unlock(&jobs_mutex);
	return 0;
	}
thread_mutex_unlock(&jobs_mutex);
return 1;
}

int all_done(void)
{
thread_mutex_lock(&jobs_mutex);
if(jobs_submitted==jobs_done){
	jobs_free=0;
	next_job_to_do=0;
	thread_mutex_unlock(&jobs_mutex);
	return 1;
	}
if(next_job_to_do==jobs_free){
	thread_cond_wait(&all_done_condition, &jobs_mutex);
	}
thread_mutex_unlock(&jobs_mutex);
return 0;
}

void wait_for_all_done(void)
{
thread_mutex_lock(&jobs_mutex);
if(jobs_submitted==jobs_done) {
	jobs_free=0;
	next_job_to_do=0;
	thread_mutex_unlock(&jobs_mutex);
	return;
	}
thread_cond_wait(&all_done_condition, &jobs_mutex);
thread_mutex_unlock(&jobs_mutex);
}

void wait_for_more_jobs(void)
{
thread_mutex_lock(&jobs_mutex);
if(jobs_submitted==jobs_done){
	thread_cond_wait(&wait_for_more_jobs_condition, &jobs_mutex); 
	}
thread_mutex_unlock(&jobs_mutex);
}

void wake_crunchers(void)
{
thread_mutex_lock(&jobs_mutex);
thread_cond_broadcast(&wait_for_more_jobs_condition);
thread_mutex_unlock(&jobs_mutex);
}

void reset_jobs_done_ratio(void)
{
jobs_reset_count=jobs_done;
}

float jobs_done_ratio(void)
{
return((1.0*(jobs_done-jobs_reset_count))/(jobs_submitted-jobs_reset_count));
}

void print_jobs_stats(void)
{
fprintf(stderr,"jobs_submitted=%ld jobs_done=%ld\n", jobs_submitted, jobs_done);
}

