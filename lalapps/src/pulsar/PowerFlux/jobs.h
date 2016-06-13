#ifndef __JOBS_H__
#define __JOBS_H__

/* Low level thread interface */

/* the following interfaces pthread library */
#ifdef LINUX_THREADS

#include <pthread.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef pthread_mutex_t MUTEX;
typedef pthread_cond_t CONDITION;

#define thread_mutex_lock pthread_mutex_lock
#define thread_mutex_unlock pthread_mutex_unlock
void thread_mutex_init(MUTEX *mutex);

#define thread_cond_wait  pthread_cond_wait
#define thread_cond_broadcast  pthread_cond_broadcast
#define thread_cond_signal pthread_cond_signal
void thread_cond_init(CONDITION *cond);

#ifdef __cplusplus
}
#endif

#else

typedef int MUTEX;
typedef int CONDITION;

#define thread_mutex_lock(mutex)
#define thread_mutex_unlock(mutex)
#define thread_mutex_init(mutex)

#define thread_cond_wait(cond, mutex)
#define thread_cond_broadcast(cond)
#define thread_cond_signal(cond)
#define thread_cond_init(cond)


#endif


#ifdef __cplusplus
extern "C" {
#endif

/* Inialize thread engine. 
   (max_threads-1) specifies how many threads to spawn *in* *addition* to the main one,
   so max_threads refers to total number of threads.
   Pass max_threads=-1 to guess number of available cores using /proc/cpuinfo */

void init_threads(int max_threads);
int get_max_threads(void);
int get_concurrency(void);
void set_concurrency(int num);

/* High level thread (job) interface */

/* By a silly convention thread id is from 0 to get_max_threads()-2, 
   with -1 passed for the main thread */

typedef void (*JOB_CRUNCHER)(int thread_id, void *data);

void init_jobs(void);
void submit_job(JOB_CRUNCHER job, void *data);
int do_single_job(int thread_id);
int all_done(void);
void wait_for_more_jobs(void);
void wake_crunchers(void);
void reset_jobs_done_ratio(void);
float jobs_done_ratio(void);
void print_jobs_stats(void);
void wait_for_all_done(void);

#ifdef __cplusplus
}
#endif


#endif
