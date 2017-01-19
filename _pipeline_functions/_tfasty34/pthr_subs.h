

/* $Name: fa_34_26_5 $ - $Id: pthr_subs.h,v 1.1.1.1 1999/10/22 20:56:01 wrp Exp $ */


#include <pthread.h>

/* error macro for thread calls */

#define check(status,string) \
   if (status != 0) {fprintf(stderr,string); \
     fprintf(stderr,"%s\n",strerror(status)); } /* error macro */

/*
#define check(status,string) \
     if (status == -1) perror(string)  */ /* error macro for thread calls */


#ifndef XTERNAL
pthread_t threads[MAX_WORKERS];

/* mutex stuff */

pthread_mutex_t reader_mutex;      /* empty buffer pointer structure lock */
pthread_mutex_t worker_mutex;      /* full buffer pointer structure lock */

/* condition variable stuff */

pthread_cond_t reader_cond_var;    /* condition variable for reader */
pthread_cond_t worker_cond_var;    /* condition variable for workers */

pthread_mutex_t start_mutex;       /* start-up synchronisation lock */
pthread_cond_t start_cond_var;     /* start-up synchronisation condition variable */

extern pthread_t threads[];

/* mutex stuff */

extern pthread_mutex_t reader_mutex;
extern pthread_mutex_t worker_mutex;

/* condition variable stuff */

extern pthread_cond_t reader_cond_var;
extern pthread_cond_t worker_cond_var;

extern pthread_mutex_t start_mutex;
extern pthread_cond_t start_cond_var;
extern int start_thread;

#endif
