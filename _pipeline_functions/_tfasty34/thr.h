
/***************************************/
/* thread global variable declarations */
/***************************************/

/* $Name: fa_34_26_5 $ - $Id: thr.h,v 1.2 1999/12/30 01:26:59 wrp Exp $ */

#ifndef MAX_WORKERS
#define MAX_WORKERS 2
#endif
#define NUM_WORK_BUF 2*MAX_WORKERS

#ifndef XTERNAL
struct buf_head *worker_buf[NUM_WORK_BUF];  /* pointers to full buffers */
struct buf_head *reader_buf[NUM_WORK_BUF];  /* pointers to empty buffers */

/* protected by worker_mutex/woker_cond_var */
int worker_buf_workp, worker_buf_readp; /* indices into full-buffers ptrs */
int num_worker_bufs;
int reader_done;

/* protected by reader_mutex/reader_cond var */
int reader_buf_workp, reader_buf_readp; /* indices into empty-buffers ptrs */
int num_reader_bufs;

/* protected by start_mutex/start_cont_var */
int start_thread=1;        /* start-up predicate, 0 starts */
#else
extern struct buf_head *worker_buf[];
extern struct buf_head *reader_buf[];
extern int num_worker_bufs, reader_done, num_reader_bufs;
extern int worker_buf_workp, worker_buf_readp;
extern int reader_buf_workp, reader_buf_readp;

extern int start_thread;
#endif

