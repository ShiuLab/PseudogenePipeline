/* copyright (c) 1996, 1997, 1998, 1999 William R. Pearson and the
   U. of Virginia */

/* $Name: fa_34_26_5 $ - $Id: work_thr.c,v 1.23 2007/04/26 18:33:20 wrp Exp $ */

/* work_thr.c - threaded worker */

/* modified 21-Oct-1998 to work with reverse complement for DNA */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include <signal.h>

#include "defs.h"		/* various constants */
#include "mw.h"			/* defines beststr */
#include "structs.h"
#include "param.h"		/* pstruct, thr_str, buf_head, rstruct */

/***************************************/
/* thread global variable declarations */
/***************************************/

#define XTERNAL
#include "thr.h"
#undef XTERNAL

void alloc_pam (int, int, struct pstruct *);
int **alloc_pam2p(int, int);
void revcomp(unsigned char *seq, int n, int *c_nt);
#ifdef WIN32
void pthread_exit(void *);
#else
void THR_EXIT(void *);
#endif

/* functions getting/sending buffers to threads (thr_sub.c) */
extern void wait_thr(void);
extern int get_wbuf(struct buf_head **cur_buf, int max_work_buf);
extern void put_wbuf(struct buf_head *cur_buf, int max_work_buf);

/* dropxx.c functions */
extern void init_work (unsigned char *aa0, int n0,
		       struct pstruct *ppst, void **f_arg);

extern void do_work (unsigned char *aa0, int n0, unsigned char *aa1, int n1,
		     int frame,
		     struct pstruct *ppst, void *f_str, int qr_flg,
		     struct rstruct *rst);

extern void close_work (unsigned char *, int, struct pstruct *, void **);

extern void irand(int);
extern int shuffle(unsigned char *, unsigned char *, int);
extern int wshuffle(unsigned char *, unsigned char *, int, int, int *);
extern void qshuffle(unsigned char *aa0, int n0, int nm0);
extern void free_pam2p(int **);

void
work_thread (struct thr_str *work_info)
{
  struct buf_head *cur_buf;
  struct buf_str *cur_buf_p;
  struct buf_str *p_rbuf;
  unsigned char *aa1s;
  int cur_cnt;
  int my_worker;
  int i, j, npam, n0, nm0;
  int ix_score, debug_lib, zsflag, zs_win, do_shuffle, ieven=0;
  int frame;

  struct rstruct rrst;
  struct pstruct my_pst, *my_ppst;
  unsigned char *aa0[6], *aa0s;
  void *f_str[6], *qf_str;

  my_worker = work_info->worker;

  wait_thr();	/* wait for start_thread predicate to drop to  0 */

  /* do init_work */

  /* let each thread have its own copy of the query */
  n0 = work_info->n0;
  nm0 = work_info->nm0;

  if ((aa0[0]=(unsigned char *)calloc((size_t)n0+2,sizeof(unsigned char)))
      ==NULL) {
    fprintf(stderr," cannot allocate aa00[%d] for worker %d\n",
	    n0, my_worker);
    exit(1);
  }
  *aa0[0]='\0';
  aa0[0]++;
  memcpy(aa0[0],work_info->aa0,n0+1);

  /* make certain that all but 0 have their own copy of pst */
  if (my_worker) {
    my_ppst = &my_pst;
    memcpy(my_ppst,work_info->ppst,sizeof(struct pstruct));

    alloc_pam(MAXSQ, MAXSQ, my_ppst);

    npam = (my_pst.ext_sq_set) ? my_pst.nsqx : my_pst.nsq;

    for (i=0; i<=npam; i++) {
      for (j=0; j<=npam; j++) {
	my_pst.pam2[0][i][j] = work_info->ppst->pam2[0][i][j];
	my_pst.pam2[1][i][j] = work_info->ppst->pam2[1][i][j];
      }
    }

    if (work_info->ppst->pam_pssm && work_info->ppst->pam2p[0]) {
      my_ppst->pam2p[0] = alloc_pam2p(n0,npam);
      my_ppst->pam2p[1] = alloc_pam2p(n0,npam);
      for (i=0; i<n0; i++) {
	for (j=0; j <= npam; j++) {
	  my_pst.pam2p[0][i][j] = work_info->ppst->pam2p[0][i][j];
	  my_pst.pam2p[1][i][j] = work_info->ppst->pam2p[1][i][j];
	}
      }
    }
  }
  else my_ppst=work_info->ppst;

  /* note that aa[5,4,3,2] are never used, but are provided so that frame
     can range from 0 .. 5; likewise for f_str[5..2] */

  aa0[5] = aa0[4] = aa0[3] = aa0[2] = aa0[1] = aa0[0];
  init_work (aa0[0], n0, my_ppst, &f_str[0]);

  f_str[5] = f_str[4] = f_str[3] = f_str[2] = f_str[1] = f_str[0];

  if (work_info->qframe == 2) {
    if ((aa0[1]=(unsigned char *)calloc((size_t)n0+2,sizeof(unsigned char)))==NULL) {
      fprintf(stderr," cannot allocate aa01[%d] for worker %d\n",
	    n0, my_worker);
    }
    *aa0[1]='\0';
    aa0[1]++;
    memcpy(aa0[1],work_info->aa0,n0+1);
    revcomp(aa0[1],n0,my_ppst->c_nt);
    init_work (aa0[1], n0, my_ppst, &f_str[1]);
  }

  if (work_info->qshuffle) {
    if ((aa0s=(unsigned char *)calloc(n0+2,sizeof(char)))==NULL) {
      fprintf(stderr,"cannot allocate aa0s[%d]\n",n0+2);
      exit(1);
    }
    *aa0s='\0';
    aa0s++;
    memcpy(aa0s,aa0[0],n0);
    qshuffle(aa0s,n0,nm0);
    init_work (aa0s, n0, my_ppst, &qf_str);
  }

  ix_score = my_ppst->score_ix;
  debug_lib = my_ppst->debug_lib;
  zsflag = my_ppst->zsflag;
  zs_win = my_ppst->zs_win;

  if (zsflag >= 10) {
    if((aa1s=calloc(work_info->max_tot+1,sizeof(char))) == NULL) {
      fprintf(stderr,"unable to allocate shuffled library sequence\n");
    }
    else {
      *aa1s=0;
      aa1s++;
      do_shuffle =1;
      irand(0);
    }
  }
  else {do_shuffle = 0;}

  /* main work loop */
  while (get_wbuf(&cur_buf,work_info->max_work_buf)) {

    cur_cnt = cur_buf->buf_cnt;
    if (cur_cnt == -1) break;
    cur_buf_p = cur_buf->buf;

    while (cur_cnt--) { /* count down the number of sequences */
      p_rbuf = cur_buf_p++;   /* step through each sequence */
      p_rbuf->rst.score[0] = p_rbuf->rst.score[1] = p_rbuf->rst.score[2] = 0;
      frame = p_rbuf->frame;

#ifdef DEBUG
      if (debug_lib) {
	if (frame >= 2) fprintf(stderr,"* frame: %d\n",frame);
	for (i=0; i<p_rbuf->n1; i++)
	  if (p_rbuf->aa1b[i]>=my_ppst->nsqx) {
	    fprintf(stderr,
		    "%s residue[%d/%d] %d range (%d)\n",
		    p_rbuf->libstr,i,p_rbuf->n1,p_rbuf->aa1b[i],my_ppst->nsqx);
	    p_rbuf->aa1b[i]=0;
	    p_rbuf->n1=i-1;
	    break;
	  }
      }
#endif

      do_work (aa0[frame], n0, p_rbuf->aa1b, p_rbuf->n1, frame,
	       my_ppst, f_str[frame], 0, &p_rbuf->rst);
      
      if (work_info->qshuffle) {
	do_work(aa0s,n0,p_rbuf->aa1b, p_rbuf->n1, frame,
		my_ppst, qf_str, 1, &rrst);
	p_rbuf->qr_score = rrst.score[ix_score];
	p_rbuf->qr_escore = rrst.escore;
      }

      if (do_shuffle) {
	if (zs_win > 0) wshuffle(p_rbuf->aa1b,aa1s,p_rbuf->n1,zs_win,&ieven);
	else shuffle(p_rbuf->aa1b,aa1s,p_rbuf->n1);

	do_work (aa0[frame], n0, aa1s, p_rbuf->n1, frame,
		 my_ppst, f_str[frame], 0, &rrst);
	p_rbuf->r_score = rrst.score[ix_score];
	p_rbuf->r_escore = rrst.escore;
      }
    }

    cur_buf->have_results = 1;

    put_wbuf(cur_buf,work_info->max_work_buf);

  } /* end main while */

  close_work(aa0[0], n0, my_ppst, &f_str[0]);
  free(aa0[0]-1);
  if (work_info->qframe == 2) {
    close_work(aa0[1], n0, my_ppst, &f_str[1]);
    free(aa0[1]-1);
  }

  if (do_shuffle) free(aa1s-1);

  if (my_worker) {
    free(my_pst.pam2[1][0]);
    free(my_pst.pam2[0][0]);
    free(my_pst.pam2[1]);
    free(my_pst.pam2[0]);
    
    if (my_pst.pam_pssm) {
      free_pam2p(my_pst.pam2p[0]);
      free_pam2p(my_pst.pam2p[1]);
    }
  }

#ifdef WIN32
  pthread_exit(&work_info->status);
#else
  THR_EXIT(&work_info->status);
#endif

}  /* end work_thread */

